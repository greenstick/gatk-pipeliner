#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

# Load ~/.bash_profile if Not Found
if [ -z $PIPELINE_HOME ]; then
    echo "Reloading ~/.bash_profile"
    source ~/.bash_profile
fi

# Assign Arguments
for i in "$@"
    do case $i in

    # Standard Arguments

        -f=*|--fileprefix=*)
        fileprefix="${i#*=}"
        shift # Access & Write Files With This Prefix
        ;;
        -s=*|--subset=*)
        subset="${i#*=}"
        shift # Access & Write Files With This Subset
        ;;
        -c=*|--condition=*)
        condition="${i#*=}"
        shift # Access & Write Files With This Condition
        ;;
        -x=*|--experiment=*)
        experiment="${i#*=}"
        shift # Access & Write Files With This Experiment
        ;;
        -p=*|--parameters=*)
        parameters="${i#*=}"
        shift # Access & Write Files With This Parameter Set
        ;;

    # Optional Arguments With Defaults

        -r=*|--reads=*)
        readsOpt="${i#*=}"
        shift # n Reads Per GB or Memory
        ;;
        -n=*|--ncores=*)
        ncoresOpt="${i#*=}"
        shift # Number of Cores to Use
        ;;
        -m=*|--memory=*)
        memoryOpt="${i#*=}"
        shift # Per Core Memory Requirement
        ;;
        -d=*|--debug=*)
        debugOpt="${i#*=}"
        shift # Trigger Debugging Available in Tools
        ;;

    # Directory Cleanup (Voids All Other Parameters)

        --clean)
        cleanOpt=true
        ;;

    # Invalid Argument Handler

        *)
        # invalid option
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

# Defaults if No Arguments Passed
ncoresDef="16"
memoryDef="6G"
cleanDef=false
readsDef=150000
debugDef=false

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
clean=${cleanOpt:-$cleanDef}
reads=${readsOpt:-$readsDef}
debug=${debugOpt:-$debugDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

# Max Reads in RAM - 200,000 per GB
maxReads=$((allocMemory * $reads))

printf "\nPARAMETERS: 
Picard Directory    = $PICARD
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Max Reads in Memory = $maxReads
Debug               = $debug
Do Cleanup          = $clean
\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
tmpDir=$PIPELINE_HOME/$subset/tmp

# Tool Specific Debugging - Picard
verbosity="INFO"

if $debug; then 
    verbosity="DEBUG"
fi

format_status "Running Mark Duplicates Script"

# If Norealignment, Get Data From Download Directory & Skip Sorting
if [ "$experiment" = "norealign" ]; then

    #
    # Mark Duplicates
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1"
    if !(has_state $state); then

        format_status "MarkDuplicates Start"
        # Define Command
        call="java -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir \
        VERBOSITY=$verbosity"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "Mark Duplicates Complete"

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2"
    if !(has_state $state); then

        format_status "Indexing BAM"
        # Define Command
        call="samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
        # Print & Call
        format_status "Command:\n$call"
        eval $call


        # Update State on Exit
        put_state $? $state
        format_status "BAM Indexing Complete"

    fi

# Otherwise, Grab Data From Merged Alignment Directory & Sort
else

    # 
    # Sort BAM
    # 

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1"
    if !(has_state $state); then

        format_status "Sorting BAM"
        # Define Command
        call="samtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "Sort BAM Complete"

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2"
    if !(has_state $state); then

        format_status "Indexing BAM Output"
        # Define Command
        call="samtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "BAM Indexing Complete"

    fi

    #
    # Mark Duplicates
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:3"
    if !(has_state $state); then

        format_status "MarkDuplicates Start"
        # Define Command
        call="java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "Mark Duplicates Complete"

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:4"
    if !(has_state $state); then

        format_status "Indexing BAM Output"
        # Define Command
        call="samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "BAM Indexing Complete"

    fi

fi

format_status "Done"

