#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

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

        -n=*|--ncores=*)
        ncoresOpt="${i#*=}"
        shift # Number of Cores to Use
        ;;
        -m=*|--memory=*)
        memoryOpt="${i#*=}"
        shift # Per Core Memory Requirement
        ;;

    # Invalid Argument Handler

        *)
        # invalid option
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

# Defaults if No Arguments Passed
ncoresDef="12"
memoryDef="6G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

# Max Reads in RAM
maxReads=$((allocMemory * 250000))

format_status "PARAMETERS:
Picard Directory    = $PICARD
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Max Reads in Memory = $maxReads"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
tmpDir=$PIPELINE_HOME/$subset/tmp

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
        format_status "Command:\njava -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir"
        java -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir

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
        format_status "Command:\nsamtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
        samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam
        
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
        format_status "Command:\nsamtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam"
        samtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam
        
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
        format_status "Command:\nsamtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam"
        samtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam

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
        format_status "Command:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir

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
        format_status "Command:\nsamtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
        samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam

        # Update State on Exit
        put_state $? $state
        format_status "BAM Indexing Complete"

    fi

fi

format_status "Done"

