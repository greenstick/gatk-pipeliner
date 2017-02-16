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
memoryDef="8G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

printf "\nPARAMETERS: 
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
tmpDir=$PIPELINE_HOME/$subset/tmp

format_status "Running BAM to FastQ Script"

#
# Shuffle & Split BAM
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition:BAMTOFASTQ:1"
if !(has_state $state); then

    format_status "Splitting Merged BAM"
    # Define Command
    call="bam splitBam -i $dataDir/downloaded/$fileprefix.$subset.$condition.bam -o $dataDir/downloaded/split/$fileprefix.$subset.$condition"
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    put_state $? $state

fi

#
# Sort BAMs
#

# State Management Delegated to Substates
format_status "Running Samtools Sort BAM"
# Retrieve Files
files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

for file in $files
    # In Parallel
    do ( 
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        substate="$fileprefix.$subset.$condition.$readgroup:BAMTOFASTQ:2"
        
        # State Check - Run Block if it Has Not Already Been Executed Successfully
        if !(has_state $substate); then
            
            # Call Bam to FastQ
            # Define Command
            call="samtools sort -n -m $memory -@ $ncores $dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam > $dataDir/downloaded/split/$fileprefix.$subset.$condition.sorted.$readgroup.bam"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Check for failed parallel call
            put_state $? $substate

        fi

    )

done

#
# Bam to FastQ
#

# State Management Delegated to Substates
format_status "Running Samtools BAM to FastQ"
# Retrieve Files
files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.sorted.*.bam))

for file in $files
    # In Parallel
    do ( 
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.sorted.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        substate="$fileprefix.$subset.$condition.$readgroup:BAMTOFASTQ:3"
        
        # State Check - Run Block if it Has Not Already Been Executed Successfully
        if !(has_state $substate); then
            
            # Call Bam to FastQ
            # Define Command
            call="samtools fastq -t -O $dataDir/downloaded/split/$fileprefix.$subset.$condition.sorted.$readgroup.bam > $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Check for failed parallel call
            put_state $? $substate

        fi

    )

done

format_status "BAM to FASTQ Complete"

format_status "Done"

