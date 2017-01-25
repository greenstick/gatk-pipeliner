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

format_status "PARAMETERS:
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
tmpDir=$PIPELINE_HOME/$subset/tmp

format_status "Running BAM to FASTQ Script"

#
# Shuffle & Split BAM
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition:BAMTOFASTQ:1"
if [ (state_registered $state) != 0 ]; then

    format_status "Shuffling & Splitting Merged BAM"
    format_status "Command:\nsamtools collate -uO $dataDir/downloaded/$fileprefix.$subset.$condition.bam $tmp | samtools split -f $dataDir/downloaded/split/$fileprefix.$subset.$condition.%%!.bam -"
    samtools collate -uO $dataDir/downloaded/$fileprefix.$subset.$condition.bam $tmpDir | bam splitBam -o $dataDir/downloaded/split/$fileprefix.$subset.$condition -v -L $PIPELINE_HOME/logs/splitbam_$fileprefix.$subset.$condition.log
    
    # Update State on Exit
    register_state $? $state
    format_status "Shuffling & Splitting Merged BAM Complete"

fi

#
# Bam to FastQ
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition:BAMTOFASTQ:2"
if [ (state_registered $state) != 0 ]; then

    format_status "Running Picard BAM to FASTQ"
    # Retrieve Files
    files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

    for file in $files
        # In Parallel
        do ( 
            # Get Read Group to Process
            suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
            readgroup=$(echo "$suffix" | sed "s|.bam$||")
            substate="$fileprefix.$subset.$condition.$readgroup:BAMTOFASTQ:2"
            
            # Run Command
            if [ (state_registered $substate) != 0 ]; then
                
                # Call Bam to FastQ
                format_status "Command:\nsamtools fastq -t $dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam > $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq"
                samtools fastq -t $dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam > $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq
            
                # Check for failed parallel call
                register_state $? $substate

            fi
        ) &

    done
    wait # Prevent Premature Exiting of Script

    # Update State on Exit
    register_state $? $state
    format_status "BAM to FASTQ Complete"

fi

format_status "Done"

