#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
set -o errexit

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

printf "\nPARAMETERS:
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
\n\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning BAM to FASTQ Script"

#
# Shuffle & Split BAM
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix.$subset.$condition:BAMTOFASTQ:1" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    printf "\n\nShuffling & Splitting Merged BAM"
    printf "\n\nCommand:\nsamtools collate -uO $dataDir/downloaded/$fileprefix.$subset.$condition.bam $tmp | samtools split -f $dataDir/downloaded/split/$fileprefix.$subset.$condition.%%!.bam -"
    samtools collate -uO $dataDir/downloaded/$fileprefix.$subset.$condition.bam $tmpDir | bam splitBam -o $dataDir/downloaded/split/$fileprefix.$subset.$condition -v -L $PIPELINE_HOME/logs/splitbam_$fileprefix.$subset.$condition.log
    
    # Update State on Exit
    statuscode=$?
    if [ $statuscode = 0 ]; then
        # Export Pipeline State
        echo "$fileprefix.$subset.$condition:BAMTOFASTQ:1" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nShuffling & Splitting Merged BAM Complete"
    else
        printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition:BAMTOFASTQ:1"
    fi

fi

#
# Bam to FastQ
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix.$subset.$condition:BAMTOFASTQ:2" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    printf "\n\nRunning Picard Bam to FastQ"
    # Retrieve Files
    files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))
    failures=0

    for file in $files
        # In Parallel
        do ( 
            # Get Read Group to Process
            suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
            readgroup=$(echo "$suffix" | sed "s|.bam$||")
            # Call Bam to FastQ
            printf "\n\nCommand:\nsamtools fastq -t $dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam > $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq\n"
            samtools fastq -t $dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam > $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq
            
            # Check for failed parallel call
            if [ $? != 0 ]; then
                failures=$((failures + 1))
            fi
        ) &

    done
    wait # Prevent Premature Exiting of Script

    # Update State on Exit
    if [ $failures = 0 ]; then
        # Export Pipeline State
        echo "$fileprefix.$subset.$condition:BAMTOFASTQ:2" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nBam to FastQ Complete"
    else
        printf "\n\n$failures Failures, Exiting - $fileprefix.$subset.$condition:BAMTOFASTQ:2"
    fi

fi

printf "\n\nDone\n"

