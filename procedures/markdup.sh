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
\n\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Mark Duplicates Script"

# If Norealignment, Get Data From Download Directory & Skip Sorting
if [ "$experiment" = "norealign" ]; then

    #
    # Mark Duplicates
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nMarkDuplicates Start"
        printf "\n\nCommand:\njava -Xmx$maxMemory \
        -jar $PICARD MarkDuplicates \
        I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        PG=null \
        TMP_DIR=$tmpDir\n"
        java -Xmx$memory \
        -jar $PICARD MarkDuplicates \
        I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        PG=null \
        TMP_DIR=$tmpDir

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nMark Duplicates Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1"
        fi

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nIndexing BAM"
        printf "\n\nCommand:\nsamtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam\n"
        samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam
        
        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBAM Indexing Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2"
        fi

    fi

# Otherwise, Grab Data From Merged Alignment Directory & Sort
else

    # 
    # Sort BAM
    # 

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nSorting BAM"
        printf "\n\nCommand:\nsamtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam\n"
        samtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam
        
        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nSort BAM Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1"
        fi

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nIndexing BAM Output"
        printf "\n\nCommand:\nsamtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam\n"
        samtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBAM Indexing Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2"
        fi

    fi

    #
    # Mark Duplicates
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:3" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nMarkDuplicates Start"
        printf "\n\nCommand:\njava -Xmx$memory \
        -jar $PICARD MarkDuplicates \
        I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        PG=null \
        TMP_DIR=$tmpDir\n"
        java -Xmx$memory \
        -jar $PICARD MarkDuplicates \
        I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        PG=null \
        TMP_DIR=$tmpDir

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:3" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nMark Duplicates Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:3"
        fi

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:4" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nIndexing BAM Output"
        printf "\n\nCommand:\nsamtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam\n"
        samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:4" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBAM Indexing Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:4"
        fi

    fi

fi

printf "\n\nDone\n"

