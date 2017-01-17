#! /usr/bin/bash

# Exit on First Error
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
        -q=*|--qualitymodel=*)
        qualitymodel="${i#*=}"
        shift # Access & Write Files With This Quality Model
        ;;

    # Additional Arguments

        -j=*|--jar=*)
        jar="${i#*=}"
        shift
        ;;
        -b=*|--sortbam=*)
        sortbamOpt="${i#*=}"
        shift # Build BAM Index?
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
sortbamDef=true

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
sortbam=${sortbamOpt:-$sortbamDef}

# Get Max Allowable Memory
allocMemory=$(echo "$memory" | sed "s|[GMKgmk]||")
allocSize=$(echo "$memory" | sed "s|[0-9]*||")
maxMemory=$(($allocMemory * $ncores))$allocSize

printf "\nPARAMETERS:
Picard Directory    = $jar
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Sort BAM?           = $sortbam
Recalibration Model = $qualitymodel
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
\n\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Mark Duplicates Script"

# If Norealignment, Get Data From Download Directory & Skip Sorting
if [ "$experiment" = "norealign" ]; then

    #
    # Mark Duplicates
    #

    printf "\n\nMarkDuplicates Start"
    printf "\n\nCommand:\njava -Xmx$maxMemory \
    -jar $jar MarkDuplicates \
    I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
    O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
    PG=null \
    TMP_DIR=$tmpDir\n"
    java -Xmx$memory \
    -jar $jar MarkDuplicates \
    I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
    O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
    PG=null \
    TMP_DIR=$tmpDir
    printf "\n\nMark Duplicates Complete"

    #
    # Create New BAM Index
    #

    printf "\n\nIndexing BAM Output"
    printf "\n\nCommand:\nsamtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam\n"
    samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam
    printf "\n\nBAM Indexing Complete"

# Otherwise, Grab Data From Merged Alignment Directory & Sort
else

    if [ "$sortbam" = "true" ]; then
    
        # 
        # Sort BAM
        # 

        printf "\n\nSorting BAM"
        printf "\n\nCommand:\nsamtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam\n"
        samtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam
        printf "\n\nSorting Complete"

        #
        # Create New BAM Index
        #

        printf "\n\nIndexing BAM Output"
        printf "\n\nCommand:\nsamtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam\n"
        samtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam
        printf "\n\nBAM Indexing Complete"

    fi

    #
    # Mark Duplicates
    #

    printf "\n\nMarkDuplicates Start"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $jar MarkDuplicates \
    I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
    PG=null \
    TMP_DIR=$tmpDir\n"
    java -Xmx$memory \
    -jar $jar MarkDuplicates \
    I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam \
    O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
    PG=null \
    TMP_DIR=$tmpDir
    printf "\n\nMark Duplicates Complete"

    #
    # Create New BAM Index
    #

    printf "\n\nIndexing BAM Output"
    printf "\n\nCommand:\nsamtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam\n"
    samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam
    printf "\n\nBAM Indexing Complete"

fi

printf "\n\nDone\n"

