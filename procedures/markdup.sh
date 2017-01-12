#! /usr/bin/bash

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

    # Optional Arguments With Defaults

        -b=*|--buildindex=*)
        indexOpt="${i#*=}"
        shift # Build BAM Index?
        ;;
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
indexDef=true

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
index=${indexOpt:-$indexDef}

printf "\nPARAMETERS:
Picard Directory    = $jar
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Build BAM Index     = $index
Recalibration Model = $qualitymodel
Memory              = $memory
\n\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Mark Duplicates Script"

cd $dataDir

# If Norealignment, Get Data From Download Directory & Skip Sorting
if [ "$experiment" = "norealign" ]; then

    #
    # Mark Duplicates
    #

    printf "\n\nMarkDuplicates Start"
    printf "\n\nCommand:\njava -Xmx$memory \
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

    # 
    # Sort BAM
    # 

    printf "\n\nSorting BAM"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $jar SortSam \ 
    I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \ 
    O=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \ 
    SORT_ORDER=coordinate \
    TMP_DIR=$tmpDir\n"
    java -Xmx$memory \
    -jar $jar SortSam \
    I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    O=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    SORT_ORDER=coordinate \
    TMP_DIR=$tmpDir
    printf "\n\nSorting Complete"

    #
    # Create New BAM Index
    #

    printf "\n\nIndexing BAM Output"
    printf "\n\nCommand:\nsamtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam\n"
    samtools index $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam
    printf "\n\nBAM Indexing Complete"

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
    I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
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

