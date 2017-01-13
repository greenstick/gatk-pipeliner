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
        -b=*|--splitbam=*)
        splitbamOpt="${i#*=}"
        shift
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
memoryDef="6G"
splitbamDef=true

# Set Optional Values
memory=${memoryOpt:-$memoryDef}
splitbam=${splitbamOpt:-$splitbamDef}

printf "\nPARAMETERS:
Picard Directory    = $jar
Reference Directory = $PIPELINE_REF
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Split BAM           = $splitbam
Memory              = $memory
\n\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning BAM to FASTQ Script"

#
# Shuffle & Split BAM
#

if [ "$splitbam" = "true" ]; then
    printf "\n\nShuffling & Splitting Merged BAM"
    printf "\n\nCommand:\nsamtools collate -uO $dataDir/downloaded/$fileprefix.$subset.$condition.bam $tmp | samtools split -f $dataDir/downloaded/split/$fileprefix.$subset.$condition.%%!.bam -"
    samtools collate -uO $dataDir/downloaded/$fileprefix.$subset.$condition.bam $tmpDir | samtools split -f $dataDir/downloaded/split/$fileprefix.$subset.$condition.%!.bam -
    printf "\n\nShuffling & Splitting Merged BAM Complete"
fi

#
# Bam to FastQ
#

# Retrieve Files
files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

printf "\n\nRunning Picard Bam to FastQ"
for file in $files
    # In Parallel
    do ( 
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        # Call Bam to FastQ
        printf "\n\nCommand:\njava -Xmx$memory \
        -jar $jar SamToFastq \
        I=$dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam \
        F=$dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        INTERLEAVE=true \
        INCLUDE_NON_PF_READS=true \
        TMP_DIR=$tmpDir\n"
        java -Xmx$memory \
        -jar $jar SamToFastq \
        I=$dataDir/downloaded/split/$fileprefix.$subset.$condition.$readgroup.bam \
        F=$dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        INTERLEAVE=true \
        INCLUDE_NON_PF_READS=true \
        TMP_DIR=$tmpDir
    ) &
done
wait # Prevent Premature Exiting of Script

printf "\n\nBam to FastQ Complete"

printf "\n\nDone\n"

