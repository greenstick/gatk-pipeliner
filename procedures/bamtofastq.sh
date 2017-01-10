#! /usr/bin/bash

# Assign Arguments
for i in "$@"
    do case $i in

    # Standard Arguments

        -r=*|--ref=*)
        reference="${i#*=}"
        shift # Reference Sequence Directory
        ;;
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
memoryDef="20G"

# Set Optional Values
memory=${memoryOpt:-$memoryDef}

printf "\nPARAMETERS:
Picard Directory    = $jar
Reference Directory = $reference
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Memory              = $memory
\n\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning BAM to FASTQ Script"

cd $dataDir

#
# Shuffle & Split BAM
#

printf "\n\nShuffling & Splitting Merged BAM"
samtools collate -uO downloaded/$fileprefix.$subset.$condition.bam $tmp | samtools split -f split/$fileprefix.$subset.$condition.%!.bam -
 
#
# Bam to FastQ
#

# Retrieve Files
files=$(echo $(ls split/$fileprefix.$subset.$condition.*.bam))

printf "\n\nRunning Picard Bam to FastQ"
for file in $files
    # In Parallel
    do ( 
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s/split\/$fileprefix\.$subset\.$condition\.//")
        readgroup=$(echo "$suffix" | sed "s/.bam$//")
        # Call Bam to FastQ
        printf "\n\nCommand:\njava -Xmx$memory \
        -jar $jar SamToFastq \
        I=$file \
        F=fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        OUTPUT_PER_RG=true \
        INTERLEAVE=true \
        INCLUDE_NON_PF_READS=true \
        TMP_DIR=$tmpDir\n"
        java -Xmx$memory \
        -jar $jar SamToFastq \
        I=$file \
        F=fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        OUTPUT_PER_RG=true \
        INTERLEAVE=true \
        INCLUDE_NON_PF_READS=true \
        TMP_DIR=tmp
    ) &
done
wait # Prevent Premature Exiting of Script

printf "\n\nBam to FastQ Complete"

cd ../procedures

printf "\n\nDone\n"

