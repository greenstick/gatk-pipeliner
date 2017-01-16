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
ncoresDef="16"
memoryDef="5G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
 
printf "\nPARAMETERS: 
Picard Directory    = $jar
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Cores               = $ncores
Memory              = $memory
\n\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Merge BAMs Script"

# 
# Inserting Read Groups
# 

printf "\n\nSamtools AddReplaceRG"
# Retrieve Files to Process
files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

for file in $files
    # In Parallel
    do (
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        # Get Read Group Arguments to Pass to Samtools
        rgArgs=$(samtools view -H $file | grep '@RG' | awk -F '\t' '{print $2,$3,$4,$5,$8}' | sed "s|[A-Z][A-Z]:[a-zA-Z0-9\.\-]*|-r '&'|g")
        printf "\n\nCommand:\nsamtools addreplacerg $rgArgs $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.bam"
        # Insert Read Groups into New BAM - WARNING USES EVAL
        eval "samtools addreplacerg ${rgArgs[@]} $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.bam"
    ) &
done
wait # Prevent Premature Exiting of Script
printf "\n\nSamtools AddReplaceRG Complete"

# 
# Merge BAMs to Single BAM
# 

# Retrieve Files to Process
files=$(echo $(ls $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.*.bam))

printf "\n\nSamtools Merge"
printf "\n\nCommand:\nsamtools merge -r $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files"
samtools merge -r $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files
printf "\n\nSamtools Merge Complete"

# printf "\n\nDone\n"