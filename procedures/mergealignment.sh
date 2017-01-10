#! /usr/bin/bash

# Assign Arguments
for i in "$@"
do
case $i in

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
ncoresDef="16"
memoryDef="8G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
index=${indexOpt:-$indexDef}
 
printf "\nPARAMETERS: 
Picard Directory    = $jar
Reference Directory = $reference
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Memory              = $memory
\n\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Merge BAMs Script"

cd $dataDir

# Retrieve Files to Process
files=$(echo $(ls $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.*.sam))

# 
# Merge BAMs to Single BAM
# 

printf "\n\nSamtools Merge"
printf "\n\nCommand:\nsamtools merge $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files"
samtools merge $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files
printf "\n\nSamtools Merge Complete"

printf "\n\nDone\n"