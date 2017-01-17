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

        -i=*|--index=*)
        indexOpt="${i#*=}"
        shift
        ;;
        -a=*|--align=*)
        alignOpt="${i#*=}"
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
        printf "\nInvalid/Unused Parameter: $i\n"
        ;;
        
    esac
done

# Defaults if No Arguments Passed
ncoresDef="20"
memoryDef="4G"
indexDef=false
alignDef="mem"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
index=${indexOpt:-$indexDef}
align=${alignOpt:-$alignDef}

# Get Max Allowable Memory
allocMemory=$(echo "$memory" | sed "s|[GMKgmk]||")
allocSize=$(echo "$memory" | sed "s|[0-9]*||")
maxMemory=$(($allocMemory * $ncores))$allocSize
 
printf "\n\nPARAMETERS: 
Picard Directory    = $jar
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Build BWA Index?    = $index
BWA Alignment       = $align
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning BWA Script\n"

#
# BWA Index
#

if [ "$index" = "true" ]; then
    printf "\n\nBWA Index\n"
    printf "\n\nCommand:\nbwa index -a bwtsw $PIPELINE_REF/Homo_sapiens_assembly19.fasta\n"
    bwa index -a bwtsw $PIPELINE_REF/Homo_sapiens_assembly19.fasta
    printf "\n\nBWA Index Complete\n"
fi

#
# BWA
#

# Retrieve Files
files=$(echo $(ls $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.*.fastq))

# Mem
if [ "$align" = "mem" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        # In Parallel
        do (
            # Extract Read Group to Pass to BWA mem
            suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
            readgroup=$(echo "$suffix" | sed "s|.fastq$||")
            # BWA mem
            printf "\n\nCommand:\nbwa mem -M -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
            bwa $align -M -t -R $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
        ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"

# BWASW
elif [ "$align" = "bwasw" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        # In Parallel
        do (
            # Extract Read Group to Pass to BWA bwasw
            suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
            readgroup=$(echo "$suffix" | sed "s|.fastq$||")
            # BWA bwasw
            printf "\n\nCommand:\nbwa $align -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
            bwa $align -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
        ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"  

# SAMPE / SAMSE (Backtrack)
elif [ "$align" = "samse" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        # In Parallel
        do (
            # Extract Read Group to Pass to BWA aln then samse or sampe
            suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
            readgroup=$(echo "$suffix" | sed "s|.fastq$||")
            # BWA aln
            printf "\n\nCommand:\nbwa aln -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai\n"
            bwa aln -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai
            # BWA samse
            printf "\n\nCommand:\nbwa $align -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
            bwa $align $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
        ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"

elif [ "$align" = "sampe" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        # In Parallel
        do (
            # Extract Read Group to Pass to BWA aln then samse or sampe
            suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
            readgroup=$(echo "$suffix" | sed "s|.fastq$||")
            # Split Paired End FastQ Files
            seqtk seq -l0 -1 $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq
            seqtk seq -l0 -2 $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq
            # BWA aln
            printf "\n\nCommand:\nbwa aln -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai\n"
            bwa aln -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai
            printf "\n\nCommand:\nbwa aln -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai\n"
            bwa aln -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sai
            # BWA sampe
            printf "\n\nCommand:\nbwa $align $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sai $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
            bwa $align $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sai $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
        ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"

else

    printf "\n\nInvalid BWA algorithm parameter: $align\n"

fi

# Backtrack

printf "\n\nDone\n"
