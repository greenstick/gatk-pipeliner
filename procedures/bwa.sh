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
ncoresDef="4"
memoryDef="8G"
indexDef=false
alignDef="sampe"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
index=${indexOpt:-$indexDef}
align=${alignOpt:-$alignDef}
 
printf "\n\nPARAMETERS: 
Picard Directory    = $jar
Reference Directory = $reference
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Run BWA Index?      = $index
BWA Alignment       = $align
Cores               = $ncores
Memory              = $memory
\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning BWA Script\n"

cd $dataDir

#
# BWA Index
#

if $index; then
    printf "\n\nBWA Index\n"
    printf "\n\nCommand:\nbwa index -a bwtsw $reference/Homo_sapiens_assembly19.fasta\n"
    bwa index -a bwtsw $reference/Homo_sapiens_assembly19.fasta
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
        do
            # In Parallel
            (
                # Extract Read Group to Pass to BWA mem
                suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # BWA mem
                printf "\n\nCommand:\nbwa mem -M -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
                bwa $align -M -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
            ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"

# BWASW
elif [ "$align" = "bwasw" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        do
            # In Parallel
            (
                # Extract Read Group to Pass to BWA bwasw
                suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # BWA bwasw
                printf "\n\nCommand:\nbwa $align -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
                bwa $align -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
            ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"  

# SAMPE / SAMSE (Backtrack)
elif [ "$align" = "samse" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        do
            # In Parallel
            (
                # Extract Read Group to Pass to BWA aln then samse or sampe
                suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # BWA aln
                printf "\n\nCommand:\nbwa aln -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai\n"
                bwa aln -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai
                # BWA samse
                printf "\n\nCommand:\nbwa $align -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam\n"
                bwa $align $reference/Homo_sapiens_assembly19.fasta $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sai $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
            ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"

elif [ "$align" = "sampe" ]; then

    printf "\n\nBWA $align\n"
    for file in $files
        do
            # In Parallel
            (
                # Extract Read Group to Pass to BWA aln then samse or sampe
                suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Split Paired End FastQ Files
                grep '@.*/1' -A 3 --no-group-separator $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq
                grep '@.*/2' -A 3 --no-group-separator $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq
                # BWA aln
                printf "\n\nCommand:\nbwa aln -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai\n"
                bwa aln -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai
                printf "\n\nCommand:\nbwa aln -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai\n"
                bwa aln -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq > $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sai
                # BWA sampe
                printf "\n\nCommand:\nbwa $align -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sam\n"
                bwa $align $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sai $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.1.sam
                printf "\n\nCommand:\nbwa $align -t $ncores $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sai $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sam\n"
                bwa $align $reference/Homo_sapiens_assembly19.fasta $paramDir/modeled/indexes/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sai $paramDir/modeled/pairs/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.2.sam
            ) &
    done
    wait # Prevent Premature Exiting of Script
    printf "\n\nBWA $align Complete\n"

else

    printf "\n\nInvalid BWA algorithm parameter: $align\n"

fi

# Backtrack

printf "\n\nDone\n"
