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

    # Additional Arguments

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
alignDef="mem"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
align=${alignOpt:-$alignDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize
 
printf "\n\nPARAMETERS: 
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
BWA Alignment       = $align
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
\n"

# Set Directories
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

printf "\n\nRunning BWA Script\n"

#
# BWA Index
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix:BWA:1" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    printf "\n\nBWA Index\n"
    printf "\n\nCommand:\nbwa index -a bwtsw $PIPELINE_REF/Homo_sapiens_assembly19.fasta\n"
    bwa index -a bwtsw $PIPELINE_REF/Homo_sapiens_assembly19.fasta
    
    # Update State on Exit
    statuscode=$?
    if [ $statuscode = 0 ]; then
        # Export Pipeline
        echo "$fileprefix:BWA:1" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nBWA Index Complete\n"
    else
        printf "\n\nUnexpected Exit $statuscode - $fileprefix:BWA:1"
    fi

fi

#
# BWA - mem or bwasw
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:BWA:2" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    # Retrieve Files
    files=$(echo $(ls $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.*.fastq))

    # Mem
    if [ "$align" = "mem" ]; then

        failures=0
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
            
                # Check for failed parallel call
                if [ $? != 0 ]; then
                    failures=$((failures + 1))
                fi
            ) &
        done
        wait # Prevent Premature Exiting of Script

        # Update State on Exit
        if [ $failures = 0 ]; then
            # Export Pipeline
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:BWA:2" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBWA $align Complete\n"
        else
            printf "\n\n$failures Failures, Exiting - $fileprefix.$subset.$condition.$experiment.$parameters:BWA:2"
        fi

    # BWASW
    elif [ "$align" = "bwasw" ]; then

        failures=0
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
                
                # Check for failed parallel call
                if [ $? != 0 ]; then
                    failures=$((failures + 1))
                fi
            ) &
        done
        wait # Prevent Premature Exiting of Script

        # Update State on Exit
        if [ $failures = 0 ]; then
            # Export Pipeline
            echo "$fileprefix.$subset.$condition.$experiment.$parameters:BWA:2" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBWA $align Complete\n"
        else
            printf "\n\n$failures Failures, Exiting - $fileprefix.$subset.$condition.$experiment.$parameters:BWA:2"
        fi

    else

        printf "\n\nInvalid BWA algorithm parameter: $align\n"

    fi

fi

# Backtrack

printf "\n\nDone\n"
