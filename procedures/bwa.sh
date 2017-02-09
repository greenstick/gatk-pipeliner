#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

# Load ~/.bash_profile if Not Found
if [ -z $PIPELINE_HOME ]; then
    source ~/.bash_profile
fi

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
        printf "Invalid/Unused Parameter: $i"
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
 
printf "\nPARAMETERS: 
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

format_status "Running BWA Script"

#
# BWA Index
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix:BWA:1"
if !(has_state $state); then

    format_status "BWA Index"
    format_status "Command:\nbwa index -a bwtsw $PIPELINE_REF/Homo_sapiens_assembly19.fasta"
    bwa index -a bwtsw $PIPELINE_REF/Homo_sapiens_assembly19.fasta
    
    # Update State on Exit
    put_state $? $state
    format_status "BWA Index Complete"

fi

#
# BWA - mem or bwasw
#

# State Management Delegated to Substates
format_status "Running BWA $align"
# Retrieve Files
files=$(echo $(ls $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.*.fastq))

# Mem
if [ "$align" = "mem" ]; then

    for file in $files
        do (
            # Extract Read Group to Pass to BWA mem
            suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
            readgroup=$(echo "$suffix" | sed "s|.fastq$||")
            substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BWA:2"
            
            # Run Command
            if !(has_state $substate); then

                # Test for Paired Ends
                head -n 1000 $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq | grep -qE "^@.*/3(\s|\t|\n)"
                end1=$?
                head -n 1000 $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq | grep -qE "^@.*/3(\s|\t|\n)"
                end2=$?

                # Is Interleaved?
                paired=""
                if [ $end1 ] && [ $end2 ]; then
                    paired="-p"
                fi

                # Call BWA mem
                format_status "Command:\nbwa mem -M -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam"
                bwa $align $paired -M -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
                
                # Check for failed parallel call
                put_state $? $substate

            fi
        )
    done

    # Update State on Exit
    put_state $? $state
    format_status "BWA $align Complete"

grep --no-filename @HWUSI-EAS100R:6:73:941:1973 *.fastq | cut -d' ' -f1 | sort | uniq -c | sort -rgk 1,1 | head

# BWASW
elif [ "$align" = "bwasw" ]; then

    for file in $files
        do (
            # Extract Read Group to Pass to BWA mem
            suffix=$(echo "$file" | sed "s|$paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.||")
            readgroup=$(echo "$suffix" | sed "s|.fastq$||")
            substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BWA:2"
            
            # Run Command
            if !(has_state $substate); then

                # BWA bwasw
                format_status "Command:\nbwa $align -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam"
                bwa $align -t $ncores $PIPELINE_REF/Homo_sapiens_assembly19.fasta $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam
            
                # Check for failed parallel call
                put_state $? $substate

            fi
        )
    done

    # Update State on Exit
    put_state $? $state
    format_status "BWA $align Complete"

else

    format_status "Invalid BWA algorithm parameter: $align"

    fi

# Backtrack

format_status "Done"
