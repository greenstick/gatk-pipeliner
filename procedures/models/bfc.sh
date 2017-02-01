#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

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
        -g=*|--readgroup=*)
        readgroup="${i#*=}"
        shift # Access & Write Files With This Read Group
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

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
allocMax=$((allocMemory * ncores))
maxMemory=$((allocMemory * ncores))$allocSize

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

# 
# Run BFC
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BFC:1"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running BFC - Default Parameters"
        # Call Error Model - Size Parameter = 3 Gigabase Human Genome (hg19)
        format_status "Command:\n
        $BFC \
        -s 3g \
        -t $ncores \
        -prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        > $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup.fastq"
        $BFC \
        -s 3g \
        -t $ncores \
        -prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        > $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup.fastq

    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #

        format_status "Running BFC - Custom Parameters"
        # Call Error Model
        format_status "Command:\n
        $BFC \
        -s 3g \
        -t $ncores \
        -prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        > $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup.fastq"
        $BFC \
        -s 3g \
        -t $ncores \
        -prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        > $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup.fastq

    fi
    
    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "BFC ($parameters $readgroup) Complete"
    return $status

fi

# ./bfc
# Usage: bfc [options] <to-count.fq> [to-correct.fq]
# Options:
#   -s FLOAT     approx genome size (k/m/g allowed; change -k and -b) [unset]
#   -k INT       k-mer length [33]
#   -t INT       number of threads [1]
#   -b INT       set Bloom filter size to pow(2,INT) bits [33]
#   -H INT       use INT hash functions for Bloom filter [4]
#   -d FILE      dump hash table to FILE [null]
#   -E           skip error correction
#   -R           refine bfc-corrected reads
#   -r FILE      restore hash table from FILE [null]
#   -w INT       no more than 5 ec or 2 highQ ec in INT-bp window [10]
#   -c INT       min k-mer coverage [3]
#   -Q           force FASTA output
#   -1           drop reads containing unique k-mers
#   -v           show version number
#   -h           show command line help
