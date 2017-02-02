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
# Run Rcorrector
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:RCORRECTOR:1"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running Rcorrector - Default Parameters"
        # Call Error Model
        format_status "Command:\n
        perl $RCORRECTOR \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -k 31 \
        -s 3g \
        -t $ncores \
        -od $paramDir/modeled
        -verbose"
        perl $RCORRECTOR \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -k 31 \
        -s 3g \
        -t $ncores \
        -od $paramDir/modeled
        -verbose

    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #

        format_status "Running Rcorrector - Custom Parameters"
        # Call Error Model
        format_status "Command:\n
        perl $RCORRECTOR \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -k 31 \
        -s 3g \
        -t $ncores \
        -od $paramDir/modeled
        -verbose"
        perl $RCORRECTOR \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -k 31 \
        -s 3g \
        -t $ncores \
        -od $paramDir/modeled
        -verbose

    fi
    
    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "Rcorrector ($parameters $readgroup) Complete"
    return $status

fi

# Usage: perl ./run_rcorrector.pl [OPTIONS]
# OPTIONS:
# Required parameters:
# 	-s seq_files: comma separated files for single-end data sets
# 	-1 seq_files_left: comma separated files for the first mate in the paried-end data sets
# 	-2 seq_files_right: comma separated files for the second mate in the paired-end data sets
# 	-i seq_files_interleaved: comma sperated files for interleaved paired-end data sets
# Other parameters:
# 	-k kmer_length (<=32, default: 23)
# 	-od output_file_directory (default: ./)
# 	-t number_of_threads (default: 1)
# 	-maxcorK INT: the maximum number of correction within k-bp window (default: 4)
# 	-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold, lower for more divergent genome (default: 0.95)
# 	-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)
# 	-stdout: output the corrected reads to stdout (default: not used)
# 	-verbose: output some correction information to stdout (default: not used)
# 	-stage INT: start from which stage (default: 0)
# 		0-start from begining(storing kmers in bloom filter);
# 		1-start from count kmers showed up in bloom filter;
# 		2-start from dumping kmer counts into a jf_dump file;
# 		3-start from error correction.