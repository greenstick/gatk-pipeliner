#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

# Load ~/.bash_profile if Not Found
if [ -z $PIPELINE_HOME ]; then
    echo "Reloading ~/.bash_profile"
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
        printf "Invalid Parameter: $i"
        ;;
    esac
done

ncoresDef="16"
memoryDef="15G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

# 
# Run Quorum
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters:QUORUM:1"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then

        #
        # Default Parameters
        #
        
        format_status "Running Quorum - $parameters Parameters"
        # Call Error Model
        format_status "Command:\n \
        quorum \
        $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        --prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        -t $ncores \
        --size 100000000000000 \
        --no-discard \
        --min-q-char 33 \
        --debug"
        quorum \
        $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        --prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        -t $ncores \
        --size 100000000000000 \
        --no-discard \
        --min-q-char 33 \
        --debug

    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #
        
        format_status "Running Quorum - $parameters Parameters"
        # Call Error Model
        format_status "Command:\n \
        quorum \
        $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        --prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        -t $ncores \
        --size 100000000000000 \
        --no-discard \
        --min-q-char 33 \
        --debug"
        quorum \
        $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        --prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
        -t $ncores \
        --size 100000000000000 \
        --no-discard \
        --min-q-char 33 \
        --debug

    fi

    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "Quorum ($parameters $readgroup) Complete"

    return $status

fi

#
# Quorum Arguments
#

# -s, --size              Mer database size (default 200M)
# -t, --threads           Number of threads (default number of cpus)
# -p, --prefix            Output prefix (default quorum_corrected)
# -k, --kmer-len          Kmer length (default 24)
# -q, --min-q-char        Minimum quality char. Usually 33 or 64 (autodetect)
# -m, --min-quality       Minimum above -q for high quality base (5)
# -w, --window            Window size for trimming
# -e, --error             Maximum number of errors in a window
#     --min-count         Minimum count for a k-mer to be good
#     --skip              Number of bases to skip to find anchor kmer
#     --anchor            Numer of good kmer in a row for anchor
#     --anchor-count      Minimum count for an anchor kmer
#     --contaminant       Contaminant sequences
#     --trim-contaminant  Trim sequences with contaminant mers
# -d, --no-discard        Do not discard reads, output a single N (false)
# -P, --paired-files      Preserve mate pairs in two files
#     --homo-trim         Trim homo-polymer on 3' end
#     --debug             Display debugging information
#     --version           Display version
# -h, --help              This message
