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
        printf "Invalid Parameter: $i"
        ;;
    esac
done

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
grep -q "$fileprefix.$subset.$condition.$experiment.$parameters:QUORUM:1" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    # Retrieve Files
    files=$(echo $(ls $dataDir/fastq/split/$fileprefix.$subset.$condition.*.bam))
    failures=0

    if [ "$parameters" = "default" ]; then
        
        printf "\n\nRunning Quorum - Default Parameters"
        for file in $files
            # In Parallel
            do (
                # Get Read Group to Process
                suffix=$(echo "$file" | sed "s/split\/$fileprefix\.$subset\.$condition\.//")
                readgroup=$(echo "$suffix" | sed "s/.bam$//")
                # Call Error Model
                printf "\n\nCommand:\n \
                $QUORUM \
                $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
                --prefix $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup \
                -t $ncores \
                --size 20G \
                --no-discard \
                --debug"
                $QUORUM \
                $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
                --prefix $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup \
                -t $ncores \
                --size 20G \
                --no-discard \
                --debug

                # Check for failed parallel call
                if [ $? != 0 ]; then
                    failures=$((failures + 1))
                fi
            ) &
        done
        wait # Prevent Premature Exiting of Script

    elif [ "$parameters" = "custom" ]; then
        
        printf "\n\nRunning Quorum - Custom Parameters"
        for file in $files
            # In Parallel
            do (
                # Get Read Group to Process
                suffix=$(echo "$file" | sed "s/split\/$fileprefix\.$subset\.$condition\.//")
                readgroup=$(echo "$suffix" | sed "s/.bam$//")
                # Call Error Model
                printf "\n\nCommand:\n \
                $QUORUM \
                $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
                --prefix $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup \
                -t $ncores \
                --size 20G \
                --no-discard \
                --debug"
                $QUORUM \
                $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
                --prefix $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup \
                -t $ncores \
                --size 20G \
                --no-discard \
                --debug

                # Check for failed parallel call
                if [ $? != 0 ]; then
                    failures=$((failures + 1))
                fi
            ) &
        done
        wait # Prevent Premature Exiting of Script

    fi

    # Update State on Exit
    if [ $failures = 0 ]; then
        # Export Pipeline State
        echo "$fileprefix.$subset.$condition.$experiment.$parameters:QUORUM:1" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nQuorum Complete"
    else
        printf "\n\n$failures Failures, Exiting - $fileprefix.$subset.$condition.$experiment.$parameters:QUORUM:1"
    fi

fi

printf "\n\nDone\n"

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
