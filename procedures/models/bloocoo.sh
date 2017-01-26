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
maxMemory=$((allocMemory * ncores))$allocSize

# Workaround to Bloocoo's Funky Output Scheme (See State BLOOCOO:2)
corrected="_corrected"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

# 
# Run Bloocoo
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters:BLOOCOO:1"
if !(state_registered $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running Bloocoo - $parameters Parameters"
        # Call Error Model & Move Outputs to Output Directory
        format_status "Command:\n
        $BLOOCOO \
        -file $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -nb-cores $ncores"
        $BLOOCOO \
        -file $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -nb-cores $ncores

    elif [ "$parameters" = "custom" ]; then
        
        format_status "Running Bloocoo - Custom Parameters"
        # Call Error Model
        format_status "Command:\n
        $BLOOCOO \
        -file $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -nb-cores $ncores \
        -slow \
        -high-precision"
        $BLOOCOO \
        -file $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -nb-cores $ncores \
        -slow \
        -high-precision

    fi

    # Update State on Exit
    status=$?
    register_state $status $state

fi

#
# Move Output & Cleanup Extra Files
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters:BLOOCOO:2"
if !(state_registered $state); then

    format_status "Moving Output & Cleaning Up"
    format_status "Command:\n
    mv $proceduresDir/$fileprefix.$subset.$condition.$readgroup$corrected.fastq $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq \
    && mv $proceduresDir/$fileprefix.$subset.$condition.$readgroup.h5 $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup.h5 \
    && rm $proceduresDir/$fileprefix.$subset.$condition.$readgroup$corrected.fastq
    && rm $proceduresDir/$fileprefix.$subset.$condition.$readgroup.h5
    "
    mv $proceduresDir/$fileprefix.$subset.$condition.$readgroup$corrected.fastq $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq \
    && mv $proceduresDir/$fileprefix.$subset.$condition.$readgroup.h5 $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup.h5 \
    && rm $proceduresDir/$fileprefix.$subset.$condition.$readgroup$corrected.fastq
    && rm $proceduresDir/$fileprefix.$subset.$condition.$readgroup.h5

    # Update State on Exit
    status=$?
    register_state $status $state

    format_status "Bloocoo ($parameters $readgroup) Complete"
    return $status

fi

#
# Bloocoo Arguments
#

# -file : the read file name, can be fasta, fastq, gzipped or not.
# -kmer-size : the k-mer size (typically ~31)
# -abundance-min : the minimal abundance threshold defining solid k-mers (typically  between 3 and 6, but depends on the read depth, you can also use 'auto' and it is automatically inferred from the data)
# -nb-cores  : number of threads used
# -high-recall  :  correct more errors but can also introduce more mistakes
# -slow : slower modes with more pass, but better correction
# -high-precision :  correct safely, correct less errors but introduce less mistakes
# -ion : (experimental) mode for correcting indels present in  ion torrent reads