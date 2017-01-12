#! /usr/bin/bash

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

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# 
# Run BayesHammer
# 

files=$(echo $(ls $dataDir/fastq/split/$fileprefix.$subset.$condition.*.bam))

if [ "$parameters" = "default" ]; then
    
    printf "\n\nRunning Bloocoo - Default Parameters"
    for file in $files
        # In Parallel
        do (
            # Get Read Group to Process
            suffix=$(echo "$file" | sed "s/split\/$fileprefix\.$subset\.$condition\.//")
            readgroup=$(echo "$suffix" | sed "s/.bam$//")
            # Call Error Model
            printf "\n\nCommand:\n
            $BLOOCOO \
            -file $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
            -nb-cores $ncores"
            $BLOOCOO \
            -file $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
            -nb-cores $ncores
            # Move File to Output Directory
            mv $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup_corrected.fastq $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq  
        ) &
    done
    wait # Prevent Premature Exiting of Script

elif [ "$parameters" = "custom" ]; then
    
    printf "\n\nRunning Bloocoo - Custom Parameters"
    for file in $files
        # In Parallel
        do (
            # Get Read Group to Process
            suffix=$(echo "$file" | sed "s/split\/$fileprefix\.$subset\.$condition\.//")
            readgroup=$(echo "$suffix" | sed "s/.bam$//")
            # Call Error Model
            printf "\n\nCommand:\n
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
            # Move File to Output Directory
            mv $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup_corrected.fastq $paramDir/modeled/$prefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq   
        ) &
    done
    wait # Prevent Premature Exiting of Script

fi

printf "\n\nDone\n"

#
# Bloocoo Arguments
#

# Usage: /home/users/cordier/packages/bayes-hammer/SPAdes-3.9.0-Linux/bin/spades.py [options] -o <output_dir>

# -file : the read file name, can be fasta, fastq, gzipped or not.
# -kmer-size : the k-mer size (typically ~31)
# -abundance-min : the minimal abundance threshold defining solid k-mers (typically  between 3 and 6, but depends on the read depth, you can also use 'auto' and it is automatically inferred from the data)
# -nb-cores  : number of threads used
# -high-recall  :  correct more errors but can also introduce more mistakes
# -slow : slower modes with more pass, but better correction
# -high-precision :  correct safely, correct less errors but introduce less mistakes
# -ion : (experimental) mode for correcting indels present in  ion torrent reads