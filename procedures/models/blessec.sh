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
# Run Bless-EC
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BLESSEC:1"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running Bless-EC - Custom Parameters"
        # Call Error Model
        format_status "Command:\n
        bless \
        -prefix $paramDir/modeled/ \
        -read $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -kmerlength 31 \
        -fpr 0.001 \
        -max_mem $memory \
        -notrim"
        bless \
        -prefix $paramDir/modeled/ \
        -read $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -kmerlength 31 \
        -fpr 0.001 \
        -max_mem $memory \
        -notrim

    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #

        format_status "Running Bless-EC - Custom Parameters"
        # Call Error Model
        format_status "Command:\n
        bless \
        -prefix $paramDir/modeled/ \
        -read $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -kmerlength 31 \
        -fpr 0.001 \
        -max_mem $memory \
        -notrim"
        bless \
        -prefix $paramDir/modeled/ \
        -read $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -kmerlength 31 \
        -fpr 0.001 \
        -max_mem $memory \
        -notrim

    fi
    
    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "Bless-EC ($parameters $readgroup) Complete"
    return $status

fi

# $BLESSEC -prefix /home/exacloud/lustre1/MRD_aml/ForBackup/pipeline/set1/model/blessec/param/default/modeled/synthetic.challenge.set1.tumor.blessec.default.C09DF.2 -read /home/exacloud/lustre1/MRD_aml/ForBackup/pipeline/set1/fastq/split/synthetic.challenge.set1.tumor.blessec.default.C09DF.2.fastq -kmerlength 31 -fpr 0.001 -max_mem 80G -notrim

# bless

# ----------------------------------------------------------------------
#            BLESS: Bloom-filter based error correction tool
# ----------------------------------------------------------------------
# VERSION: 0.24
# DATE   : Jan. 29, 2015
# ----------------------------------------------------------------------

# 1. USAGE
#      bless <OPTIONS>

# 2. OPTIONS
#      1) REQUIRED
#      -read <file>: Unpaired input fastq file. One of "-read" or
#           "-read1/-read2" should be used.
#      -read1 <file>: Forward fastq file of paired-end reads. It should
#           be used with "-read2".
#      -read2 <file>: Reverse fastq file of paired-end reads. It should
#           be used with "-read1".
#      -kmerlength <number>: Length of k-mers.
#      -prefix <string>: Output file prefix.

#      2) OPTIONAL
#      -count <integer>: Minimum occurrences for solid k-mers.
#           This number is automatically determined when it is not
#           given.
#      -extend <integer>: Read extension amount. Default: 5.
#      -fpr <float>: Target false positive probability for the Bloom
#           filter. Default: 0.001.
#      -smpthread <integer>: Number of threads used in a SMP node.
#           Default: number of cores in each SMP node.
#      -load <prefix>: Skip the solid k-mer finding step and load
#           pre-built Bloom filter data. When BLESS is executed with
#           "-prefix <prefix>" option <prefix>.bf.data and
#           <prefix>.bf.size are generated. If "-load <prefix>" option
#           is used, BLESS does not try to fiile solid k-mers in inputs
#           and load the existing Bloom filter data files.
#      -max_mem <4-1024>: Set maximum memory usage for KMC. Default: 4.
#      -notrim: Turn trimming off.
#      -seed <integer>: Set a seed for random number generation.
#           Default: 0.

# 3. EXAMPLES
#      1) PAIRED INPUT READS
#      bless -read1 in1.fastq -read2 in2.fastq -prefix \
#           directory/prefix -kmerlength 31

#      2) UNPAIRED INPUT READ
#      bless -read in.fastq -prefix directory/prefix -kmerlength 31

# 4. CONTACT
#      Yun Heo (yunheo1@illinois.edu)