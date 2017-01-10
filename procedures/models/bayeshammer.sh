#! /usr/bin/bash

# Assign Arguments
for i in "$@"
    do case $i in

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
        -g=*|--readgroup=*)
        readgroup="${i#*=}"
        shift # Access & Write Files With This Read Group
        ;;

    # Additional Arguments

        -a=*|--arguments=*)
        arguments="${i#*=}"
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

# Set Directories
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# 
# Run BQSR
# 

cd $dataDir

printf "\n\nRunning BayesHammer"

# Quorum

printf "\n\nCommand:\n \
$BAYESHAMMER --prefix experiments/bayeshammer/$prefix.$subset.$condition.bayeshammer --no-discard --debug"

$BAYESHAMMER \
-o $paramDir \
--12 fastq/$fileprefix.$subset.$condition.$readgroup.fastq \
--only-error-correction \
--threads 4 \
--memory 92

cd ../procedures

printf "\n\nDone\n"

#
# BayesHammer Arguments
#

# Usage: /home/users/cordier/packages/bayes-hammer/SPAdes-3.9.0-Linux/bin/spades.py [options] -o <output_dir>

# Basic options:
# -o  <output_dir>    directory to store all the resulting files (required)
# --sc            this flag is required for MDA (single-cell) data
# --meta          this flag is required for metagenomic sample data
# --rna           this flag is required for RNA-Seq data
# --plasmid       runs plasmidSPAdes pipeline for plasmid detection
# --iontorrent        this flag is required for IonTorrent data
# --test          runs SPAdes on toy dataset
# -h/--help       prints this usage message
# -v/--version        prints version

# Input data:
# --12    <filename>  file with interlaced forward and reverse paired-end reads
# -1  <filename>  file with forward paired-end reads
# -2  <filename>  file with reverse paired-end reads
# -s  <filename>  file with unpaired reads
# --pe<#>-12  <filename>  file with interlaced reads for paired-end library number <#> (<#> = 1,2,..,9)
# --pe<#>-1   <filename>  file with forward reads for paired-end library number <#> (<#> = 1,2,..,9)
# --pe<#>-2   <filename>  file with reverse reads for paired-end library number <#> (<#> = 1,2,..,9)
# --pe<#>-s   <filename>  file with unpaired reads for paired-end library number <#> (<#> = 1,2,..,9)
# --pe<#>-<or>    orientation of reads for paired-end library number <#> (<#> = 1,2,..,9; <or> = fr, rf, ff)
# --s<#>      <filename>  file with unpaired reads for single reads library number <#> (<#> = 1,2,..,9)
# --mp<#>-12  <filename>  file with interlaced reads for mate-pair library number <#> (<#> = 1,2,..,9)
# --mp<#>-1   <filename>  file with forward reads for mate-pair library number <#> (<#> = 1,2,..,9)
# --mp<#>-2   <filename>  file with reverse reads for mate-pair library number <#> (<#> = 1,2,..,9)
# --mp<#>-s   <filename>  file with unpaired reads for mate-pair library number <#> (<#> = 1,2,..,9)
# --mp<#>-<or>    orientation of reads for mate-pair library number <#> (<#> = 1,2,..,9; <or> = fr, rf, ff)
# --hqmp<#>-12    <filename>  file with interlaced reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)
# --hqmp<#>-1 <filename>  file with forward reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)
# --hqmp<#>-2 <filename>  file with reverse reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)
# --hqmp<#>-s <filename>  file with unpaired reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)
# --hqmp<#>-<or>  orientation of reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9; <or> = fr, rf, ff)
# --nxmate<#>-1   <filename>  file with forward reads for Lucigen NxMate library number <#> (<#> = 1,2,..,9)
# --nxmate<#>-2   <filename>  file with reverse reads for Lucigen NxMate library number <#> (<#> = 1,2,..,9)
# --sanger    <filename>  file with Sanger reads
# --pacbio    <filename>  file with PacBio reads
# --nanopore  <filename>  file with Nanopore reads
# --trusted-contigs   <filename>  file with trusted contigs
# --untrusted-contigs <filename>  file with untrusted contigs

# Pipeline options:
# --only-error-correction runs only read error correction (without assembling)
# --only-assembler    runs only assembling (without read error correction)
# --careful       tries to reduce number of mismatches and short indels
# --continue      continue run from the last available check-point
# --restart-from  <cp>    restart run with updated options and from the specified check-point ('ec', 'as', 'k<int>', 'mc')
# --disable-gzip-output   forces error correction not to compress the corrected reads
# --disable-rr        disables repeat resolution stage of assembling

# Advanced options:
# --dataset   <filename>  file with dataset description in YAML format
# -t/--threads    <int>       number of threads
#                 [default: 16]
# -m/--memory <int>       RAM limit for SPAdes in Gb (terminates if exceeded)
#                 [default: 250]
# --tmp-dir   <dirname>   directory for temporary files
#                 [default: <output_dir>/tmp]
# -k      <int,int,...>   comma-separated list of k-mer sizes (must be odd and
#                 less than 128) [default: 'auto']
# --cov-cutoff    <float>     coverage cutoff value (a positive float number, or 'auto', or 'off') [default: 'off']
# --phred-offset  <33 or 64>  PHRED quality offset in the input reads (33 or 64)
#                 [default: auto-detect]