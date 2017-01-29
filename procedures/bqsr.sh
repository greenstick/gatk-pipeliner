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
        -q=*|--qualitymodel=*)
        qualitymodel="${i#*=}"
        shift # Access & Write Files With This Quality Model
        ;;

    # Optional Arguments With Defaults

        -r=*|--reads=*)
        readsOpt="${i#*=}"
        shift # n Reads Per GB or Memory
        ;;
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
ncoresDef="10"
memoryDef="8G"
readsDef=150

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
reads=${readsOpt:-$readsDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

# Max Reads in RAM - 200,000 per GB
maxReads=$((allocMemory * $reads))

printf "\nPARAMETERS: 
GATK Directory      = $GATK
Reference Directory = $PIPELINE_REF
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Max Reads in Memory = $maxReads
\n"

# Set Directories
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# 
# Run BQSR
# 

format_status "Running BQSR Script"

if [ "$qualitymodel" = "nobqsr" ]; then

    #
    # No BQSR - Step 1
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition:NOBQSR:1"
    if !(has_state $state); then

        # Copy Files & Rename to Maintain Consistency
        format_status "Copying BAM & BAI Files to Model Directory..."
        format_status "Command:\ncp $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.ba* $recalDir \
        && mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
        && mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam.bai $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam.bai"
        cp $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.ba* $recalDir \
        && mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
        && mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam.bai $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam.bai
        
        # Update State on Exit
        put_state $? $state
        format_status "BAM & BAI Copy Complete"

    fi

fi

if [ "$qualitymodel" = "bqsr" ]; then

    #
    # BQSR - Step 1
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:1"
    if !(has_state $state); then

        format_status "BQSR - Step 1 Start"
        format_status "Command:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -o $recalDir/logs/bqsr/recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal1.txt \
        -nct $ncores \
        --lowMemoryMode \
        --read_buffer_size $maxReads"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -o $recalDir/logs/bqsr/recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal1.txt \
        -nct $ncores \
        --lowMemoryMode \
        --read_buffer_size $maxReads

        # Update State on Exit
        put_state $? $state
        format_status "BQSR - Step 1 Complete"

    fi

    #
    # BQSR - Step 2
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:2"
    if !(has_state $state); then

        format_status "BQSR - Step 2 Start"
        format_status "Command:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/logs/bqsr/post_recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal2.txt \
        -nct $ncores \
        --lowMemoryMode \
        --read_buffer_size $maxReads"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/logs/bqsr/post_recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal2.txt \
        -nct $ncores \
        --lowMemoryMode \
        --read_buffer_size $maxReads

        # Update State on Exit
        put_state $? $state
        format_status "BQSR - Step 2 Complete"

    fi

    #
    # BQSR - Step 3
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:3"
    if !(has_state $state); then

        format_status "BQSR - Step 3 Start"
        format_status "Command:\njava -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T AnalyzeCovariates \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -before $recalDir/logs/bqsr/recal_data_$condition.table \
        -after $recalDir/logs/bqsr/post_recal_data_$condition.table \
        -plots $recalDir/logs/bqsr/recalibration_plots_$condition.pdf \
        --log_to_file $recalDir/logs/bqsr/log_$condition-generateplots.txt"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T AnalyzeCovariates \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -before $recalDir/logs/bqsr/recal_data_$condition.table \
        -after $recalDir/logs/bqsr/post_recal_data_$condition.table \
        -plots $recalDir/logs/bqsr/recalibration_plots_$condition.pdf \
        --log_to_file $recalDir/logs/bqsr/log_$condition-generateplots.txt

        # Update State on Exit
        put_state $? $state
        format_status "BQSR - Step 3 Complete"

    fi

    #
    # BQSR - Step 4
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:4"
    if !(has_state $state); then

        format_status "BQSR - Step 4 Start"
        format_status "Command:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T PrintReads \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
        --log_to_file $recalDir/logs/bqsr/log_$condition-printreads.txt \
        -nct $ncores \
        --read_buffer_size $maxReads"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T PrintReads \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
        --log_to_file $recalDir/logs/bqsr/log_$condition-printreads.txt \
        -nct $ncores \
        --read_buffer_size $maxReads

        # Update State on Exit
        put_state $? $state
        format_status "BQSR - Step 4 Complete"

    fi

    #
    # Create New BAM Index (Not Strictly Neccessary, Could Rename the Old One, but Does Appear to Prevent Downstream Errors)
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:5"
    if !(has_state $state); then

        format_status "Indexing BAM Output"
        format_status "Command:\nsamtools index $recalDir/$fileprefix.$subset.$condition.$experiment.$qualitymodel.bam"
        samtools index $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam

        # Update State on Exit
        put_state $? $state
        format_status "BAM Indexing Complete"

    fi

fi

format_status "Done"
