#! /usr/bin/bash

# Exit on First Error
set -o errexit

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

    # Additional Arguments

        -j=*|--jar=*)
        jar="${i#*=}"
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
ncoresDef="12"
memoryDef="6G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
ncores=${ncoresOpt:-$ncoresDef}

# Get Max Allowable Memory
allocMemory=$(echo "$memory" | sed "s|[GMKgmk]||")
allocSize=$(echo "$memory" | sed "s|[0-9]*||")
maxMemory=$(($allocMemory * $ncores))$allocSize

printf "\nPARAMETERS:
GATK Directory      = $jar
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
\n\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# 
# Run BQSR
# 

printf "\n\nRunning BQSR Script"

if [ "$qualitymodel" = "nobqsr" ]; then

    printf "\n\nCopying BAM & BAI Files to Model Directory...\n"
    cp $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.ba* $recalDir
    # Rename to Maintain Consistency
    mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam
    mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam.bai $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam.bai
    printf "\n\nDone"

fi

if [ "$qualitymodel" = "bqsr" ]; then

    #
    # BQSR - Step 1
    #

    printf "\n\nBQSR - Step 1 Start"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $jar -T BaseRecalibrator \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
    -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
    -o $recalDir/logs/bqsr/recal_data_$condition.table \
    --log_to_file $recalDir/logs/bqsr/log_$condition-recal1.txt \
    -nct $ncores\n"
    java -Xmx$memory \
    -jar $jar -T BaseRecalibrator \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
    -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
    -o $recalDir/logs/bqsr/recal_data_$condition.table \
    --log_to_file $recalDir/logs/bqsr/log_$condition-recal1.txt \
    -nct $ncores
    printf "\n\nBQSR - Step 1 Complete"

    #
    # BQSR - Step 2
    #

    printf "\n\nBQSR - Step 2 Start"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $jar -T BaseRecalibrator \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
    -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
    -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
    -o $recalDir/logs/bqsr/post_recal_data_$condition.table \
    --log_to_file $recalDir/logs/bqsr/log_$condition-recal2.txt \
    -nct $ncores\n"
    java -Xmx$memory \
    -jar $jar -T BaseRecalibrator \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
    -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
    -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
    -o $recalDir/logs/bqsr/post_recal_data_$condition.table \
    --log_to_file $recalDir/logs/bqsr/log_$condition-recal2.txt \
    -nct $ncores
    printf "\n\nBQSR - Step 2 Complete"

    #
    # BQSR - Step 3
    #

    printf "\n\nBQSR - Step 3 Start"
    printf "\n\nCommand:\njava -Xmx$maxMemory \
    -jar $jar -T AnalyzeCovariates \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -before $recalDir/logs/bqsr/recal_data_$condition.table \
    -after $recalDir/logs/bqsr/post_recal_data_$condition.table \
    -plots $recalDir/logs/bqsr/recalibration_plots_$condition.pdf \
    --log_to_file $recalDir/logs/bqsr/log_$condition-generateplots.txt\n"
    java -Xmx$memory \
    -jar $jar -T AnalyzeCovariates \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -before $recalDir/logs/bqsr/recal_data_$condition.table \
    -after $recalDir/logs/bqsr/post_recal_data_$condition.table \
    -plots $recalDir/logs/bqsr/recalibration_plots_$condition.pdf \
    --log_to_file $recalDir/logs/bqsr/log_$condition-generateplots.txt
    printf "\n\nBQSR - Step 3 Complete"

    #
    # BQSR - Step 4
    #

    printf "\n\nBQSR - Step 4 Start"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $jar -T PrintReads \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
    -o $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
    --log_to_file $recalDir/logs/bqsr/log_$condition-printreads.txt \
    -nct $ncores\n"
    java -Xmx$memory \
    -jar $jar -T PrintReads \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
    -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
    -o $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
    --log_to_file $recalDir/logs/bqsr/log_$condition-printreads.txt \
    -nct $ncores
    printf "\n\nBQSR - Step 4 Complete"

    #
    # Create New BAM Index (Not Strictly Neccessary, Could Rename the Old One, but Does Appear to Prevent Downstream Errors)
    #

    printf "\n\nIndexing BAM Output"
    printf "\n\nCommand:\nsamtools index $recalDir/$fileprefix.$subset.$condition.$experiment.$qualitymodel.bam\n"
    samtools index $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam
    rm $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam.bai
    printf "\n\nBAM Indexing Complete"

fi

printf "\n\nDone\n"
