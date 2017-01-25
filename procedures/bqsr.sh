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

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

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
\n\n"

# Set Directories
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# 
# Run BQSR
# 

printf "\n\nRunning BQSR Script"

if [ "$qualitymodel" = "nobqsr" ]; then

    #
    # No BQSR - Step 1
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:NOBQSR:1" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        failures=0
        printf "\n\nCopying BAM & BAI Files to Model Directory...\n"
        cp $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.ba* $recalDir
        # Check for failed command
        if [ $? != 0 ]; then
            failures=$((failures + 1))
        fi
        # Rename to Maintain Consistency
        mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam
        # Check for failed command
        if [ $? != 0 ]; then
            failures=$((failures + 1))
        fi
        mv $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.bam.bai $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam.bai
        # Check for failed command
        if [ $? != 0 ]; then
            failures=$((failures + 1))
        fi
        
        # Update State on Exit
        if [ $failures = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:NOBQSR:1" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nDone"
        else
            printf "\n\n$failures Failures, Exiting - $fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:NOBQSR:1"
            exit 1
        fi

    fi

fi

if [ "$qualitymodel" = "bqsr" ]; then

    #
    # BQSR - Step 1
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:1" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nBQSR - Step 1 Start"
        printf "\n\nCommand:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -o $recalDir/logs/bqsr/recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal1.txt \
        -nct $ncores\n"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -o $recalDir/logs/bqsr/recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal1.txt \
        -nct $ncores

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:1" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBQSR - Step 1 Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:1"
            exit $statuscode
        fi

    fi

    #
    # BQSR - Step 2
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:2" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nBQSR - Step 2 Start"
        printf "\n\nCommand:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T BaseRecalibrator \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -knownSites $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        -knownSites $PIPELINE_REF/Mills_and_1000G_gold_standard.indels.hg19.sites_modified.vcf \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/logs/bqsr/post_recal_data_$condition.table \
        --log_to_file $recalDir/logs/bqsr/log_$condition-recal2.txt \
        -nct $ncores\n"
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
        -nct $ncores

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:2" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBQSR - Step 2 Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:2"
            exit $statuscode
        fi

    fi

    #
    # BQSR - Step 3
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:3" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nBQSR - Step 3 Start"
        printf "\n\nCommand:\njava -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T AnalyzeCovariates \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -before $recalDir/logs/bqsr/recal_data_$condition.table \
        -after $recalDir/logs/bqsr/post_recal_data_$condition.table \
        -plots $recalDir/logs/bqsr/recalibration_plots_$condition.pdf \
        --log_to_file $recalDir/logs/bqsr/log_$condition-generateplots.txt\n"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T AnalyzeCovariates \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -before $recalDir/logs/bqsr/recal_data_$condition.table \
        -after $recalDir/logs/bqsr/post_recal_data_$condition.table \
        -plots $recalDir/logs/bqsr/recalibration_plots_$condition.pdf \
        --log_to_file $recalDir/logs/bqsr/log_$condition-generateplots.txt

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:3" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBQSR - Step 3 Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:3"
            exit $statuscode
        fi

    fi

    #
    # BQSR - Step 4
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:4" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nBQSR - Step 4 Start"
        printf "\n\nCommand:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T PrintReads \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
        --log_to_file $recalDir/logs/bqsr/log_$condition-printreads.txt \
        -nct $ncores\n"
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T PrintReads \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        -BQSR $recalDir/logs/bqsr/recal_data_$condition.table \
        -o $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam \
        --log_to_file $recalDir/logs/bqsr/log_$condition-printreads.txt \
        -nct $ncores

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:4" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBQSR - Step 4 Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:4"
            exit $statuscode
        fi

    fi

    #
    # Create New BAM Index (Not Strictly Neccessary, Could Rename the Old One, but Does Appear to Prevent Downstream Errors)
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    grep -q "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:5" $PIPELINE_HOME/pipeline.state
    if [ $? != 0 ]; then

        printf "\n\nIndexing BAM Output"
        printf "\n\nCommand:\nsamtools index $recalDir/$fileprefix.$subset.$condition.$experiment.$qualitymodel.bam\n"
        samtools index $recalDir/$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel.bam

        # Update State on Exit
        statuscode=$?
        if [ $statuscode = 0 ]; then
            # Export Pipeline State
            echo "$fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:5" >> $PIPELINE_HOME/pipeline.state
            printf "\n\nBAM Indexing Complete"
        else
            printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$condition.$experiment.$parameters.$qualitymodel:BQSR:5"
            exit $statuscode
        fi

    fi

fi

printf "\n\nDone\n"
