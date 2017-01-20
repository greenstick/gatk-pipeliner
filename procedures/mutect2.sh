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
ncoresDef="10"
memoryDef="8G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize
 
printf "\nPARAMETERS:
GATK Directory      = $GATK
Data File Prefix    = $fileprefix
Data Subset         = $subset
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
\n\n"

# Set Directories
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel

printf "\n\nRunning Contamination Estimation & Mutect2 Script"

#
# ContEst
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:1" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    printf "\n\nContEst Start"
    printf "\n\nCommand:\njava -Xmx$maxMemory \
    -jar $GATK -T ContEst \
    --precision 0.001 \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I:eval $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
    -I:genotype $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
    -pf $PIPELINE_REF/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
    -isr INTERSECTION \
    --population ALL \
    --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
    -o $recalDir/logs/contest/cont_est_recal_$experiment.txt\n"
    java -Xmx$memory \
    -jar $GATK -T ContEst \
    --precision 0.001 \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I:eval $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
    -I:genotype $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
    -pf $PIPELINE_REF/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
    -isr INTERSECTION \
    --population ALL \
    --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
    -o $recalDir/logs/contest/cont_est_recal_$experiment.txt
    
    # Update State on Exit
    statuscode=$?
    if [ $statuscode = 0 ]; then
        # Export Pipeline State
        echo "$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:1" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nContEst Complete"
    else
        printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:1"
    fi

fi

#
# Mutect2
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    printf "\n\nMuTect2 Start"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $GATK -T MuTect2 \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
    -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
    --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
    --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
    --tumor_lod 10.0 \
    --contamination_fraction_to_filter 0.01 \
    -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.raw.snps.indels.vcf \
    --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
    -nct $ncores\n"
    java -Xmx$memory \
    -jar $GATK -T MuTect2 \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
    -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
    --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
    --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
    --tumor_lod 10.0 \
    --contamination_fraction_to_filter 0.01 \
    -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
    --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
    --graphOutput $recalDir/logs/mutect2/assembly_graph_info.txt \
    -nct $ncores

    # Update State on Exit
    statuscode=$?
    if [ $statuscode = 0 ]; then
        # Export Pipeline State
        echo "$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nMuTect2 Complete"
    else
        printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2"
    fi

fi

# 
# Copy VCFS to User I/O Directory
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
grep -q "$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:3" $PIPELINE_HOME/pipeline.state
if [ $? != 0 ]; then

    printf "\n\nCopying VCFs to I/O Directory..."
    cp $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf /home/users/$USER/io/

    # Update State on Exit
    statuscode=$?
    if [ $statuscode = 0 ]; then
        # Export Pipeline State
        echo "$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:3" >> $PIPELINE_HOME/pipeline.state
        printf "\n\nVCFs Copied to I/O Directory"
    else
        printf "\n\nUnexpected Exit $statuscode - $fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:3"
    fi

fi

printf "\n\nDone\n"

