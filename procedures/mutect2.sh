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

        -C=*|--contamination=*)
        contaminationOpt="${i#*=}"
        shift # Sample Contamination
        ;;
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
        -d=*|--debug=*)
        debugOpt="${i#*=}"
        shift # Trigger Debugging Available in Tools
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
readsDef="0"
debugDef=false
contaminationDef=true

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
reads=${readsOpt:-$readsDef}
debug=${debugOpt:-$debugDef}
contamination=${contaminationOpt:-$contaminationDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize
 
# Max Reads in RAM - 200,000 per GB
maxReads=$((allocMemory * $reads))

# Use Default n Reads or User Defines
if [ "$maxReads" = "0" ]; then
    maxReads="default"
    readbuffersize=""
else
    readbuffersize="--read_buffer_size "$maxReads
fi

printf "\nPARAMETERS:
GATK Directory      = $GATK
Data File Prefix    = $fileprefix
Data Subset         = $subset
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Contamination       = $contamination
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Max Reads in Memory = $maxReads
Debug               = $debug
\n"

# Set Directories
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# Tool Specific Debugging - GATK
monitorThreads=""
performanceLog=false
loggingLevel="INFO"

if $debug; then 
    monitorThreads="--monitorThreadEfficiency"
    performanceLog=true
    loggingLevel="DEBUG"
fi

format_status "Running Contamination Estimation & Mutect2 Script"

#
# If Contamination Not Specified, Run ContEst
#

if [ "$contamination" = "false" ]; then

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:1"
    if !(has_state $state); then

        format_status "ContEst Start"
        format_status "Command:\njava -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T ContEst \
        --precision 0.001 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:eval $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:genotype $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        -pf $PIPELINE_REF/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
        -isr INTERSECTION \
        --population ALL \
        --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
        -o $recalDir/logs/contest/cont_est_recal_$experiment.txt \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        java -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T ContEst \
        --precision 0.001 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:eval $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:genotype $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        -pf $PIPELINE_REF/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
        -isr INTERSECTION \
        --population ALL \
        --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
        -o $recalDir/logs/contest/cont_est_recal_$experiment.txt \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize # Additional Optional Args
    
        # Update State on Exit
        put_state $? $state
        format_status "ContEst Complete"

    fi

    #
    # Mutect2
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2"
    if !(has_state $state); then

        # Get Contamination
        contamination=$(awk -F '\t' 'NR >=2 {print $4}'  $recalDir/logs/contest/cont_est_recal_$experiment.txt)
        format_status "Proportion Contamination: $contamination"

        format_status "MuTect2 Start"
        format_status "Command:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect2 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect2 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
        --graphOutput $recalDir/logs/mutect2/assembly_graph_info.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize


        # Update State on Exit
        put_state $? $state
        format_status "MuTect2 Complete"

    fi

else

    #
    # Mutect2
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2"
    if !(has_state $state); then

        format_status "MuTect2 Start"
        format_status "Command:\njava -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect2 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect2 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
        --graphOutput $recalDir/logs/mutect2/assembly_graph_info.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize # Additional Optional Args


        # Update State on Exit
        put_state $? $state
        format_status "MuTect2 Complete"

    fi

fi

# 
# Copy VCFS to User I/O Directory
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:3"
if !(has_state $state); then

    format_status "Copying VCFs to I/O Directory"
    cp $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf /home/users/$USER/io/

    # Update State on Exit
    put_state $? $state
    format_status "VCFs Copied to I/O Directory"

fi

format_status "Done"

