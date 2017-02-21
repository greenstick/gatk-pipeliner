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

    # Directory Cleanup (Voids All Other Parameters)

        --clean)
        cleanOpt=true
        ;;

    # Invalid Argument Handler

        *)
        # invalid option
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

# Defaults if No Arguments Passed
ncoresDef="16"
memoryDef="6G"
cleanDef=false
readsDef="0"
debugDef=false
contaminationDef=true

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
clean=${cleanOpt:-$cleanDef}
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
Do Cleanup          = $clean
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

format_status "Running Contamination Estimation & Mutect Script"

#
# If Contamination Not Specified, Run ContEst
#

if $contamination; then

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT:1"
    if !(has_state $state); then

        format_status "ContEst Start"
        # Define Command
        call="java -Xmx$maxMemory \
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
        # Print & Call
        format_status "Command:\n$call"
        eval $call


        # Update State on Exit
        put_state $? $state
        format_status "ContEst Complete"

    fi

    #
    # Mutect
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT:2"
    if !(has_state $state); then

        # Get Contamination
        contaminationPercent=$(awk -F '\t' 'NR >=2 {print $4}'  $recalDir/logs/contest/cont_est_recal_$experiment.txt)
        contamination=$(python -c "print($contaminationPercent/100)")
        format_status "Proportion Contamination: $contamination"

        format_status "MuTect Start"
        # Define Command

        call="$JAVA7 -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect/log_mutect_$experiment.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "MuTect Complete"

    fi

else

    #
    # Mutect
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT:2"
    if !(has_state $state); then

        format_status "MuTect Start"
        # Define Command
        call="$JAVA7 -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect/log_mutect_$experiment.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "MuTect Complete"

    fi

fi

# 
# Copy VCFS to User I/O Directory
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT:3"
if !(has_state $state); then

    format_status "Copying VCFs to I/O Directory"
    # Define Command
    call="cp $recalDir/logs/mutect/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf /home/users/$USER/io/"
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    put_state $? $state
    format_status "VCFs Copied to I/O Directory"

fi

format_status "Done"

