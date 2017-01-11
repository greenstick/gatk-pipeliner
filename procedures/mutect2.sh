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

    # Additional Arguments

        -j=*|--jar=*)
        jar="${i#*=}"
        shift
        ;;
        -e=*|--contest=*)
        contestOpt="${i#*=}"
        shift # Run ContEst (0 or 1)
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
ncoresDef="8"
memoryDef="8G"
contestDef="true"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
contest=${contestOpt:-$contestDef}
 
printf "\nPARAMETERS:
GATK Directory      = $jar
Run ContEst         = $contest
Reference Directory = $reference
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Memory              = $memory
Cores               = $ncores
\n\n"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Contamination Estimation & Mutect2 Script"

echo "$recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam"

cd $dataDir

echo $(pwd)
#
# ContEst
#

if [ "$contest" = "true" ]; then

    printf "\n\nContEst Start"
    printf "\n\nCommand:\njava -Xmx$memory \
    -jar $jar -T ContEst \
    --precision 0.001 \
    -R $reference/Homo_sapiens_assembly19.fasta \
    -I:eval $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
    -I:genotype $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
    -pf $reference/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
    -isr INTERSECTION \
    --population ALL \
    --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
    -o $recalDir/logs/contest/cont_est_recal_$experiment.txt\n"
    java -Xmx$memory \
    -jar $jar -T ContEst \
    --precision 0.001 \
    -R $reference/Homo_sapiens_assembly19.fasta \
    -I:eval $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
    -I:genotype $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
    -pf $reference/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
    -isr INTERSECTION \
    --population ALL \
    --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
    -o $recalDir/logs/contest/cont_est_recal_$experiment.txt
    printf "\n\nContEst Complete"

fi

#
# Mutect2
#

printf "\n\nMuTect2 Start"
printf "\n\nCommand:\njava -Xmx$memory \
-jar $jar -T MuTect2 \
-R $reference/Homo_sapiens_assembly19.fasta \
-I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
-I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
--dbsnp $reference/dbsnp_138.hg19_modified.vcf \
--cosmic $reference/b37_cosmic_v54_120711_modified.vcf \
--tumor_lod 10.0 \
--contamination_fraction_to_filter 0.01 \
-o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.raw.snps.indels.vcf \
--log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
-nct $ncores\n"
java -Xmx$memory \
-jar $jar -T MuTect2 \
-R $reference/Homo_sapiens_assembly19.fasta \
-I:tumor $recalDir/$fileprefix.$subset.tumor.$experiment.$parameters.$qualitymodel.bam \
-I:normal $recalDir/$fileprefix.$subset.normal.$experiment.$parameters.$qualitymodel.bam \
--dbsnp $reference/dbsnp_138.hg19_modified.vcf \
--cosmic $reference/b37_cosmic_v54_120711_modified.vcf \
--tumor_lod 10.0 \
--contamination_fraction_to_filter 0.01 \
-o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
--log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
--graphOutput $recalDir/logs/mutect2/assembly_graph_info.txt \
-nct $ncores
printf "\n\nMuTect2 Complete"

# 
# Copy VCFS to User I/O Directory
# 

printf "\n\nCopying VCFs to I/O Directory..."
cp $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf /home/users/$USER/io/
printf "\n\nDone"

printf "\n\nDone\n"

