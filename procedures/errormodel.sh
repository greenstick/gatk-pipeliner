#! /usr/bin/bash

# Assign Arguments
for i in "$@"
do
case $i in

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
ncoresDef="4"
memoryDef="8G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

printf "\nPARAMETERS:
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
dataDir=$PIPELINE_HOME/$subset
modelDir=$PIPELINE_HOME/$subset/model/$experiment
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

printf "\n\nRunning Error Model Script"

#
# Delegate Args & Call Experiment Script
#

# Model scripts output FASTQ to: 
# $paramDir/pre-align/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq

# And output logs to:
# $paramDir/logs

case "$experiment" in

    "bayeshammer")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Bayes Hammer"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    printf "\n\nRunning Picard Bam to FastQ"
                    # python $BAYESHAMMER
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/bayeshammer.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/bayeshammer.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
    ;;

    "blessec")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Bless-EC"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # $BLESSEC
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/blessec.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/blessec.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "bloocoo")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Bloocoo"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # $BLOOCOO
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/bloocoo.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/bloocoo.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;
    "decgpu")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Dec-GPU"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # $DECGPU
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/decgpu.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/decgpu.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "karect")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Karect"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # $KARECT
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/karect.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/karect.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "kgem")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning KGEM"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # java Xmx$memory -jar $ERIF
                    # java Xmx$memory -jar $KGEM
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/kgem.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/kgem.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "musket")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Musket"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # $DECGPU
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/musket.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/musket.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "quorum")
        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Quorum"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # $QUORUM
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/quorum.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/quorum.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "rcorrector")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Rcorrector"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # perl $RCORRECTOR
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/rcorrector.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/rcorrector.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "seecer")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Seecer"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # bash $SEECER
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/seecer.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/seecer.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "shorah")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning SHoRAH"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    # python $SHORAH
                    arguments=""
                    printf "\n\nCommand:\n \
                    source models/shorah.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments\n"
                    source models/shorah.sh -r=$reference -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory -a=$arguments
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "nomodel")
        printf "\n\nNo Model Selected"
        printf "\n\nMoving FastQ to Modeled Directory...\n"
        files=$(echo $(ls $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do
                # In Parallel
                (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$dataDir/fastq/split/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    cp $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "norealign")

        printf "\n\nNo Realignment Selected"
        printf "\n\nMoving FastQ to Merged Directory...\n"
        files=$(echo $(ls $dataDir/downloaded/$fileprefix.$subset.$condition.ba*))
        for file in $files
            do
                # In Parallel
                (
                    suffix=$(echo "$file" | sed "s|$dataDir/downloaded/$fileprefix.$subset.$condition.||")
                    cp $file $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameter.$suffix
                ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    *)

    printf "\n\nInvalid Experiment (Error Model) Parameter: $experiment"
    ;;
esac

printf "\n\nDone\n"
