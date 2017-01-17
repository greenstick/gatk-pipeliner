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
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=$(echo "$memory" | sed "s|[GMKgmk]||")
allocSize=$(echo "$memory" | sed "s|[0-9]*||")
maxMemory=$(($allocMemory * $ncores))$allocSize

printf "\nPARAMETERS:
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

printf "\n\nRunning Error Model Script"

#
# Delegate Args & Call Experiment Script
#

case "$experiment" in

    "bayeshammer")

        printf "\n\nCopying Raw FastQ Files to Error Model Directory\n"
        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        printf "\n\nRunning Bayes Hammer"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                printf "\n\nRunning Picard Bam to FastQ"
                # python $BAYESHAMMER
                printf "\n\nCommand:\n \
                source $proceduresDir/models/bayeshammer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/bayeshammer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # $BLESSEC
                printf "\n\nCommand:\n \
                source $proceduresDir/models/blessec.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/blessec.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # $BLOOCOO
                printf "\n\nCommand:\n \
                source $proceduresDir/models/bloocoo.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/bloocoo.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # $DECGPU
                printf "\n\nCommand:\n \
                source $proceduresDir/models/decgpu.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/decgpu.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # $KARECT
                printf "\n\nCommand:\n \
                source $proceduresDir/models/karect.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/karect.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # java Xmx$memory -jar $ERIF
                # java Xmx$memory -jar $KGEM
                printf "\n\nCommand:\n \
                source $proceduresDir/models/kgem.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/kgem.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # $DECGPU
                printf "\n\nCommand:\n \
                source $proceduresDir/models/musket.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/musket.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # $QUORUM
                printf "\n\nCommand:\n \
                source $proceduresDir/models/quorum.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/quorum.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # perl $RCORRECTOR
                printf "\n\nCommand:\n \
                source $proceduresDir/models/rcorrector.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/rcorrector.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # bash $SEECER
                printf "\n\nCommand:\n \
                source $proceduresDir/models/seecer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/seecer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
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
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # python $SHORAH
                printf "\n\nCommand:\n \
                source $proceduresDir/models/shorah.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory\n"
                source $proceduresDir/models/shorah.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -q=$qualitymodel -n=$ncores -m=$memory
            ) &
        done
        wait # Prevent Premature Exiting of Script
    ;;

    "nomodel")
        printf "\n\nNo Model Selected"
        printf "\n\nMoving FastQ to Modeled Directory...\n"
        files=$(echo $(ls $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            # In Parallel
            do (
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
        printf "\n\nNothing to do...Exiting"
    ;;

    # Catch Any Invalid Error Models & Output Error
    *)

    printf "\n\nInvalid Experiment (Error Model) Parameter: $experiment"
    ;;

esac

printf "\n\nDone\n"
