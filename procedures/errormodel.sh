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
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

format_status "PARAMETERS:
Reference Directory = $PIPELINE_REF
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

format_status "Running Error Model Script"

# A Really Long & Dirty Conditional
if  [ "$experiment" = "bayeshammer" ] || \
    [ "$experiment" = "blessec" ]     || \
    [ "$experiment" = "bloocoo" ]     || \
    [ "$experiment" = "decgpu" ]      || \
    [ "$experiment" = "karect" ]      || \
    [ "$experiment" = "kgem" ]        || \
    [ "$experiment" = "musket" ]      || \
    [ "$experiment" = "quorum" ]      || \
    [ "$experiment" = "rcorrector" ]  || \
    [ "$experiment" = "seecer" ]      || \
    [ "$experiment" = "nomodel" ]; then 
    echo "meow"

    state="$fileprefix.$subset.$condition.$experiment.$parameters:ERRORMODEL:1"
    if state_registered $state; then
        echo "moo"

        if [ "$experiment" = "nomodel" ]; then
            echo "bark"

            format_status "Copying Read FASTQ Files to Modeled Directory..."
            files=$(echo $(ls $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq))
            
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$dataDir/fastq/split/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:1"
                    
                    # Run Command
                    if state_registered $substate; then

                        format_status "Command:\ncp $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq"
                        cp $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq
                    
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script

        else

            echo "oink"
            format_status "Copying Read FASTQ Files to Pre-Alignment Directory..."
            format_status "Command:\ncp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/"
            cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/
        
        fi

        # Update State on Exit
        register_state $? $state
        format_status "FASTQ Copy Complete"
    fi

fi

#
# Delegate Args & Call Experiment Script
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters:ERRORMODEL:2"
if state_registered $state; then

    case "$experiment" in

        "bayeshammer")

            format_status "Running Bayes Hammer"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then

                        format_status "Command:\nsource $proceduresDir/models/bayeshammer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/bayeshammer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                    
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &

            done
            wait # Prevent Premature Exiting of Script
        ;;

        "blessec")

            format_status "Running Bless-EC"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        format_status "Command:\nsource $proceduresDir/models/blessec.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/blessec.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                    
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "bloocoo")

            format_status "Running Bloocoo"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        format_status "Command:\nsource $proceduresDir/models/bloocoo.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/bloocoo.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory

                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;
        "decgpu")

            format_status "Running Dec-GPU"
             # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        format_status "Command:\nsource $proceduresDir/models/decgpu.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/decgpu.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "karect")

            format_status "Running Karect"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        format_status "Command:\nsource $proceduresDir/models/karect.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/karect.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "kgem")

            format_status "Running KGEM"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        # java Xmx$memory -jar $ERIF
                        # java Xmx$memory -jar $KGEM
                        format_status "Command:\nsource $proceduresDir/models/kgem.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/kgem.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "musket")

            format_status "Running Musket"
             # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        format_status "Command:\nsource $proceduresDir/models/musket.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/musket.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "quorum")
            
            format_status "Running Quorum"
             # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        format_status "Command:\nsource $proceduresDir/models/quorum.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/quorum.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate
                        
                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "rcorrector")

            format_status "Running Rcorrector"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        # perl $RCORRECTOR
                        format_status "Command:\nsource $proceduresDir/models/rcorrector.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/rcorrector.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "seecer")

            format_status "Running Seecer"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        # source $SEECER
                        format_status "Command:\nsource $proceduresDir/models/seecer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/seecer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "shorah")

            format_status "Running SHoRAH"
            # Retrieve Files to Process
            files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
            for file in $files
                # In Parallel
                do (
                    # Extract Read Group to Pass Through
                    suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                    readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                    substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:2"
                    
                    # Run Command
                    if state_registered $substate; then
                        # python $SHORAH
                        format_status "Command:\n \
                        source $proceduresDir/models/shorah.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                        source $proceduresDir/models/shorah.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory
                        
                        # Check for failed parallel call
                        register_state $? $substate

                    fi
                ) &
            done
            wait # Prevent Premature Exiting of Script
        ;;

        "nomodel")
            format_status "No Model Selected"
        ;;

        "norealign")
            format_status "No Realignment Selected"
            format_status "Nothing to do...Exiting"
        ;;

        # Catch Any Invalid Error Models & Output Error
        *)

        format_status "Invalid Experiment (Error Model) Parameter: $experiment"
        ;;

    esac

    # Update State on Exit
    register_state $? $state
    format_status "Error Model Complete"

fi

format_status "Done"
