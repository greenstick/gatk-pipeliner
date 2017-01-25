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
ncoresDef="12"
memoryDef="6G"

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize
 
format_status "PARAMETERS: 
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

format_status "Running Merge BAMs Script"

# 
# Inserting Read Groups
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state"$fileprefix.$subset.$condition.$experiment.$parameters:MERGEALIGNMENT:1"
if [ (state_registered $state) != 0 ]; then

    format_status "Samtools AddReplaceRG"
    # Retrieve Files to Process
    files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

    for file in $files
        # In Parallel
        do (
            # Get Read Group to Process
            suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
            readgroup=$(echo "$suffix" | sed "s|.bam$||")
            substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:MERGEALIGNMENT:1"
            
            # Run Command
            if [ (state_registered $substate) != 0 ]; then
                # Get Read Group Arguments to Pass to Samtools
                rgArgs=$(samtools view -H $file | grep '@RG' | awk -F '\t' '{print $2,$3,$4,$5,$6,$7,$8}' | sed "s|[A-Z][A-Z]:[a-zA-Z0-9\.\-\:]*|-r &|g")
                
                # Check for failed parallel call
                register_state $? $substate

                format_status "Command:\nsamtools addreplacerg $rgArgs $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.bam" 
                # Insert Read Groups into New BAM - WARNING USES EVAL
                eval "samtools addreplacerg ${rgArgs[@]} $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.bam"
            
                # Check for failed parallel call
                register_state $? $substate

            fi
        ) &
    done
    wait # Prevent Premature Exiting of Script

    # Update State on Exit
    register_state $? $state
    format_status "Samtools AddReplaceRG Complete"

fi

# 
# Merge BAMs to Single BAM
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state"$fileprefix.$subset.$condition.$experiment.$parameters:MERGEALIGNMENT:2"
if [ (state_registered $state) != 0 ]; then

    # Retrieve Files to Process
    files=$(echo $(ls $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.*.bam))

    format_status "Samtools Merge"
    format_status "Command:\nsamtools merge -r $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files"
    samtools merge -r -f $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files
    
    # Update State on Exit
    register_state $? $state
    format_status "Samtools Merge Complete"

fi

format_status "Done"
