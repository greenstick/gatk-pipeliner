# Dependencies
#
# Java          1.8.0_112
# GATK          3.6+
# Picard Tools  2.6.0
# Samtools      1.3.1
# BWA           7.15
# Seqtk         1.2-r94

# Processing Steps:
# 
# Step      File        Function        Output Directory
# 1.        bam         -- bamtofastq   -> fastq/split
# 2.        fastq       -- fastqc       -> fastq/split
# 3.        fastq       -- errormodel*  -> model/{}/param/{}/modeled
# 4.        fastq       -- fastqc*      -> model/{}/param/{}/pre-align
# 5.        fastq       -- bwa          -> model/{}/param/{}/post-align
# 6.        fastq       -- mergealign   -> model/{}/param/{}/merged
# 7.        bam         -- fastqc       -> model/{}/param/{}/post-align
# 9.        bam         -- markdup      -> model/{}/param/{}/post-align
# 9.        bam         -- recal*       -> model/{}/param/{}/recal/{}/log/bqsr
# 10.       bam         -- contest      -> model/{}/param/{}/recal/{}/log/contest
# 11.       bam         -- mutect2      -> model/{}/param/{}/recal/{}/log/mutect2


#
# References / Locations
#

printf "\nSetting Pipeline Variables...\n"

# Set File System Root Directory
root_dir=/home/exacloud/lustre1

# Set User Directory
user_dir=/home/users/cordier

# Set Pipeline Directory
pipeline_dir=$root_dir/MRD_aml/ForBackup/pipeline

# Set Reference Genome Directory/File
reference_dir=$root_dir/MRD_aml/ForBackup/pipeline/reference/hg19

# Set Output Directory
output_dir=$user_dir/io # Set to home directory off lustre FS for FTP access

# Set Logging Directory
logging_dir=$pipeline_dir/logs

#
# Exports
#

while true; do
    IFS= read -p "\nWrite pipeline exports to ~/.bash_profile (default is yes)? " yn
    case $yn in
        [Yy]* ) (
                printf "Exporting Directory Locations...\n"

                # Header for Bash Profile
                echo -e '\n#\n# Pipeline Management\n#' >> ~/.bash_profile

                # Write Pipeline Home Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_HOME=$pipeline_dir" >> ~/.bash_profile
                # Set Pipeline Home Directory
                export PIPELINE_HOME=$pipeline_dir

                # Write Pipeline Output Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_OUT=$output_dir" >> ~/.bash_profile
                # Set Pipeline Output Directory
                export PIPELINE_OUT=$output_dir
                
                # Write Pipeline Reference Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_REF=$reference_dir" >> ~/.bash_profile
                # Set Pipeline Reference Directory
                export PIPELINE_REF=$reference_dir
                
                # Write Pipeline Logging Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_LOG=$logging_dir" >> ~/.bash_profile
                # Set Pipeline Logging Directory
                export PIPELINE_LOG=$logging_dir

            )
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Please answer yes [y] or no [n].";;
    esac
done

#
# Logging Management
#

while true; do
    IFS= read -p "\nEnforce TMUX session logging via ~/.bash_profile (default is yes)? " yn
    case $yn in
        [Yy]* )  
            # Inject Logging Automation Script to ~/.bash_profile
            echo -e '\n# Pipeline Logging Management\nif [[ $TERM = "screen" ]] && [[ $(ps -p $PPID -o comm=) = "tmux" ]] && [[ $PWD/ = $PIPELINE_HOME/* ]]; then\n\t# Prompt User / Read Input\n\tIFS= read -p "Enter job logging name: " name\n\tmkdir $PIPELINE_HOME/logs 2> /dev/null\n\tlogname="tmux.$name.$(date "+%d%m%Y%H%M%S").log"\n\tscript -f $PIPELINE_LOG/${logname}\n\texit\nfi' >> ~/.bash_profile
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Please answer yes [y] or no [n].";;
    esac
done

#
# Scaffolding
#

printf "\nScaffolding Directories...\n"

# Scaffold Pipeline Directory & Change Directory
cd $root_dir && mkdir -p $pipeline_dir && cd $pipeline_dir

# Pipeline Scaffolding
mkdir -p set{1..6}/model/{bayeshammer,bless,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{logs,post-align,pre-align,markdup,recal,merged,modeled}
mkdir -p set{1..6}/model/{bayeshammer,bless,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/modeled/{pairs,indexes}
mkdir -p set{1..6}/model/{bayeshammer,bless,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{post-align,pre-align}/fastq
mkdir -p set{1..6}/model/{bayeshammer,bless,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/{bqsr,nobqsr}/logs/{contest,mutect2}
mkdir -p set{1..6}/model/{bayeshammer,bless,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/bqsr/logs/bqsr
mkdir -p set{1..6}/{fastq,downloaded,tmp}
mkdir -p set{1..6}/downloaded/{intervals,metrics,original,split}
mkdir -p set{1..6}/fastq/{fastqc,split}
mkdir -p tmp

printf "\nDone\n"
