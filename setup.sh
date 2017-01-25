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
# Step    Procedure           State               Input       Output      
# 1       bamtofastq.sh       BAMTOFASTQ:1        bam         bam
# 2       bamtofastq.sh       BAMTOFASTQ:2        bam         fastq
# 3       errormodel.sh       ERRORMODEL:1        fastq       fastq
# 4       errormodel.sh       BAYESHAMMER:1       fastq       fastq
# 5       errormodel.sh       BLOOCOO:1           fastq       fastq
# 6       errormodel.sh       QUORUM:1            fastq       fastq
# 7       bwa.sh              BWA:1               fastq       fai
# 8       bwa.sh              BWA:2               fastq       bam
# 9       mergealignment.sh   MERGEALIGNMENT:1    bam         bam
# 10      mergealignment.sh   MERGEALIGNMENT:2    bam         bam
# 11      markduplicates.sh   MARKDUPLICATES:1    bam         bam
# 12      markduplicates.sh   MARKDUPLICATES:2    bam         bam
# 13      markduplicates.sh   MARKDUPLICATES:3    bam         bam
# 14      markduplicates.sh   MARKDUPLICATES:4    bam         bam
# 15      bqsr.sh             NOBQSR:1            bam         bam
# 16      bqsr.sh             BQSR:1              bam         table
# 17      bqsr.sh             BQSR:2              bam         table
# 18      bqsr.sh             BQSR:3              bam         pdf
# 19      bqsr.sh             BQSR:4              bam         bam
# 20      bqsr.sh             BQSR:5              bam         bai
# 21      mutect2.sh          MUTECT2:1           bam         txt
# 22      mutect2.sh          MUTECT2:2           bam         vcf
# 23      mutect2.sh          MUTECT2:3           vcf         vcf

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
# Scaffolding
#

printf "\nScaffolding Directories...\n"

# Scaffold Pipeline Directory & Change Directory
cd $root_dir && mkdir -p $pipeline_dir && cd $pipeline_dir

# Pipeline Scaffolding
mkdir -p set{1..6}/model/{bayeshammer,bless,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{logs,post-align,pre-align,markdup,recal,merged,modeled}
mkdir -p set{1..6}/model/{bayeshammer,bless,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/modeled/{pairs,indexes}
mkdir -p set{1..6}/model/{bayeshammer,bless,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{post-align,pre-align}/fastq
mkdir -p set{1..6}/model/{bayeshammer,bless,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/{bqsr,nobqsr}/logs/{contest,mutect2}
mkdir -p set{1..6}/model/{bayeshammer,bless,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/bqsr/logs/bqsr
mkdir -p set{1..6}/{fastq,downloaded,tmp}
mkdir -p set{1..6}/downloaded/{intervals,metrics,original,split}
mkdir -p set{1..6}/fastq/{fastqc,split}

#
# Exports
#

while true; do
    echo
    IFS= read -p "Write pipeline exports to ~/.bash_profile & ~/.bashrc? " yn
    case $yn in
        [Yy]* ) (

                printf "Exporting Directory Locations...\n"

                # Header for Bash Profile
                echo -e '\n#\n# Pipeline Management\n#' >> ~/.bash_profile
                echo -e '\n#\n# Pipeline Management\n#' >> ~/.bashrc

                # Write Pipeline Home Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_HOME=$pipeline_dir" >> ~/.bash_profile
                echo -e "\nexport PIPELINE_HOME=$pipeline_dir" >> ~/.bashrc
                # Set Pipeline Home Directory
                export PIPELINE_HOME=$pipeline_dir

                # Write Pipeline Output Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_OUT=$output_dir" >> ~/.bash_profile
                echo -e "\nexport PIPELINE_OUT=$output_dir" >> ~/.bashrc
                # Set Pipeline Output Directory
                export PIPELINE_OUT=$output_dir
                
                # Write Pipeline Reference Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_REF=$reference_dir" >> ~/.bash_profile
                echo -e "\nexport PIPELINE_REF=$reference_dir" >> ~/.bashrc
                # Set Pipeline Reference Directory
                export PIPELINE_REF=$reference_dir
                
                # Write Pipeline Logging Directory to ~/.bash_profile
                echo -e "\nexport PIPELINE_LOG=$logging_dir" >> ~/.bash_profile
                echo -e "\nexport PIPELINE_LOG=$logging_dir" >> ~/.bashrc
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
    echo
    IFS= read -p "Enforce TMUX session logging via ~/.bash_profile & ~/.bashrc? " yn
    case $yn in
        [Yy]* )  
            # Inject Logging Automation Script to ~/.bash_profile
            echo -e '\n# Pipeline Logging Management\nif [[ $TERM = "screen" ]] && [[ $(ps -p $PPID -o comm=) = "tmux" ]] && [[ $PWD/ = $PIPELINE_HOME/* ]]; then\n\t# Prompt User / Read Input\n\tcd $PIPELINE_HOME/procedures\n\tIFS= read -p "Enter job logging name: " name\n\tmkdir $PIPELINE_HOME/logs 2> /dev/null\n\tlogname="tmux.$name.$(date "+%d%m%Y%H%M%S").log"\n\tscript -f $PIPELINE_LOG/${logname}\n\texit\nfi' >> ~/.bash_profile
            echo -e '\n# Pipeline Logging Management\nif [[ $TERM = "screen" ]] && [[ $(ps -p $PPID -o comm=) = "tmux" ]] && [[ $PWD/ = $PIPELINE_HOME/* ]]; then\n\t# Prompt User / Read Input\n\tcd $PIPELINE_HOME/procedures\n\tIFS= read -p "Enter job logging name: " name\n\tmkdir $PIPELINE_HOME/logs 2> /dev/null\n\tlogname="tmux.$name.$(date "+%d%m%Y%H%M%S").log"\n\tscript -f $PIPELINE_LOG/${logname}\n\texit\nfi' >> ~/.bashrc
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Please answer yes [y] or no [n].";;
    esac
done

#
# State Management
#

while true; do
    echo
    IFS= read -p "Setup pipeline state manager? " yn
    case $yn in
        [Yy]* )  
            # Inject State Manager Link to ~/.bash_profile
            echo -e '\n# Pipeline State Manger Link\nif [ -f $PIPELINE_HOME/core/pipeline.core ]; then\n\tsource $PIPELINE_HOME/core/pipeline.core\nfi' >> ~/.bash_profile
            echo -e '\n# Pipeline State Manger Link\nif [ -f $PIPELINE_HOME/core/pipeline.core ]; then\n\tsource $PIPELINE_HOME/core/pipeline.core\nfi' >> ~/.bashrc
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Please answer yes [y] or no [n].";;
    esac
done
