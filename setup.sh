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
# References
#

printf "\nSetting References...\n"

# Set Root Directory Reference
root_dir=/home/exacloud/lustre1

# Set Pipeline Directory Reference
pipeline_dir=MRD_aml/data/dream-synthetic/pipeline

#
# Exports
#

while true; do
    read -p "Write pipeline home directory export to ~/.bash_profile (default is yes)?" yn
    case $yn in
        [Yy]* ) (printf "\nExporting Locations...\n"
            # Write Pipeline Home Directory to ~/.bash_profile
            echo "export PIPELINE_HOME=$root_dir/$pipeline_dir" >> ~/.bash_profile
            # Set Pipeline Home Directory
            export PIPELINE_HOME=$root_dir/$pipeline_dir
            ); break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
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
mkdir -p set{1..6}/downloaded/{intervals,metrics,original}
mkdir -p set{1..6}/fastq/{fastqc,split}
mkdir -p tmp

printf "\nDone\n"
