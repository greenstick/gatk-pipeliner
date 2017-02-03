# Thesis Pipeline

A system to run experimental conditions for my master's thesis. 

## Overview & Limitations

Rather than a pipeline, this is more of a non-scalable pipeline management system that has been specifically designed to implement GATK based genomics pipelines on OHSU's ExaCloud cluster. It's design goals are:

* To maintain a consistent API for calling disparate bioinformatics tools (e.g. GATK, Picard, Samtools)
* To enforce requirements for the experimental design of the thesis
* To enable complete logging of processes

Given the scoping constraints within the shell environment, the system relies almost entirely on namespacing for scope management. Further, if logging procedures (e.g. variant calling with Mutect2) is required, each module must be called within the context of a TMUX session that has been initialized within the scope of the systems' `$PIPELINE_HOME` directory. 

As of yet, the system does not:

* Guarantee each module will succeed in its computation
* Expose all parameters of the algorithms being run

Regarding the last point, unexposed parameters are hard-coded into the scripts and can be updated to a user's needs. For example, the reference genome used is hard coded into all scripts (gross, I know). 

Finally, there is a basic state management system that, while imperfect, helps to reduce the opportunity cost of script failures on the cluster. The state management system makes heavy use of error codes from the software modules. This is a big caveat – if the software module does not pass relevant error codes, the state management system may indicate a procedure has run when it actually failed. 

## Installation

To install:

    git clone https://github.com/greenstick/thesis-pipeline.git
    source setup.sh

## How to Run a Procedure

To run a procedure, the general syntax is:

    source [some procedure] \
    -f=$prefix \
    -s=$dataset \
    -c=$condition \
    -x=$experiment \
    -p=$parameters \
    -q=$bqsr \
    -m=$memory_per_thread \
    -n=$threads
    # Note: not all parameters are accepted by all procedure scripts
    
Where:

* `$prefix` is the file prefix of a given dataset; this can be an empty string or something like 'data-from-arnold'
* `$dataset` is the dataset name, for example 'dataset1' or 'dna-seq-experiment-1'
* `$condition` is whether the condition is a 'control' or 'experimental'
* `$experiment` is the main experiment being run. This could be anything, but you'll likely have to write your own procedure
* `$parameters` is whether the parameters are default for the experimental condition or custom
* `$bqsr` is whether BQSR is being used; again this may not be relevant

A more concrete example:

    source errormodel.sh \
    -f="synthetic.challenge" \
    -s="set1" \
    -c="tumor" \
    -x="bayeshammer" \
    -p="default" \
    -m="12G" \
    -n=8
    # Note that for the -s argument, this should be both the name of a data folder in the pipeline directory AND
    # the part of the BAM / FASTQ files to be processed. For Example: synthetic.challenge.set1.tumor.bam
    
## Directory Structure

```
pipeline
├── datasetX
|   └── ...
└── datasetY
    ├── tmp
    ├── downloaded
    │   ├── intervals
    │   ├── metrics
    │   ├── original
    |   ├── split
    │   └── *data_to_process_goes_here*
    ├── fastq
    │   ├── fastqc
    │   └── split
    │       ├── unmerged
    │       └── unpaired
    ├── model
    │   └── errormodelX                                                                                                                                                                               
    │       └── param
    │           ├── custom
    │           │   ├── logs
    │           │   ├── markdup
    │           │   ├── merged
    │           │   ├── modeled
    │           │   │   ├── indexes
    │           │   │   └── pairs
    │           │   ├── post-align
    │           │   │   └── fastq
    │           │   ├── pre-align
    │           │   │   └── fastq
    │           │   └── recal
    │           │       ├── bqsr
    │           │       │   └── logs
    │           │       │       ├── bqsr
    │           │       │       ├── contest
    │           │       │       └── mutect2
    │           │       └── nobqsr
    │           │           └── logs
    │           │               ├── contest
    │           │               └── mutect2
    │           └── default
    │               ├── logs
    │               ├── markdup
    │               ├── merged
    │               ├── modeled
    │               │   ├── indexes
    │               │   └── pairs
    │               ├── post-align
    │               │   └── fastq
    │               ├── pre-align
    │               │   └── fastq
    │               └── recal
    │                   ├── bqsr
    │                   │   └── logs
    │                   │       ├── bqsr
    │                   │       ├── contest
    │                   │       └── mutect2
    │                   └── nobqsr
    │                       └── logs
    │                           ├── contest
    │                           └── mutect
    └── tmp 
```

   
