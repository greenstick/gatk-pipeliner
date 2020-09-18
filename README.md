# GATK Pipeliner

This is a system that was developed to run experimental conditions for my master's thesis. It is based on the GATK pipeline and designed to allow for the inclusion of arbitrary error modeling modules into GATK. The ultimate aim of this research (completedd circa 2017) was to benchmark the GATK pipeline with andd without the BQSR module for the specific use case of detecting minimal residual disease (MRD) burden in acute myeloid leukemia (AML) – see the [greenstick/thesis-analyses](https://github.com/greenstick/thesis-analyses)) for one approach to doing this kind of benchmarking. If you're interested in using this as a tool to generate bioinformatics pipelines, I would suggest that a much more robust solution is [snakemake](https://snakemake.github.io/). 

## Overview & Limitations

Rather than a pipeline, this is more of a non-scalable pipeline management system that has been specifically designed to implement GATK based genomics pipelines on OHSU's ExaCloud cluster. It's design goals are:

* To maintain a consistent API for calling disparate bioinformatics tools (e.g. GATK, Picard, Samtools)
* To enforce requirements for the experimental design of the thesis
* To enable complete logging and auditing of pipeline processes

Given the scoping constraints within the shell environment, the system relies almost entirely on namespacing for scope management. Further, if logging procedures (e.g. variant calling with Mutect2) are required, each module must be called within the context of a TMUX session that has been initialized within the scope of the systems' `$PIPELINE_HOME` directory. 

Note, the system does not:

* Guarantee each module will succeed in its computation
* Expose system-level variables for supplemental data files used in processing (i.e. reference genomes, call sets, etc.)
* Expose all parameters of the algorithms being run

Regarding the last point, unexposed parameters are hard-coded into the scripts and can be updated to as needed. For example, the HG19 reference genome used is hard coded into all scripts (gross, I know). 

Finally, there is a basic state management system that, while imperfect, helps to reduce the opportunity cost of script failures on the cluster – HPC systems with distributed file systems are sometimes unstable, after all. The state management system makes heavy use of error codes from the software modules. This is a big caveat – if the software module does not pass relevant error codes, the state management system may indicate a procedure has run when it actually failed. If you intend to include novel modules into your pipelines managed by this tool, be mindful of the error codes they return. 

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

Note that the `dataset`, `tmp`, and `logs` directories are scaffolded when you call `source setup.py`. The reference directory is not and must be managed entirely by the user. For development, that mean setting up an `alias` in my `~/.bash_profile` & `~/.bashrc` to point to the directory. From there, scripts within the directory directly refer to the reference genomes, any VCF callsets, etc. (e.g. the version of the hg19 reference genome as a FASTA or a 1000 genomes VCF). Finally, the data that will be processed is placed in the `downloaded` child directory of the `dataset` directory. A future update will likely provide a variable defined during the setup procedure that will store a reference to a user-defined home data directory. The same *may* go for any reference genomes.

```
pipeline
├── setup.sh
├── submit.sh
├── README.md
├── LICENSE
├── core
│   ├── pipeline.config
│   ├── pipeline.core
│   ├── pipeline.state
│   └── states.manifest
├── procedures
│   ├── bamtofastq.sh
│   ├── bqsr.sh
│   ├── bwa.sh
│   ├── errormodel.sh
│   ├── markduplicates.sh
│   ├── mergealignment.sh
│   ├── models
│   ├── mutect2.sh
|   └── utils
|       ├── emit-error.py
|       ├── fasta-qual-to-fastq.py
|       ├── fastq-to-fasta-qual.py
|       ├── paired-end-to-single-ends.py
|       └── single-ends-to-paired-end.py
├── logs
├── tmp
├── reference
│   ├── extra
│   └── hg19
├── dataset
|   └── ...
└── dataset
    ├── downloaded
    │   ├── intervals
    │   ├── metrics
    │   ├── original
    │   └── split
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

   
