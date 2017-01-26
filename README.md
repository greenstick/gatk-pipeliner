# Thesis Pipeline

A system to run experimental conditions for my master's thesis. 

Rather than a pipeline, this is more of a non-scalable pipeline management system that has been specifically designed to implement GATK based genomics pipelines on OHSU's ExaCloud cluster. It's design goals are:

* To maintain a consistent API for calling disparate bioinformatics tools (e.g. GATK, Picard, Samtools)
* To enforce requirements for the experimental design of the thesis
* To enable complete logging of processes

For logging, the system requires that each module is called within the context of a TMUX session. 

The system does not:

* Guarantee each module will succeed in its computation
* Expose all parameters of the algorithms being run

Regarding the last point, unexposed parameters are hard-coded into the scripts and can be updated to a user's needs. For example, the reference genome used is hard coded into all scripts (gross, I know). 

Finally, there is a basic state management system that, while imperfect, helps to reduce the opportunity cost of script failures on the cluster. The state management system makes heavy use of error codes from the software modules. This is a big caveat – if the software module does not pass relevant error codes, the state management system may indicate a procedure has run when it actually failed. 

To install:

    git clone https://github.com/greenstick/thesis-pipeline.git
    source setup.sh

To run a procedure, the general syntax is:

	  source [some procedure] \
	  -f=$prefix \
	  -s=$dataset \
	  -c=$condition \
	  -x=$experiment \
	  -p=$parameters \
	  -q=$bqsr \
	  -m=$memory_per_thread
	  -n=$threads
    # Note not all parameters are accepted by all procedure scripts
    
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
    
