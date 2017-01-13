# Thesis Pipeline

A system to run experimental conditions for my master's thesis. 

Rather than a pipeline, this is more of a non-scalable pipeline management system that has been specifically designed to implement GATK based genomics pipelines on OHSU's ExaCloud cluster. It's design goals are:

* To maintain a consistent API for calling disparate bioinformatics tools (e.g. GATK, Picard, Samtools)
* To enforce requirements for the experimental design of the thesis
* To enable complete logging of processes

For logging, the system requires that each module is called within the context of a TMUX session. 

The system does not:

* Guarantee each module will succeed in its computation
* Manage error codes whatsoever (big caveat!)
* Expose all parameters of the algorithms being run

Regarding the last point, unexposed parameters are hard-coded into the scripts and can be updated to a user's needs. For example, the reference genome used is hard coded into all scripts (gross, I know).

