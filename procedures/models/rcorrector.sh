# Usage: perl ./run_rcorrector.pl [OPTIONS]
# OPTIONS:
# Required parameters:
# 	-s seq_files: comma separated files for single-end data sets
# 	-1 seq_files_left: comma separated files for the first mate in the paried-end data sets
# 	-2 seq_files_right: comma separated files for the second mate in the paired-end data sets
# 	-i seq_files_interleaved: comma sperated files for interleaved paired-end data sets
# Other parameters:
# 	-k kmer_length (<=32, default: 23)
# 	-od output_file_directory (default: ./)
# 	-t number_of_threads (default: 1)
# 	-maxcorK INT: the maximum number of correction within k-bp window (default: 4)
# 	-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold, lower for more divergent genome (default: 0.95)
# 	-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)
# 	-stdout: output the corrected reads to stdout (default: not used)
# 	-verbose: output some correction information to stdout (default: not used)
# 	-stage INT: start from which stage (default: 0)
# 		0-start from begining(storing kmers in bloom filter);
# 		1-start from count kmers showed up in bloom filter;
# 		2-start from dumping kmer counts into a jf_dump file;
# 		3-start from error correction.