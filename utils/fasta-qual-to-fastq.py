#! /home/users/cordier/.linuxbrew/bin/python3

if __name__ == "__main__":


	# Imports
	import sys, getopt
	# Library Import
	from Bio import SeqIO
	from Bio.SeqIO.QualityIO import PairedFastaQualIterator

	# Define References
	fasta = ""
	qual = ""
	prefix = ""

	# 
	# Parse Arguments
	# 

	try:
		opts, args = getopt.getopt(sys.argv, "hf:q:o", ["--inputfasta", "--inputqual", "--outputprefix", "--help"])
	except getopt.GetoptError:
		print("fasta-qual to fastq.py -f <input fasta> -q <input qual> -o <output prefix>")
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print("Merge Fasta and Qual files into a FastQ file.\nExample usage:\nfasta-qual-to-fastq.py -f <input fasta> -q <input qual> -o <output prefix>")
		elif opt in ("-f", "--inputfasta"):
			fasta = arg
		elif opt in ("-q", "--inputqual"):
			qual = arg
		elif opt in ("-o", "--outputprefix"):
			prefix = arg

	# Assertions for Required Input
	assert (len(fasta) > 0), "No Fasta input provided!"
	assert (len(qual) > 0), "No Qual input provided!"

	# 
	# Conversion
	#

	#  If No Prefix, Use Same as FASTQ
	if len(prefix) < 1: 
		prefix = fasta.split(".")[0:-1]

	# Merge Fasta & Qual into FastQ
	records = PairedFastaQualIterator(open(fasta), open(qual))
	SeqIO.write(records, prefix + ".fastq", "fastq")


else:

	pass
