#! /home/users/cordier/.linuxbrew/bin/python3

if __name__ == "__main__":


	# Imports
	import sys, getopt
	# Library Import
	from Bio import SeqIO

	# Define References
	fastq = ""
	prefix = ""

	# 
	# Parse Arguments
	# 

	try:
		opts, args = getopt.getopt(sys.argv, "hi:o", ["--input", "--outputprefix", "--help"])
	except getopt.GetoptError:
		print("fastq-to-fasta-qual.py -i <input fastq> -o <output prefix>")
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print("Split a FASTQ into Seperate FASTA and QUAL files.\nExample usage:\nfastq-to-fasta-qual.py -i <input fastq> -o <output prefix>")
		elif opt in ("-i", "--input"):
			fastq = arg
		elif opt in ("-o", "--outputprefix"):
			prefix = arg

	print(fastq)

	# Assertions for Required Input
	assert (len(fastq) > 0), "No FASTQ input provided!"

	# 
	# Conversion
	#

	#  If No Prefix, Use Same as FastQ
	if len(prefix) < 1: 
		prefix = fastq.split(".")[0:-1]

	# Split FastQ into Fasta & Qual
	SeqIO.convert(fastq, "fastq", prefix + ".fasta", "fasta")
	print("Fasta file written to %s.fasta" % prefix)
	SeqIO.convert(fastq, "fastq", prefix + ".qual", "qual")
	print("Qual file written to %s.fasta" % prefix)
	

else:

	pass
