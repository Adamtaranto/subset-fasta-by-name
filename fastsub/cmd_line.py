#!/usr/bin/env python
#python 2.7.5 requires biopython
#fastsub
#Version 1.0.0 Adam Taranto, Sep 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

#################################################################################################
# Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta. #
#################################################################################################

import argparse
import fastsub as fs


def main():
	# Get cmd line args
	args = getArgs()
	# Check for output directories
	outdir = fs.getOut(args.outdir)
	# Read in names or ortho groups file
	readNames = fs.getTargetNames(args.nameFile, args)
	# Create dictionary of all sequences
	SeqMaster = fs.readFasta(args.inFasta)
	# Write to cluster files if OrthoMCL mode.
	if args.OrthoMCLMode:
		fs.orthomode(readNames, SeqMaster,outdir, args.nameFile)
	# Write each seq from namefile to individual fasta.
	# If no namefile, write all input seqs to own file.
	elif args.splitMode:
		fs.splitmode(SeqMaster,outdir)
	# Write selected seqs to single fasta
	# If no namefile write all seqs back out to single file. Useful for format cleanup.
	else:
		fs.filtermode(readNames, SeqMaster,outdir, args.outName, args.nameFile)

def getArgs():
	###Argument handling.
	parser = argparse.ArgumentParser(
						description='Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta. Or splits multi-fasta into single files. ',
						prog='fastsub')
	parser.add_argument("-i", "--inFasta",
						required=True,
						type=str,
						default= None,
						help="Multi fasta to extract subset from")
	parser.add_argument("-n", "--nameFile",
						type=str,
						default= None,
						help="Comma delimited file with target seq names in column one. If none given, all sequences will be returned.")
	parser.add_argument("-o", "--outName",
						type=str,
						default= "filtered_output.fa", 
						help="File for filtered sequence file to be written to.")
	parser.add_argument("-d", "--outdir",
						type=str,
						default= None, 
						help="Directory for new sequence files to be written to.")
	parser.add_argument('--splitMode',
						action='store_true',
						default=False,
						help='If set split each sequence into new fasta file.')
	parser.add_argument('--OrthoMCLMode',
						action='store_true',
						default=False,
						help='If set treat nameFile as "groups" output from OrthoMCL clustering. Write member sequences to Cluster output files as multifasta.')
	args = parser.parse_args()
	return args
