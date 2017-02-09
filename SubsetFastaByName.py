#!/usr/bin/env python
#python 2.7.5 requires biopython
#SubsetFastaByName.py
#Version 1. Adam Taranto, June 2015
#Contact, Adam Taranto, adam.taranto@anu.edu.au

#################################################################################################
# Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta. #
#################################################################################################

import csv
import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def tempPathCheck(args):
    absOutDir = os.path.abspath(args.outDir)
    if not os.path.isdir(absOutDir):
        os.makedirs(absOutDir)

def chunkstring(string, length=80):
	return (string[0+i:length+i] for i in range(0, len(string), length))

def splitmode(SeqMaster, args):
	#Write records for seqs in name file to new fasta
	for key in SeqMaster.keys():
		if args.outDir:
			outName	= os.path.join(args.outDir, key + ".fa")
		else:
			outName	= key + ".fa"
		fasta_name	= ">%s" % (key)
		sequence	= "%s" %(SeqMaster[key])
		#Open new file
		fasta_file	= open(outName,'w')
		#Write seq name line
		fasta_file.write(fasta_name + "\n")
		#Write sequence with line wrapping
		for line in chunkstring(sequence):
			fasta_file.write(line + "\n")
		fasta_file.close()

def filtermode(reader, SeqMaster, outName, args):
	#Open output fasta
	if args.outDir:
			outName	= os.path.join(args.outDir, outName)
	fasta_file=open(outName,'w')
	#Open log file for names not found in master set
	error_list=open(str('NotFound_' + args.nameFile),'w')

	#Write records for seqs in name file to new fasta
	for row in reader:
		name=row[0]
		try:
			SeqMaster[name]
		except:
			print('bad: ' + name)
			error_list.write(name+"\n")
		else:
			fasta_name	= ">%s" % (name)
			sequence 	= "%s" %(SeqMaster[name])
			fasta_file.write(fasta_name+"\n")
			for line in chunkstring(sequence):
				fasta_file.write(line+"\n")

	fasta_file.close()
	error_list.close()

def main(args):

	if args.inFasta is None:
		sys.exit('No input fasta provided')	

	if args.outDir:
		tempPathCheck(args)
	
	if args.splitMode: #If running in splitmode
		outName = None
	else:
		if args.outName is None:
			outName = "filtered_output.fa"
		else:
			outName = args.outName
		if args.nameFile is None:
			sys.exit('No name list provided')
		else: #Read in name list
			f = open(args.nameFile, 'rt')
			reader = csv.reader(f,dialect='excel')

	#Create empty dictionary
	SeqMaster={}

	#Populate dictionary with master set of fasta records
	for seq_record in SeqIO.parse(args.inFasta, "fasta"):
		SeqMaster[seq_record.id]=str(seq_record.seq)

	if args.splitMode:
		splitmode(SeqMaster, args)
	else:
		filtermode(reader, SeqMaster, outName, args)

if __name__== '__main__':
	###Argument handling.
	parser = argparse.ArgumentParser(
		description='Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta. Or splits multi-fasta into single files. ',
		prog='SubsetFastaByName')
	parser.add_argument("-i", "--inFasta",
		type=str,
		default= None,
		help="Multi fasta to extract subset from")
	parser.add_argument("-n", "--nameFile",
		type=str,
		default= None,
		help="Comma delimited file with target seq names in column one")
	parser.add_argument("-o", "--outName",
		type=str,
		default= None, 
		help="File for filtered sequence file to be written to.")
	parser.add_argument("-d", "--outDir",
		type=str,
		default= None, 
		help="Directory for new sequence files to be written to.")
	parser.add_argument('--splitMode',
						action='store_true',
						default=False,
						help='If set split each sequence into new fasta file.')

	args = parser.parse_args()

	main(args);