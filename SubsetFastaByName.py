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
import argparse;
from Bio import SeqIO
from Bio.Seq import Seq

def main(inFasta=None, nameFile=None, outName='filtered_seqs.fa'):

	if inFasta is None:
		sys.exit('No input fasta provided')

	if nameFile is None:
		sys.exit('No name list provided')

	#Create empty dictionary
	SeqMaster={}

	#populate dictionary with master set of fasta records
	for seq_record in SeqIO.parse(inFasta, "fasta"):
		SeqMaster[seq_record.id]=str(seq_record.seq)

	f = open(nameFile, 'rt')
	reader = csv.reader(f,dialect='excel')

	#Open output fasta
	fasta_file=open(outName,'w')

	#Open log file for names not found in master set
	error_list=open(str('NotFound_'+nameFile),'w')

	#Write records for seqs in name file to new fasta
	for row in reader:
		name=row[0]
		try:
			SeqMaster[name]
		except:
			print 'bad: ' + name
			error_list.write(name+"\n")
		else:
			fasta_name= ">%s" % (name)
			seq_denovo= "%s" %(SeqMaster[name])
			fasta_file.write(fasta_name+"\n")
			fasta_file.write(seq_denovo+"\n")

	fasta_file.close()
	error_list.close()

if __name__== '__main__':
	###Argument handling.
	arg_parser = argparse.ArgumentParser(description='Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta');
	arg_parser.add_argument("-i","--inFasta", help="Multi fasta to extract subset from");
	arg_parser.add_argument("-n","--nameFile", help="Comma delimited file with target seq names in column one");
	arg_parser.add_argument("-o","--outName", default='filtered_seqs.fa', help="Directory for new sequence file to be written to.");
	args = arg_parser.parse_args();

	###Variable Definitions
	inFasta=args.inFasta;
	nameFile=args.nameFile;
	outName=args.outName;

	main(inFasta, nameFile, outName);