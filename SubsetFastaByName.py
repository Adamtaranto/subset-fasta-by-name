#!/usr/bin/env python
#python 2.7.5 requires biopython
#SubsetFastaByName.py
#Version 1.1.0 Adam Taranto, May 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

#################################################################################################
# Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta. #
#################################################################################################

import csv
import sys
import os
import argparse
from Bio import SeqIO
from _version import __version__

__version__ = '1.1.0'

def tempPathCheck(args):
	if args.outDir:
	    absOutDir = os.path.abspath(args.outDir)
	    if not os.path.isdir(absOutDir):
	        os.makedirs(absOutDir)

def readFasta(args):
	#Populate empty dictionary with master set of fasta records
	SeqMaster= dict()
	for seq_record in SeqIO.parse(args.inFasta, "fasta"):
		SeqMaster[seq_record.id]=(str(seq_record.description),str(seq_record.seq))
	return SeqMaster

def getTargetNames(args):
	'''If namefile or othoMCL groups file provided, return generator object'''
	if args.OrthoMCLMode and args.nameFile is None:
		print("Sequence names required if running OrthoMCLMode. Provide groups.txt file.")
		sys.exit(1)
	if args.nameFile is None:
		readNames = None
		print("No namefile provided. Returning all sequences.")
	else: #Read in name list
		f = open(args.nameFile, 'rt')
		readNames = csv.reader(f,dialect='excel')
	return readNames

def chunkstring(string, length=80):
	return (string[0+i:length+i] for i in range(0, len(string), length))

def ortho2dict(reader):
	groupDict = dict()
	for line in reader:
		line = line[0].split()
		clusterID = line[0].strip(':')
		if clusterID not in groupDict.keys():
			groupDict[clusterID] = list()
		for name in line[1:]:
			groupDict[clusterID].append(name)
	return groupDict

def orthomode(reader, SeqMaster, args):
	# Load reader into seq by cluster dict
	groupNames = ortho2dict(reader)
	# Open log file for names not found in master set
	error_list = open(str('NotFound_' + args.nameFile),'w')
	# Process each cluster
	for group in groupNames.keys():
		# Set cluster output path
		if args.outDir:
			outName	= os.path.join(args.outDir, group + ".fa")
		else:
			outName	= group + ".fa"
		# Open cluster output file
		fasta_file=open(outName,'w')
		# Retrive member sequences
		for name in groupNames[group]:
			try:
				SeqMaster[name]
			except:
				print('Sequence not found: ' + name)
				error_list.write(group + "\t" + name + "\n")
			else:
				fasta_name	= ">%s" % (SeqMaster[name][0])
				sequence 	= "%s" %(SeqMaster[name][1])
				fasta_file.write(fasta_name + "\n")
				for line in chunkstring(sequence):
					fasta_file.write(line + "\n")
		# Close Cluster output
		fasta_file.close
	# Close error log
	error_list.close()
	
def splitmode(SeqMaster, args):
	#Write records for seqs in name file to new fasta
	for key in SeqMaster.keys():
		if args.outDir:
			outName	= os.path.join(args.outDir, key + ".fa")
		else:
			outName	= key + ".fa"
		fasta_name	= ">%s" % (SeqMaster[key][0])
		sequence	= "%s" %(SeqMaster[key][1])
		#Open new file
		fasta_file	= open(outName,'w')
		#Write seq name line
		fasta_file.write(fasta_name + "\n")
		#Write sequence with line wrapping
		for line in chunkstring(sequence):
			fasta_file.write(line + "\n")
		fasta_file.close()

def filtermode(reader, SeqMaster, args):
	#Open output fasta
	if args.outDir:
			outName = os.path.join(args.outDir, args.outName)
	fasta_file=open(outName,'w')
	#Write records for seqs in name file to single new fasta file
	if reader:
		#Open log file for names not found in master set
		error_list=open(str('NotFound_' + args.nameFile),'w')
		for row in reader:
			name=row[0]
			try:
				SeqMaster[name]
			except:
				print('bad: ' + name)
				error_list.write(name+"\n")
			else:
				fasta_name	= ">%s" % (SeqMaster[name][0])
				sequence 	= "%s" %(SeqMaster[name][1])
				fasta_file.write(fasta_name + "\n")
				for line in chunkstring(sequence):
					fasta_file.write(line+"\n")
		error_list.close()
	else:
		# If no names provided read all sequences to new file
		for key in SeqMaster.keys():
			fasta_name	= ">%s" % (SeqMaster[key][0])
			sequence	= "%s" %(SeqMaster[key][1])
			#Write seq name line
			fasta_file.write(fasta_name + "\n")
			#Write sequence with line wrapping
			for line in chunkstring(sequence):
				fasta_file.write(line + "\n")
	fasta_file.close()

def main(args):
	# Check for output directories
	tempPathCheck(args)
	# Read in names or ortho groups file
	readNames = getTargetNames(args)
	# Create dictionary of all sequences
	SeqMaster = readFasta(args)
	# Write to cluster files if OrthoMCL mode.
	if args.OrthoMCLMode:
		orthomode(readNames, SeqMaster, args)
	# Write each seq from namefile to individual fasta.
	# If no namefile, write all input seqs to own file.
	elif args.splitMode:
		splitmode(SeqMaster, args)
	# Write selected seqs to single fasta
	# If no namefile write all seqs back out to single file. Useful for format cleanup.
	else:
		filtermode(readNames, SeqMaster, args)

if __name__== '__main__':
	###Argument handling.
	parser = argparse.ArgumentParser(
						description='Takes a multi-fasta file and a list of sequence names, prints named sequences to a new fasta. Or splits multi-fasta into single files. ',
						prog='SubsetFastaByName')
	parser.add_argument('-v', '--version', 
						action='version', 
						version='%(prog)s {version}'.format(version=__version__))
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
	parser.add_argument("-d", "--outDir",
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

	main(args);