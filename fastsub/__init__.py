#!/usr/bin/env python

import csv
import sys
import os
from Bio import SeqIO

def getOut(outdir):
	if outdir:
	    absOutDir = os.path.abspath(outdir)
	    if not os.path.isdir(absOutDir):
	        os.makedirs(absOutDir)
	    return absOutDir
	else: 
		return os.getcwd()

def readFasta(infile):
	#Populate empty dictionary with master set of fasta records
	SeqMaster= dict()
	for seq_record in SeqIO.parse(infile, "fasta"):
		SeqMaster[seq_record.id]=(str(seq_record.description),str(seq_record.seq))
	return SeqMaster

def getTargetNames(nameFile, args):
	'''If namefile or othoMCL groups file provided, return generator object'''
	if args.OrthoMCLMode and nameFile is None:
		print("Sequence names required if running OrthoMCLMode. Provide groups.txt file.")
		sys.exit(1)
	if nameFile is None:
		readNames = None
		print("No namefile provided. Returning all sequences.")
	else: #Read in name list
		f = open(nameFile, 'rt')
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

def orthomode(reader, SeqMaster, outdir, nameFile):
	# Load reader into seq by cluster dict
	groupNames = ortho2dict(reader)
	# Open log file for names not found in master set
	error_list = open(str('NotFound_' + nameFile),'w')
	# Process each cluster
	for group in groupNames.keys():
		# Set cluster output path
		if outdir:
			outName	= os.path.join(outdir, group + ".fa")
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
	
def splitmode(SeqMaster, outdir):
	#Write records for seqs in name file to new fasta
	for key in SeqMaster.keys():
		if outdir:
			outName	= os.path.join(outdir, key + ".fa")
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

def filtermode(reader, SeqMaster, outdir,outName, nameFile):
	#Open output fasta
	if outdir:
			outName = os.path.join(outdir, outName)
	fasta_file=open(outName,'w')
	#Write records for seqs in name file to single new fasta file
	if reader:
		#Open log file for names not found in master set
		error_list=open(str('NotFound_' + nameFile),'w')
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