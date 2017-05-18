# subset-fasta-by-name

A collection of functions for subsetting and reformatting sequences from fasta files.

## Uses  
  - Extract a subset of sequences to a new multifasta file (--nameFile list.txt)
  - Extract a subset of sequences to individual fasta files (--nameFile list.txt, --splitMode)
  - Write all sequences to new multifasta (--nameFile *None*) [Use to clean up formatting in multifasta.]
  - Split all sequences to individual fasta files (--splitMode)
  - Write sequences to multifasta by OrthoMCL Cluster membership (--OrthoMCLmode)

Settings:

**-h, --help**:  
  - Show this help message and exit
**-v, --version**:  
  - Show program's version number and exit
**-i INFASTA, --inFasta INFASTA**:  
  - Multi fasta to extract subset from
**-n NAMEFILE, --nameFile NAMEFILE**:  
  - Comma delimited file with target seq names in column one. If none given, all sequences will be returned.
**-o OUTNAME, --outName OUTNAME**:  
  - File for filtered sequence file to be written to.
**-d OUTDIR, --outDir OUTDIR**:  
  - Directory for new sequence files to be written to.
**--splitMode**:  
  - If set split each sequence into new fasta file.
**--OrthoMCLMode**:  
  - If set treat nameFile as "groups" output from OrthoMCL clustering. Write member sequences to Cluster output files as multifasta.

