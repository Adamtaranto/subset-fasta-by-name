# fastsub

A collection of functions for subsetting and reformatting sequences from fasta files.

# Table of contents

* [Installing Fastsub](#installing-fastsub)
* [Example usage](#example-usage)
* [Standard options](#standard-options)
* [License](#license)

# Installing Fastsub

Clone and install from this repository:
```bash
git clone https://github.com/Adamtaranto/subset-fasta-by-name.git && cd subset-fasta-by-name && pip install -e .
```

# Example usage

Extract a subset of sequences to a new multifasta file

```bash
fastsub -i input.fa -o subset.fa --nameFile list.txt
```

Extract a subset of sequences to individual fasta files

```bash
fastsub -i input.fa -o subset.fa --nameFile list.txt --splitMode
```

Write all sequences to new multifasta [Use to clean up formatting in multifasta.]

```bash
fastsub -i input.fa -o clean_input.fa
```

Split all sequences to individual fasta files

```bash
fastsub -i input.fa -d split_files --splitMode
```

Write sequences to multifasta by OrthoMCL Cluster membership

```bash
fastsub -i input.fa -d ortho_clusters --OrthoMCLmode
```

# Standard options

```
fastsub [-h] -i INFASTA [-n NAMEFILE] [-o OUTNAME] [-d OUTDIR]
               [--splitMode] [--OrthoMCLMode]

Takes a multi-fasta file and a list of sequence names, prints named sequences
to a new fasta. Or splits multi-fasta into single files.

optional arguments:
  -h, --help            show this help message and exit
  -i, --inFasta
                        Multi fasta to extract subset from
  -n, --nameFile
                        Comma delimited file with target seq names in column
                        one. If none given, all sequences will be returned.
  -o, --outName
                        File for filtered sequence file to be written to.
  -d, --outdir
                        Directory for new sequence files to be written to.
  --splitMode           If set split each sequence into new fasta file.
  --OrthoMCLMode        If set treat nameFile as "groups" output from OrthoMCL
                        clustering. Write member sequences to Cluster output
                        files as multifasta.
```

# License

Software provided under MIT license.
