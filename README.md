# Genome assembly using De Bruijn graphs

*debruijn.py* is a script for performing genome assembly relying on De Bruijn graphs.

## Requirements

- Python >= 3.6
- NetworkX

## Quickstart

The script takes a [fastq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/) file
and outputs a file with the contigs built.

### Options

- `-h, --help`:
Provide a help page

- `-i, --input`:
Name of the fastq file

- `-o, --contig`:
Give the name of the output file

- `-k, --kmer`:
Size of the kmers (default: 21)

### Examples

Asumming the requirements are satisfied, using the following command in *debruijn/*
```shell
$ python3 debruijn.py -i file.fastq -o sequences.contigs
```
will create `sequences.contigs`, a file containing the contigs assembled and their size.
