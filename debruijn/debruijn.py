#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:14:42 2019

@author: lyang_crosson
"""

import argparse
from pathlib import Path


"""
-i fichier fastq single end
-k taille des kmer (optionnel - default 21)
-o fichier config
"""

def read_fastq(fastq_file):
    """Reads a fastq file and returns a sequence generator.
    """
    with open(fastq_file, "r") as filin:
        for line_number, line in enumerate(filin):
            print(line_number, line.strip())
    


def cut_kmer():
    """Cuts and returns k-mer iterator.
    """
    pass


def build_kmer_dict():
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Read single-end fastq file and returns.")

    # arguments
    parser.add_argument("-i", "--input", required=True,
                        help=("name of the fastq file."))
    parser.add_argument("-k", "--kmer", type=int, const=21, nargs="?",
                        help="length of kmers.")
    parser.add_argument("-o", "--config", type=str,
                        help="name of config file.")

    # get all arguments
    options = parser.parse_args()

    return options
    


if __name__ == "__main__":
    options = main()
    print(options)
#    data_dir = Path(__name__).resolve().parent.parent.joinpath("data/")
#    read_fastq(data_dir.joinpath("eva71_two_reads.fq"))
