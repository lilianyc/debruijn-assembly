#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:14:42 2019

@author: lyang_crosson
"""

import argparse
from pathlib import Path

import networkx as nx

# =============================================================================
# part a
# =============================================================================

def read_fastq(fastq_file):
    """Reads a fastq file and returns a sequence generator.
    """
    with open(fastq_file, "r") as filin:
        for line_number, line in enumerate(filin):
            if line_number % 4 == 1:
                yield line.strip()
    

def cut_kmer(sequence, kmer_size):
    """Cuts and returns k-mer iterator.
    """
    sequence_length = len(sequence)
    offset = 0
    while offset + kmer_size <= sequence_length:
        yield sequence[offset:(offset + kmer_size)]
        offset += 1


def build_kmer_dict(fastq_file, kmer_size):
    """Returns a dict of kmer counts.
    """
    kmer_count = {}
    sequences = read_fastq(fastq_file)
    for sequence in sequences:
        kmers = cut_kmer(sequence, kmer_size)
        for kmer in kmers:
            if not kmer in kmer_count:
                kmer_count[kmer] = 1
            else:
                kmer_count[kmer] += 1
    return kmer_count

# =============================================================================
# part b
# =============================================================================

def build_seq(kmer_count):
    graph = nx.DiGraph()
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

    data_dir = Path(__name__).resolve().parent.joinpath("data/")
    sequences = read_fastq(data_dir.joinpath("eva71_two_reads.fq"))
    seq1 = next(sequences)
    kmers = cut_kmer(seq1, 2)

    kmer_count = build_kmer_dict(data_dir.joinpath("eva71_two_reads.fq"), 2)
    kmer_count
