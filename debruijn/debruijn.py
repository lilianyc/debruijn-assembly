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
# Identify unique kmers.
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
# Build the De Bruijn graph.
# =============================================================================

def build_graph(kmer_count):
    """Build a networkx.DiGraph from a dict of kmer counts.
    """
    graph = nx.DiGraph()
    for kmer in kmer_count.keys():
        node_1 = kmer[:-1]
        node_2 = kmer[1:]
        graph.add_edge(node_1, node_2, weight=kmer_count[kmer])
    return graph

# =============================================================================
# Graph analysis.
# =============================================================================

def get_starting_nodes(graph):
    """Returns a list of starting nodes.
    """
    starting_node_list = []
    for node in graph.nodes():
        if not list(graph.predecessors(node)):
            starting_node_list.append(node)
    return starting_node_list


def get_sink_nodes(graph):
    """Returns a list of sink nodes.
    """
    sink_node_list = []
    for node in graph.nodes():
        if not list(graph.successors(node)):
            sink_node_list.append(node)
    return sink_node_list


def get_contigs(graph, starting_nodes, sink_nodes):
    pass


def save_contigs():
    pass


def std():
    pass


def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Read single-end fastq file and returns.")

    # arguments
    parser.add_argument("-i", "--input", required=True,
                        help=("name of the fastq file."))
    parser.add_argument("-k", "--kmer", type=int, default=21,
                        help="length of kmers (default: 21).")
    parser.add_argument("-o", "--config",
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

    kmer_count = build_kmer_dict(data_dir.joinpath("eva71_two_reads.fq"), 3)
    print(kmer_count)

    graph = build_graph(kmer_count)
