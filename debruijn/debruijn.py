#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:14:42 2019

@author: lyang_crosson
"""

import argparse
import os
from pathlib import Path
import random
import statistics

import networkx as nx

# =============================================================================
# 1. Create the De Bruijn graph.
# =============================================================================

# =============================================================================
# a. Identify unique kmers.
# =============================================================================

def read_fastq(fastq_file):
    """Reads a fastq file and returns a generator from sequences."""
    with open(fastq_file, "r") as filin:
        for line_number, line in enumerate(filin):
            if line_number % 4 == 1:
                yield line.strip()


def cut_kmer(sequence, kmer_size):
    """Cuts a sequence into a k-mer iterator."""
    sequence_length = len(sequence)
    offset = 0
    while offset + kmer_size <= sequence_length:
        yield sequence[offset:(offset + kmer_size)]
        offset += 1


def build_kmer_dict(fastq_file, kmer_size):
    """Returns a dict of kmer counts."""
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
# b. Build the De Bruijn graph.
# =============================================================================

def build_graph(kmer_count):
    """Build a networkx.DiGraph from a dict of kmer counts."""
    graph = nx.DiGraph()
    for kmer in kmer_count.keys():
        node_1 = kmer[:-1]
        node_2 = kmer[1:]
        graph.add_edge(node_1, node_2, weight=kmer_count[kmer])
    return graph

# =============================================================================
# 2. Graph analysis.
# =============================================================================

def get_starting_nodes(graph):
    """Get the list of starting nodes in a networkx graph."""
    starting_node_list = []
    for node in graph.nodes():
        if not list(graph.predecessors(node)):
            starting_node_list.append(node)
    return starting_node_list


def get_sink_nodes(graph):
    """Get the list of sink nodes in a networkx graph."""
    sink_node_list = []
    for node in graph.nodes():
        if not list(graph.successors(node)):
            sink_node_list.append(node)
    return sink_node_list


def get_contigs(graph, starting_nodes, sink_nodes):
    """Returns a list of tuple (contigs, len(contigs))."""
    contigs = []
    for start_node in starting_nodes:
        for sink_node in sink_nodes:
            try:
                path = nx.shortest_path(graph, start_node, sink_node)
            except:
                continue
            contig = "".join([node[0] for node in path[:-1]] + [path[-1]])
            contigs.append((contig, len(contig)))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contig_tuples, output_filename):
    """Save a list of (contig, len(contig)) tuples in a file."""
    annotation_template = ">contig_{} len={}\n"
    with open(output_filename, "w") as filout:
        for i, contig_tuple in enumerate(contig_tuples):
            filout.write(annotation_template.format(i, contig_tuple[1]))
            filout.write(fill(contig_tuple[0])+"\n")

# =============================================================================
# 3. Graph simplification.
# =============================================================================

# =============================================================================
# a. Bubble resolution
# =============================================================================

def std(values):
    """Compute the standard deviation of a list of values."""
    return statistics.stdev(values)


def path_average_weight(graph, path):
    """Compute the average weight on a path."""
    total = 0
    for node_1, node_2 in zip(path[:-1], path[1:]):
        total += graph[node_1][node_2]["weight"]
    return total/(len(path)-1)


def remove_paths(graph, paths, delete_entry_node, delete_sink_node):
    """Return graph with paths removed.

    Not elegant and potentially dangerous.  Try using remove_nodes_from().
    """
    for path in paths:
        for node_1 in path[1:-1]:
            try:
                graph.remove_node(node_1)
            except:
                pass
        if delete_entry_node:
            try:
                graph.remove_node(path[0])
            except:
                pass
        if delete_sink_node:
            try:
                graph.remove_node(path[-1])
            except:
                pass
    return graph



def select_best_path(graph, paths, path_lengths, avg_path_weights,
                     delete_entry_node=False, delete_sink_node=False):
    """Return a cleaned graph with the supposedly best path.

    Ugly.
    """
    # We put a random seed over 9000.
    random.seed(9001)

    # Sort by weight then by length
    best_weight_indexes = [i for i, weight in enumerate(avg_path_weights)
                           if weight == max(avg_path_weights)]
    best_length_and_weights = [length for i, length in enumerate(path_lengths)
                               if i in best_weight_indexes]
    # Do on length
    best_path_indexes = [i for i in best_weight_indexes
                         if path_lengths[i] == max(best_length_and_weights)]
    #print(best_path_indexes)
    best_path_index = random.choice(best_path_indexes)
    graph = remove_paths(graph, paths[:best_path_index]+paths[(best_path_index+1):],
                         delete_entry_node, delete_sink_node)
    return graph


def solve_bubble(graph, start_node, sink_node):
    """Remove bubbles between start_node and sink_node.
    """
    paths = list(nx.all_simple_paths(graph, start_node, sink_node))
    path_lengths = [len(path) for path in paths]
    avg_path_weights = [path_average_weight(graph, path) for path in paths]

    return select_best_path(graph, paths, path_lengths, avg_path_weights)


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
