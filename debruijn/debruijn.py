#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:14:42 2019

@author: lyang_crosson
"""
# pylint: disable=too-many-arguments

import argparse
import os
import random
import statistics

import networkx as nx

# =============================================================================
# 1. Create the De Bruijn graph
# =============================================================================

# =============================================================================
# a. Identify unique kmers
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
# b. Build the De Bruijn graph
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
# 2. Graph analysis
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
            except (nx.NodeNotFound, nx.NetworkXNoPath):
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
# 3. Graph simplification
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
        try:
            total += graph[node_1][node_2]["weight"]
        # No path between the nodes.
        except KeyError:
            pass
    return total/(len(path)-1) if total else 0


def remove_paths(graph, paths, delete_entry_node, delete_sink_node):
    """Return graph with paths removed."""
    for path in paths:
        # Use the boolean properties to save some lines.
        graph.remove_nodes_from(path[(not delete_entry_node):
                                     (None if delete_sink_node else -1)])
    return graph


def select_best_path(graph, paths, path_lengths, avg_path_weights,
                     delete_entry_node=False, delete_sink_node=False):
    """Return a cleaned graph with the supposedly best path.

    The function assumes paths not empty.

    """
    # We put a random seed over 9000.
    random.seed(9001)

    # Indexes of paths with best weight.
    best_weight_indexes = [i for i, weight in enumerate(avg_path_weights)
                           if weight == max(avg_path_weights)]
    # Paths with best lengths for paths with best weights.
    best_length_and_weights = [length for i, length in enumerate(path_lengths)
                               if i in best_weight_indexes]
    # Indexes of paths with best weights and length.
    best_path_indexes = [i for i in best_weight_indexes
                         if path_lengths[i] == max(best_length_and_weights)]

    best_path_index = random.choice(best_path_indexes)
    graph = remove_paths(graph, paths[:best_path_index]+paths[(best_path_index+1):],
                         delete_entry_node, delete_sink_node)
    return graph


def solve_bubble(graph, start_node, sink_node):
    """Remove bubbles between start_node and sink_node."""
    paths = list(nx.all_simple_paths(graph, start_node, sink_node))
    path_lengths = [len(path) for path in paths]
    avg_path_weights = [path_average_weight(graph, path) for path in paths]

    return select_best_path(graph, paths, path_lengths, avg_path_weights)


def simplify_bubbles(graph):
    """Remove all bubbles from a graph."""
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    # Find all bubbles.
    for start_node in starting_nodes:
        for sink_node in sink_nodes:
            # Find the smallest englobing bubble.
            current_entry = start_node
            current_exit = sink_node
            successors = list(graph.successors(current_entry))
            predecessors = list(graph.predecessors(current_exit))
            while len(successors) < 2 and successors:
                current_entry = successors[0]
                successors = list(graph.successors(current_entry))
            while len(predecessors) < 2 and predecessors:
                current_exit = predecessors[0]
                predecessors = list(graph.predecessors(current_exit))
            # A path exists between the nodes.
            if list(nx.all_simple_paths(graph, current_entry, current_exit)):
                graph = solve_bubble(graph, current_entry, current_exit)
    return graph

# =============================================================================
# b. Tips detection
# =============================================================================

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips from an oriented graph."""
    # Sanitize the graph.
    graph = simplify_bubbles(graph)
    # All existing tips.
    tips = []
    # Resolve all tips.
    for start_node in starting_nodes:
        current_node = start_node
        path = [current_node]
        successors = list(graph.successors(current_node))
        predecessors = list(graph.predecessors(current_node))
        # Iterate till the node has 2 predecessors/successors, or is sink node.
        while len(successors) < 2 and len(predecessors) < 2 and successors:
            current_node = successors[0]
            path.append(current_node)
            successors = list(graph.successors(current_node))
            predecessors = list(graph.predecessors(current_node))
        tips.append(path)

#    if tips:
    path_lengths = [len(path) for path in tips]
    avg_path_weights = [path_average_weight(graph, path) for path in tips]

    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_entry_node=True)

    return graph


def solve_out_tips(graph, sink_nodes):
    """Remove out tips from an oriented graph."""
    graph = simplify_bubbles(graph)

    tips = []
    # Resolve all tips.
    for sink_node in sink_nodes:
        current_node = sink_node
        path = [current_node]
        successors = list(graph.successors(current_node))
        predecessors = list(graph.predecessors(current_node))
        # Iterate till node has 2 predecessors/successors, or is start node.
        while len(successors) < 2 and len(predecessors) < 2 and predecessors:
            current_node = predecessors[0]
            path.append(current_node)
            successors = list(graph.successors(current_node))
            predecessors = list(graph.predecessors(current_node))
        # Reverse the path for path_average_weight.
        tips.append(path[::-1])

    path_lengths = [len(path) for path in tips]
    avg_path_weights = [path_average_weight(graph, path) for path in tips]

    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_sink_node=True)
    return graph


def user_input():
    """Parse the submitted command line."""
    parser = argparse.ArgumentParser(
        description="Read single-end fastq file and returns.")

    # arguments
    parser.add_argument("-i", "--input", required=True,
                        help=("name of the fastq file."))
    parser.add_argument("-k", "--kmer", type=int, default=21,
                        help="length of kmers (default: 21).")
    parser.add_argument("-o", "--contig", required=True,
                        help="name of contig file.")

    # get all arguments
    options = parser.parse_args()

    return options


def main():
    """Entry point for debruijn's CLI."""
    options = user_input()
    print(f"kmer size: {options.kmer}")
    kmer_count = build_kmer_dict(options.input, options.kmer)
    print(f"number of kmer found: {len(kmer_count)}")
    graph = build_graph(kmer_count)

    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    # Simplify the graph.
    graph = simplify_bubbles(graph)
    print(f"after bubbles {len(graph.edges)}")
    graph = solve_entry_tips(graph, starting_nodes)
    print(f"after entry tips {len(graph.edges)}")
    graph = solve_out_tips(graph, sink_nodes)
    print(f"after out tips {len(graph.edges)}")

    new_starting_nodes = get_starting_nodes(graph)
    new_sink_nodes = get_sink_nodes(graph)

    print("Computing contigs...")
    contig_tuples = get_contigs(graph, new_starting_nodes, new_sink_nodes)

    print("Saving the contigs into a file...")
    save_contigs(contig_tuples, options.contig)


if __name__ == "__main__":
    main()
