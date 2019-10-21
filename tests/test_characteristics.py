"""Tests for graph characteristic"""
import pytest
import os
import networkx as nx
import hashlib
from .context import debruijn
#from .context import debruijn_comp
from debruijn import get_starting_nodes
from debruijn import get_sink_nodes
from debruijn import get_contigs
from debruijn import save_contigs


def test_get_starting_nodes():
    graph = nx.DiGraph()
    graph.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    nodes = get_starting_nodes(graph)    
    assert len(nodes) == 2
    assert 1 in nodes
    assert 3 in nodes

def test_get_sink_nodes():
    graph = nx.DiGraph()
    graph.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    nodes = get_sink_nodes(graph)
    assert len(nodes) == 2
    assert 6 in nodes
    assert 7 in nodes

def test_get_contigs():
    graph = nx.DiGraph()
    graph.add_edges_from([("TC", "CA"), ("AC", "CA"), ("CA", "AG"), ("AG", "GC"), ("GC", "CG"), ("CG", "GA"), ("GA", "AT"), ("GA", "AA")])
    contig_list = get_contigs(graph, ["TC", "AC"], ["AT" , "AA"])
    results = ["TCAGCGAT", "TCAGCGAA", "ACAGCGAT", "ACAGCGAA"]
    assert len(contig_list) == 4
    for contig in contig_list:
        assert contig[0] in results
        assert contig[1] == 8


# def test_get_contigs_comp():
#     graph = nx.DiGraph()
#     graph.add_edges_from([(("AG", "TC"), ("CA", "GT")), (("AC", "TG"), ("CA", "GT")), (("CA", "GT"), ("AG", "TC")), 
#         (("AG", "TC"), ("CG", "GC")), (("CG", "GC"), ("CG", "GC")), (("CG", "GC"), ("CT", "GA")), (("CT", "GA"), ("AT", "TC")), 
#         (("CT", "GA"), ("AA", "TT"))])
#     contig_list = get_contigs(graph, ["TC", "AC"], ["AT" , "AA"])
#     results = ["TCAGCGAT", "TCAGCGAA", "ACAGCGAT", "ACAGCGAA"]
#     assert len(contig_list) == 4
#     for contig in contig_list:
#         assert contig[0] in results
#         assert contig[1] == 8


def test_save_contigs():
    test_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "test.fna"))
    contig = [("TCAGCGAT", 8), ("TCAGCGAA",8), ("ACAGCGAT", 8), ("ACAGCGAA", 8)]
    save_contigs(contig, test_file)
    with open(test_file, 'rb') as contig_test:
        assert hashlib.md5(contig_test.read()).hexdigest() == "ca84dfeb5d58eca107e34de09b3cc997"