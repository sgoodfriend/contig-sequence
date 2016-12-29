#!/usr/bin/env python

import argparse
from collections import defaultdict

from Bio import SeqIO

"""
sequence_joiner outputs the reconstructed sequence of DNA fragments.
"""
__author__ = 'goodfriend-scott'


class SeqHash(object):
    """
    SeqHash is a class object that contains hashing constants and helper method used in the
    SeqPrefixHashMap.
    """
    Q = 36028797018963913  # Next prime below 2^55.
    R = 128  # Size of alphabet
    RQ = [1]

    @classmethod
    def rq_for_length(cls, length):
        """
        Returns the multiplier modded by Q of the 0th-index of a string of characters of length
        `length`.

        :param length: The length of the string of characters.
        :return: The multiplier of the 0-th index of a string of `length` characters.
        """
        while len(cls.RQ) < length:
            cls.RQ.append((cls.RQ[-1] * cls.R) % cls.Q)
        return cls.RQ[length - 1]


class SeqNode(object):
    """
    Graph node representing a SeqRecord.

    .. warnings:: Hash h must be set before you can use __hash_ or __eq__.
    """

    def __init__(self, seq_record):
        self.seq_record = seq_record
        self.h = None

    def seq(self):
        """
        Helper function that returns the seq_record's Seq object.

        :return: The seq_record's Seq object.
        """
        return self.seq_record.seq

    def seq_hash(self, idx=0, length=None):
        """
        Calculates the substring hash of the sequence.

        :param idx: The index to start the substring from. Defaults to 0, the beginning.
        :param length: The length of the substring. Defaults to the remaining length of the seq
            from idx.
        :return: The int hash of the subsequence.
        """
        _seq = self.seq()
        if length is None:
            length = len(_seq) - idx
        h = 0
        for i in xrange(idx, idx + length):
            h = (h * SeqHash.R + ord(_seq[i])) % SeqHash.Q
        return h

    def __hash__(self):
        assert self.h is not None, 'Hash h must be set before hash function can be used.'
        return (self.h << 8) ^ hash(self.seq_record.name)

    def __eq__(self, other):
        return (self.h, self.seq_record.name) == (other.h, other.seq_record.name)

    def __repr__(self):
        return '{}'.format(self.seq_record.name)


class SeqJoinEdge(object):
    """
    Graph edge for SeqJoinGraph.

    SeqJoinEdge represents the source and destination nodes and the index in the source sequence
    where overlap begins.
    """

    def __init__(self, source, source_idx, destination):
        self.source = source
        self.source_idx = source_idx
        self.destination = destination

    def __repr__(self):
        return '{0} at {1} to {2}'.format(
            self.source,
            self.source_idx,
            self.destination
        )


class SeqPrefixHashMap(object):
    """
    SeqPrefixHashMap maps subsequence hashes to SeqNodes that have matching prefixes.

    SeqPrefixHashMap uses a multiple pattern Rabin-Karp string searching algorithm by first
    building up a 2-level dictionary (key: subsequence length, subsequence hash, value:
    list of SeqNodes that have matching prefix length and hash). node_to_edges() next goes through
    all the sequences checking all possible subsequence suffixes that can match the prefix using
    the map.
    """
    def __init__(self, _seq_records):
        """
        Builds up the SeqPrefixHashMap from an iterable of SeqRecords.
        :param _seq_records: Iterable of SeqRecords
        .. complexity:: Time: O(N*l)=O(n), Space: O(N*l)=O(n), where N is number of SeqRecords,
            l is the average length of fragments, and n is the input size (N*l, effectively)
        """
        self.min_join_length = (min([len(s.seq) for s in _seq_records]) + 1) / 2
        self.length_to_hash_to_seqs_map = defaultdict(self.list_generating_defaultdict)
        self.nodes = []
        for s_rec in _seq_records:
            node = SeqNode(s_rec)
            _seq = s_rec.seq
            h = node.seq_hash(length=self.min_join_length)
            self.length_to_hash_to_seqs_map[self.min_join_length][h].append(node)
            for l in xrange(self.min_join_length + 1, len(_seq)):
                h = (h * SeqHash.R + ord(_seq[l - 1])) % SeqHash.Q
                self.length_to_hash_to_seqs_map[l][h].append(node)
            node.h = h
            self.nodes.append(node)

    @classmethod
    def list_generating_defaultdict(cls):
        return defaultdict(list)

    def node_to_edges(self):
        """
        Returns the mapping of nodes to SeqJoinEdges using the internal SeqPrefixHashMap.

        For each SeqNode, calculate a subsequence hash for each possible suffix length from
        min_join_length to length and see if there are candidate SeqNodes with matching prefix
        length and hash in the map. For each candidate explicitly check to make sure the suffix
        matches the candidate's prefix character by character (Las Vegas variant of Rabin-Karp).

        :return: Dictionary of nodes to SeqJoinEdges.
        .. complexity:: Typical Time: O(N*l)=O(n), Space: O(N*l)=O(n), where N is number of
            SeqRecords, l is the average length of fragments, and n is the input size (N*l,
            effectively).
            Worst-case O(l^2*N^2)=O(n^2) time complexity can occur if all the hashes at each length
            collide during building of SeqPrefixHashMap. The additional N term comes from having
            to check each sequence, the additional l term comes from the character by character
            check.
        """
        node_to_edges = defaultdict(list)
        for n in self.nodes:
            seq = n.seq_record.seq
            idx = self.min_join_length
            h = n.seq_hash(idx)
            while True:
                length = len(seq) - idx
                if h in self.length_to_hash_to_seqs_map[length]:
                    for candidate_node in self.length_to_hash_to_seqs_map[length][h]:
                        if seq[idx:] == candidate_node.seq()[:length]:
                            node_to_edges[n].append(SeqJoinEdge(n, idx, candidate_node))
                idx -= 1
                if idx < 0:
                    break
                rq = SeqHash.rq_for_length(len(seq) - idx)
                h = (h + rq * ord(seq[idx])) % SeqHash.Q
        return node_to_edges


class SeqJoinGraph(object):
    """
    SeqJoinGraph outputs the reconstructed sequence by finding the maximum path of SeqJoinEdges.

    SeqJoinGraph's nodes are SeqNodes, which represent SeqRecords with precomputed hashing data.
    The edges are SeqJoinEdges, which represent a possible way to combine the source and
    destination sequences and the index at where the overlap begins. SeqJoinGraph's primary method
    sequence() finds the non-branching maximum path of edges that traverses all nodes and then
    uses the edges to output the reconstructed sequence as a string.
    """
    def __init__(self, nodes, node_to_edges):
        self.nodes = nodes
        self.node_to_edges = node_to_edges

    @classmethod
    def graph_from_seq_records(cls, _seq_records):
        """
        Primary constructor method that generates a SeqPrefixHashMap to generate nodes and edges.

        :param _seq_records: Iterable of SeqRecords.
        :return: SeqJoinGraph object populated with nodes and edges.
        ..complexity:: Typical: Time & Space O(n)=O(N*l), Worst-case: Time O(n^2)=O(N^2*l^2), where
            N is number of SeqRecords, l is the average length of fragments, and n is the input
            size (N*l, effectively).
            This is simply the SeqPrefixHashMap map generation followed by edge finding.
        """
        prefix_map = SeqPrefixHashMap(_seq_records)
        return cls(prefix_map.nodes, prefix_map.node_to_edges())

    def __repr__(self):
        lines = ['{0} nodes, {1} edges'.format(
            len(self.nodes), sum([len(edges) for n, edges in self.node_to_edges.iteritems()])
        )]
        for n in self.nodes:
            edges = self.node_to_edges[n]
            lines.append(' {0}: {1}'.format(
                n, ', '.join([str(e) for e in edges])
            ))
        return '\n'.join(lines)

    def root_node(self):
        """
        Returns the node that is not a destination of any edge, if it is the only one.

        :return: SeqNode that is not a destination of any edge, if only one. Otherwise, returns
            None.
        ..complexity:: Time O(N) where N is the number of SeqNodes.
        """
        unvisitable_nodes = set(self.nodes)
        for n, edges in self.node_to_edges.iteritems():
            for e in edges:
                unvisitable_nodes.remove(e.destination)
        if len(unvisitable_nodes) == 1:
            return unvisitable_nodes.pop()
        else:
            return None

    def end_node(self):
        """
        Returns the first node that has no outgoing edges.

        :return: First SeqNode with no outgoing edges.
        .. complexity:: Time O(N) where N is the number of SeqNodes.
        """
        for n, edges in self.node_to_edges.iteritems():
            if len(edges) == 0:
                return n
        return None

    def maximum_path(self):
        """
        Returns the maximum path of edges starting from the root node to the end node.

        :return: List of SeqJoinEdges representing a path from the root to the end node, passing
            through every node.
        .. warning:: Current implementation assumes a single unbranched path exists between root
            and end nodes.
        .. complexity:: Time O(N) where N is the number of SeqNodes thanks to the unbranched
            path assumption.
        """
        assert (sum([len(edges) == 1 for n, edges in self.node_to_edges.iteritems()]) ==
                len(self.nodes) - 1), \
            'maximum_path only supports unbranched graphs.'
        path = []
        node = self.root_node()
        while len(self.node_to_edges[node]):
            edge = self.node_to_edges[node][0]
            path.append(edge)
            node = edge.destination
        return path

    def sequence(self):
        """
        Returns the reconstructed sequence.

        :return: String reconstructed sequence.
        .. warning:: Current implementation assumes a single unbranched path exists between root
            and end nodes.
        .. note:: Iterate through each SeqJoinEdge of the maximum path. At each edge, append the
            non-overlapping sequence from the origin into a list of subsequences that are joined
            together. Since we have not been recording the destination node's sequence, we need
            to cap off the list with the end_node's sequence.
        .. complexity:: Time O(N), where N is the number of SeqNodes.
        """
        maximum_path = self.maximum_path()
        seq_list = []
        for e in maximum_path:
            seq_list.append(e.source.seq()[:e.source_idx])
        end_node = self.end_node()
        seq_list.append(end_node.seq())
        return ''.join([str(seq) for seq in seq_list])

if __name__ == '__main__':
    # Parse filename argument.
    parser = argparse.ArgumentParser(
        description='Outputs the reconstructed sequence of DNA fragments.')
    parser.add_argument('filename', help='FASTA filename')
    args = parser.parse_args()
    # Parse SeqRecords from FASTA file.
    seq_records = [s_record for s_record in SeqIO.parse(args.filename, 'fasta')]
    # Generate SeqJoinGraph, which is able to return the reconstructed sequence.
    graph = SeqJoinGraph.graph_from_seq_records(seq_records)
    rec_seq = graph.sequence()
    print rec_seq
