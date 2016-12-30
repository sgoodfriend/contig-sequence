#!/usr/bin/env python

import argparse
import sys

from Bio import SeqIO

from seq_join_graph_components import SeqPrefixHashMap

"""
sequence_joiner outputs the reconstructed sequence of DNA fragments. If no sequence was found,
outputs nothing and returns error code -1.
"""
__author__ = 'goodfriend-scott'


class SeqPath(object):
    """
    SeqPath is a path container useful when you want to be able to update the path during a
    recursive call.
    """
    def __init__(self):
        self.path = []

    def __len__(self):
        return len(self.path)


class SeqJoinLongestPath(object):
    """
    SeqJoinLongestPath finds the longest path of edges in a given SeqJoinGraph.

    SeqJoinLongestPath finds the longest path with two different algorithms with vastly different
    performance guarantees. If the graph is a single-component linear path, then we can traverse
    the nodes from the root to the end in O(N) time. If there is branching, then we fall back
    to a depth first search longest path algorithm that can handle cycles. The worst case runtime
    for this algorithm is O(2^E) if every node is connected to every other node for every possible
    length (think all fragments are just A's). E is the number of edges, which is N*L in the
    worst possible case.
    """
    def __init__(self, seq_graph):
        self.graph = seq_graph
        self.root_nodes = self.graph.root_nodes()
        self.end_nodes = self.graph.end_nodes()
        if self.is_linear_graph():
            self.longest_path = self.linear_longest_path()
        else:
            self.longest_path = self.dfs_longest_path()

    def is_linear_graph(self):
        """
        Returns if the graph is a linear single-component graph.

        :return: Boolean of if the graph is a linear single component graph.
        .. complexity:: O(N) time.
        """
        if len(self.root_nodes) != 1 or len(self.end_nodes) != 1:
            return False
        num_degree_one_nodes = sum([len(edges) == 1
                                    for n, edges in self.graph.node_to_edges.iteritems()])
        return num_degree_one_nodes == len(self.graph.nodes) - 1

    def linear_longest_path(self):
        """
        Computes the longest path of edges assuming a linear single-component graph.
        :return: List representing the longest path of edges in the graph.
        .. complexity:: O(N) time.
        """
        path = []
        node = self.root_nodes[0]
        while len(self.graph.node_to_edges[node]):
            edge = self.graph.node_to_edges[node][0]
            path.append(edge)
            node = edge.destination
        return path

    def dfs_longest_path(self):
        """
        Computes the longest path of edges using the general longest path solving algorithm.
        :return: List representing the longest path of edges in the graph.
        .. complexity:: Worst-case O(2^E) given a graph where every possible length of node
            overlap connects to every other node. Typical performance is likely closer to O(N) if
            you assume that branching occurs sparsely. Performance has an added N multiplier for
            each inability to determine start and end node (O(N^3) if you cannot determine either
            start or end node before starting).
        .. warning:: Uses recursion where depth of the stack can get up to the number of nodes.
            The default sys.recursionlimit is likely 1000, so if you are going above this number
            of fragments, you'll should increase the recursionlimit with sys.setrecursionlimit.
        .. todo:: Switch recursion to iterative solution to support graphs above size ~1000.
        """
        start_candidates = [self.root_nodes[0]] if len(self.root_nodes) == 1 else self.graph.nodes
        end_candidates = [self.end_nodes[0]] if len(self.end_nodes) == 1 else self.graph.nodes
        longest_path = SeqPath()
        for s_node in start_candidates:
            for e_node in end_candidates:
                self.__dfs(s_node, e_node, longest_path)
        return longest_path.path

    def __dfs(self, s_node, e_node, longest_path, current_path=None, visiting_set=None):
        if current_path is None:
            current_path = []
        if visiting_set is None:
            visiting_set = set()
        if s_node == e_node:
            if len(current_path) > len(longest_path):
                longest_path.path = current_path[:]
            return
        visiting_set.add(s_node)
        for edge in self.graph.node_to_edges[s_node]:
            if edge.destination not in visiting_set:
                current_path.append(edge)
                self.__dfs(edge.destination, e_node, longest_path, current_path, visiting_set)
                current_path.pop()
        visiting_set.remove(s_node)


class NoFragmentsException(Exception):
    """
    Exception raised when graph has no nodes (no SeqRecords inputted) when non-zero nodes is
    assumed.
    """
    pass


class SingleSequenceAssertionFailure(Exception):
    """
    Exception raised if graph doesn't form a singly connected structure.
    """
    pass


class LinearSequenceAssertionFailure(Exception):
    """
    Exception raised if the singly connected structure path would form a cycle because the last
    node can join with the first node.
    """
    pass


class SeqJoinGraph(object):
    """
    SeqJoinGraph outputs the reconstructed sequence by finding the longest path of SeqJoinEdges.

    SeqJoinGraph's nodes are SeqNodes, which represent SeqRecords with precomputed hashing data.
    The edges are SeqJoinEdges, which represent a possible way to combine the source and
    destination sequences and the index at where the overlap begins. SeqJoinGraph's primary method
    sequence() finds the non-branching longest path of edges that traverses all nodes and then
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

    @classmethod
    def graph_from_filename(cls, filename):
        return cls.graph_from_seq_records([s for s in SeqIO.parse(filename, 'fasta')])

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

    def root_nodes(self):
        """
        Returns the nodes that are not a destination of any edge.

        :return: List of SeqNodes that are not a destination of any edge.
        ..complexity:: Time O(N) where N is the number of SeqNodes.
        """
        unvisitable_nodes = set(self.nodes)
        for n, edges in self.node_to_edges.iteritems():
            for e in edges:
                if e.destination in unvisitable_nodes:
                    unvisitable_nodes.remove(e.destination)
        return list(unvisitable_nodes)

    def end_nodes(self):
        """
        Returns nodes with no outgoing edges.

        :return: List of SeqNodes with no outgoing edges.
        .. complexity:: Time O(N) where N is the number of SeqNodes.
        """
        nodes = [n for n in self.nodes if len(self.node_to_edges[n]) == 0]
        return nodes

    def __longest_path(self):
        """
        Returns the longest path of edges starting from the root node to the end node.

        :return: List of SeqJoinEdges representing a path from the root to the end node, passing
            through every node.
        .. complexity:: If the graph forms a single, linearly connected structure: O(N).
            Otherwise, the solver falls back to a dfs longest path solver that accounts for cycles.
            Worst case for dfs: O(2^E), where E is the number of edges. In the worst-case the
            number of edges is N*L if every possible length also matches (think of fragments of
            all A's).
        """
        solver = SeqJoinLongestPath(self)
        return solver.longest_path

    def __assert_linearity(self, longest_path):
        """
        Checks to make sure the destination does not have an edge leading to the source. Raises a
        LinearSequenceAssertionFailure if an edge connects the end to the beginning.

        :param longest_path: The longest path of edges existing in the graph.
        :raises: LinearSequenceAssertionFailure if an edge connects the end to the beginning.
        """
        source = longest_path[0].source
        destination = longest_path[-1].destination
        for edges in self.node_to_edges[destination]:
            if edges.destination == source:
                raise LinearSequenceAssertionFailure(self)

    def sequence(self):
        """
        Returns the reconstructed sequence.

        :return: String reconstructed sequence.
        :raises: NoFragmentsException if graph has no nodes (no SeqRecords inputted).
            SingleSequenceAssertionFailure if graph doesn't form a singly connected structure.
            LinearSequenceAssertionFailure if the singly connected structure path would form a
            cycle because the last node can join with the first node.
        .. complexity:: Time O(N) for a single, linearly connected structure. Worst-case: O(2^E),
            where E can be N*L if every possible suffix length matches every node (think of
            fragments of all A's).
        """
        if len(self.nodes) == 0:
            raise NoFragmentsException(self)
        longest_path = self.__longest_path()
        if len(longest_path) != len(self.nodes) - 1:
            raise SingleSequenceAssertionFailure(self)
        if len(longest_path) == 0:
            return str(self.nodes[0].seq())
        else:
            self.__assert_linearity(longest_path)
            seq_list = []
            for e in longest_path:
                seq_list.append(e.source.seq()[:e.source_idx])
            seq_list.append(longest_path[-1].destination.seq())
            return ''.join([str(seq) for seq in seq_list])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Outputs the single reconstructed sequence of DNA fragments, otherwise '
                    'exits with error code and outputs nothing. '
                    'Error codes: -1 no inputted fragments, -2 could not find a single sequence '
                    'that accounts for all fragments, -3 sequence is found to be circular.')
    parser.add_argument('filename', help='FASTA filename')
    args = parser.parse_args()
    graph = SeqJoinGraph.graph_from_filename(args.filename)
    try:
        rec_seq = graph.sequence()
        print rec_seq
    except NoFragmentsException:
        sys.exit(-1)
    except SingleSequenceAssertionFailure:
        sys.exit(-2)
    except LinearSequenceAssertionFailure:
        sys.exit(-3)
