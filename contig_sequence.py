#!/usr/bin/env python

import argparse
import sys

try:
    from Bio import SeqIO
except ImportError:
    print 'Biopython must be installed. If you have pip, run "pip install biopython"'
    sys.exit(1)

from seq_join_graph_components import SeqPrefixHashMap

"""
contig_sequence outputs the reconstructed sequence of DNA fragments. If no sequence was found,
outputs nothing and returns error code -1.
"""
__author__ = 'goodfriend-scott'


class NoFragmentsException(Exception):
    """
    Exception raised when graph has no nodes (no SeqRecords inputted); however, non-zero nodes are
    assumed.
    """
    pass


class SingleSequenceAssertionFailure(Exception):
    """
    Exception raised if longest path doesn't traverse all nodes through the graph.
    """
    pass


class LinearSequenceAssertionFailure(Exception):
    """
    Exception raised if the longest path would form a cycle because the last node can join with
    the first node.
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
    def graph_from_seq_records(cls, seq_records):
        """
        Primary constructor method that generates a SeqPrefixHashMap to generate nodes and edges.

        :param seq_records: Iterable of SeqRecords.
        :return: SeqJoinGraph object populated with nodes and edges.
        ..complexity:: Typical: Time & Space O(n)=O(N*l),
            Worst-case: Time & Space O(n^2)=O(N^2*l^2), where N is number of SeqRecords, l is the
            average length of fragments, and n is the input size (N*l, effectively). This is simply
            the SeqPrefixHashMap map generation followed by edge finding.
        """
        prefix_map = SeqPrefixHashMap(seq_records)
        return cls(prefix_map.nodes, prefix_map.node_to_edges())

    @classmethod
    def graph_from_filename(cls, filename):
        """
        Construct SeqJoinGraph from a FASTA file.
        """
        return cls.graph_from_seq_records([s for s in SeqIO.parse(filename, 'fasta')])

    def __repr__(self):
        lines = ['{0} nodes, {1} edges'.format(
            len(self.nodes), sum([len(edges) for edges in self.node_to_edges.values()])
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
        .. complexity:: Time O(N+E) where N is the number of SeqNodes, E is the number
            SeqJoinEdges. Typical: O(n). Worst-case: O(n^2).
        """
        unvisitable_nodes = set(self.nodes)
        for edges in self.node_to_edges.values():
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
            number of edges is N^2*L if every possible length also matches (think of fragments of
            all A's).
        """
        solver = SeqJoinLongestPath(self)
        return solver.longest_path

    def __assert_linearity(self, longest_path):
        """
        Checks to make sure the destination does not have an edge leading to the source. Raises a
        LinearSequenceAssertionFailure if an edge connects the end to the beginning.

        :param longest_path: The non-zero length longest path of edges existing in the graph.
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
            SingleSequenceAssertionFailure if longest path doesn't traverse all nodes through the
                graph.
            LinearSequenceAssertionFailure if the longest path would form a cycle because the last
            node can join with the first node.
        .. complexity:: Time O(N) for a single, linearly connected structure. Worst-case: O(2^E),
            where E can be N^2*L if every possible suffix length matches every node (think of
            fragments of all A's).
        """
        if len(self.nodes) == 0:
            raise NoFragmentsException(self)
        longest_path = self.__longest_path()
        if len(longest_path) != len(self.nodes) - 1:
            raise SingleSequenceAssertionFailure(self)
        if len(longest_path) == 0:
            return str(self.nodes[0].get_seq())
        else:
            self.__assert_linearity(longest_path)
            seq_list = []
            for e in longest_path:
                seq_list.append(e.source.get_seq()[:e.source_idx])
            seq_list.append(longest_path[-1].destination.get_seq())
            return ''.join([str(seq) for seq in seq_list])


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
    length (think all fragments are just A's). E is the number of edges, which is N^2*L in the
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
        Returns if the graph is likely a linear single-component graph.

        is_linear_graph checks if there is 1 root node, 1 end node, and all nodes but the end
        node have 1 edge. This does not guarantee a single-component linear graph since there
        could be one linear component and a number of simple cycles.
        :return: Boolean of if the graph is likely a linear single component graph.
        .. complexity:: O(N) time.
        """
        if len(self.root_nodes) != 1 or len(self.end_nodes) != 1:
            return False
        num_degree_one_nodes = sum([len(edges) == 1
                                    for edges in self.graph.node_to_edges.values()])
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
            start or end node before starting). Each additional edge above degree 1 doubles the
            time taken.
        .. warning:: Uses recursion where maximum depth of the stack is the number of nodes.
            The default sys.recursionlimit is likely 1000. So if you are going getting near this
            number of fragments, you should increase the recursionlimit with sys.setrecursionlimit.
        .. todo:: Switch recursion to iterative solution to support graphs above size ~1000.
        """
        start_candidates = [self.root_nodes[0]] if len(self.root_nodes) == 1 else self.graph.nodes
        end_candidates = [self.end_nodes[0]] if len(self.end_nodes) == 1 else self.graph.nodes
        cur_longest_path = SeqPath()
        for s_node in start_candidates:
            for e_node in end_candidates:
                self.__dfs(s_node, e_node, cur_longest_path)
        return cur_longest_path.path

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Outputs the single reconstructed sequence of DNA fragments, otherwise '
                    'exits with error code and outputs nothing. '
                    'Error codes: -1 no inputted fragments, -2 could not find a single sequence '
                    'that accounts for all fragments, -3 sequence is found to be circular, '
                    '1 Python cannot find the Biopython package.')
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
