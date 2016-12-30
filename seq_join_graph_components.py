from collections import defaultdict

"""
seq_join_graph_components contains helper implementations for SeqJoinGraph creations.
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

    def __hash__(self):
        return hash(self.source) ^ hash(self.destination) << 1 ^ self.source_idx

    def __signature(self):
        return self.source, self.source_idx, self.destination

    def __eq__(self, other):
        return self.__signature() == other.__signature()


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
        self.length_to_hash_to_seqs_map = defaultdict(self.list_generating_defaultdict)
        self.nodes = []
        for s_rec in _seq_records:
            node = SeqNode(s_rec)
            _seq = s_rec.seq
            min_join_length = (len(_seq) + 1) / 2
            h = node.seq_hash(length=min_join_length)
            self.length_to_hash_to_seqs_map[min_join_length][h].append(node)
            for l in xrange(min_join_length + 1, len(_seq) + 1):
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
        A fragment joining with itself is explicitly not allowed.

        :return: Dictionary of nodes to SeqJoinEdges.
        .. complexity:: Typical Time: O(N*l)=O(n), Space: O(N*l)=O(n), where N is number of
            SeqRecords, l is the average length of fragments, and n is the input size (N*l,
            effectively).
            Worst-case O(l^2*N^2)=O(n^2) time complexity can occur if all the hashes at each length
            collide during building of SeqPrefixHashMap. The additional N term comes from having
            to check each sequence, the additional l term comes from the character by character
            check.
        .. note:: Self joining is explicitly not allowed.
        """
        node_to_edges = defaultdict(list)
        for n in self.nodes:
            seq = n.seq_record.seq
            idx = len(seq) - (len(seq) + 1) / 2
            h = n.seq_hash(idx)
            while True:
                length = len(seq) - idx
                if h in self.length_to_hash_to_seqs_map[length]:
                    for candidate_node in self.length_to_hash_to_seqs_map[length][h]:
                        if candidate_node == n:
                            continue
                        if seq[idx:] == candidate_node.seq()[:length]:
                            node_to_edges[n].append(SeqJoinEdge(n, idx, candidate_node))
                idx -= 1
                if idx < 0:
                    break
                rq = SeqHash.rq_for_length(len(seq) - idx)
                h = (h + rq * ord(seq[idx])) % SeqHash.Q
        return node_to_edges
