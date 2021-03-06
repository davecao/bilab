from __future__ import print_function
# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-26 14:34:10
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-07-29 11:03:03

__all__ = ['DeBruijnGraph']


class DeBruijnGraph:
    ''' A de Bruijn directed multigraph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes of
        the de Bruijn graph are k-1-mers and edges correspond to the
        k-mer that joins a left k-1-mer to a right k-1-mer. '''

    @staticmethod
    def chop(st, k):
        ''' Chop a string up into k mers of given length '''
        for i in xrange(0, len(st) - (k - 1)):
            yield (st[i:i + k], st[i:i + k - 1], st[i + 1:i + k])

    class Node:
        ''' Node in a de Bruijn graph, representing a k-1 mer.  We keep
            track of # of incoming/outgoing edges so it's easy to check
            for balanced, semi-balanced. '''
        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0

        def isSemiBalanced(self):
            return abs(self.nin - self.nout) == 1

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash(self.km1mer)

        def __str__(self):
            return self.km1mer

    def __init__(self, strIter, k, circularize=False):
        ''' Build de Bruijn multigraph given string iterator and k-mer
            length k '''
        self.G = {}     # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        for st in strIter:
            if circularize:
                st += st[:k - 1]
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)
                nodeL.nout += 1
                nodeR.nin += 1
                self.G.setdefault(nodeL, []).append(nodeR)
        # Iterate through nodes and tally how many are balanced,
        # semi-balanced, or neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian walk (not cycle)
        self.head, self.tail = None, None
        for node in self.nodes.itervalues():
            if node.isBalanced():
                self.nbal += 1
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def nnodes(self):
        ''' Return # nodes '''
        return len(self.nodes)

    def nedges(self):
        ''' Return # edges '''
        return len(self.G)

    def hasEulerianWalk(self):
        ''' Return true iff graph has Eulerian walk. '''
        return self.nneither == 0 and self.nsemi == 2

    def hasEulerianCycle(self):
        ''' Return true iff graph has Eulerian cycle. '''
        return self.nneither == 0 and self.nsemi == 0

    def isEulerian(self):
        ''' Return true iff graph has Eulerian walk or cycle '''
        # technically, if it has an Eulerian walk
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    def eulerianWalkOrCycle(self):
        ''' Find and return Eulerian walk or cycle (as appropriate) '''
        assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            assert self.head is not None
            assert self.tail is not None
            g.setdefault(self.tail, []).append(self.head)
        # graph g has an Eulerian cycle
        tour = []
        src = g.iterkeys().next()  # pick arbitrary starting node

        def __visit(n):
            while len(g[n]) > 0:
                dst = g[n].pop()
                __visit(dst)
            tour.append(n)

        __visit(src)
        tour = tour[::-1][:-1]  # reverse and then take all but last node

        if self.hasEulerianWalk():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        # Return node list
        return map(str, tour)
