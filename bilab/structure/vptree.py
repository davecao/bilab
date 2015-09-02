# -*- coding: utf-8 -*-
"""
    This module defines :class:`VPTree` 
    To do list:
    1. HeapItem 
    2. DistanceComparator
"""
from __future__ import print_function

import uuid
import heapq
import random
import sys
import types
import numpy as np
from bilab import LOGGER
from bilab.structure.Heap import Heap
from bilab.utilities import nth_element
from itertools import permutations

__all__ = ['VPTree']

class PriorityQueue(object):

    """ Implementation of Priority Queue """
    def __init__(self, sequence=None, key=None, inplace=True):
        self._queue = Heap(sequence=None, key=None, inplace=inplace)

    def push(self, item):
        """
            Push item into priority queue

        Args:
            item (object): an object stored in heap.
        """
        #heapq.heappush(self._queue, item)
        #self._index += 1
        self._queue.push(item)

    def pop(self):
        """
            Removes item with the smallest value off the priority queue.

        Return:
            an object which has the smallest item off the heap,
            maintaining the heap invariant.
        """
        return self._queue.pop()

    def popm(self):
        """
            Removes item with the largest value off the priority queue.

        Return:
            an object which has the largest item off the heap,
            maintaining the heap invariant.
        """
        return self._queue.popm()

    def size(self):
        """
        Return:
            the size of current queue
        """
        return len(self._queue)

    def empty(self):
        """
        Return:
            return whether the priority queue is empty: i.e., whether its size
            is zero.
        """
        return self.size() == 0

    def top(self):
        """
            return the item with the largest value 
        """
        largest_ = self._queue.largest()
        if len(largest_) == 1:
            largest_ = largest_[0]
        return largest_

class Node(object):
    def __init__(self, identifier=None, index=0, value=0.0, name=None, threshold=0.0):
        super(Node, self).__init__()
        self.__identifier = (str(uuid.uuid1()) if identifier is None else self.__sanitize_id(str(identifier)))
        # index of self._items but pay attention to that self._items
        #  is rearranged while building up the tree
        self.index = index 
        # 
        self.value = value
        self.name = name
        self.threshold = threshold
        self.left  = None
        self.right = None

    def __iter__(self):
        """ return the iterator that iterates through the elements in the BST 
        rooted at this node in an inorder sequence SLOWER SLOWER"""

        if self.left:
            for elem in self.left:
                yield elem
        yield (self.index, self.threshold)
        if self.right:
            for elem in self.right:
                yield elem

    def __str__(self):
        return "Node:{}, index={}, threshold={}".format(self.__identifier, 
                  self.index, self.threshold)

    def __sanitize_id(id):
        return id.strip().replace(" ", "")

    def is_leaf(self):
        return (self.left is None) and (self.right is None)

class Visit:
    def __init__(self, node):
        self.node = node 

""" Visitor """
class NodeVisitor(object):
    """ depth first traverse """
    def visit(self, node, callback=None):
        stack = [Visit(node)]
        last_result = None
        while stack:
            try:
                last = stack[-1]
                if isinstance(last, types.GeneratorType):
                    stack.append(last.send(last_result))
                    last_result = None
                elif isinstance(last, Visit):
                    stack.append(self._visit(stack.pop().node))
                else:
                    last_result = stack.pop()
            except StopIteration:
                stack.pop()
        return last_result
    
    def _visit(self, node):
        methname = 'visit_' + type(node).__name__
        meth = getattr(self, methname, None)
        if meth is None:
            meth = self.generic_visit
        return meth(node)

    def visit_Node(self, node):
        print(node)
        if node.left:
            return Visit(node.left)
        if node.right:
            return Visit(node.right)
        return

    def generic_visit(self, node):
        raise RuntimeError('No {} method'.format('visit_' + type(node).__name__))

class VPTree(object):
    """ 
        Implementation of vantage point tree

    reference: 
         This code was adopted with from Steve Hanov's great tutorial 
         at http://stevehanov.ca/blog/index.php?id=130
         and 
         vptree.h is used in t-SNE 
    """

    def __init__(self, items, itemComparator=lambda x,y: x<y, 
                            distance=lambda x, y: np.linalg.norm(x-y),
                            strategy="random"):
        """
        Args: 
            items (list) : a list of objects
            strategy (function) : used to select a node when building up the tree.
            itemComparator (function) : 'random' or 'median', 
                                    compare two objects in the list.
            distance (function) : compute the distance between two objects.
        """
        super(VPTree, self).__init__()
        self._root = None
        self._tau = 0.0
        # used in nth_element when building vp tree. 
        self._itemComparator = itemComparator 
        # used for meauring the distance between any two elements of items
        self._distance = distance 
        # used to select a node when building up the tree
        if strategy == 'random':
            self._strategy = lambda l, u: int(np.random.rand() * (u - l - 1)) + l
        elif strategy == 'median':
            self._strategy = lambda l, u: (u + l) >> 1
        else:
            print("Unknown the strategy for node seletion: {}".format(strategy))
            sys.exit(1)

        # a list of object. and sort items
        self._items = items
        self._items.sort() 
        # Start to build a vptree. zero-based 
        self._root = self.__build_tree(0, len(items)-1) 

    def __eq__(self, other):
        if other is self:
            if isinstance(other, self.__class__):
                return self.__dict__ == other.__dict__
            return NotImplemented
            # more strictly check
            #if type(other) is type(self):
            #    return self.__dict__ == other.__dict__
            #return False
        else:
            return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def __build_tree(self, lower, upper):
        if upper == lower:
            return None
        # lower index is center of current node
        node = Node(index=lower,name=self._items[lower].name, value=self._items[lower].value)
        if upper - lower > 1:
            #if we did not arrive at leaf yet
            #Choose an arbitrary point(vantage point) and move it to the start
            #i = int(np.random.rand() * (upper - lower - 1)) + lower
            i = self._strategy(lower, upper)
            # swap: _items[lower], _items[i]
            self._items[i], self._items[lower] = self._items[lower], self._items[i]
            # Partition around the median distance
            median = (upper + lower) >> 1
            # sort elements in two range:
            #  1. lower+1 and median 
            #  2. median and upper
            #  pivot with lower
            # in C++ std::nth_element(_items.begin() + lower + 1,
            #                 _items.begin() + median,
            #                 _items.begin() + upper,
            #                 DistanceComparator(_items[lower]))
            nth_element(self._items, lower+1, median, upper, comp=self._itemComparator)
            # Threshold of the new node will be the distance to the median
            node.threshold = self._distance(self._items[lower], self._items[median])
            # Recursively build tree
            node.left = self.__build_tree(lower + 1, median)
            node.right = self.__build_tree(median, upper)
        return node
    
    def __search(self, node, target, k, heapQueue):
        """ 
            Helper function that searches the tree
        Args:
            node (class):
            target :
            k :
            heapQueue(list) : a priority heap queue 
        """
        if node is None:
            return 

        # Compute distance between target and current node
        dist = self._distance(self._items[node.index], target)
        # Current node within radius tau
        if (dist < self._tau):
            if heapQueue.size() == k:
                # remove furthest node from result list (if we already have k results)
                heapQueue.popm()
            # add current node to result list
            heapQueue.push((dist, node.index))
            if heapQueue.size() == k:
                # update tau (farthest point in the result list)
                self._tau = heapQueue.top()[0]

        # Return if we arrived at a leaf
        if node.is_leaf():
            return

        if dist < node.threshold:
            # If the target lies within the radius of ball
            if (dist - self._tau) <= node.threshold:
                # if there can still be neighbors inside the ball, 
                # recursively search left child first
                self.__search(node.left, target, k, heapQueue)
            if (dist + self._tau) >= node.threshold:
                # if there can still be neighbors inside the ball, 
                # recursively search right child first
                self.__search(node.right, target, k, heapQueue)
        else:
            # If the target lies outsize the radius of the ball
            if (dist + self._tau) >= node.threshold:
                self.__search(node.right, target, k, heapQueue)
            if (dist - self._tau) <= node.threshold:
                self.__search(node.left, target, k, heapQueue)

    def __nw_format(self, node, format=1):
        """ output format for newick """
        if format == 1:
            return "{}:{}".format(node.name, node.threshold)
        elif format == 7 or format == 8:
            return "{}".format(node.name)
        elif format == 9 or format == 100:
            return ","

    def __toNewickFormat(self, node, format=1):
        """ 
            Return a newick format with internal node name 
        """
        #if (node.left is None) and (node.right is None) :
        if node.is_leaf():
            return self.__nw_format(node, format)

        if (node.left is not None) and (node.right is not None):
            s = ""
            s += "("
            s += self.__toNewickFormat(node.left, format=format)
            s += ","
            s += self.__toNewickFormat(node.right, format=format)
            s += ")"
            s += self.__nw_format(node, format=format)
            return s
        elif (node.left is not None) or (node.right is not None):
            if node.left is not None:
                s = ""
                s += "("
                s += self.__toNewickFormat(node.left, format=format)
                s += ")"
                s += self.__nw_format(node, format=format)
                return s
            elif node.right is not None:
                s = ""
                s += "("
                s += self.__toNewickFormat(node.right, format=format)
                s += ")"
                s += self.__nw_format(node, format=format)
                return s

    def search(self, target, k, results, distances):
        """ 
        Function that uses the tree to find the k nearest neighbors of target
        """
        # Use a priority queue to store intermediate results on
        heap = PriorityQueue()

        # Variable that tracks the distance to the farthest point in our results
        self._tau = sys.float_info.max
        
        # Perform the search
        self._search(self._root, target, k, heap)
        
        # Gather final results
        #results = []
        #distances = []
      
        while not heap.empty():
            # stored in decreased order
            results.append(self._items[heap.top()[1]])
            distances.append(heap.top())
            heap.popm()
        
        # Results are in reverse order
        # results.reverse()
        # distances.reverse()

    def bfs_traverse(self, callback=lambda l,n:print("Level {}: {}".format(l,n))):
        """ 
            Bread first traverse 

        Kwargs:
            callback (function) : hold two parameters in order,
                        1. current level (int)
                        2. current node (Node)

        """
        list = [self._root]
        level = 0;
        while len(list) > 0:
            for n in list:
                callback(level, n)
            list = [n.left for n in list if n.left] + \
                   [n.right for n in list if n.right]
            level += 1

    def __preOrder(self, callback=lambda n:print(n)):
        """
            Pre-order traverse: root, left, right 
        """
        stack = []
        current = self._root
        while True:
            while current is not None:
                callback(current)
                stack.append(current)
                current = current.left
            if not stack:
                return
            current = stack.pop()
            while current.right is None and stack:
                current = stack.pop()
            current = current.right

    def __inOrder(self, callback=lambda n:print(n)):
        """
            In-order traverse: left, root, right 
        """
        current = self._root
        stack = []
        while True:
            while current is not None:
                stack.append(current)
                current = current.left;
            if not stack:
                return
            current = stack.pop()
            callback(current)
            while current.right is None and stack:
                current = stack.pop()
                callback(current)
            current = current.right 

    def __postOrder(self, callback=lambda n:print(n)):
        """
            Post-order traverse: left, right, root 
        """
        stack = []
        current = self._root
        while True:
            while current is not None:
                stack.append((current,False))
                current = current.left;
            if stack:
                current, visited = stack.pop()
            else:
                return 
            while((current.right is None) or (visited is True)) and stack:
                callback(current)
                current, visited = stack.pop()
            else:
                if not stack and visited:
                    callback(current)
                    return
                stack.append((current,True))
                current = current.right

    def dfs_traverse(self, strategy='inorder', callback=lambda n:print(n)):
        """ 
            Depth first traverse 
                1. Pre-order
                2. In-order
                3. Post-order
        Kwargs:
            strategy (str) : 'preoder','inorder' or 'postorder', default is 'inorder'.
            callback (function) : hold one parameter, current visiting node.
        """
        if strategy == 'preoder':
            self.__preOrder(callback=callback)
        elif strategy == 'inorder':
            self.__inOrder(callback=callback)
        elif strategy == 'postorder':
            self.__postOrder(callback=callback)
        else:
            print("Input argument strategy Unknown: {}".format(strategy))
            sys.exit(1)

    def compare(self, other):
        """
            Tree comparison
        """
        # check equality of class/object of inputs
        if not self == other:
            return NotImplemented



    def save(self, format=1):
        """ 
            Save a tree to a specified format 

        .. table::
            ======  ==============================================
            FORMAT  DESCRIPTION
            ======  ==============================================
            0        flexible with support values
            1        flexible with internal node names
            2        all branches + leaf names + internal supports
            3        all branches + all names
            4        leaf branches + leaf names
            5        internal and leaf branches + leaf names
            6        internal branches + leaf names
            7        leaf branches + all names
            8        all names
            9        leaf names
            100      topology only
            ======  ==============================================
            Format 1 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);
            Format 7 = (A,(B,(D,G)E)C);
            Format 8 = (A,(B,(D,G)));
            Format 9 = (,(,(,)));
        """
        FORMAT = [1, 7, 8, 9, 100]
        if not format in FORMAT:
            print("Error: unsurpport specified format for newick format now.")
            sys.exit(1)
        s = "(" + self.__toNewickFormat(self._root, format=format) + ");"
        return s

    def draw(self, treeGraph, width=183, height=None, units='mm', dpi=300):
        try:
            from ete3 import Tree, TreeStyle
        except ImportError:
            raise ImportError('ete3 is a required package')
        # Creates an empty tree
        s = "(" + self.__toNewickFormat(self._root) + ");"
        t = Tree(self.save(), format=1)
        print(t.get_ascii(show_internal=True))
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        ts.show_branch_support = True
        t.render(treeGraph, w=width, h=height, units=units, dpi=dpi, tree_style=ts)


