from __future__ import print_function
# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-26 14:34:10
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-08-08 15:48:11
import uuid
import sys
from six import iteritems
from bilab import string_types
from bilab.graph import VertexProperty
import pprint


class Vertex(VertexProperty):
    # def __new__(cls, *args):
    #     """ create an object """
    #
    def __init__(self, name, vtype="basic", desp=None, species=None):
        super(Vertex, self).__init__()
        self.id = uuid.uuid4()
        self.name = name
        self.desp = desp
        self.species = species
        self.adjacent = {}
        self.vtype = vtype
        self.weight = None  # vertex weight
        self.degree = 0  # for undirected graph
        self.in_degree = 0  # for directed graph
        self.out_degree = 0  # for directed graph

    def __str__(self):
        return "{}".format(self.name)

    def __repr__(self):
        return '{}'.format(self.__class__)

    def __eq__(self, other):
        if isinstance(other, Vertex):
            return self.name == other.name

    def __hash__(self):
        return hash(self.name)
#        # return hash((self.__class__,) + (self.__getinitargs__()))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __iter__(self):
        # loop over adjacent
        for neighbor, cost in iteritems(self.adjacent):
            yield (neighbor, cost)

    def __getitem__(self, neighbor):

        if isinstance(neighbor, Vertex):
            if neighbor in self.adjacent:
                return self.adjacent[neighbor]

        if isinstance(neighbor, string_types):
            for x in self.adjacent:
                if x.name == neighbor:
                    return x
        return None

    def __getstate__(self):
        try:
            # state = vars(self).copy()
            state = self.__dict__
        except TypeError as e:
            print("Vertex: TypeError {}".format(e))
            state = {}
        return state

    def __setstate__(self, state):
        for k in state:
            setattr(self, k, state[k])

    def __reduce__(self):
        from copy_reg import __newobj__
        if hasattr(self, '__getnewargs__'):
            args = self.__getnewargs__()
        else:
            args = (self.name, )

        if hasattr(self, '__getstate__'):
            state = self.__getstate__()
        elif hasattr(type(self), '__slots__'):
            state = self.__dict__,
            {
                k: getattr(self, k) for k in type(self).__slots__
            }
        else:
            state = self.__dict__

        if isinstance(self, list):
            listitems = self
        else:
            listitems = None

        if isinstance(self, dict):
            dictitems = self.items()
        else:
            dictitems = None

        return __newobj__, (type(self),) + args, state, listitems, dictitems

    def __getnewargs__(self):
        """ parameters for creating an object """
        return (Vertex.__str__(self), )

    # def __getinitargs__(self):
    #     """ init args when unpickle an object """
    #     return (self.name, )

    def add_neighbor(self, neighbor, weight=0):
        self.adjacent[neighbor] = weight

    def get_connections(self):
        return self.adjacent.keys()

    def get_num_neighbors(self):
        return len(self.adjacent.keys())

    def get_id(self):
        return self.id

    def get_degree(self):
        # for undirected graph
        return len(self.adjacent)

    def get_weight(self):
        # if None then use degree
        if self.weight is None:
            return self.degree
        else:
            return self.weight

    def get_neighbor_cost(self, n):
        # find cost in adjacency list
        if n in self.adjacent.keys():
            return self.adjacent[n]
        return None

    def remove_neighbor(self, v):
        status = False
        if v in self.adjacent.keys():
            del self.adjacent[v]
            status = True
        return status

    def toGraphviz(self, dest=sys.stdout):
        graphviz_str = '    {} [label="{}",{}];'.format(
                self.name, self.name, self.serialize()
            )
        print(graphviz_str, file=dest)

    def toJson(self, dest=sys.stdout):
        json_str = """
{{
  "classes": "",
  "selected": false,
  "group": "nodes",
  "grabbable": true,
  "locked": false,
  "selectable": true,
  "removed": false,
  "data": {{
      "occ": 10.0,
      "type": {},
      "id": {},
      "weight": 0.02429838415745353,
  }},
  "grabbed": false
}}""".format(self.vtype, self.name)
        print(json_str, file=dest)
