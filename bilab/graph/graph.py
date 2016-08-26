from __future__ import print_function

# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-26 14:35:11
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-08-09 14:49:07
import types
import sys
from copy import deepcopy
# from copy_reg import __newobj__ as reduce_newobj
# from six import iteritems
from bilab import string_types, get_timestamp
from bilab.graph import GraphProperty, Vertex, Edge


__all__ = ['Graph']


class Graph(GraphProperty):

    def __init__(self, name="unamed", directed=0, func=None):
        super(Graph, self).__init__()
        self.vertices_d = {}
        self.edges_l = []
        # 0: undirected, 1: directed not implemented yet.
        self.directed = directed
        self.num_vertices = 0
        self.num_edges = 0
        self.label = name
        self.execute = func
        if func is not None:
            self.execute = types.MethodType(func, self)

    @property
    def __deepcopy_keyattrs(self):
        return {
            'vertices_d': self.vertices_d,
            'edges_l': self.edges_l,
            'directed': self.directed,
            'num_vertices': self.num_vertices,
            'num_edges': self.num_edges,
            'label': self.label
        }

    def __deepcopy__(self, memo):
        """ Copy object """
        kwds = self.__deepcopy_keyattrs
        cls = self.__class__
        result = cls.__new__(cls, **kwds)
        memo[id(self)] = result
        result.__init__(**kwds)
        return result

    def __str__(self):
        """ informal of class """
        return 'Graph label: {}'.format(self.label)

    def __repr__(self):
        return '{}'.format(self.__class__)

    def __getstate__(self):
        """ state """
        try:
            # state = vars(self).copy()
            state = dict(self.__dict__)
            del state["execute"]
        except TypeError as e:
            print("TypeError {}".format(e))
            state = {}

        return state

    def __setstate__(self, state):
        for k in state:
            setattr(self, k, state[k])

    def __iter__(self):
        return iter(self.vertices_d.values())

    def __contains__(self, v):
        # the vertex object in or not
        return v in self.vertices_d

    def __getitem__(self, k):
        """ return a vertex """
        if isinstance(k, Vertex) and k in self.vertices_d:
            return self.vertices_d[k]

        if isinstance(k, string_types):
            for x in self.vertices_d:
                if x.name == k:
                    return x

        if isinstance(k, (int, long)):
            return self.vertices_d[self.vertices_d.keys()[0]]

    def __reduce__(self):
        from copy_reg import __newobj__

        if hasattr(self, '__getnewargs__'):
            args = self.__getnewargs__()
        else:
            args = (self.label, )

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

    def add_vertex(self, node):
        self.num_vertices += 1
        # new_vertex = Vertex(node)
        # self.vertices_d[node] = new_vertex
        if not isinstance(node, Vertex):
            raise TypeError("Inappropriate argument type for add_vertex")
        if node not in self.vertices_d:
            self.vertices_d[node] = node
        # return new_vertex

    def get_vertex(self, n):
        if n in self.vertices_d:
            return self.vertices_d[n]
        else:
            return None

    def add_edge(self, frm, to, cost=1, nonloop=True):
        if frm not in self.vertices_d:
            self.add_vertex(frm)
        if to not in self.vertices_d:
            self.add_vertex(to)
        if nonloop and frm == to:
            return
        self.vertices_d[frm].add_neighbor(self.vertices_d[to], cost)
        self.vertices_d[to].add_neighbor(self.vertices_d[frm], cost)
        edge = Edge(frm, to)
        if edge not in self.edges_l:
            self.edges_l.append(edge)
            self.num_edges += 1

    def get_vertices(self):
        return self.vertices_d.values()

    def get_edges(self):
        #  convert adjacency list to edge list
        return self.edges_l

    def get_edge(self, frm, to):
        # find the cost
        hashId = hash(frm.name) + hash(to.name)
        for inx, e in enumerate(self.edges_l):
            if hashId == e.hashId:
                return e
        return None

    def get_num_vertices(self):
        return self.num_vertices

    def get_num_edges(self):
        return self.num_edges

    def get_most_connected_vertex(self):
        return max(self.vertices_d.values(), key=lambda x: x.get_degree())

    def check_integrity(self):
        """ Check edges in edges_l vs adjacency """
        for v in self.vertices_d:
            for n in v.get_connections():
                found = self.get_edge(v, n)
                if found is None:
                    print("Not found {} - {}".format(v.name, n.name))
        for e in self.edges_l:
            f1 = self.vertices_d[e.frm][e.to]
            f2 = self.vertices_d[e.to][e.frm]
            if not(f1 and f2):
                print("Not found in vertices list: {} - {}".format(
                    e.frm.name,
                    e.to.name))

    def remove_vertices(self, vlist):
        # dict comprehesion
        # new_vertices_d = {k: v for k, v
        #                   in iteritems(self.vertices_d)
        #                   if k in vlist}
        # self.vertices_d = new_vertices_d
        for v in vlist:
            # remove v in the root
            if v in self.vertices_d.keys():
                for n in v.get_connections():
                    self.remove_edge(v, n)
                del self.vertices_d[v]
                self.num_vertices -= 1
        # remove edges that connected to vertices
        for v in self.vertices_d:
            for n in v.get_connections():
                if n in vlist:
                    v.remove_neighbor(n)
                    self.remove_edge(v, n)

    def remove_edge(self, frm, to):
        found = False
        # if frm not in self.vertices_d:
        #     print("Vertex:{} - not found".format(frm.name))
        # if to not in self.vertices_d:
        #     print("Vertex:{} - not found".format(to.name))
        hashId = hash(frm.name) + hash(to.name)
        for inx, e in enumerate(self.edges_l):
            if hashId == e.hashId:
                found = True
                del self.edges_l[inx]
                self.num_edges -= 1
        return found

    def copy(self):
        return deepcopy(self)

    def subgraph(self, vlist, name="unamed", directed=0, func=None):
        """ given a list of nodes, return a subgraph """
        # create an empty graph
        sgraph = self.__class__(name=name, directed=directed,
                                func=func)
        # add vertices and edges
        for v in vlist:
            # sgraph.vertices_d[v] = v
            # edges
            for neighbor, cost in v:
                    sgraph.add_edge(v, neighbor, cost=cost)
        return sgraph

    def toMetisCSR(self, v_weight=False, e_weight=False):
        """ prepare metis's arguments
        xadj : vertices
        adjncy: adjacency list in one dimension
        e.g. adjacency list =
          [                                xadj[0] = 0
           [1, 2, 3, 4],        vertex 0 - xadj[1] = 4
           [0],                 vertex 1 - xadj[2] = 5
           [0],                 vertex 2 - xadj[3] = 6
           [0],                 vertex 3 - xadj[4] = 7
           [0, 5],              vertex 4 - xadj[5] = 9
           [4, 6],              vertex 5 - xadj[6] = 11
           [13, 5, 7],          vertex 6 - xadj[7] = 14
           [8, 6],              vertex 7 - xadj[8] = 16
           [9, 10, 11, 12, 7],  vertex 8 - xadj[9] = 21
           [8],                 vertex 9 - xadj[10] = 22
           [8],                 vertex 10 - xadj[11] = 23
           [8],                 vertex 11 - xadj[12] = 24
           [8],                 vertex 12 - xadj[13] = 25
           [14, 6],             vertex 13 - xadj[14] = 27
           [13, 15],            vertex 14 - xadj[15] = 29
           [16, 17, 18, 14],    vertex 15 - xadj[16] = 33
           [15],                vertex 16 - xadj[17] = 34
           [15],                vertex 17 - xadj[18] = 35
           [15]                 vertex 18 - xadj[19] = 36
          ]

        metis:
        xadj = [0,4,5,6,7,9,11,14,16,21,22,23,24,25,27,29,33,34,35,36]
        adjncy = [
            1,2,3,4,
            0,
            0,
            0,
            0,5,
            4,6,
            13,5,7,
            8,6,
            9,10,11,12,7,
            8,
            8,
            8,
            8,
            14,6,
            13,15,
            16,17,18,14,
            15,
            15,
            15]
        """
        xadj_inx = {}

        xadj = [0]
        xadj_append = xadj.append

        vwgt = []
        vwgt_append = vwgt.append

        adjncy = []
        adjncy_append = adjncy.append

        adjwgt = []
        adjwgt_append = adjwgt.append

        i = 1
        for v in self.vertices_d:
            if v_weight:
                vwgt_append(v.get_weight())
            xadj_inx[v.name] = (i, v.get_num_neighbors())
            i += 1
        # k: v.name
        # val: (number, num_of_neighbors) a tuple
        sorted_inx = sorted([(val, k) for (k, val) in xadj_inx.items()])
        neighbor_counter = 0
        for val, n_inx in sorted_inx:
            neighbors = []
            frm = self[n_inx]
            for n in frm.get_connections():
                if self.get_edge(frm, n):
                    neighbors.append(xadj_inx[n.name][0])
                    neighbor_counter += 1
                else:
                    print("Edge: {} - {} not found.".format(frm.name, n.name))
                if e_weight:
                    adjwgt_append(frm.get_neighbor_cost(n))
            adjncy_append(neighbors)
            xadj_append(neighbor_counter)

#        neighbor_counter = 0
#        for v in self.vertices_d:
#            neighbors = []
#            for n in v.get_connections():
#                if self.get_edge(v, n):
#                    neighbors.append(xadj_inx[n.name][0])
#                    neighbor_counter += 1
#                else:
#                    print("Edge: {} - {} not found.".format(v.name, n.name))
#                if e_weight:
#                    adjwgt_append(v.get_neighbor_cost(n))
#            adjncy_extend(neighbors)
#            xadj_append(neighbor_counter)
        metis_args = (self.num_vertices, self.num_edges,
                      xadj_inx, sorted_inx, xadj, adjncy)
        if v_weight:
            metis_args += (vwgt,)
        else:
            metis_args += (None,)
        if e_weight:
            metis_args += (adjwgt,)
        else:
            metis_args += (None,)

        return metis_args

    def save2metis(self, v_weight=False, e_weight=False,
                   message=None, srcfile=None, file=sys.stdout):
        """ Save graph to metis file format """
        timestamp = get_timestamp()
        metis_file_desp = ''
        if not (srcfile is None):
            metis_file_desp += "% {}\n".format(srcfile)
        metis_file_desp += "% Metis format created on {}\n".format(timestamp)
        num_vertices, num_edges, xadj_inx, sorted_inx, xadj, adjncy,\
            vwgt, adjwgt = self.toMetisCSR(
                v_weight=v_weight, e_weight=e_weight)

        if not(message is None):
            metis_file_desp += "% {}\n".format(message)
        for vname, key in sorted_inx:
            metis_file_desp += "% {} : {} {}\n".format(key, vname[0], vname[1])
        g_metis = ""
        for inx, v in enumerate(adjncy):
            n = "  ".join(map(str, v))
            # neighbors
            v_neighbor_list = []
            for v_neighbor in v:
                v_neighbor_list.append(sorted_inx[v_neighbor - 1][1])
            # vertex
            test_num = len(v_neighbor_list) == sorted_inx[inx][0][1]
            metis_file_desp += "% vertex {}: {} {} - {}\n".format(
                sorted_inx[inx][1],
                sorted_inx[inx][0][0],
                sorted_inx[inx][0][1],
                test_num)
            metis_file_desp += "% neighbors {}\n".format(
                ",".join(v_neighbor_list))
            g_metis += " {}\n".format(n)
        metis_file_desp += "% #vertices # edges\n"
        metis_file_desp += " {} {}\n".format(num_vertices, num_edges)
        metis_file_desp += "%\n"
        metis_file_desp += \
            "% Here begins the list of vertex neighbors for each vertex\n"
        metis_file_desp += "%\n"
        metis_file_desp += g_metis
        print("{}".format(metis_file_desp), file=file)
