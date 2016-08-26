from __future__ import print_function

# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-26 14:33:28
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-07-26 14:38:22
import sys

__all__ = ['toJson', 'toGraphviz']


def toJson(self, dest=sys.stdout):
    #  create a json format
    # json format:
    # {"nodes":[
    #  { "classes": "",
    #    "selected": false,
    #    "group": "nodes",
    #    "grabbable": true,
    #    "locked": false,
    #    "selectable": true,
    #    "removed": false,
    #    "data": {
    #        "occ": 10.0,
    #        "type": "basic",
    #        "id": "GO:0002028",
    #        "weight": 0.02429838415745353,
    #        "name": "regulation of sodium ion transport"
    #     },
    #    "grabbed": false
    #  }, ...
    # ],
    # "edges":[
    #
    # ]}
    nodes = self.get_vertices()
    edges = self.get_edges()
    print("{\"nodes\":[", file=dest)
    for node in nodes:
        node.toJson(dest=dest)
        print(",", file=dest)
    print("],", file=dest)
    print("\"edges:[", file=dest)
    for edge in edges:
        edge.toJson(dest=dest)
        print(",", file=dest)
    print("]}", file=dest)


def toGraphviz(self, dest=sys.stdout):
    #  vertexes = self.get_vertices()
    #  graph definition
    g_desp = """
graph {} {{
    graph [
    {}
    ];
    // node define
""".format(self.label, self.serialize())
    print(g_desp, file=dest)
    #  nodes
    vertexes = self.get_vertices()
    for vertex in vertexes:
        vertex.toGraphviz(dest=dest)
    #  edges
    edges = self.get_edges()
    edge_desp = ""
    for edge in edges:
        edge_desp += "    " + str(edge) + ";\n"
    print(edge_desp, file=dest)
    print("}\n", file=dest)
