# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-08-01 23:49:35
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-08-01 23:52:06
from bilab.graph import Graph
from bilab.graph.visitor import bfs_visitor

__all__ = ["connected_components",
           "connected_component_subgraphs"]


def connected_components(G):
    if not isinstance(G, Graph):
        raise TypeError("Inappropriate argument type for bfs_visitor")
    visited = set()
    for v in G:
        if v not in visited:
            c = set(bfs_visitor(G, v))
            yield c
            visited.update(c)


def connected_component_subgraphs(G, copy=True):
    for c in connected_components(G):
        if copy:
            yield G.subgraph(c).copy()
        else:
            yield G.subgraph(c)
