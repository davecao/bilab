# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-29 12:10:01
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-08-02 00:44:58

from bilab.graph import Graph

__all__ = ["dfs_visitor",
           "bfs_visitor"]


def dfs_visitor(G, source=None):
    """ depth first search - iteative"""
    stack = []
    visited = []
    if source is None:
        stack.append(G[0])
        # visited.append(G[0])
    else:
        stack.append(source)
        # visited.append(source)

    while len(stack):
        v = stack.pop()
        if not (v in visited):
            stack_aux = []
            visited.append(v)
            yield v
            for w in v.get_connections():
                if not (w in visited):
                    stack_aux.append(w)
            while len(stack_aux):
                stack.append(stack_aux.pop())
    # return visited


def bfs_visitor(G, source):
    """ bread first search """
    if not isinstance(G, Graph):
        raise TypeError("Inappropriate argument type for bfs_visitor")
    visited = set()
    # start = G.get_most_connected_vertex()
    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in visited:
                yield v
                visited.add(v)
                nextlevel.update(v.get_connections())
