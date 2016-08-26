# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-26 14:31:57
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-07-26 15:45:52

__all__ = ['VertexProperty', 'GraphProperty']


class VertexProperty(object):
    """ Define vertex property """
    def __init__(self, **kwargs):
        super(VertexProperty, self).__init__()
        # Schema
        self.schema = kwargs.get('schema', "default")
        # color and stroke
        self.color = kwargs.get('color', "black")
        self.fillcolor = kwargs.get('fillcolor', "white")
        # font
        self.fontcolor = kwargs.get('fontcolor', "black")
        self.fontsize = kwargs.get('fontsize', 10)
        self.fontname = kwargs.get("fontname", "Helvetica")
        # shape
        self.shape = kwargs.get('shape', "circle")
        self.style = kwargs.get('style', "filled")
        # size
        self.width = kwargs.get("width", 1.0)
        self.height = kwargs.get("height", 0.8)

    def __str__(self):
        return self.schema

    def serialize(self):
        return 'color="{}", shape="{}", fillcolor="{}", fontcolor="{}", style="{}"'.format(
            self.color, self.shape, self.fillcolor,
            self.fontcolor, self.style
        )


class EdgeProperty(object):
    """ Define graph properties"""
    def __init__(self, **kwargs):
        super(EdgeProperty, self).__init__()
        # Schema
        self.schema = kwargs.get('schema', "default")


class GraphProperty(object):
    """ Define graph properties"""
    def __init__(self, **kwargs):
        super(GraphProperty, self).__init__()
        # Schema
        self.schema = kwargs.get('schema', "default")
        self.charset = kwargs.get('charset', "UTF-8")
        self.label = kwargs.get('label', "unnamed graph")
        self.labelloc = kwargs.get('labelloc', "b")
        self.labeljust = kwargs.get('labeljust', "r")
        self.bgcolor = kwargs.get('bgcolor', "white")
        self.fontcolor = kwargs.get('fontcolor', "black")
        self.fontsize = kwargs.get('fontsize', 18)
        self.style = kwargs.get('style', "filled")
        self.rankdir = kwargs.get('rankdir', "TB")
        self.margin = kwargs.get('margin', 0.5)
        self.splines = kwargs.get('splines', "spline")
        self.ranksep = kwargs.get('ranksep', 1.0)
        self.nodesep = kwargs.get('nodesep', 0.9)

    def __str__(self):
        return self.schema

    def serialize(self):
        return \
"""  charset = "{}";
    label = "{}",
    labelloc = "{}",
    labeljust = "{}",
    bgcolor = {},
    fontcolor = {},
    fontsize = {},
    style = "{}",
    rankdir = {},
    margin = {},
    splines = {},
    ranksep = {},
    nodesep = {}""".format(
        self.charset, self.label, self.labelloc, self.labeljust, self.bgcolor,
        self.fontcolor, self.fontsize, self.style, self.rankdir, self.margin,
        self.splines, self.ranksep, self.nodesep )
