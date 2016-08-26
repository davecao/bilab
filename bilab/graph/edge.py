from __future__ import print_function

# -*- coding: utf-8 -*-
# @Author: Wei Cao
# @Date:   2016-07-26 14:34:45
# @Last Modified by:   Wei Cao
# @Last Modified time: 2016-08-08 15:22:20
import sys
# from copy_reg import __newobj__ as reduce_newobj

__all__ = ["Edge"]


class Edge(object):

    def __init__(self, frm, to, etype=0, cost=1, desp=None):
        super(Edge, self).__init__()
        self.frm = frm
        self.to = to
        self.etype = etype  # 0: undirect edge in a graph 1: directed edge
        self.desp = desp
        self.cost = cost
        self.hashId = self.__hash__()

    def __str__(self):
        return self.frm.name + ' -- ' + self.to.name

    def __repr__(self):
        return "{}".format(self.__class__)

    def __eq__(self, other):
        # in 3.0 __cmp__
        if isinstance(other, Edge):
            return self.hashId == other.hashId

    def __hash__(self):
        return hash(self.frm.name) + hash(self.to.name)

    def __getstate__(self):
        """ state """
        try:
            # state = vars(self).copy()
            state = dict(self.__dict__)
        except TypeError as e:
            print("TypeError {}".format(e))
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
            args = (self.frm, self.to, )

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

    def toJson(self, dest=sys.stdout):
        json_str = """
{{
    "classes": "multi-unbundled-bezier",
    "selected": false,
    "locked": false,
    "grabbable": true,
    "grabbed": false,
    "selectable": true,
    "removed": false,
    "data": {{
        "networkId": {},
        "source": {},
        "group": "edges",
        "target": {},
        "weight": {},
        "networkGroupId": {},
        "intn": true,
        "id": {},
        "rIntnId": 15
    }}
}}""".format(self.frm.name, self.frm.name, self.to.name,
             self.to.name, self.to.name, self.cost)
        print(json_str, file=dest)
