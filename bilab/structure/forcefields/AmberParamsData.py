# -*- coding: utf-8 -*-
# from __future__ import absolute_import

__all__ = ['AmberAtomType',
           'AmberBondParameters',
           'AmberBondAngleParameters',
           'AmberDihedralParameters',
           'AmberImproperParameters',
           'AmberHbondParameters',
           'AmberLJParameterSet']


class AmberAtomType(object):
    """ Atom Type """
    def __init__(self, name, mass):
        self.name = name
        self.mass = mass
        self.hydrophylic = False


class AmberBondParameters(object):
    """ Bond  """
    def __init__(self, a1, a2, k, l):
        self.a1 = a1
        self.a2 = a2
        self.k = k
        self.l = l


#
# Bond angle parameter class
#
class AmberBondAngleParameters(object):

    def __init__(self, a1, a2, a3, k, angle):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.k = k
        self.angle = angle


#
# Dihedral angle parameter class
#
class AmberDihedralParameters(object):

    def __init__(self, a1, a2, a3, a4, divf, k, delta, n):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.terms = [(k/divf, delta, int(abs(n)))]

    def addTerm(self, divf, k, delta, n):
        self.terms.append((k/divf, delta, int(abs(n))))

    def __repr__(self):
        if self.a1 is None and self.a4 is None:
            return 'X-' + self.a2.name + '-' + self.a3.name + '-X'
        else:
            return self.a1.name + '-' + self.a2.name + '-' + \
                   self.a3.name + '-' + self.a4.name
    __str__ = __repr__


#
# Improper dihedral angle parameter class
#
class AmberImproperParameters(object):

    def __init__(self, a1, a2, a3, a4, divf, k, delta, n):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.terms = [(0.5*k, delta, int(abs(n)))]

    def addTerm(self, divf, k, delta, n):
        self.terms.append((0.5*k, delta, int(abs(n))))

    def __repr__(self):
        if self.a4 is None:
            if self.a3 is None:
                return self.a1.name + '-' + self.a2.name + '-X-X'
            else:
                return self.a1.name + '-' + self.a2.name + '-' + \
                       self.a3.name + '-X'
        else:
            return self.a1.name + '-' + self.a2.name + '-' + \
                   self.a3.name + '-' + self.a4.name
    __str__ = __repr__


#
# H-bond parameter class
#
class AmberHbondParameters(object):

    def __init__(self, a1, a2, a, b):
        self.a1 = a1
        self.a2 = a2
        self.a = a
        self.b = b


#
# Lennard-Jones parameter sets
#
class AmberLJParameterSet(object):

    def __init__(self, name, type):
        self.name = name
        self.type = type
        self.entries = {}

    def addEntry(self, name, p1, p2, p3):
        if self.type == 'SK':
            self.entries[name] = (p1, p2, p3)
        else:
            self.entries[name] = (p1, p2)

    def getEntry(self, name):
        return self.entries[name]
