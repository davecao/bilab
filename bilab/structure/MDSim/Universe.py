# -*- coding: utf-8 -*-

import sys
import numpy as np

from bilab.chemicals import ElementData
from bilab.structure.atomic import AtomGroup
from bilab.structure.forcefields import ForceField


class Universe(object):

    __slots__ = ['_masses', '_bonds', '_angles', 'dihedrals', 'impropers']

    """docstring for Universe"""
    def __init__(self, ag, forcefield):
        super(Universe, self).__init__()
        # self._MAXCOLVALENT = 2.6
        # self._LEEWAY = 1.1
        if isinstance(ag, AtomGroup):
            self._ag = ag
            self._ag.foundCovalentBonds()
        else:
            print("TypeError: input should be an object of AtomGroup")
            sys.exit(1)
        if isinstance(forcefield, ForceField):
            self._ff = forcefield
        else:
            print("TypeError: input should be an object of ForceField")
            sys.exit(1)

    @property
    def _setMass(self):
        for atom in self._ag.iterAtoms():
            ele = ElementData.get(atom.getName()[0].lower())
        #return np.array([ElementData.])

    def _setAngles(self):
        pass
