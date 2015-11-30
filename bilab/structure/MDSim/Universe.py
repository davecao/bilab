# -*- coding: utf-8 -*-

import sys
import numpy as np

from bilab.chemicals import ElementData
from bilab.structure.atomic import AtomGroup
from bilab.structure.forcefields import ForceField
from bilab import reduce_newobj


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
            self._setMass()
        else:
            print("TypeError: input should be an object of AtomGroup")
            sys.exit(1)
        if isinstance(forcefield, ForceField):
            self._forcefield = forcefield
        else:
            print("TypeError: input should be an object of ForceField")
            sys.exit(1)
        self._evaluator = {}
        self._atom_properties = {}
        self._objects = []
        self._configuration = None

    @property
    def __deepcopy_keyattrs(self):
        return {
            'masses': self._masses,
            'bonds': self._bonds,
            'angles': self._angles,
            'dihedrals': self.dihedrals,
            'impropers': self.impropers
        }

    @property
    def _setMass(self):
        masses = []
        append = masses.append
        for atom in self._ag.iterAtoms():
            ele = ElementData.get(atom.getName()[0].lower())
            append(float(ele.getAtomicMass()[-1]))
        self._ag.setData('Mass', np.array(masses))

    @property
    def _setAngles(self):
        pass

    @property
    def _setCharges(self):
        """ Set charges defined in forcefield """
        return self.__setCharges

    @property
    def __deepcopy__(self, memo):
        """ Copy object """
        kwds = self.__deepcopy_keyattrs
        args = self.__deepcopy_args
        cls = self.__class__
        result = cls.__new__(cls, *args, **kwds)
        memo[id(self)] = result
        result.__init__(*args, **kwds)
        return result

    def _get_slots(self):
        all_slots = (
            getattr(cls, '__slots__', ()) for cls in self.__class__.__mro__)
        r = set(slot for slots in all_slots for slot in slots)
        return r

    def __getstate__(self):
        # subclass do not must be slots
        try:
            state = vars(self).copy()
        except TypeError:
            state = {}

        for slot in self._get_slots():
            try:
                val = getattr(self, slot)
                state[slot] = val
            except AttributeError:
                pass
        return state

    def __setstate__(self, state):
        for k in state:
            setattr(self, k, state[k])

    def __reduce__(self):
        args = (self._ag, self._ff)
        state = self.__getstate__()

        return reduce_newobj, (type(self),)+args, state, None, None
