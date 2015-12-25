# -*- coding: utf-8 -*-

import sys
import numpy as np

from bilab.chemicals import ElementData
from bilab.structure.atomic import AtomGroup
from bilab.structure.forcefields import ForceField
from bilab import reduce_newobj


class Universe(object):

    __slots__ = [
        '_masses', '_bonds',
        '_angles', '_dihedrals',
        '_impropers', '_charges']

    """docstring for Universe"""
    def __init__(self, ag, forcefield):
        super(Universe, self).__init__()

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
        self._masses = []
        self._bonds = []
        self._angles = []
        self._dihedrals = []
        self._impropers = []
        self._charges = []

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
            'dihedrals': self._dihedrals,
            'impropers': self._impropers,
            'charges': self._charges
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
    def angles(self):
        """
        take settings from forcefield
        """
        return self._angles

    @property
    def dihedrals(self):
        return self._dihedrals

    @property
    def impropers(self):
        return self._impropers

    @property
    def charges(self):
        """ get charges defined in forcefield """
        return self._charges

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

    def __add__(self, other):
        if isinstance(other, AtomGroup):
            self.ag += other
        elif isinstance(other, Universe):
            self.ag += other.ag
        else:
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

    def description(self):
        """
        Description of the Universe

        lammps:
        LAMMPS Description

                2004  atoms
                1365  bonds
                 786  angles
                 207  dihedrals
                  12  impropers
                  14  atom types
                  18  bond types
                  31  angle types
                  21  dihedral types
                   2  improper types

         36.840194 64.211560 xlo xhi
         41.013691 68.385058 ylo yhi
         29.768095 57.139462 zlo zhi
         """
        desp = """
LAMMPS Description
        {}  atoms
        {}  bonds
        {}  angles
        {}  dihedrals
        {}  impropers
        {}  atom types
        {}  bond types
        {}  angle types
        {}  dihedral types
        {}  improper types

 {} {} xlo xhi
 {} {} ylo yhi
 {} {} zlo zhi
""".format(len(self.ag), len(self._bonds), len(self._angles),
           len(self._dihedrals), len(self._impropers),
           )
