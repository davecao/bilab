# -*- coding: utf-8 -*-

import sys
from bilab.chemicals import ElementData
from bilab.structure.atom import AtomGroup
from bilab.geometry.distance import euclidean
from bilab.structure.forcefields import ForceField
from bilab.structure import findNeighbors


class Universe(object):

    """docstring for Universe"""
    def __init__(self, ag, forcefield):
        super(Universe, self).__init__()
        self._MAXCOLVALENT = 2.6
        self._LEEWAY = 1.1
        if isinstance(ag, AtomGroup):
            self._ag = ag
        else:
            print("TypeError: input should be an object of AtomGroup")
            sys.exit(1)
        if isinstance(forcefield, ForceField):
            self._ff = forcefield
        else:
            print("TypeError: input should be an object of ForceField")
            sys.exit(1)
        bonds = self._foundCovalentBonds()
        self._ag.setBonds(bonds)

    def _foundCovalentBonds(self):
        """
        find covalent bonds since the connet info in pdb sometimes is not
        reliable.
        """
        LEEWAY = self._LEEWAY
        MAXCOLVALENT = self._MAXCOLVALENT
        ag = self._ag
        # covalent_bonds = {}
        covalent_bonds_inx = []
        append = covalent_bonds_inx.append
        for atom in ag:
            ele = ElementData.get(atom.getName()[0].lower())
            covalent_radius = ele.getCovalentRadius()[0]/100
            # r_ext = (covalent_r + 2) * PointsPerAngstrom * LEEWAY
            r_ext = (covalent_radius + MAXCOLVALENT) * LEEWAY
            neighs = findNeighbors(atom, r_ext, ag)
            for n in neighs:
                n1 = n[1]
                ele1 = ElementData.get(n1.getName()[0].lower())
                covalent_r1 = ele1.getCovalentRadius()[0]/100
                if n1 == atom:
                    continue
                if euclidean(atom.getCoords(), n1.getCoords()) <\
                        (covalent_radius + covalent_r1) * LEEWAY:
                    # covalent_bonds[tuple(sorted([atom, n1]))] = 0
                    append[tuple(atom.getIndex(), n1.getIndex())]
        # return covalent_bonds, covalent_bonds_inx
        return covalent_bonds_inx
