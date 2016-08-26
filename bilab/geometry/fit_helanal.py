#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2016-01-12 12:37:56
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-12 14:47:17
# 1. Generate definition of the module
#   f2py-2.7 helanal.f95 -m helanal -h helanal.pyf --overwrite-signature
# 2. Compile
#   f2py-2.7 -c helanal.pyf helanal.f95 -o _helanal.so
#
from __future__ import division, print_function, absolute_import

import numpy as np

try:
    from bilab.geometry._helanal import helanal
except ImportError:
    raise ImportError('helanal module could not be imported. '
                      'Check the installation of bilab. or compile '
                      'helanal.f95 with f2py again to solve the problem.')

def helanal_helix_params(mol):
    """
        Calculate the bending angles along the helix axis using a sliding
        window of 9 residues

    Args:
        mol(class): bilab.structure.
    """
    ca_data = mol.select("name CA")

    # get coordinates
    ca_coords = ca_data.getCoords()

    helanal.fit(ca_coords)
