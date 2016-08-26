#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2016-01-01 15:42:39
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-21 02:32:58
from __future__ import division, print_function, absolute_import

import numpy as np

"""
 spline

                           [2, -2, 1,   1][x1, y1, z1]
[x, y, z] =[u^3, u^2, u, 1][-3, 3, -2, -1][x2, y2, z2]
 point       parameter     [0, 0,    1, 0][x1',y1',z1']
  on          vector       [1, 0,    0, 0][x2',y2',z2']
 the spline                    basis        control matrix
"""


def b_splines_3Dbasis(s=1.0/6):
    """ B spline basis """
    bb = [
        [-1, 3, -3, 1],
        [3, -6, 3, 0],
        [-3, 0, 3, 0],
        [1, 4, 1, 0]
    ]
    return np.asmatrix(bb)*s


def cubic_hermite_3Dbasis():
    """ Cubic Hermite spline basis
        [u*u*u, u*u, u, 1] parameter vector
    """
    bb = [
        [2, -2, 1, 1],
        [-3, 3, -2, -1],
        [0, 0, 1, 0],
        [1, 0, 0, 0]
    ]
    return np.asmatrix(bb)


def catmull_rom_3Dbasis(s=0.5):
    """ Return catmull rom basis
    Args:
        s (float): tension, typically 1/2
                   =1
[ 2, -2, 1,  1][ 0, 1, 0, 0]
[-3,  3,-2, -1][ 0, 0, 1, 0]
[ 0,  0, 1,  0][-s, 0, s, 0]
[ 1,  0, 0,  0][ 0,-s, 0, s]
    """
    cr = [
        [-s, 2-s, s-2, s],
        [2*s, s-3, 3-2*s, -s],
        [-s, 0, s, 0],
        [0, 1, 0, 0]
    ]
    return np.asmatrix(cr)


def bezier_3Dbasis():
    """ Also known as Bernstein polynomial
    bb = [ 2, -2,  1,  1] [1, 0, 0, 0]
         [-3,  3, -2, -1] [0, 0, 0, 1] =Bezier basis
         [ 0,  0,  1,  0] [-3,3, 0, 0]
         [ 1,  0,  0,  0] [0, 0,-3, 3]
           Hermite basis   Bezier to Hermite
    """

    bb = [
        [-1, 3, -3, 1],
        [3, -6, 3, 0],
        [-3, 3, 0, 0],
        [1, 0, 0, 0]
    ]
    return np.asmatrix(bb)


def interp4points(x0, basis, nPoints=100):
    """
    Interpolation for 4 points
    """
    x0 = np.asmatrix(x0)
    u_range = np.arange(0.0, 1.01, 1.0/nPoints)  # 1.01 to include the last
    d = np.zeros((len(u_range), 3))
    for i, u in enumerate(u_range):
        coef = np.asmatrix([u*u*u, u*u, u, 1])
        p = coef * basis * x0
        d[i] = p
    return d.tolist()


def interp3d(x0, basis, nPoints=100):
    """ 3D spline fitting

    Args:
        x0: ndarray, control vectors, 4x3
        u: float, [0, 1]
    """
    x0 = np.asmatrix(x0)
    sz = x0.shape[0]
    # The curve C will contain an array of (x,y) points.
    C = []
    extend = C.extend
    # process first segment: P0-P0-P1-P2
    first_seg = x0[[0, 0, 1, 2], :]
    c = interp4points(first_seg, basis, nPoints=nPoints)
    extend(c)
    for i in range(sz-3):
        c = interp4points(x0[i:i+4], basis, nPoints=nPoints)
        extend(c)
    # process last segment: Pn-2, Pn-1, Pn, Pn
    last_seg = x0[[sz-3, sz-2, sz-1, sz-1], :]
    c = interp4points(last_seg, basis, nPoints=nPoints)
    extend(c)
    return np.asmatrix(C)
