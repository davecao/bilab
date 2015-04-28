# -*- coding: utf-8 -*-
import math
from bilab.geometry.tess import *
from bilab.utilities import checkCoords
from bilab.chemicals import ElementData
from bilab.structure.contacts import findNeighbors
from timeit import default_timer as timer

__all__ = [ 'calcASA' ]

def pos_distance(p1, p2):
  return math.sqrt(pos_distance_sq(p2, p1))

def pos_distance_sq(p1, p2):
  x = p1.x - p2.x
  y = p1.y - p2.y
  z = p1.z - p2.z
  return x*x + y*y + z*z;

def getElementVDWRadius(elem_str):
    """ 
    Args:
        elem_str (str): a string of the name of an chemical element
    Returns:
        radius (float64) : the vdw radius of the chemical element
    """
    elem = ElementData.get(elem_str.lower())
    return elem.getVanDerWaalsRadius()[0]/100

def find_neighbor_indices(atoms, probe, k, verbose=False):
    """
    Args:
        atoms (class): bilab.structure.AtomGroup
        probe (float): the radius of a probe
        k (int): the source atom/point

    Returns:
        list of indices of atoms within probe distance to atom k.

    Note: AtomGroup.select much faster than findNeighbors
    """
    neighbors = []
    atom_k = atoms[k]
    serial_no = atom_k.getSerial()
    radius = getElementVDWRadius(atom_k.getName()) + probe + probe
    #start =timer()
    #neighbors_pair = findNeighbors(atom_k, radius, atoms2=atoms)
    #end=timer()
    #print ("{}".format(end-start))
    #start =timer()
    neighbors=atoms.select("within {} of center"
        .format(radius, serial_no), center=atom_k.getCoords()).copy()
    #end=timer()
    neighbors_ex = neighbors.select("not (serial {})".format(serial_no)).copy()
    #print ("{}".format(end-start))
    if verbose:
        for at in neighbors_ex.iterAtoms():
            print ("Name:{}, Res:{}, ch:{}".format(at.getName(),
                at.getResname(),at.getChid()))
    return neighbors_ex

def calcASA(atoms, probe, n_sphere_point=960):
    """
    Calculate accessible surface areas of the atoms, using the probe
    and atom radius to define the surface.

    Args:
        atoms (class): bilab.structure.AtomGroup
        probe (float): the radius of a probe

    Kwargs:
        n_sphere_point (int) : the number of point on the surface of
                               a unit sphere
    Returns:
        a tuple of two elements, the first one is summation of ASA and
        the second one is a list of accessible surface areas of the atoms,
        using the probe and atom radius to define the surface.

    """
    # check atoms
    try:
        coords = atoms.getCoords()
    except AttributeError:
        coords = atoms
        try:
            checkCoords(coords, csets=True, dtype=None, name='atoms')
        except TypeError:
            raise TypeError('atoms must be an Atomic instance')

    sphere_points = tesselate_by_sprial(n_sphere_point)

    const = 4.0 * math.pi / len(sphere_points)
    #test_point = Vector3d()
    test_point = None
    areas = []
    for i, atom_i in enumerate(atoms):
        neighbors = find_neighbor_indices(atoms, probe, i)
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0

        radius = probe + getElementVDWRadius(atom_i.getName())
        atom_i_coords = atom_i.getCoords()
        
        test_point = sphere_points * radius + atom_i_coords
        
        n_accessible_point = 0

#        for point in test_points:
#            is_accessible = True
#            test_point.x = point[0]*radius + atom_i.pos.x
#            test_point.y = point[1]*radius + atom_i.pos.y
#            test_point.z = point[2]*radius + atom_i.pos.z
#
#            cycled_indices = range(j_closest_neighbor, n_neighbor)
#            cycled_indices.extend(range(j_closest_neighbor))
#
#            for j in cycled_indices:
#                atom_j = atoms[neighbor_indices[j]]
#                r = atom_j.radius + probe
#                diff_sq = pos_distance_sq(atom_j.pos, test_point)
#                if diff_sq < r*r:
#                    j_closest_neighbor = j
#                    is_accessible = False
#                    break
#            if is_accessible:
#                n_accessible_point += 1

        area = const*n_accessible_point*radius*radius
        areas.append(area)
    print "%.1f angstrom squared." % sum(areas)
    return areas
