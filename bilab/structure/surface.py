# -*- coding: utf-8 -*-
import math

from bilab.geometry.distance import euclidean
from bilab.geometry.tess import *
from bilab.utilities import checkCoords
from bilab.chemicals import ElementData
from bilab.structure.contacts import findNeighbors
from bilab.structure.atomic import AtomGroup
from timeit import default_timer as timer

__all__ = [ 'calcASA' ]

class Point(object):
    """ 3D points """
    def __init__(self, x, y, z, **kwargs):
        """Initialization.

        Args:
           x (float): x coordinate
           y (float): y coordinate
           z (float): z coordinate
           is_accessible (bool) : is it accessible by a probe 

        Kwargs:
        """
        super(Point, self).__init__()
        self.x = x
        self.y = y
        self.z = z
        self.coords = [x, y, z]
        self.is_accessible = True

    def setAccessibility(self, probe_coords, probe_radius):
        """ set accessibility for a probe"""
        dist = euclidean(probe_coords, self.coords)
        if dist < probe_radius:
            self.is_accessible = False

class NumericSurface(object):
    """ 
        represent an atom with a numerical sphere 

    """
    def __init__(self, center, radius, n_sphere_point=960):
        sphere_points = tesselate_by_sprial(n_sphere_point)
        test_points = sphere_points * radius + center
        self.points = []
        self.numOfAccessiblePoints = 0
        for p in test_points:
            self.points.append(Point(p[0], p[1], p[2]))

    def setAccessibility(self, c, r):
        for p in self.points:
            p.setAccessibility(c, r)
            if p.is_accessible:
                self.numOfAccessiblePoints += 1

    def getNumOfAccessiblePoint(self):
        return self.numOfAccessiblePoints

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
    neighbors = AtomGroup(title="neighbors")
    atom_k = atoms[k]
    serial_no = atom_k.getSerial()
    radius = getElementVDWRadius(atom_k.getElement()) + probe + probe

    #start =timer()
    sel_center = "within {} of center".format(radius)
    found_neighbors=atoms.select(sel_center, center=atom_k.getCoords())

    if found_neighbors:
        neighbors = found_neighbors.copy()
    else:
        return neighbors
    #end=timer()
    neighbors_ex = neighbors.select("not (serial {})".format(serial_no)).copy()
    #print ("{}".format(end-start))
    if verbose:
        for at in neighbors_ex.iterAtoms():
            print ("{}_{} - Name:{}, sn:{}, Res:{}, ch:{}".format(
                atom_k, atom_k.getSerial(), at.getName(), at.getSerial(),
                at.getResname(), at.getChid()))
    return neighbors_ex

def calcASA(atoms, probe, n_sphere_point=960, verbose=False):
    """
    Calculate accessible surface areas of the atoms, using the probe
    and atom radius to define the surface. (Shrake-Rupley algorithm)
    
    Reference: A. Shrake & J. A. Rupley. "Environment and Exposure to
    Solvent of Protein Atoms. Lysozyme and Insulin." 
    J Mol Biol. 79(1973) 351- 371. 
    
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


    const = 4.0 * math.pi / n_sphere_point
    #test_point = Vector3d()
    test_point = None
    areas = []
    for i, atom_i in enumerate(atoms):
        elem_i = atom_i.getElement()
        center = atom_i.getCoords()
        radius = getElementVDWRadius(elem_i)
        # generate surface
        surface = NumericSurface(center, radius + probe, 
                                 n_sphere_point=n_sphere_point)
        neighbors = find_neighbor_indices(atoms, probe, i, verbose=verbose)
        num_neighbors = len(neighbors)
        
        for atom_j in neighbors.iterAtoms():
            elem_j = atom_j.getElement()
            atom_j_coord = atom_j.getCoords()
            r = getElementVDWRadius(elem_j)
            surface.setAccessibility(atom_j_coord, r + probe)
        # count accessible-points 
        n_accessible_point = surface.getNumOfAccessiblePoint()
#        n_accessible_point = 0
#        area = 0
        #loop over neighbours
#        for atom_j in neighbors.iterAtoms():
#            atom_j_coord = atom_j.getCoords()
#            r = getElementVDWRadius(atom_j.getName()) + probe
#            for inx, test_p in enumerate(test_points):
                
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
    #print "%.1f angstrom squared." % sum(areas)
    return areas
