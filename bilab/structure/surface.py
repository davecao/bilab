# -*- coding: utf-8 -*-
import math
import numpy as np

# from bilab import PY3K
# from bilab.structure.atomic import ATOMIC_FIELDS
# from bilab.geometry.distance import euclidean
from bilab.geometry.tess import *
from bilab.utilities import checkCoords
from bilab.chemicals import ElementData
# from bilab.structure.contacts import findNeighbors
from bilab.structure.atomic import AtomGroup
# from timeit import default_timer as timer
from bilab.utilities import profile

# from numbapro import jit, cuda, autojit, vectorize, float32, int8
# from joblib import Parallel, delayed


__all__ = ['calcASA']

# @vectorize("int8(float32[:], float32[:], float32)", target='gpu')
#@jit('void(float32[:,:,:], float32[:,:,:], float32, float32[:])', target='gpu')
def setGPUAccessibility(coords, probe_coords, probe_radius, res):
        """ set accessibility for a probe"""
        # dist = euclidean(probe_coords, self.coords)

        vx = coords[0] - probe_coords[0]
        vy = coords[1] - probe_coords[1]
        vz = coords[2] - probe_coords[2]
        dist = vx * vx + vy * vy + vz * vz
        thre_r2 = probe_radius * probe_radius
        if dist <= thre_r2:
            # self.is_accessible = False
            return 0
        return 1


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
        # dist = euclidean(probe_coords, self.coords)
        vx = self.x - probe_coords[0]
        vy = self.y - probe_coords[1]
        vz = self.z - probe_coords[2]
        dist = vx * vx + vy * vy + vz * vz
        thre_r2 = probe_radius * probe_radius
        if dist <= thre_r2:
            self.is_accessible = False


class NumericSurface(object):
    """
        represent an atom with a numerical sphere
    """
    def __init__(self, center, radius, sphere_points):
        # sphere_points = tesselate_by_sprial(n_sphere_point)
        self.test_points = np.array(sphere_points) * radius + center
        # self.points = []  # test_points
        self.numOfpoints = len(sphere_points)
        self.numOfAccessiblePoints = np.ones(
                (self.numOfpoints, 1),
                dtype=np.int8)
        # self.numOfAccessiblePoints = len(sphere_points)
        # if len(self.points):
        #    self.points = []
        # for p in self.test_points:
        #    self.points.append(Point(p[0], p[1], p[2]))

    def setAccessibility(self, c, r):
        # res = np.ones((self.numOfpoints,1), dtype=np.int8)
        test_inx = np.where(self.numOfAccessiblePoints == 1)[0]
        test_ = self.test_points[test_inx]
        c_ = np.asarray(c, dtype=np.float32)
        nc = np.tile(c_, (len(test_), 1))
        # p = np.asarray(self.points[], dtype=np.float32)
        # setGPUAccessibility(p, nc, r, res)
        res = np.where(np.sum((test_ - nc)**2, axis=1) < r*r)[0]
        # print(res)
        self.numOfAccessiblePoints[test_inx[res]] = 0
        # for inx in np.where(res==0)[0]:
        #    self.points[inx].is_accessible = False
        #for p in self.points:
        #    if p.is_accessible:
        #        p.setAccessibility(c, r)

    def getNumOfAccessiblePoint(self):
        # numOfAccessiblePoints = 0
        # for p in self.points:
        #    if p.is_accessible:
        #        numOfAccessiblePoints += 1
        numOfAccessiblePoints = len(self.numOfAccessiblePoints[
                                    self.numOfAccessiblePoints == 1])
        return numOfAccessiblePoints


def getElementVDWRadius(elem_str):
    """
    Args:
        elem_str (str): a string of the name of an chemical element
    Returns:
        radius (float64) : the vdw radius of the chemical element
    """
    elem = ElementData.get(elem_str.lower())
    return elem.getVanDerWaalsRadius()[0]/100


def find_neighbor_2NLR(atoms):
    pass


def find_neighbor_3NLR(atoms):
    pass


def find_neighbor_4NLR(atoms):
    pass


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
    r_atom_k = getElementVDWRadius(atom_k.getElement())
    serial_no = atom_k.getSerial()
    # radius = getElementVDWRadius(atom_k.getElement()) + probe + probe
    radius_offset = 1.8 * 2
    radius = r_atom_k + probe + probe + radius_offset
    # start =timer()
    sel_center = "within {} of center".format(radius)
    found_neighbors = atoms.select(sel_center, center=atom_k.getCoords())

    if found_neighbors:
        neighbors = found_neighbors.copy()
    else:
        return neighbors
    # end=timer()
    neighbors_ex = neighbors.select("not (serial {})".format(serial_no)).copy()
    # print ("{}".format(end-start))
    if verbose:
        for at in neighbors_ex.iterAtoms():
            print ("{}_{} - Name:{}, SN:{}, Res:{}, Ch:{}".format(
                atom_k, atom_k.getSerial(), at.getName(), at.getSerial(),
                at.getResname(), at.getChid()))
    return neighbors_ex


# @profile
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
                               a unit sphere. 258 ~ 1026 seem best
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

    # test_point = Vector3d()
    # test_point = None
    areas = []
    surfaces = []
    append = areas.append
    append_surf = surfaces.append

    sphere_points = tesselate_by_sprial(n_sphere_point)
    points_density = len(sphere_points)
    const = 4.0 * math.pi / points_density
    for i, atom_i in enumerate(atoms):
        elem_i = atom_i.getElement()
        center = atom_i.getCoords()
        radius = getElementVDWRadius(elem_i)
        # generate surface
        tot_radius = radius + probe
        surface = NumericSurface(center, tot_radius, sphere_points)
        neighbors = find_neighbor_indices(atoms, probe, i, verbose=False)
        # num_neighbors = len(neighbors)

        for atom_j in neighbors.iterAtoms():
            elem_j = atom_j.getElement()
            atom_j_coord = atom_j.getCoords()
            r = getElementVDWRadius(elem_j)
            surface.setAccessibility(atom_j_coord, r + probe)
        # count accessible-points
        n_accessible_point = surface.getNumOfAccessiblePoint()
        area = const*n_accessible_point*tot_radius*tot_radius
        # areas.append(area)
        append(area)
        append_surf(surface)
        if verbose:
            print(
                "{} -- Res:{}_{}_{}, atom: {}, const:{}, #access:{}/{},"
                "radius:{:.3f},"
                "Probe:{}, area:{:.3f}".format(
                    i, atom_i.getResname(),
                    atom_i.getResnum(), atom_i.getChid(),
                    atom_i.getName(), const, n_accessible_point,
                    points_density,
                    probe, radius, area))

    atoms.setData('surface_area', np.array(areas))
    atoms.setData('surface_points', np.array(surfaces))
    return areas
