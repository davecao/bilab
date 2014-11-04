# -*- coding: utf-8 -*-

from bilab.geometry.tess import *
from bilab.utilities import checkCoords


__all__ = [ 'calcASA' ]

def find_neighbor_indices(atoms, probe, k):
    """
    Args:
        atoms (class): bilab.structure.AtomGroup
        probe (float): the radius of a probe
        k (int): the source atom/point

    Returns:
        list of indices of atoms within probe distance to atom k.

    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + probe + probe
    indices = range(k)
    indices.extend(range(k+1, len(atoms)))
    for i in indices:
        atom_i = atoms[i]
        dist = pos_distance(atom_k.pos, atom_i.pos)
        if dist < radius + atom_i.radius:
            neighbor_indices.append(i)
    return neighbor_indices

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
        neighbor_indices = find_neighbor_indices(atoms, probe, i)
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0
        radius = probe + atom_i.radius

        n_accessible_point = 0
        for point in sphere_points:
            is_accessible = True

            test_point.x = point[0]*radius + atom_i.pos.x
            test_point.y = point[1]*radius + atom_i.pos.y
            test_point.z = point[2]*radius + atom_i.pos.z

            cycled_indices = range(j_closest_neighbor, n_neighbor)
            cycled_indices.extend(range(j_closest_neighbor))

            for j in cycled_indices:
                atom_j = atoms[neighbor_indices[j]]
                r = atom_j.radius + probe
                diff_sq = pos_distance_sq(atom_j.pos, test_point)
                if diff_sq < r*r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1

        area = const*n_accessible_point*radius*radius
        areas.append(area)
    print "%.1f angstrom squared." % sum(areas)
    return areas
