#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2015-12-08 14:20:11
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-02 01:04:51

# In scipy v0.15.1 suggested to use basinhopping other than annealing
# ------------------
# scipy.optimize.anneal(*args, **kwds)[source]
# anneal is deprecated! Deprecated in scipy 0.14.0, use basinhopping instead
#
# Minimize a function using simulated annealing.
#
# Uses simulated annealing, a random algorithm that uses no derivative
# information from the function being optimized. Other names for this family
# of approaches include: “Monte Carlo”, “Metropolis”, “Metropolis-Hastings”,
# etc. They all involve
#   (a) evaluating the objective function on a random set of points.
#   (b) keeping those that pass their randomized evaluation critera.
#   (c) cooling (i.e., tightening) the evaluation critera.
#   (d) repeating until their termination critera are met.
# In practice they have been used mainly in discrete rather than in continuous
# optimization.
# Available annealing schedules are ‘fast', ‘cauchy' and ‘boltzmann'.

import sys
import numpy as np
import scipy
import bilab
from bilab.structure.atomic import AtomGroup, ATOMIC_FIELDS
# from bilab.structure.measure import calcCenter
from bilab.structure import moveAtoms, calcCenter, Transformation
from bilab.geometry import unit
from bilab.geometry.distance import euclidean
from bilab.geometry.transform.rotation_matrix import euler_matrix

mem_pot_data = np.loadtxt(
        "{}/data/mem_kb_potential.dat".format(bilab.data),
        delimiter=',')
betamem_pot_data = np.loadtxt(
        "{}/data/betamem_kb_potential.dat".format(bilab.data),
        delimiter=',')

aa_order = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
aa_order_map = dict(zip(aa_order, range(len(aa_order))))


class SearchBounds(object):
    def __init__(self,
                 xmax=[2*np.pi, 2*np.pi, 15.0],
                 xmin=[0.0, 0.0, -15.0],
                 steps=[1, 1, 1]):
        """
            xmax = [alpha, beta, ztrans]
            xmin = [alpha, beta, ztrans]
        """
        # should check the length - !!!
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
        self.nparams = self.xmax.shape[0]
        self.steps = steps
        self.current = []
        self.grids = self._generate_grids()

    def __call__(self, **kwargs):
        """ Called by basinhopping
            to determine the run is accepted or rejected
            if True, accepted, otherwise rejected.
        """
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin

    def _generate_grids(self):
        grids = [np.arange(s, e, t) for s, e, t in
                 zip(self.xmin, self.xmax, self.steps)]
        return grids

    def get_random_params(self):
        """Return random selected params from xmax"""
        params = [np.random.choice(self.grids[i], 1)[0] for i in
                  range(self.nparams)]
        return params


def generate_dum_atoms(best_e, best_c,
                       xyplane=None,
                       offset=8):
    """ Create dummy atoms for the membrane

    Args:
        natoms (int): the number of dummy atoms for each layer
        xyplane(array): [xmin, xmax, ymin, ymax]
        best_e (float): the shift value along z-axis for extracellular side
        best_c (float): the shift value along z-axis for cytoplasma side
        offset (integer): offset value to dispose the boundary of xy-plane
        polar (bool):
    """
    if len(xyplane) != 4:
        print("xyplane parameter is incorrect. "
              "It should be [xmin, xmax, ymin, ymax]")
        sys.exit(1)

    dum_ag = AtomGroup(title="DUM O upper membrane layer")
    # alternative array [-1, 1, -1, 1]
    offval = np.empty((len(xyplane),))
    offval[::2] = -1
    offval[1::2] = 1
    # xyplane for the dummy atoms
    # [-1, 1, -1, 1] * offset =
    #    [xmin-offset, xmax+offset, ymin-offset, ymax+offset]
    xyplane = np.asarray(xyplane)
    xyplane = xyplane + offval*offset

    # generate coords
    # np.ceil(np.sqrt(natoms))
    num_x = np.int(np.ceil((xyplane[1] - xyplane[0])/2.0))
    num_y = np.int(np.ceil((xyplane[3] - xyplane[2])/2.0))

    # x = np.linspace(xyplane[0], xyplane[1], num=num_x, endpoint=True)
    # y = np.linspace(xyplane[2], xyplane[3], num=num_y, endpoint=True)
    # m_range = np.around((xyplane[2:] - xyplane[:1])/2.0)
    # coords_origin_O = np.append(xyplane[[0, 2]], 15.00+best_e)  # O atom
    # coords_origin_N = np.append(xyplane[[0, 2]], -15.00-best_c)  # N atom
    # print("xy-plane({}):[{},{},{},{}], num={}".format(
    #      natoms, xyplane[0], xyplane[1], xyplane[2], xyplane[3], num))
    x = range(num_x)
    # num_y = len(y)  # range(1, m_range[1]+1)
    y = range(num_y)
    natoms = num_x * num_y

    # Prepare atoms
    atomnames_O = np.zeros(natoms, dtype=ATOMIC_FIELDS['name'].dtype)
    atomnames_N = np.zeros(natoms, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(natoms, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(natoms, dtype=ATOMIC_FIELDS['resnum'].dtype)
    termini = np.zeros(natoms, dtype=bool)
    altlocs = np.zeros(natoms, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(natoms, dtype=ATOMIC_FIELDS['icode'].dtype)
    # coordinates = np.zeros((natoms, 3), dtype=float)
    chainids = np.zeros(natoms, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(natoms, dtype=bool)
    serials = np.zeros(natoms, dtype=ATOMIC_FIELDS['serial'].dtype)
    bfactors = np.zeros(natoms, dtype=ATOMIC_FIELDS['beta'].dtype)
    occupancies = np.zeros(natoms, ATOMIC_FIELDS['occupancy'].dtype)
    segnames = np.zeros(natoms, dtype=ATOMIC_FIELDS['segment'].dtype)
    elements = np.zeros(natoms, dtype=ATOMIC_FIELDS['element'].dtype)

    coordinates = np.vstack((
                    np.repeat(x, num_y)*2,
                    np.tile(y, num_x)*2,
                    np.ones(natoms))).T
    # set z-slice
    coordinates = np.append(xyplane[[0, 2]], 1) + coordinates
    # half_atoms = natoms/2
    # if polar:
    #     coordinates[:half_atoms, :] = 24.00+best_e
    #     coordinates[half_atoms:, :] = -24.00-best_c
    # else:
    coordinates[:, 2] = 15.00+best_e
    # coordinates[half_atoms:, 2] = -15.00-best_c
    # print("xy-plane:[{}]".format(",".join(
    #     ["{:.4f}".format(s) for s in xyplane])))
    # print("x:{} y:{} count:{}".format(num_x, num_y, natoms))
    # print(coordinates.shape)
    # print(atomnames_O.shape)
    # print(resnames.shape)
    # print(hetero.shape)
    # print(serials.shape)
    # print(occupancies.shape)
    # set atom names
    atomnames_O[:] = "O"

    # set residue names
    resnames[:] = "DUM"
    # set residue num
    # resnums[:] = range(start)
    # set hetero flags
    hetero[:] = True
    # set serials
    serials[:] = range(1000, natoms+1000)
    occupancies[:] = 1.0
    bfactors[:] = 0.0
    dum_ag.setCoords(coordinates)
    dum_ag.setNames(atomnames_O)
    dum_ag.setResnames(resnames)
    dum_ag.setResnums(resnums)
    dum_ag.setFlags('hetatm', hetero)
    dum_ag.setSerials(serials)
    dum_ag.setOccupancies(occupancies)
    dum_ag.setFlags('pdbter', termini)
    dum_ag.setAltlocs(altlocs)
    dum_ag.setIcodes(np.char.strip(icodes))
    dum_ag.setChids(chainids)
    dum_ag.setBetas(bfactors)
    dum_ag.setSegnames(np.char.strip(segnames))
    dum_ag.setElements(np.char.strip(elements))
    # Create N dum
    dum_ag_N = dum_ag.copy()
    dum_ag_N.setTitle("DUM N lower layer")
    atomnames_N[:] = "N"
    coordinates[:, 2] = -30.00 - best_e - best_c
    serials[:] += 1000
    dum_ag_N.setNames(atomnames_N)
    dum_ag_N.setCoords(coordinates)
    dum_ag_N.setSerials(serials)
    return dum_ag + dum_ag_N


def get_slice_index(z, slice_start=-24.0, slice_stop=24.0, interval=1.5):
    """
    48  Angstrom thick membrane
    32 * 1.5 Angstrom thick slices within the membrane (-24 to 24)
    2 slices accounting for cytoplasm/lumen (< -24, > 24)
    34 slices
    """
    offset = 0.5
    if z < slice_start:
        return 0
    if z >= slice_stop:
        return 33
    slices = np.arange(slice_start, slice_stop+offset, interval)
    inx = np.searchsorted(slices, z, side='right')
    return inx


def origin_shift(ag, center='geometric'):
    """
    Shift the atoms of atomgroup to the origin

    Arg:
        ag (object):
            the object of bilab.structure.atomic.atomgroup
        center (str):
            value is 'geometric' or 'mass'
            if mass is selected, the atomic mass info should be stored in
            the object of class atomgroup preliminarily.
    """
    if center == 'mass':
        wt = ag.getMasses()
    elif center == 'geometric':
        wt = None
    # ct = calcCenter(ag, weights=wt)
    moveAtoms(ag, to=np.zeros(3), weights=wt)


def rotation_matrix(angle, direction):
    """ Generate rotation matrix

    Args:
        angle (float): in radian
        direction (array):
            x-axis: [1, 0, 0]
            y-axis: [0, 1, 0]
            z-axis: [0, 0, 1]
    """
    sina = np.sin(angle)
    cosa = np.cos(angle)
    direction = unit(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[0.0,         -direction[2],  direction[1]],
                   [direction[2], 0.0,          -direction[0]],
                   [-direction[1], direction[0],  0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    return M


def generate_trans(alpha, beta, gamma, delta_x, delta_y, delta_z, axes='sxyz'):
    """
      Rotate X-axis --> Rotate Y-axis --> Z trans

    Args:
        alpha (float): radian, rotate x-axis
                        range in [0, 2*pi]
        beta  (float): radian, rotate y-axis
                        range in [0, 2*pi]
        z_trans (float): translate on z-axis
    """
    # gamma = 0.0
    # sina = np.sin(alpha)
    # sinb = np.sin(beta)
    # sinr = np.sin(gamma)
    # cosa = np.cos(alpha)
    # cosb = np.cos(beta)
    # cosr = np.cos(gamma)
    # Rotation Z -> Y -> X
    # Rot = [
    #     [cosb*cosr, cosr*sina*sinb-cosb*sinr, cosa*cosr*sinb+sina*sinr],
    #     [cosbsinr, cosa*cosr + sina*sinb*sinr,
    #      -cosr*sina+cosa*sinb*sinr],
    #     [-sinb, cosb*sina, cosa*cosb]]
    # Rotation X -> Y -> Z
    # Rot = [
    #     [cosb*cosr, -cosb*sinr, sinb],
    #     [cosa*sinr+sina*sinb*cosr, cosa*cosr-sina*sinb*sinr, -sina*cosb],
    #     [sina*sinr-cosa*sinb*cosr, sina*cosr+cosa*sinb*sinr, cosa*cosb]
    # ]
    # Rot = np.asarray(Rot)
    quaternions = euler_matrix(alpha, beta, gamma, axes=axes)
    # Trans = np.array([0, 0, delta_z, 1])
    quaternions[0, -1] = delta_x
    quaternions[1, -1] = delta_y
    quaternions[2, -1] = delta_z
    transformation = Transformation(quaternions)
    return transformation


class TakeStep(object):
    """ Control the behavier of stepsize of basinhopping
    see class RandomDisplacement of basinhopping.
    """
    def __init__(self, stepsize=0.5, distfun="uniform"):
        """
        Args:
            stepsize (float):
            distfun (str): uniform or gaussian
        """
        self.f = {"uniform": self._uniform,
                  "gaussian": self._gaussian}
        self.stepsize = stepsize
        if distfun not in self.f.keys():
            print("ValueError: wrong argument for parameter distfun, uniform"
                  " or gaussian")
            sys.exit(1)
        self.generator = self.f[distfun]

    def _uniform(self, x):
        """ Uniform distribution
        A random number is selected between -displacement and + displacement
        This is similar from the example, but the last coordinate to
        take larger steps then the rest of the coordinates.

        """
        s = self.stepsize
        # head is bigger
        # x[0] += np.random.uniform(-2.*s, 2.*s)
        # x[1:] += np.random.uniform(-s, s, x[1:].shape)
        # tail is bigger
        # x[:-1] += np.random.uniform(-s, s, x[1:].shape)
        # x[-1] += np.random.uniform(-2.*s, 2.*s)
        # all is the same
        x = np.random.uniform(-s, s, np.shape(x))
        return x

    def _gaussian(self, x, mu=0.0, sigma=1.0):
        """ Gaussian distribution
        A random number is selected between -displacement and + displacement
        """
        s = self.stepsize
        perturb = np.random.normal(mu, sigma, x.shape)
        n_perturb = (perturb - np.amin(perturb)) / \
                    (np.amax(perturb) - np.amin(perturb)) - s
        x += n_perturb
        # x[-1] *= 2
        # x[0] += np.random.uniform(-2.*s, 2.*s)
        # x[1:] += np.random.uniform(-s, s, x[1:].shape)
        return x

    def __call__(self, x):
        return self.generator(x)


def local_minimizer(fun, x0, args, **kwargs):
    """
    Local minimizer - extend the local minimizer for
                      scipy.optimize.basinhopping
    """
    opt = scipy.optimize.OptimizeResult(
            x=x0,
            fun=fun(x0),
            success=True,
            nfev=1)
    return opt


class Orient(object):
    """ Rotate the atoms

    Args:
        ag (object):
            the object of bilab.structure.atomic.atomgroup
        rotMat (matrix):
            rotation matrix
    Return:
        statistical energy
    """
    def __init__(self, ag, bounds, mem_pot=mem_pot_data,
                 beta_pot=betamem_pot_data,
                 origin_type='geometric',
                 forcespan=False,
                 verbose=False):
        """
        Class for problem setting

        Args:
            ag:
            mem_pot:
            beta_pot:
            origin_type: 'geometric' or 'mass'
            forcespan (bool):
        Return:

        """
        super(Orient, self).__init__()
        if not isinstance(ag, AtomGroup):
            raise ValueError("Argument error: AtomGroup should be"
                             " an object of "
                             "bilab.structure.atomic.AtomGroup")
        self.ag = ag.copy()
        self.mem_pot = mem_pot
        self.beta_pot = beta_pot
        self.bounds = bounds
        self.forcespan = forcespan
        self.nslices = mem_pot.shape[1]
        self.verbose = verbose
        # 34 layers
        # preallocate arrays to store energies for each layer
        self.energies = np.zeros(self.nslices, dtype=np.float64, order='C')
        # save total energy, summation of all layers
        self.total_energy = np.Inf
        # move to origin
        # origin_shift(self.ag, center=origin_type)

    def __sub__(self, other):
        """ Define substract operation"""
        if isinstance(other, Orient):
            return self.total_energy - other.total_energy
        return NotImplemented

    def __call__(self, params):
        # if not isinstance(transformation, Transformation):
        #    raise ValueError("Argument error: transformation should be"
        #                     " an object of "
        #                     "bilab.structure.transform.Transformation")

        top1 = 0
        top2 = 0
        # if not len(params):
        #    params = self.bounds.get_random_params()
        energies = np.zeros(self.nslices, dtype=np.float64, order='C')
        # energies = np.zeros(self.nslices, dtype=np.float64, order='C')
        alpha, beta, ztrans = params
        gamma = 0.0
        delta_x = 0.0
        delta_y = 0.0
        transformation = generate_trans(alpha, beta, gamma,
                                        delta_x, delta_y, ztrans)
        # rotate all atoms
        ag_copy = self.ag.copy()
        atoms = transformation.apply(ag_copy)
        # loop over all atoms
        for at in atoms:
            x, y, z = at.getCoords()
            if z > 17.5:
                top1 = 1
            elif z < -17.5:
                top2 = 1
            resn = at.getResname()
            inx = get_slice_index(z, slice_start=-24.0,
                                  slice_stop=24.0, interval=1.5)
            energies[inx] += self.mem_pot[aa_order_map[resn]][inx]
        # get summation of energies of each layer
        total_energy = np.sum(energies)
        if self.forcespan and ((not top1) or (not top2)):
            total_energy += 10000
        # if total_energy < self.total_energy:
        #    self.total_energy = total_energy
        #    self.energies = energies
        # if self.verbose:
        #    print("total = {}, 1-34 = {}".format(
        #          total_energy, ",".join(str(x) for x in energies.tolist())))
        # should return total energy and its gradient.
        # Not explicit formula to get gradient
        # next_params = self.bounds.get_random_params()
        return total_energy  # , next_params


def print_cb(x, f, accepted):
    test = "rejected"
    if accepted:
        test = "accepted"
        print("alpha={:.4f} beta={:.4f} ztrans={:.4f} f(x0) = {:.4f} {}"
              .format(x[0], x[1], x[2], f, test))


def print_optResults(opt, meth):
    """Print OptimizeResult"""
    x = opt['x']
    print("Summary of optimization by basinhopping")
    print("Local minimizer:{}".format(meth))
    print("Iteration:{}".format(opt['nit']))
    print("Minimization failure:{}".format(opt['minimization_failures']))
    print("Number of evaluations of the objective functions:{}"
          .format(opt['nfev']))
    print(
        "global minimum: {} at alpha={:.4f} beta={:.4f} ztrans={:.4f}"
        .format(opt['fun'], x[0], x[1], x[2]))
    print(
        "Degree: alpha={:.4f} beta={:.4f} ztrans={:.4f}"
        .format(x[0]*180/np.pi, x[1]*180/np.pi, x[2]))


def optimize_lipid_head(ag, mem_pot=mem_pot_data,
                        beta_pot=betamem_pot_data,
                        stepsize=0.25,
                        betaBarrel=False):
    """
        Optimize lipid head positions
    Args:
        beta (bool): atomgroup of a beta-barrel protein, default is False.
        lipid_head (array): an array of length 2,
        core_region(array): define the regions along z axis
    """
    # select representative points from ag

    # Define lipid head region
    lipid_top_head = [10.0, 20.0]
    lipid_bot_head = [-20.0, -10.0]
    # Define the core regrion of membrane
    core_region = [-7.5, 10.0]
    core_region2 = [-7.5, 10.0]
    stepsize = 0.25
    if betaBarrel:
        core_region = [-4.0, -10.0]
        core_region2 = [-2.0, -10.0]
    # allocate average potentials/per layer in core region
    # sum(mem_pot[start, stop])/nslices
    print("Optimize lipid head location and estimate thickness by grid search"
          " [{:.4f}, {:.4f}] and [{:.4f}, {:.4f}]"
          .format(lipid_top_head[0], lipid_top_head[1],
                  lipid_bot_head[0], lipid_bot_head[1]))
    core_slice_index = [get_slice_index(core_region[0]),
                        get_slice_index(core_region[1])]
    # core_slice_index2 = [get_slice_index(core_region2[0]),
    #                     get_slice_index(core_region2[1])]
    # nslice_core = core_slice_index[1] - core_slice_index[0] + 1
    # nslice2_core = core_slice_index2[1] - core_slice_index2[0] + 1
    # ave = np.sum(mem_pot[:, core_slice_index])/nslice_core
    # mean by column
    ave_core = np.mean(mem_pot[:, core_slice_index[0]:core_slice_index[1]],
                       axis=1)
    # allocate average potentials/per layer outside core region
    #
    ave_out = np.mean(
        np.hstack((mem_pot[:, :core_slice_index[0]],
                  mem_pot[:, core_slice_index[1]:])),
        axis=1)
    # Z axis: max and min
    cmax = np.amax(ag.getCoords(), axis=0)
    cmin = np.amin(ag.getCoords(), axis=0)
    zmax = cmax[2]
    zmin = cmin[2]
    # grid search for lipid head in top and bottom of the membrane
    lowest_energy = np.inf
    best_e = 0.0  # extracellular side
    best_c = 0.0  # cytoplasma side
    for e in np.arange(core_region[0], core_region[1], stepsize):
        energy_e = 0.0
        for atom in ag.iterAtoms():
            x, y, z = atom.getCoords()
            resn = atom.getResname()
            z -= e
            if lipid_top_head[0] <= z <= lipid_top_head[1]:
                # lipid head region
                inx = get_slice_index(z)
                energy_e += mem_pot[aa_order_map[resn]][inx]
            elif z < lipid_top_head[0]:
                # core region of the membrane
                energy_e += ave_core[aa_order_map[resn]]
            else:
                # outside of the membrane
                energy_e += ave_out[aa_order_map[resn]]

        for c in np.arange(core_region2[0], core_region2[1], stepsize):
            energy_c = 0.0
            for atom in ag.iterAtoms():
                x, y, zx = atom.getCoords()
                resn = atom.getResname()
                zx += c
                # lipid_bot_head = [-20.0, -10.0]
                if lipid_bot_head[0] <= zx <= lipid_bot_head[1]:
                    # lipid head region
                    inx = get_slice_index(zx)
                    energy_c += mem_pot[aa_order_map[resn]][inx]
                elif zx > lipid_bot_head[0]:
                    # core region of the membrane
                    energy_c += ave_core[aa_order_map[resn]]
                else:
                    # outside of the membrane
                    energy_c += ave_out[aa_order_map[resn]]
            total_energy = energy_e + energy_c
            # judgement
            if (total_energy < lowest_energy) and\
               ((15.0 + best_e) < zmax and (-15.0 - best_c) > zmin):
                print("thickness = {:.3f} - cyto:{:.4f} extra:{:.4f} "
                      "total:{}".format(
                         (30.0+best_e+best_c), best_c, best_e, total_energy))
                print("{:.4f} 0===== =====0 {:.4f}"
                      .format(15+best_e, -15-best_c))
                print("{:.4f} <===========> {:.4f}".format(zmax, zmin))
                lowest_energy = total_energy
                best_e = e
                best_c = c
    thickness = 30.0 + best_c + best_e
    print("Total energy: {:.4f}".format(lowest_energy))
    print("Extracellular side shift: {:.4f}".format(best_e))
    print("cytoplasma side shift: {:.4f}".format(best_c))
    print("Estimated thickness of the membrane: {:.4f}".format(thickness))

    return thickness, best_e, best_c


def fit_to_membrane(ag,
                    mem_pot=mem_pot_data,
                    beta_pot=betamem_pot_data,
                    global_minimizer="minimahopping",
                    minimizer='SLSQP',
                    Temperature=10.0,
                    stepsize=0.5,
                    Niter=1000,
                    distfun="uniform",
                    forcespan=False,
                    verbose=False):
    """
    Find the best fit with membrane layers

    Args:
        ag (object):
            the object of bilab.structure.atomic.atomgroup

    Statistical Potential parameters:
        mem_pot (array): 20 x nsclices, statistical potential parameters
        beta_pot (array): 20 x nsclices, statistical potential parameters
                         for beta-barrel membrane proteins
    Minimization parameters for basinhopping algorithm:
        Temperature (float):
                Temperature factor used in the metropolis criterion
                    func(xnew) - func(xold)
          exp( -1* ------------------------  )
                               T
        stepsize (float): Standard monte carlo algorithm
             Using a customize take_step function which using
             uniform distribution
        Niter (integer): Iteration times.
        distfun (str): uniform or gaussian to generate a displacement
    Miscellanous:
        forcespan(bool): default is False
        verbose (bool): default is False
    """
    # Select representative points for residues
    #  Gly: CA, others: CB
    reps = ag.select(
        "(name CB and protein) or (name CA and resname GLY)").copy()
    # temp1 = ag.select("protein and name CA CB").copy()
    print("{} atoms selected for orientation".format(reps.numAtoms()))
    # move to center
    origin_shift(reps, center="geometric")
    cb = None
    if verbose:
        cb = print_cb
    coords = reps.getCoords()
    dist_ca2ca = euclidean(np.amax(coords, axis=0), np.amin(coords, axis=0))
    # check center
    center = calcCenter(reps).round(3)
    if np.count_nonzero(center) != 0:
        print("Warning: Failed to move to the origin")

    # ---- parameter's grid ------------
    alpha_stepsize = 0.0174532925
    beta_stepsize = 0.0174532925
    z_stepsize = 0.25
    zmin = -15.0
    zmax = dist_ca2ca + 15.0 + 0.001
    #
    # alpha = np.arange(0, 2*np.pi, alpha_stepsize)
    # beta = np.arange(0, 2*np.pi, beta_stepsize)
    # zrange = np.arange(zmin, zmax, z_stepsize)
    # PI = np.pi
    DB_PI = 2*np.pi
    x_max = [DB_PI, DB_PI, zmax]
    x_min = [0.0, 0.0, zmin]
    # x_max = [PI, PI, zmax]
    # x_min = [-PI, -PI, zmin]
    print("Search range:")
    print("1. max (alpha, beta, ztrans): [{}, {}, {}]".format(
          x_max[0], x_max[1], x_max[2]))
    print("2. min (alpha, beta, ztrans): [{}, {}, {}]".format(
          x_min[0], x_min[1], x_min[2]))
    bounds = SearchBounds(xmax=x_max,
                          xmin=x_min,
                          steps=[alpha_stepsize, beta_stepsize, z_stepsize])
    initial_guess = [0.0, 0.0, 0.0]  # bounds.get_random_params()
    print("Initial guess: {}, {}, {}".format(
          initial_guess[0], initial_guess[1], initial_guess[2]))
    # [np.random.choice(alpha, 1)[0],
    #  np.random.choice(beta, 1)[0],
    #  np.random.choice(zrange, 1)[0]]
    orient = Orient(reps, bounds, mem_pot=mem_pot,
                    beta_pot=beta_pot,
                    forcespan=forcespan,
                    verbose=verbose)
    # methods
    # method : str or callable, optional
    #    Type of solver. Should be one of
    #    ‘Nelder-Mead' (see here)
    #    ‘Powell' (see here)
    #    ‘CG' (see here)
    #    ‘BFGS' (see here)
    #    ‘Newton-CG' (see here)
    #    ‘L-BFGS-B' (see here)
    #    ‘TNC' (see here)
    #    ‘COBYLA' (see here)
    #    ‘SLSQP' (see here)
    #    ‘dogleg' (see here)
    #    ‘trust-ncg' (see here)
    # custom - a callable object (added in version 0.14.0),
    # see below for description. If not given, chosen to be one of
    # BFGS, L-BFGS-B, SLSQP, depending if the problem has
    # constraints or bounds.
    # Constrained minimization:
    # 1. 'L-BFGS-B'
    # minimizer_kwargs = {'method': 'L-BFGS-B', "jac": True}
    # 2. 'SLSQP'
    # minimizer_kwargs = {'method': 'SLSQP'}
    # 3. Customize
    # minimizer_kwargs = {'method': local_minimizer}
    minimizer_kwargs = None
    if minimizer == 'local':
        minimizer_kwargs = {'method': local_minimizer}
    else:
        minimizer_kwargs = {'method': minimizer}

    takeStep = TakeStep(stepsize=stepsize, distfun=distfun)
    if global_minimizer == "minimahopping":
        ret = bilab.optimization.minimahopping(
            orient,
            initial_guess,
            niter=Niter,
            T=Temperature,
            stepsize=stepsize,
            minimizer_kwargs=minimizer_kwargs,
            take_step=takeStep,
            accept_test=bounds,
            callback=cb,
            interval=50,
            disp=False,
            niter_success=None
        )
    elif global_minimizer == "basinhopping":
        ret = scipy.optimize.basinhopping(
                orient,
                initial_guess,
                niter=Niter,
                T=Temperature,
                stepsize=stepsize,
                minimizer_kwargs=minimizer_kwargs,
                take_step=takeStep,
                accept_test=bounds,
                callback=cb,
                interval=50,
                disp=False,
                niter_success=None)
    else:
        print("Unknown the name of global minimizer: {}".format(
            global_minimizer))
        sys.exit(1)

    f_alpha, f_beta, f_ztrans = ret['x']
    f_gamma = 0.0
    f_delta_x = 0.0
    f_delta_y = 0.0
    transformation = generate_trans(f_alpha, f_beta, f_gamma,
                                    f_delta_x, f_delta_y, f_ztrans)
    if verbose:
        print_optResults(ret, minimizer)
    # optimize the polar heads of the membrane and estimate the thickness
    # center = calcCenter(reps).round(3)
    # print("center before transformation:")
    # print(center)
    reps_at = transformation.apply(reps)
    # center_after = calcCenter(reps_at).round(3)
    # print("center after transformation:")
    # print(center_after)
    # coords_max = np.amax(reps_at.getCoords(), axis=0)
    # coords_min = np.amin(reps_at.getCoords(), axis=0)
    # xyplane = [coords_min[0], coords_max[0], coords_min[1], coords_max[1]]
    # center = calcCenter(reps_at).round(3)

    thickness, best_e, best_c = optimize_lipid_head(
        reps_at, mem_pot=mem_pot, beta_pot=beta_pot, stepsize=0.25,
        betaBarrel=False)
    # rotate all atoms
    atoms = transformation.apply(ag)
    atoms_xy_max = np.amax(atoms.getCoords(), axis=0)
    atoms_xy_min = np.amin(atoms.getCoords(), axis=0)
    at_xy_plane = [atoms_xy_min[0], atoms_xy_max[0],
                   atoms_xy_min[1], atoms_xy_max[1]]
    # Add membrane
    mem_atoms = generate_dum_atoms(
        best_e, best_c,
        xyplane=at_xy_plane, offset=8)
    return atoms, mem_atoms
