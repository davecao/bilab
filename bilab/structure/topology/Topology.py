# -*- coding: utf-8 -*-
import numpy as np
from functools import wraps
from bilab.structure.topology import AMBER_TOPO_FIELDS
from bilab.structure.topology.TopoFields import wrapGetMethod, wrapSetMethod

__all__ = ['Topology', 'AmberTopology']


# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    line = method.func_code.co_firstlineno
    filename = method.func_code.co_filename

    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError(
            'Abstract method %s at File "%s", line %s'
            'should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented

AMBER_TOPO_FLAGS = (
    'TITLE',
    'POINTERS',
    'ATOM_NAME',
    'CHARGE',
    'ATOMIC_NUMBER',
    'MASS',
    'ATOM_TYPE_INDEX',
    'NUMBER_EXCLUDED_ATOMS',
    'NONBONDED_PARM_INDEX',
    'RESIDUE_LABEL',
    'RESIDUE_POINTER',
    'BOND_FORCE_CONSTANT',
    'BOND_EQUIL_VALUE',
    'ANGLE_FORCE_CONSTANT',
    'ANGLE_EQUIL_VALUE',
    'DIHEDRAL_FORCE_CONSTANT',
    'DIHEDRAL_PERIODICITY',
    'DIHEDRAL_PHASE',
    'SCEE_SCALE_FACTOR',
    'SCNB_SCALE_FACTOR',
    'SOLTY',
    'LENNARD_JONES_ACOEF',
    'LENNARD_JONES_BCOEF',
    'BONDS_INC_HYDROGEN',
    'BONDS_WITHOUT_HYDROGEN',
    'ANGLES_INC_HYDROGEN',
    'ANGLES_WITHOUT_HYDROGEN',
    'DIHEDRALS_INC_HYDROGEN',
    'DIHEDRALS_WITHOUT_HYDROGEN',
    'EXCLUDED_ATOMS_LIST',
    'HBOND_ACOEF',
    'HBOND_BCOEF',
    'HBCUT',
    'AMBER_ATOM_TYPE',
    'TREE_CHAIN_CLASSIFICATION',
    'JOIN_ARRAY',
    'IROTAT',
    'RADIUS_SET',
    'RADII',
    'SCREEN'
)

AMBER_IFBOX_FLAGS = (
    'SOLVENT_POINTERS',
    'ATOMS_PER_MOLECULE',
    'BOX_DIMENSIONS'
)

AMBER_IFCAP_FLAGS = (
    'CAP_INFO',
    'CAP_INFO2'
)

AMBER_IFPERT_FLAGS = set()

AMBER_IPOL_FLAGS = ('POLARIZABILITY')

AMBER_IPOL_IFPERT_FLAGS = set()


class Topology(object):
    """Base class """

    __slots__ = []

    def __getattribute__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'.format(
                self.__class__.__name__, name))

    def __getstate__(self):

        return dict([(slot, getattr(self, slot))
                     for slot in self.__class__.__slots__])

    def __setstate__(self, state):

        for slot in self.__class__.__slots__:
            try:
                value = state[slot]
            except KeyError:
                pass
            else:
                setattr(self, slot, value)


class AmberTopology(Topology):
    __slots__ = ['_title', '_n_atoms', '_coords',
                 '_timestamps', '_acsi',
                 '_data',
                 '_flags',
                 '_subsets']

    def __init__(self, title='Unnamed'):

        self._title = str(title)
        self._n_atoms = 0
        self._coords = None
        self._timestamps = None
        self._data = dict()
        self._acsi = None

        self._flags = None
        self._subsets = None

    def __str__(self):
        return 'AmberTopolgy ' + self._title

    def __len__(self):
        return self._n_atoms

#    def __contains__(self, item):
#        try:
#            item.getACSIndex()
#        except AttributeError:
#            return False
#        else:
#            try:
#                ag = item.getAtomGroup()
#            except AttributeError:
#                return item == self
#            else:
#                return ag == self


for fname, field in AMBER_TOPO_FIELDS.items():

    meth = field.meth_pl
    getMeth = 'get' + meth
    setMeth = 'set' + meth
    if field.call:
        # Define public method for retrieving a copy of data array
        if not field.private:
            def getData(self, var=fname, call=field.call):
                try:
                    return self._data[var].copy()
                except KeyError:
                    [getattr(self, meth)() for meth in call]
                    return self._data[var].copy()

        # Define private method for retrieving actual data array
        def _getData(self, var=fname, call=field.call):
            try:
                return self._data[var]
            except KeyError:
                [getattr(self, meth)() for meth in call]
                return self._data[var].copy()
    else:
        if not field.private:
            def getData(self, var=fname):
                try:
                    return self._data[var].copy()
                except KeyError:
                    pass

        def _getData(self, var=fname):
            return self._data.get(var)

    if not field.private:
        getData = wrapGetMethod(getData)
        getData.__name__ = getMeth
        getData.__doc__ = field.getDocstr('get')
        setattr(AmberTopology, getMeth, getData)

    _getData = wrapGetMethod(_getData)
    _getData.__name__ = '_' + getMeth
    _getData.__doc__ = field.getDocstr('_get')
    setattr(AmberTopology, '_' + getMeth, _getData)

    if field.readonly or field.private:
        continue

    # Define public method for setting values in data array
    def setData(self, array, var=fname, dtype=field.dtype,
                ndim=field.ndim, none=field.none, flags=field.flags):
        if array is None:
            self._data.pop(var, None)
        else:
            if self._n_atoms == 0:
                self._n_atoms = len(array)
            elif len(array) != self._n_atoms:
                raise ValueError('length of array must match number '
                                 'of atoms')

            if isinstance(array, list):
                array = np.array(array, dtype)
            elif not isinstance(array, np.ndarray):
                raise TypeError('array must be an ndarray or a list')
            elif array.ndim != ndim:
                    raise ValueError('array must be {0} '
                                     'dimensional'.format(ndim))
            elif array.dtype != dtype:
                try:
                    array = array.astype(dtype)
                except ValueError:
                    raise ValueError('array cannot be assigned type '
                                     '{0}'.format(dtype))
            self._data[var] = array
            if none:
                self._none(none)
            if flags and self._flags:
                self._resetFlags(var)

    setData = wrapSetMethod(setData)
    setData.__name__ = setMeth
    setData.__doc__ = field.getDocstr('set')
    setattr(AmberTopology, setMeth, setData)

del getData
del _getData
del setData

# class AmberTopology(Topology):
#     __slots__ = [
#        '_ititl',
#        '_natom', '_ntypes', '_nbonh', '_mbona', '_ntheth', '_mtheta',
#        '_nphih', '_mphia', '_nhparm', '_nparm', '_nnb', '_nres',
#        '_nbona', '_ntheta', '_nphia', '_numbnd', '_numang', '_nptra',
#        '_natyp', '_nphb', '_ifpert', '_nbper', '_ngper', '_ndper',
#        '_mbper', '_mgper', '_mdper', '_ifbox', '_nmxrs', '_ifcap',
#        '_numextra', '_ncopy',
#        '_igraph',
#        '_charge',
#        '_atnum',
#        '_amass',
#        '_iac',
#        '_numex',
#        '_ico',
#        '_lbres',
#        '_ipres',
#        '_rk',
#        '_req',
#        '_tk',
#        '_teq',
#        '_pk',
#        '_pn',
#        '_phase',
#        '_one_scee',
#        '_one_scnb',
#        '_solty',
#        '_cn1',
#        '_cn2',
#        '_ibh', '_jbh', '_icbh',
#        '_ib', '_jb', '_icb',
#        '_ith', '_jth', '_kth', '_icth',
#        '_it', '_jt', '_kt', '_ict',
#        '_iph', '_jph', '_kph', '_lph', '_icph',
#        '_ip', '_jp', '_kp', '_lp', '_icp',
#        '_inb',
#        '_asol',
#        '_bsol',
#        '_hcut',
#        '_isymbl',
#        '_itree',
#        '_join',
#        '_irotat',
#        '_type',
#        '_rborn',
#        '_fs',
#        '_iptres', '_nspm', '_nspsol',
#        '_nsp',
#        '_oldbeta', '_box1', '_box2', '_box3',
#        '_natcap',
#        '_cutcap', '_xcap', '_ycap', '_zcap',
#        '_ibper', '_jbper',
#        '_icper',
#        '_ipter', '_jpter', '_ktper',
#        '_ictper'
#        '_ipper', '_jpper', '_kpper', '_lpper',
#        '_icpper',
#        '_labres',
#        '_igrper',
#        '_ismper',
#        '_almper',
#        '_iaper',
#        '_iacper',
#        '_cgper',
#        '_atpol',
#        '_atpol1'
#     ]
#     """docstring for AmberTopology"""
#     def __init__(self):
#         super(AmberTopology, self).__init__()
#         self._ititl = None
#         # %FLAG POINTERS
#         self._natom = 0, self._ntypes = 0, self._nbonh = 0,
#         self._mbona = 0, self._ntheth = 0, self._mtheta = 0,
#         self._nphih = 0, self._mphia = 0, self._nhparm = 0,
#         self._nparm = 0, self._nnb = 0, self._nres = 0,
#         self._nbona = 0, self._ntheta = 0, self._nphia = 0,
#         self._numbnd = 0, self._numang = 0, self._nptra = 0,
#         self._natyp = 0, self._nphb = 0, self._ifpert = 0,
#         self._nbper = 0, self._ngper = 0, self._ndper = 0,
#         self._mbper = 0, self._mgper = 0, self._mdper = 0,
#         self._ifbox = 0, self._nmxrs = 0, self._ifcap = 0,
#         self._numextra = 0, self._ncopy = 0
#         # %FLAG ATOM_NAME
#         self._igraph = []
#         # %FLAGS CHARGE
#         self._charge = []
#         # %FLAG ATOMIC_NUMBER
#         self._atnum = []
#         # %FLAG MASS
#         self._amass = []
#         # %FLAG ATOM_TYPE_INDEX
#         self._iac = []
#         # %FLAG NUMBER_EXCLUDED_ATOMS
#         self._numex = []
#         # %FLAG NONBONDED_PARM_INDEX
#         self._ico = []
#         # %FLAG RESIDUE_LABEL
#         self._lbres = [],
#         self._ipres = [],
#         self._rk = [],
#         self._req = [],
#         self._tk',
#         self._teq',
#         self._pk',
#         self._pn',
#         self._phase',
#         self._one_scee',
#         self._one_scnb',
#         self._solty',
#         self._cn1',
#         self._cn2',
#         self._ibh', '_jbh', '_icbh',
#         self._ib', '_jb', '_icb',
#         self._ith', '_jth', '_kth', '_icth',
#         self._it', '_jt', '_kt', '_ict',
#         self._iph', '_jph', '_kph', '_lph', '_icph',
#         self._ip', '_jp', '_kp', '_lp', '_icp',
#         self._inb',
#         self._asol',
#         self._bsol',
#         self._hcut',
#         self._isymbl',
#         self._itree',
#         self._join',
#         self._irotat',
#         self._type',
#         self._rborn',
#         self._fs',
#         self._iptres', '_nspm', '_nspsol',
#         self._nsp',
#         self._oldbeta', '_box1', '_box2', '_box3',
#         self._natcap',
#         self._cutcap', '_xcap', '_ycap', '_zcap',
#         self._ibper', '_jbper',
#         self._icper',
#         self._ipter', '_jpter', '_ktper',
#         self._ictper'
#         self._ipper', '_jpper', '_kpper', '_lpper',
#         self._icpper',
#         self._labres',
#         self._igrper',
#         self._ismper',
#         self._almper',
#         self._iaper',
#         self._iacper',
#         self._cgper',
#         self._atpol',
#         self._atpol1'
# 
