# -*- coding: utf-8 -*-

from numpy import array

# Literature:
#   Beatriz et al. Covalent radii revisited, Dalton Trans(21):2832-2838.
#      doi:10.1039/b801115j Table 2 1-96
__all__ = ['ElementData']

DTYPE = array(['a']).dtype.char  # 'S' for PY2K and 'U' for PY3K

class ElementProperties(object):
    """Element data field."""

    __slots__ = ['name', 'dtype', 'doc', 'doc_pl', 
                 'meth', 'meth_pl','ndim',
                 'depr', 'depr_pl', 
                 'desc', 'long_desc','private', 'call',
                 'value'
                 ]

    def __init__(self, name, dtype, **kwargs):

        #: data field name used in atom selections
        self.name = name
        #: data type (primitive Python types)
        self.dtype = dtype
        #: internal variable name used as key for :class:`.AtomGroup` ``_data``
        self.doc = kwargs.get('doc', name)
        #: plural form for documentation
        self.doc_pl = kwargs.get('doc_pl', self.doc + 's')
        #: description of data field, used in documentation
        self.desc = kwargs.get('desc')
        #  long description of data field
        self.long_desc = kwargs.get('long_desc')
        #: expected dimension of the data array
        self.ndim = kwargs.get('ndim', 1)
        #: atomic get/set method name
        self.meth = kwargs.get('meth', name.capitalize())
        #: get/set method name in plural form
        self.meth_pl = kwargs.get('meth_pl', self.meth + 's')
        #: get/set 
        self.private = kwargs.get('private', False)
        #: deprecated method name
        self.depr = kwargs.get('depr')
        #: deprecated method name in plural form
        self.call = kwargs.get('call', None)
        self.depr_pl = None
        if self.depr is not None:
            self.depr_pl = kwargs.get('depr_pl', self.depr + 's')

    def getDocstr(self, meth, plural=True):
        """ Return docstring for the field """
        assert meth in ('set', 'get', '_get'), "meth must be 'set' or 'get'"
        assert isinstance(plural, bool), 'plural must be a boolean'

# Periodic table
#         Group 1, 2, ... , 18
# Period  H
#         Li

# Categories:
#  Actinide  
#  AlkaliMetal
#  Gas         standard temprature and pressure
#  Lanthanide 
#  Liquid      standard temprature and pressure
#  Ferromagnetic  standard temprature and pressure
#  Metal
#  Metalloid
#  Natural        natural material
#  Nonmetal  
#  NobleGas
#  Radioactive
#  Solid       standard temprature and pressure
#  Stable 
#  Synthetic

# Shell configuration:
#   K(1), L(2), M(3), N(4), O(5), P(6), Q(7)
# Subshell 
#  s: l=0, max=2,  every shell,  sharp 
#  p: l=1  max=6,  2nd shell and higher, principle
#  d: l=2, max=10, 3rd shell and higher, diffuse
#  f: l=3, max=14, 4th shell and higher, fundamental
#  g: l=4, max=18, 5th shell and higher,
# e.g.,
#  Atomic number, Element, No. electrons/shell, group
#   Z=1, Hydrogen, 1,   1
#   Z=2, Helium,   2,   18
#   Z=3, Lithium,  2,1, 1
#    ...

ELEMENT_FIELDS = {
    'AtomicNumber': ElementProperties('AtomicNumber', int, 
                    doc='atomic number Z'),

    'AtomicWeight': ElementProperties('AtomicWeight', float, 
                    doc='atomic weight'),

    'AtomicMass' : ElementProperties('AtomicMass', float, 
                   doc='Atomic mass of the element' ),

    'Symbol' : ElementProperties('Symbol', DTYPE+'1', 
                    doc=' the symbol of the element'),

    'Color' : ElementProperties('Color', float, ndim=3, 
                    doc=' the color of the element'),

    'Block' : ElementProperties('Block', DTYPE+'1',
                    doc='the block of the element in the periodic table'),

    'Group' : ElementProperties('Group', int,
                    doc='the group of the element in the periodic table'),

    'Period' : ElementProperties('Period', int,
                    doc='the period of the element in the periodic table'),

    'Series' : ElementProperties('Series', int,
                    doc='the series of the element in the periodic table'),

    'CrystalStructure': ElementProperties('CrystalStructure', DTYPE+'1',
                    doc='the type of crystal structure'),

    'LatticeAngles' : ElementProperties('LatticeAngles', float, ndim=3,
                    doc='lattic angle in crystal structure'),

    'SpaceGroupNumber' : ElementProperties('SpaceGroupNumber',int,
                    doc='space group number.'),

    'SpaceGroup' : ElementProperties('SpaceGroup', DTYPE+'1', 
                    doc='space group., i.e., P63/mmc'),

    'AtomicRadius' : ElementProperties('AtomicRadius', float, 
                    doc='atmoic radius of the element'),

    'CovalentRadius': ElementProperties('CovalentRadius', float, 
                    doc='covalent raidus of the element'),

    'VanDerWaalsRadius': ElementProperties('VanDerWaalsRadius', float,
                    doc='van der Waal radius of the element.'),

    'ElectronConfiguration' : ElementProperties('ElectronConfiguration', DTYPE+'1',
                    doc='The electron configuration., i.e., H for 1s1'),

    'ElectronConfigurationString' : ElementProperties('ElectronConfigurationString', DTYPE+'1',
                    doc='A string of the electron configuration'),

    'ElectronShellConfiguration': ElementProperties('ElectronShellConfiguration', DTYPE+'1',
                    doc='The electron shell configuration., i.e., '),

    'IonizationEnergies' : ElementProperties('IonizationEnergies', float,
                    doc='Ionization energy of the element. i.e., H for 1312kj/mol'),

    'QuantumNumbers' : ElementProperties('QuantumNumbers', DTYPE+'1', 
                    doc='Quantum number of the element'),

    'Abbreviation' : ElementProperties('Abbreviation', DTYPE+'1', 
                    doc='Abbreviation of the element'),

    'CASNumber' : ElementProperties('CASNumber', DTYPE+'1',
                    doc='CAS number of the element'),

    'PubChemCID' : ElementProperties('PubChemCID', DTYPE+'1', 
                    doc='PubChem CID of the element'),

    'StandardName' : ElementProperties('StandardName', DTYPE+'1', 
                    doc='Name of the element, e.g., Hydrogen'),

    'Category' : ElementProperties('Category', DTYPE+'1', 
                    doc='Category of the element, e.g., Metal, Nonmetal')
}

ATOM_COVALENT_RADIUS = {
    'h' :0.31,
    'he':0.28,
    'li':1.28,
    'be':0.96,
    'b' :0.84,
    'c' :0.76,
    'n' :0.71,
    'o' :0.66,
    'f' :0.57,
    'ne':0.58,
    'na':1.66,
    'mg':1.41,
    'al':1.21,
    'si':1.11,
    'p' :1.07,
    's' :1.05,
    'cl':1.02,
    'ar':1.06,
    'k' :2.03,
    'ca':1.76,
    'sc':1.70,
    'ti':1.60,
    'v' :1.53,
    'cr':1.39,
    'mn':1.39,
    'fe':1.32,
    'co':1.26,
    'ni':1.24,
    'cu':1.32,
    'zn':1.22,
    'ga':1.22,
    'ge':1.20,
    'as':1.19,
    'se':1.20,
    'br':1.20,
    'kr':1.16,
    'rb':2.20,
    'sr':1.95,
    'y' :1.90,
    'zr':1.75,
    'nb':1.64,
    'mo':1.54,
    'tc':1.47,
    'ru':1.46,
    'rh':1.42,
    'pd':1.39,
    'ag':1.45,
    'cd':1.44,
    'in':1.42,
    'sn':1.39,
    'sb':1.39,
    'te':1.38,
    'i' :1.39,
    'xe':1.40,
    'cs':2.44,
    'ba':2.15,
    'la':2.07,
    'ce':2.04,
    'pr':2.03,
    'nd':2.01,
    'pm':1.99,
    'sm':1.98,
    'eu':1.98,
    'gd':1.96,
    'tb':1.94,
    'dy':1.92,
    'ho':1.92,
    'er':1.89,
    'tm':1.90,
    'yb':1.87,
    'lu':1.87,
    'hf':1.75,
    'ta':1.70,
    'w' :1.62,
    're':1.51,
    'os':1.44,
    'ir':1.41,
    'pt':1.36,
    'au':1.36,
    'hg':1.32,
    'tl':1.45,
    'pb':1.46,
    'bi':1.48,
    'po':1.40,
    'at':1.50,
    'rn':1.50,
    'fr':2.60,
    'ra':2.21,
    'ac':2.15,
    'th':2.06,
    'pa':2.00,
    'u' :1.96,
    'np':1.90,
    'pu':1.87,
    'am':1.80,
    'cm':1.69,
    'bk':0.0,
    'cf':0.0,
    'es':0.0,
    'fm':0.0,
    'md':0.0,
    'no':0.0,
    'lr':0.0,
    'rf':0.0,
    'db':0.0,
    'sg':0.0,
    'bh':0.0,
    'hs':0.0,
    'mt':0.0,
    'ds':0.0,
    'rg':0.0,
    'cn':0.0,
    'uut':0.0,
    'uuq':0.0,
    'uup':0.0,
    'uuh':0.0,
    'uus':0.0,
    'uuo':0.0
}

def get_covalent_radius(element):
    """ find the covalent radius of a chemical element in periodic table.

    Args:
        element (str): the name of a chemical element

    Returns:
        covalent radius (float)

    Example:

    .. ipython:: python

        import bilab
        bilab.chemicals.get_covalent_radius('h')

    """

    if element.lower() in ATOM_COVALENT_RADIUS:
        return ATOM_COVALENT_RADIUS[element.lower()]
    return None

class ElementData(object):
    """ Checmical elements data in periodic table """
    __slots__ = []

    def __init__(self, title='Periodic Table'):
        self._title = str(title)


def wrapGetMethod(fn):
    def getMethod(self):
        return fn(self)
    return getMethod


def wrapSetMethod(fn):
    def setMethod(self, data):
        return fn(self, data)
    return setMethod

for fname, eprop in ELEMENT_FIELDS.items():
    meth = eprop.meth_pl
    getMeth = 'get' + meth
    setMeth = 'set' + meth
    if eprop.call:
        # define public method 
        if eprop.private:
            def getData(self, var=fname, call=eprop.call):
                try:
                    return self._data[var].copy()
                except KeyError:
                    [ getattr(self, meth)() for meth in call ]
                return self._data[var].copy()
        # Define private method for retrieving actual data array
        def _getData(self, var=fname, call=eprop.call):
            try:
                return self._data[var]
            except KeyError:
                [getattr(self, meth)() for meth in call]
                return self._data[var].copy()
    else:
        if not eprop.private:
            def getData(self, var=fname):
                try:
                    return self._data[var].copy()
                except KeyError:
                    pass

        def _getData(self, var=fname):
            return self._data.get(var)

    if not eprop.private:
        getData = wrapGetMethod(getData)
        getData.__name__ = getMeth
        getData.__doc__ = eprop.getDocstr('get')
        setattr(ElementData, getMeth, getData)

    _getData = wrapGetMethod(_getData)
    _getData.__name__ = '_' + getMeth
    _getData.__doc__ = eprop.getDocstr('_get')
    setattr(ElementData, '_' + getMeth, _getData)


    # Define public method for setting values in data array
    def setData(self, array, var=fname, dtype=eprop.dtype,
                ndim=eprop.ndim, none=eprop.none, flags=eprop.flags):
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

    setData = wrapSetMethod(setData)
    setData.__name__ = setMeth
    setData.__doc__ = eprop.getDocstr('set')
    setattr(eprop, setMeth, setData)

del getData
del _getData
del setData
