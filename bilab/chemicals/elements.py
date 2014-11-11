# -*- coding: utf-8 -*-

import sys
import numpy as np

# Literature:
#   Beatriz et al. Covalent radii revisited, Dalton Trans(21):2832-2838.
#      doi:10.1039/b801115j Table 2 1-96
__all__ = [ 'ElementData', 
            'Element',
            'ElementCollection',
            'PeriodicTable',
            'get_covalent_radius']

#DTYPE = np.array(['a']).dtype.char  # 'S' for PY2K and 'U' for PY3K
DTYPE = np.str
PI = np.pi

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
        #self.meth = kwargs.get('meth', name.capitalize())
        self.meth = kwargs.get('meth', name)
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

#  Block
#  g1                                                                       g18
# H 1s  | g2                                               g13,14,15,16,17 1s->|
# Li <- 2s ->|                                            |<---- 2p          ->|
# Na <- 3s ->|            g3 g4 g5 g6 g7 g8 g9 g10 g11 g12|<---- 3p          ->|
# K  <- 4s ->|           |<--        3d                -->|<---- 4p          ->|
# Rb <- 5s ->|           |<--        4d                -->|<---- 5p          ->|
# Cs <- 6s ->|<-- 4f  -->|<--        5d                -->|<---- 6p          ->|
# Fr <- 7s ->|<-- 5f  -->|<--        6d                -->|<---- 7p          ->|
# <-s-block->|<-f-block->|<--       d-block            -->|<---- p-block     ->|

# Categories/Series:
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

    'StandardName' : ElementProperties('StandardName', np.str, 
                    doc='Name of the element, e.g., Hydrogen'),

    'Symbol' : ElementProperties('Symbol', DTYPE, 
                    doc=' the symbol of the element'),

    'Abbreviation' : ElementProperties('Abbreviation', DTYPE, 
                    doc='Abbreviation of the element'),

    'AtomicNumber': ElementProperties('AtomicNumber', int, 
                    doc='atomic number Z', 
                    desc="The atomic number of a chemical element (also known"\
                    +" as its proton number) is the number of protons found"\
                    +" in nucleus of an atom of the element."),

    'AtomicWeight': ElementProperties('AtomicWeight', float, 
                    doc='atomic weight'),

    'AtomicMass' : ElementProperties('AtomicMass', DTYPE, 
                   doc='Atomic mass of the element' ),


    'Color' : ElementProperties('Color', float, ndim=3, 
                    doc=' the color of the element'),

    'Block' : ElementProperties('Block', DTYPE,
                    doc='the block of the element in the periodic table'),

    'Group' : ElementProperties('Group', int,
                    doc='the group of the element in the periodic table'),

    'Period' : ElementProperties('Period', int,
                    doc='the period of the element in the periodic table'),

    'CrystalStructure': ElementProperties('CrystalStructure', DTYPE,
                    doc='the type of crystal structure'),

    'LatticeAngles' : ElementProperties('LatticeAngles', float, ndim=3,
                    doc='lattic angle in crystal structure'),
    
    'LatticeConstants' : ElementProperties('LatticeAngles', float, ndim=3,
                    doc='lattic constants(x, y, z) in crystal structure'),

    'SpaceGroupNumber' : ElementProperties('SpaceGroupNumber',int,
                    doc='space group number.'),

    'SpaceGroup' : ElementProperties('SpaceGroup', DTYPE, 
                    doc='space group., i.e., P63/mmc'),

    'AtomicRadius' : ElementProperties('AtomicRadius', float, 
                    doc='atmoic radius of the element in pm 10x10^-12'),

    'CovalentRadius': ElementProperties('CovalentRadius', float, 
                    doc='covalent raidus of the element'),

    'VanDerWaalsRadius': ElementProperties('VanDerWaalsRadius', float,
                    doc='van der Waal radius of the element.'),

    'ElectronConfiguration' : ElementProperties('ElectronConfiguration', DTYPE,
                    doc='The electron configuration., i.e., H for 1s1'),

    'ElectronConfigurationString' : ElementProperties('ElectronConfigurationString', DTYPE,
                    doc='A string of the electron configuration'),

    'ElectronShellConfiguration': ElementProperties('ElectronShellConfiguration', DTYPE,
                    doc='The electron shell configuration., i.e., '),

    'IonizationEnergies' : ElementProperties('IonizationEnergies', float,
                    doc='First ionization energy of the element. i.e., H for 1312kj/mol'),

    'QuantumNumbers' : ElementProperties('QuantumNumbers', DTYPE, 
                    doc='Quantum number of the element'),


    'CASNumber' : ElementProperties('CASNumber', DTYPE,
                    doc='CAS number of the element'),

    'PubChemCID' : ElementProperties('PubChemCID', DTYPE, 
                    doc='PubChem CID of the element'),

    'Category' : ElementProperties('Category', DTYPE, 
                    doc='Category of the element, e.g., Metal, Nonmetal')
}

class Element(object):
    """ Checmical element's data in periodic table """
    __slots__ = ['_el', '_index', '_data']

    def __init__(self):
        super(Element, self).__init__()
        self._data = dict()
        self._el = []
        self._index = []

    def getData(self, label):
        """ 
            Return a copy of data associated with *label*, if it is present. 
        """
        try:
            data = self._el._getData(label)
        except KeyError:
            pass
        else:
            if data.ndim > 1:
                return data[self._index]
            else:
                return data[self._index].copy()
    _getData = getData

    def setData(self, label, data):
        """ Update *data* associated with *label*.

        :raise AttributeError: when *label* is not in use or read-only"""
        if label in ELEMENT_FIELDS:
            getattr(self, 'set' + ELEMENT_FIELDS[label].meth)(data)
        else:
            try:
                self._el._data[label][self._index] = data
            except KeyError:
                raise AttributeError('data with label {0} must be set for'
                                ' Element first'.format(repr(label)))
        return self

def wrapGetMethod(fn):
    def getMethod(self):
        return fn(self)
    return getMethod


def wrapSetMethod(fn):
    def setMethod(self, data):
        return fn(self, data)
    return setMethod

for fname, eprop in ELEMENT_FIELDS.items():
    #meth = eprop.meth_pl
    meth = eprop.meth
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
        setattr(Element, getMeth, getData)

    _getData = wrapGetMethod(_getData)
    _getData.__name__ = '_' + getMeth
    _getData.__doc__ = eprop.getDocstr('_get')
    setattr(Element, '_' + getMeth, _getData)


    # Define public method for setting values in data array
    def setData(self, array, var=fname, dtype=eprop.dtype, ndim=eprop.ndim):
        if array is None:
            self._data.pop(var, None)
        else:

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
    setattr(Element, setMeth, setData)

del getData
del _getData
del setData

class ElementCollection(object):
    """
        Periodic table
    """
    __slots__ = ['_elename', '_z', '_el']

    def __init__(self):
        """ Initialization """
        super(ElementCollection, self).__init__()
        self._elename = []
        self._z = []
        self._el = []

    def add(self, el):
        """ Add each element """
        if isinstance(el, Element):
            self._add_element(el)
        elif isinstance(el, dict):
            element_instance = Element()
            for key, val in el.iteritems():
                print("key:{} val:{}".format(key, val))
                if val.isdigit():
                    element_instance.setData(key, np.array(val))
                else:
                    element_instance.setData(key, val)

            self._add_element(element_instance)
        else:
            raise TypeError('el must be an instance of Element')
    def _add_element(self, el):
        self._elename.append(el.getAbbreviation())
        self._z.append(el.getAtomicNumber())
        self._el.append(el)

    def __iter__(self):
        return iter(self._el)

    def __getitem__(self, i):
        return self._el[i]

    def __len__(self):
        return len(self._el)

    def get(self, name_or_z):
        """ 
            Access the element object by Abbreviation name or atomic number

        Args:
            name_or_z (str or int): Abbreviation name of an element
                    e.g., 'Au' or 79
        
        Kwargs:

        Returns:
            An object of an element

        .. ipython:: python

        import bilab
        bilab.chemicals.ElementData.get('au')

        """
        string_types = str
        PY3 = sys.version[0] == 3
        if PY3:
            string_types = basestring
        
        if isinstance(name_or_z, string_types):
            # the input is abbreviation name of an element
            key_attr = '_elename'
        elif isinstance(name_or_z, int):
            key_attr = '_z'

        try:
            inx = getattr(self, key_attr).index(name_or_z)
            return self._el[inx]
        except ValueError:
            return None

    def getElementList(self):
        return self._elename;

# 1 angstrom = 100 pm
# atomic radius: mainly use empirical data, 
#      calculated data will be used if empirical ones are not available
#      even calculated ones are not available. see Wikioedia.org
#   http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)

PeriodicTable =[{
    'StandardName':'Hydrogen',
    'Symbol': 'H',
    'Abbreviation' : 'h',
    'AtomicNumber': 1,
    'AtomicWeight': 1.00794,
    'AtomicMass'  : '1.00794u',
    'Color' : [],
    'Block' : 's',
    'Group' : 1,
    'Period': 1,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2, PI/2, 2*PI/3],
    'LatticeConstants': [470, 470, 340],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 25,
    'CovalentRadius': 31,
    'VanDerWaalsRadius': 120,
    'ElectronConfiguration' : '1s1',
    'ElectronConfigurationString' : '1s1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1312,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '1333-74-0',
    'PubChemCID' : '783',
    'Category' : ''
    },
{
    'StandardName':'Helium',
    'Symbol': 'He',
    'Abbreviation' : 'he',
    'AtomicNumber': 2,
    'AtomicWeight': 4.002602,
    'AtomicMass'  : '4.002602u',
    'Color' : [],
    'Block' : 's',
    'Group' : 18,
    'Period': 1,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [424.2, 424.2, 424.2],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 31,
    'CovalentRadius': 28,
    'VanDerWaalsRadius': 140,
    'ElectronConfiguration' : '1s2',
    'ElectronConfigurationString' : '1s2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 2372.3,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7440-59-7',
    'PubChemCID' : '23987',
    'Category' : ''
},
{
    'StandardName':'Lithium',
    'Symbol': 'Li',
    'Abbreviation' : 'li',
    'AtomicNumber': 3,
    'AtomicWeight': 6.941,
    'AtomicMass'  : '6.941u',
    'Color' : [192, 192, 192],
    'Block' : 's',
    'Group' : 1,
    'Period': 2,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [351, 351, 351],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 145,
    'CovalentRadius': 128,
    'VanDerWaalsRadius': 182,
    'ElectronConfiguration' : '[He]2s1',
    'ElectronConfigurationString' : '[He]2s1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 520.2,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7439-93-2',
    'PubChemCID' : '3028194',
    'Category' : ''
},
{
    'StandardName':'Beryllium',
    'Symbol': 'Be',
    'Abbreviation' : 'be',
    'AtomicNumber': 4,
    'AtomicWeight': 9.012182,
    'AtomicMass'  : '9.012182u',
    'Color' : [112, 128, 144],
    'Block' : 's',
    'Group' : 2,
    'Period': 2,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2, PI/2, 2*PI/3],
    'LatticeConstants': [228.58, 228.58, 358.43],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 105,
    'CovalentRadius': 96,
    'VanDerWaalsRadius': 153,
    'ElectronConfiguration' : '[He]2s2',
    'ElectronConfigurationString' : '[He]2s2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 889.5,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7440-41-7',
    'PubChemCID' : '5460467',
    'Category' : ''
},
{
    'StandardName':'Boron',
    'Symbol': 'B',
    'Abbreviation' : 'b',
    'AtomicNumber': 5,
    'AtomicWeight': 10.811,
    'AtomicMass'  : '10.811u',
    'Color' : [0, 0, 0],
    'Block' : 'p',
    'Group' : 13,
    'Period': 2,
    'CrystalStructure': 'Simple Trigonal',
    'LatticeAngles' : [1.01334, 1.01334, 1.01334],
    'LatticeConstants': [506, 506, 506],
    'SpaceGroupNumber' : 166,
    'SpaceGroup' : 'R_3m',
    'AtomicRadius' : 85,
    'CovalentRadius': 84,
    'VanDerWaalsRadius': None,
    'ElectronConfiguration' : '[He]2s1',
    'ElectronConfigurationString' : '[He]2s1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 889.5,
    'QuantumNumbers' : '2P1/2',
    'CASNumber' : '7440-42-8',
    'PubChemCID' : '5462311',
    'Category' : ''
},

{
    'StandardName':'Carbon',
    'Symbol': 'C',
    'Abbreviation' : 'c',
    'AtomicNumber': 6,
    'AtomicWeight': 12.0107,
    'AtomicMass'  : '12.0107u',
    'Color' : [0, 0, 0],
    'Block' : 'p',
    'Group' : 14,
    'Period': 2,
    'CrystalStructure': 'Simple Hexagonale',
    'LatticeAngles' : [PI/2, PI/2, 2*PI/3],
    'LatticeConstants': [246.4, 246.4, 671.1],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 70,
    'CovalentRadius': 76,
    'VanDerWaalsRadius': 170,
    'ElectronConfiguration' : '[He]2s2,2p2',
    'ElectronConfigurationString' : '[He]2s2,2p2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1086.5,
    'QuantumNumbers' : '3P0',
    'CASNumber' : '7440-44-0',
    'PubChemCID' : '5462310',
    'Category' : ''
},
{
    'StandardName':'Nitrogen',
    'Symbol': 'N',
    'Abbreviation' : 'n',
    'AtomicNumber': 7,
    'AtomicWeight': 14.0067,
    'AtomicMass'  : '14.0067u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 15,
    'Period': 2,
    'CrystalStructure': 'Simple Hexagonale',
    'LatticeAngles' : [PI/2, PI/2, 2*PI/3],
    'LatticeConstants': [386.1, 386.1, 626.5],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 65,
    'CovalentRadius': 71,
    'VanDerWaalsRadius': 155,
    'ElectronConfiguration' : '[He]2s2,2p3',
    'ElectronConfigurationString' : '[He]2s2,2p3',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1402.3,
    'QuantumNumbers' : '4S3/2',
    'CASNumber' : '7727-37-9',
    'PubChemCID' : '947',
    'Category' : ''
},
{
    'StandardName':'Oxygen',
    'Symbol': 'O',
    'Abbreviation' : 'o',
    'AtomicNumber': 8,
    'AtomicWeight': 15.9994,
    'AtomicMass'  : '15.9994u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 16,
    'Period': 2,
    'CrystalStructure': 'Base Centered Monoclinic',
    'LatticeAngles' : [PI/2, 2.313085, PI/2],
    'LatticeConstants': [540.3, 342.9, 508.6],
    'SpaceGroupNumber' : 12,
    'SpaceGroup' : 'C12/m1',
    'AtomicRadius' : 60,
    'CovalentRadius': 66,
    'VanDerWaalsRadius': 152,
    'ElectronConfiguration' : '[He]2s2,2p4',
    'ElectronConfigurationString' : '[He]2s2,2p4',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1313.9,
    'QuantumNumbers' : '3P2',
    'CASNumber' : '7782-44-7',
    'PubChemCID' : '977',
    'Category' : ''
},
{
    'StandardName':'Fluorine',
    'Symbol': 'F',
    'Abbreviation' : 'f',
    'AtomicNumber': 9,
    'AtomicWeight': 18.9984032,
    'AtomicMass'  : '18.9984032u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 17,
    'Period': 2,
    'CrystalStructure': 'Base Centered Monoclinic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [550, 328, 728],
    'SpaceGroupNumber' : 15,
    'SpaceGroup' : 'C12/c1',
    'AtomicRadius' : 50,
    'CovalentRadius': 57,
    'VanDerWaalsRadius': 147,
    'ElectronConfiguration' : '[He]2s2,2p5',
    'ElectronConfigurationString' : '[He]2s2,2p5',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1681,
    'QuantumNumbers' : '2P3/2',
    'CASNumber' : '7782-41-4',
    'PubChemCID' : '24524',
    'Category' : ''
},
{
    'StandardName':'Neon',
    'Symbol': 'Ne',
    'Abbreviation' : 'ne',
    'AtomicNumber': 10,
    'AtomicWeight': 20.1797,
    'AtomicMass'  : '20.1797u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 18,
    'Period': 2,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [442.9, 442.9, 442.9],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 50,
    'CovalentRadius': 57,
    'VanDerWaalsRadius': 147,
    'ElectronConfiguration' : '[He]2s2,2p6',
    'ElectronConfigurationString' : '[He]2s2,2p6',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1681,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7782-41-4',
    'PubChemCID' : '24524',
    'Category' : ''
},
{
    'StandardName':'Sodium',
    'Symbol': 'Na',
    'Abbreviation' : 'na',
    'AtomicNumber': 11,
    'AtomicWeight': 22.98977,
    'AtomicMass'  : '22.98977u',
    'Color' : [192, 192, 192],
    'Block' : 's',
    'Group' : 1,
    'Period': 3,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [429.06, 429.06, 429.06],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 180,
    'CovalentRadius': 166,
    'VanDerWaalsRadius': 227,
    'ElectronConfiguration' : '[Ne]3s1',
    'ElectronConfigurationString' : '[Ne]3s1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 495.8,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7440-23-5',
    'PubChemCID' : '5360545',
    'Category' : ''
},
{
    'StandardName':'Magnesium',
    'Symbol': 'Mg',
    'Abbreviation' : 'mg',
    'AtomicNumber': 12,
    'AtomicWeight': 24.305,
    'AtomicMass'  : '24.305u',
    'Color' : [192, 192, 192],
    'Block' : 's',
    'Group' : 2,
    'Period': 3,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2, PI/2, PI/3],
    'LatticeConstants': [320.94, 320.94, 521.08],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 150,
    'CovalentRadius': 141,
    'VanDerWaalsRadius': 173,
    'ElectronConfiguration' : '[Ne]3s2',
    'ElectronConfigurationString' : '[Ne]3s2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 737.7,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7439-95-4',
    'PubChemCID' : '5462224',
    'Category' : ''
},
{
    'StandardName':'Aluminum',
    'Symbol': 'Al',
    'Abbreviation' : 'al',
    'AtomicNumber': 13,
    'AtomicWeight': 26.981538,
    'AtomicMass'  : '26.981538u',
    'Color' : [192, 192, 192],
    'Block' : 'p',
    'Group' : 13,
    'Period': 3,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [404.95, 404.95, 404.95],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 125,
    'CovalentRadius': 121,
    'VanDerWaalsRadius': 173 ,
    'ElectronConfiguration' : '[Ne]3s2,3p1',
    'ElectronConfigurationString' : '[Ne]3s2,3p1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 577.5,
    'QuantumNumbers' : '2P1/2',
    'CASNumber' : '7429-90-5',
    'PubChemCID' : '5359268',
    'Category' : ''
},
{
    'StandardName':'Silicon',
    'Symbol': 'Si',
    'Abbreviation' : 'si',
    'AtomicNumber': 14,
    'AtomicWeight': 28.0855,
    'AtomicMass'  : '28.0855u',
    'Color' : [128, 128, 128],
    'Block' : 'p',
    'Group' : 14,
    'Period': 3,
    'CrystalStructure': 'Tetrahedral Packing',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [543.09, 543.09, 543.09],
    'SpaceGroupNumber' : 227,
    'SpaceGroup' : 'Fd_3m',
    'AtomicRadius' : 110,
    'CovalentRadius': 111,
    'VanDerWaalsRadius': 210 ,
    'ElectronConfiguration' : '[Ne]3s2,3p2',
    'ElectronConfigurationString' : '[Ne]3s2,3p2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 786.5,
    'QuantumNumbers' : '3P0',
    'CASNumber' : '7440-21-3',
    'PubChemCID' : '5461123',
    'Category' : ''
},
{
    'StandardName':'Phosphorus',
    'Symbol': 'P',
    'Abbreviation' : 'p',
    'AtomicNumber': 15,
    'AtomicWeight': 30.973761,
    'AtomicMass'  : '30.973761u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 15,
    'Period': 3,
    'CrystalStructure': 'Simple Triclinic',
    'LatticeAngles' : [1.25384, 1.57725, 1.24896],
    'LatticeConstants': [1145, 550.3, 1126.1],
    'SpaceGroupNumber' : 2,
    'SpaceGroup' : 'P-1',
    'AtomicRadius' : 100,
    'CovalentRadius': 107,
    'VanDerWaalsRadius': 180 ,
    'ElectronConfiguration' : '[Ne]3s2,3p3',
    'ElectronConfigurationString' : '[Ne]3s2,3p3',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1011.8,
    'QuantumNumbers' : '4S3/2',
    'CASNumber' : '7723-14-0',
    'PubChemCID' : '',
    'Category' : ''
},
{
    'StandardName':'Sulfur',
    'Symbol': 'S',
    'Abbreviation' : 's',
    'AtomicNumber': 16,
    'AtomicWeight': 32.065,
    'AtomicMass'  : '32.065u',
    'Color' : [255,255, 0],
    'Block' : 'p',
    'Group' : 16,
    'Period': 3,
    'CrystalStructure': 'Face Centered Orthorhombic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [1043.7, 1284.5, 2436.9],
    'SpaceGroupNumber' : 70,
    'SpaceGroup' : 'Fddd',
    'AtomicRadius' : 100,
    'CovalentRadius': 105,
    'VanDerWaalsRadius': 180 ,
    'ElectronConfiguration' : '[Ne]3s2,3p4',
    'ElectronConfigurationString' : '[Ne]3s2,3p4',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 999.6,
    'QuantumNumbers' : '3P2',
    'CASNumber' : '7704-34-9',
    'PubChemCID' : '',
    'Category' : ''
},
{
    'StandardName':'Chlorine',
    'Symbol': 'Cl',
    'Abbreviation' : 'cl',
    'AtomicNumber': 17,
    'AtomicWeight': 35.453,
    'AtomicMass'  : '35.453u',
    'Color' : [255,255, 0],
    'Block' : 'p',
    'Group' : 17,
    'Period': 3,
    'CrystalStructure': 'Basic Centered Orthorhombic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [622.35, 445.61, 817.85],
    'SpaceGroupNumber' : 64,
    'SpaceGroup' : 'Cmca',
    'AtomicRadius' : 100,
    'CovalentRadius': 102,
    'VanDerWaalsRadius': 175 ,
    'ElectronConfiguration' : '[Ne]3s2,3p5',
    'ElectronConfigurationString' : '[Ne]3s2,3p5',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1251.2,
    'QuantumNumbers' : '2P3/2',
    'CASNumber' : '7782-50-5',
    'PubChemCID' : '24526',
    'Category' : ''
},
{
    'StandardName':'Argon',
    'Symbol': 'Ar',
    'Abbreviation' : 'ar',
    'AtomicNumber': 18,
    'AtomicWeight': 39.948,
    'AtomicMass'  : '39.948u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 18,
    'Period': 3,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [525.6, 525.6, 525.6],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 71,
    'CovalentRadius': 106,
    'VanDerWaalsRadius': 188 ,
    'ElectronConfiguration' : '[Ne]3s2,3p6',
    'ElectronConfigurationString' : '[Ne]3s2,3p6',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1520.6,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7440-37-1',
    'PubChemCID' : '23968',
    'Category' : ''
},
{
    'StandardName':'Potassium',
    'Symbol': 'K',
    'Abbreviation' : 'k',
    'AtomicNumber': 19,
    'AtomicWeight': 39.0983,
    'AtomicMass'  : '39.0983u',
    'Color' : [192,192,192],
    'Block' : 's',
    'Group' : 1,
    'Period': 4,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [532.8, 532.8, 532.8],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 220,
    'CovalentRadius': 203,
    'VanDerWaalsRadius': 275 ,
    'ElectronConfiguration' : '[Ar]4s1',
    'ElectronConfigurationString' : '[Ar]4s1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 418.8,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7440-09-7',
    'PubChemCID' : '5462222',
    'Category' : ''
},
{
    'StandardName':'Calcium',
    'Symbol': 'Ca',
    'Abbreviation' : 'ca',
    'AtomicNumber': 20,
    'AtomicWeight': 40.078,
    'AtomicMass'  : '40.078u',
    'Color' : [192,192,192],
    'Block' : 's',
    'Group' : 2,
    'Period': 4,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [558.84, 558.84, 558.84],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 180,
    'CovalentRadius': 176,
    'VanDerWaalsRadius': 231 ,
    'ElectronConfiguration' : '[Ar]4s2',
    'ElectronConfigurationString' : '[Ar]4s2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 589.8,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7440-70-2',
    'PubChemCID' : '5460341',
    'Category' : ''
},
{
    'StandardName':'Scandium',
    'Symbol': 'Sc',
    'Abbreviation' : 'sc',
    'AtomicNumber': 21,
    'AtomicWeight': 44.95591,
    'AtomicMass'  : '44.95591u',
    'Color' : [192,192,192],
    'Block' : 'd',
    'Group' : 3,
    'Period': 4,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2,PI/2,2*PI/3],
    'LatticeConstants': [330.9, 330.9, 527.33],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 160,
    'CovalentRadius': 170,
    'VanDerWaalsRadius': 211 ,
    'ElectronConfiguration' : '[Ar]4s2,3d1',
    'ElectronConfigurationString' : '[Ar]4s2,3d1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 633.1,
    'QuantumNumbers' : '2D3/2',
    'CASNumber' : '7440-20-2',
    'PubChemCID' : '23952',
    'Category' : ''
},
{
    'StandardName':'Titanium',
    'Symbol': 'Ti',
    'Abbreviation' : 'ti',
    'AtomicNumber': 22,
    'AtomicWeight': 47.867,
    'AtomicMass'  : '47.867u',
    'Color' : [192,192,192],
    'Block' : 'd',
    'Group' : 4,
    'Period': 4,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2,PI/2,2*PI/3],
    'LatticeConstants': [295.08, 295.08, 468.55],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 140,
    'CovalentRadius': 160,
    'VanDerWaalsRadius': None ,
    'ElectronConfiguration' : '[Ar]4s2,3d2',
    'ElectronConfigurationString' : '[Ar]4s2,3d2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 658.8,
    'QuantumNumbers' : '3F2',
    'CASNumber' : '7440-32-6',
    'PubChemCID' : '23963',
    'Category' : ''
},
{
    'StandardName':'Vanadium',
    'Symbol': 'V',
    'Abbreviation' : 'v',
    'AtomicNumber': 23,
    'AtomicWeight': 50.9415,
    'AtomicMass'  : '50.9415u',
    'Color' : [192,192,192],
    'Block' : 'd',
    'Group' : 5,
    'Period': 4,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [303, 303, 303],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 135,
    'CovalentRadius': 153,
    'VanDerWaalsRadius':None ,
    'ElectronConfiguration' : '[Ar]4s2,3d3',
    'ElectronConfigurationString' : '[Ar]4s2,3d3',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 650.9,
    'QuantumNumbers' : '4F3/2',
    'CASNumber' : '7440-62-2',
    'PubChemCID' : '23990',
    'Category' : ''
},
{
    'StandardName':'Chromium',
    'Symbol': 'Cr',
    'Abbreviation' : 'cr',
    'AtomicNumber': 24,
    'AtomicWeight': 51.9961,
    'AtomicMass'  : '51.9961u',
    'Color' : [192,192,192],
    'Block' : 'd',
    'Group' : 6,
    'Period': 4,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [291, 291, 291],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 140,
    'CovalentRadius': 166,
    'VanDerWaalsRadius':None ,
    'ElectronConfiguration' : '[Ar]4s2,3d5',
    'ElectronConfigurationString' : '[Ar]4s2,3d5',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 652.9,
    'QuantumNumbers' : '7S3',
    'CASNumber' : '7440-47-3',
    'PubChemCID' : '23976',
    'Category' : ''
},
{
    'StandardName':'Manganese',
    'Symbol': 'Mn',
    'Abbreviation' : 'mn',
    'AtomicNumber': 25,
    'AtomicWeight': 54.938045,
    'AtomicMass'  : '54.938045u',
    'Color' : [192,192,192],
    'Block' : 'd',
    'Group' : 7,
    'Period': 4,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [891.25, 891.25, 891.25],
    'SpaceGroupNumber' : 217,
    'SpaceGroup' : 'l_43m',
    'AtomicRadius' : 140,
    'CovalentRadius': 139,
    'VanDerWaalsRadius':None ,
    'ElectronConfiguration' : '[Ar]4s2,3d5',
    'ElectronConfigurationString' : '[Ar]4s2,3d5',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 717.3,
    'QuantumNumbers' : '6S5/2',
    'CASNumber' : '7439-96-5',
    'PubChemCID' : '23930',
    'Category' : ''
},
{
    'StandardName':'Iron',
    'Symbol': 'Fe',
    'Abbreviation' : 'fe',
    'AtomicNumber': 26,
    'AtomicWeight': 55.845,
    'AtomicMass'  : '55.845u',
    'Color' : [128, 128, 128],
    'Block' : 'd',
    'Group' : 8,
    'Period': 4,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [286.65, 286.65, 286.65],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 140,
    'CovalentRadius': 156,
    'VanDerWaalsRadius':None ,
    'ElectronConfiguration' : '[Ar]4s2,3d6',
    'ElectronConfigurationString' : '[Ar]4s2,3d6',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 762.5,
    'QuantumNumbers' : '5D4',
    'CASNumber' : '7439-89-6',
    'PubChemCID' : '23925',
    'Category' : ''
},
{
    'StandardName':'Cobalt',
    'Symbol': 'Co',
    'Abbreviation' : 'co',
    'AtomicNumber': 27,
    'AtomicWeight': 58.933195,
    'AtomicMass'  : '58.933195u',
    'Color' : [128, 128, 128],
    'Block' : 'd',
    'Group' : 9,
    'Period': 4,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2,PI/2,2*PI/3],
    'LatticeConstants': [250.71, 250.71, 406.95],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 135,
    'CovalentRadius': 126,
    'VanDerWaalsRadius':None ,
    'ElectronConfiguration' : '[Ar]4s2,3d7',
    'ElectronConfigurationString' : '[Ar]4s2,3d7',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 760.4,
    'QuantumNumbers' : '4F9/2',
    'CASNumber' : '7440-48-4',
    'PubChemCID' : '104730',
    'Category' : ''
},
{
    'StandardName':'Nickel',
    'Symbol': 'Ni',
    'Abbreviation' : 'ni',
    'AtomicNumber': 28,
    'AtomicWeight': 58.6934,
    'AtomicMass'  : '58.6934u',
    'Color' : [128, 128, 128],
    'Block' : 'd',
    'Group' : 10,
    'Period': 4,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [352.4, 352.4, 352.4 ],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 135,
    'CovalentRadius': 124,
    'VanDerWaalsRadius': 163 ,
    'ElectronConfiguration' : '[Ar]4s2,3d8',
    'ElectronConfigurationString' : '[Ar]4s2,3d8',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 737.1,
    'QuantumNumbers' : '3F4',
    'CASNumber' : '7440-02-0',
    'PubChemCID' : '935',
    'Category' : ''
},
{
    'StandardName':'Copper',
    'Symbol': 'Cu',
    'Abbreviation' : 'cu',
    'AtomicNumber': 29,
    'AtomicWeight': 63.546,
    'AtomicMass'  : '63.546u',
    'Color' : [184, 115, 51],
    'Block' : 'd',
    'Group' : 11,
    'Period': 4,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [361.49, 361.49, 361.49],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 135,
    'CovalentRadius': 138,
    'VanDerWaalsRadius': 140 ,
    'ElectronConfiguration' : '[Ar]4s1,3d10',
    'ElectronConfigurationString' : '[Ar]4s1,3d10',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 745.5,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7440-50-8',
    'PubChemCID' : '23978',
    'Category' : ''
},
{
    'StandardName':'Zinc',
    'Symbol': 'Zn',
    'Abbreviation' : 'zn',
    'AtomicNumber': 30,
    'AtomicWeight': 65.38,
    'AtomicMass'  : '65.38u',
    'Color' : [112, 128, 144],
    'Block' : 'd',
    'Group' : 12,
    'Period': 4,
    'CrystalStructure': 'Simple Hexagonal',
    'LatticeAngles' : [PI/2,PI/2,2*PI/3],
    'LatticeConstants': [266.49, 266.49, 494.68],
    'SpaceGroupNumber' : 194,
    'SpaceGroup' : 'P63/mmc',
    'AtomicRadius' : 135,
    'CovalentRadius': 122,
    'VanDerWaalsRadius': 139 ,
    'ElectronConfiguration' : '[Ar]4s2,3d10',
    'ElectronConfigurationString' : '[Ar]4s2,3d10',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 906.4,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7440-66-6',
    'PubChemCID' : '23994',
    'Category' : ''
},
{
    'StandardName':'Gallium',
    'Symbol': 'Ga',
    'Abbreviation' : 'ga',
    'AtomicNumber': 31,
    'AtomicWeight': 69.723,
    'AtomicMass'  : '69.723u',
    'Color' : [192, 192, 192],
    'Block' : 'p',
    'Group' : 13,
    'Period': 4,
    'CrystalStructure': 'Base Centered Orthorhombic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [451.97, 766.33, 452.6],
    'SpaceGroupNumber' : 64,
    'SpaceGroup' : 'Cmca',
    'AtomicRadius' : 130,
    'CovalentRadius': 122,
    'VanDerWaalsRadius': 187 ,
    'ElectronConfiguration' : '[Ar]4s2,3d10,4p1',
    'ElectronConfigurationString' : '[Ar]4s2,3d10,4p1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 578.8,
    'QuantumNumbers' : '2P1/2',
    'CASNumber' : '7440-55-3',
    'PubChemCID' : '5360835',
    'Category' : ''
},
{
    'StandardName':'Germanium',
    'Symbol': 'Ge',
    'Abbreviation' : 'ge',
    'AtomicNumber': 32,
    'AtomicWeight': 72.64,
    'AtomicMass'  : '72.64u',
    'Color' : [128, 128, 128],
    'Block' : 'p',
    'Group' : 14,
    'Period': 4,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2,PI/2,PI/2],
    'LatticeConstants': [565.75, 565.75, 565.75 ],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 125,
    'CovalentRadius': 120,
    'VanDerWaalsRadius': 211,
    'ElectronConfiguration' : '[Ar]4s2,3d10,4p2',
    'ElectronConfigurationString' : '[Ar]4s2,3d10,4p2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 762.0,
    'QuantumNumbers' : '3P0',
    'CASNumber' : '7440-56-4',
    'PubChemCID' : '6326954',
    'Category' : ''
},
{
    'StandardName':'Arsenic',
    'Symbol': 'As',
    'Abbreviation' : 'as',
    'AtomicNumber': 33,
    'AtomicWeight': 74.9216,
    'AtomicMass'  : '74.9216u',
    'Color' : [192, 192, 192],
    'Block' : 'p',
    'Group' : 15,
    'Period': 4,
    'CrystalStructure': 'Simple Trigonal',
    'LatticeAngles' : [PI/2,PI/2,2*PI/3],
    'LatticeConstants': [375.98, 375.98, 1054.75 ],
    'SpaceGroupNumber' : 166,
    'SpaceGroup' : 'R_3m',
    'AtomicRadius' : 115,
    'CovalentRadius': 119,
    'VanDerWaalsRadius': 185,
    'ElectronConfiguration' : '[Ar]4s2,3d10,4p3',
    'ElectronConfigurationString' : '[Ar]4s2,3d10,4p3',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 947.0,
    'QuantumNumbers' : '4S3/2',
    'CASNumber' : '7440-38-2',
    'PubChemCID' : '5359596',
    'Category' : ''
},
{
    'StandardName':'Selenium',
    'Symbol': 'Se',
    'Abbreviation' : 'se',
    'AtomicNumber': 34,
    'AtomicWeight': 78.96,
    'AtomicMass'  : '78.96u',
    'Color' : [128, 128, 128],
    'Block' : 'p',
    'Group' : 16,
    'Period': 4,
    'CrystalStructure': 'Simple Monoclinic',
    'LatticeAngles' : [PI/2, 1.58493, PI/2],
    'LatticeConstants': [905.4, 908.3, 1160.1],
    'SpaceGroupNumber' : 14,
    'SpaceGroup' : 'P121/c1',
    'AtomicRadius' : 115,
    'CovalentRadius': 120,
    'VanDerWaalsRadius': 190,
    'ElectronConfiguration' : '[Ar]4s2,3d10,4p4',
    'ElectronConfigurationString' : '[Ar]4s2,3d10,4p4',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 941.0,
    'QuantumNumbers' : '3P2',
    'CASNumber' : '7782-49-2',
    'PubChemCID' : '6326970',
    'Category' : ''
},
{
    'StandardName':'Bromine',
    'Symbol': 'Br',
    'Abbreviation' : 'br',
    'AtomicNumber': 35,
    'AtomicWeight': 79.904,
    'AtomicMass'  : '79.904u',
    'Color' : [255, 0, 0],
    'Block' : 'p',
    'Group' : 17,
    'Period': 4,
    'CrystalStructure': 'Base Centered Orthorhombic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [672.65, 464.51, 870.23],
    'SpaceGroupNumber' : 64,
    'SpaceGroup' : 'Cmca',
    'AtomicRadius' : 115,
    'CovalentRadius': 120,
    'VanDerWaalsRadius': 185,
    'ElectronConfiguration' : '[Ar]4s2,3d10,4p5',
    'ElectronConfigurationString' : '[Ar]4s2,3d10,4p5',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1139.9,
    'QuantumNumbers' : '2P3/2',
    'CASNumber' : '7726-95-6',
    'PubChemCID' : '24408',
    'Category' : ''
},
{
    'StandardName':'Krypton',
    'Symbol': 'Kr',
    'Abbreviation' : 'kr',
    'AtomicNumber': 36,
    'AtomicWeight': 83.798,
    'AtomicMass'  : '83.798u',
    'Color' : [],
    'Block' : 'p',
    'Group' : 18,
    'Period': 4,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [570.6, 570.6, 570.6],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 88,
    'CovalentRadius': 116,
    'VanDerWaalsRadius': 202,
    'ElectronConfiguration' : '[Ar]4s2,3d10,4p6',
    'ElectronConfigurationString' : '[Ar]4s2,3d10,4p6',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 1350.8,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7726-95-6',
    'PubChemCID' : '24408',
    'Category' : ''
},
{
    'StandardName':'Rubidium',
    'Symbol': 'Rb',
    'Abbreviation' : 'rb',
    'AtomicNumber': 37,
    'AtomicWeight': 85.4678,
    'AtomicMass'  : '85.4678u',
    'Color' : [192,192, 192],
    'Block' : 's',
    'Group' : 1,
    'Period': 5,
    'CrystalStructure': 'Body Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [558.5, 558.5, 558.5],
    'SpaceGroupNumber' : 229,
    'SpaceGroup' : 'lm_3m',
    'AtomicRadius' : 235,
    'CovalentRadius': 220,
    'VanDerWaalsRadius': 303,
    'ElectronConfiguration' : '[Kr]5s1',
    'ElectronConfigurationString' : '[Kr]5s1',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 403.0,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7440-17-7',
    'PubChemCID' : '5357696',
    'Category' : ''
},
{
    'StandardName':'Strontium',
    'Symbol': 'Sr',
    'Abbreviation' : 'sr',
    'AtomicNumber': 38,
    'AtomicWeight': 87.62,
    'AtomicMass'  : '87.62u',
    'Color' : [192,192, 192],
    'Block' : 's',
    'Group' : 2,
    'Period': 5,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [608.49, 608.49, 608.49],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 200,
    'CovalentRadius': 195,
    'VanDerWaalsRadius': 249,
    'ElectronConfiguration' : '[Kr]5s2',
    'ElectronConfigurationString' : '[Kr]5s2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 549.5,
    'QuantumNumbers' : '1S0',
    'CASNumber' : '7440-24-6',
    'PubChemCID' : '5359327',
    'Category' : ''
},
{
    'StandardName':'Silver',
    'Symbol': 'Ag',
    'Abbreviation' : 'ag',
    'AtomicNumber': 47,
    'AtomicWeight': 107.8682,
    'AtomicMass'  : '107.8682u',
    'Color' : [192,192, 192],
    'Block' : 'd',
    'Group' : 11,
    'Period': 5,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [408.53, 408.53, 408.53],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 160,
    'CovalentRadius': 145,
    'VanDerWaalsRadius': 172,
    'ElectronConfiguration' : '[Kr]5s1,4d10',
    'ElectronConfigurationString' : '[Kr]5s1,4d10',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 731.0,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7440-22-4',
    'PubChemCID' : '23954',
    'Category' : ''
},

{
    'StandardName':'Platinum',
    'Symbol': 'Pt',
    'Abbreviation' : 'pt',
    'AtomicNumber': 78,
    'AtomicWeight': 195.084,
    'AtomicMass'  : '195.084u',
    'Color' : [128, 128, 128],
    'Block' : 'd',
    'Group' : 10,
    'Period': 6,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [392.42, 392.42, 392.42],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 135,
    'CovalentRadius': 136,
    'VanDerWaalsRadius': 175,
    'ElectronConfiguration' : '[Xe]6s1,4f14,5d9',
    'ElectronConfigurationString' : '[Xe]6s1,4f14,5d9',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 870.0,
    'QuantumNumbers' : '3D3',
    'CASNumber' : '7440-06-4',
    'PubChemCID' : '23939',
    'Category' : ''
},

{
    'StandardName':'Gold',
    'Symbol': 'Au',
    'Abbreviation' : 'au',
    'AtomicNumber': 79,
    'AtomicWeight': 196.966569,
    'AtomicMass'  : '196.966569u',
    'Color' : [255, 215, 0],
    'Block' : 'd',
    'Group' : 11,
    'Period': 6,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [407.82, 407.82, 407.82],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 135,
    'CovalentRadius': 136,
    'VanDerWaalsRadius': 166,
    'ElectronConfiguration' : '[Xe]6s1,4f14,5d10',
    'ElectronConfigurationString' : '[Xe]6s1,4f14,5d10',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 890.1,
    'QuantumNumbers' : '2S1/2',
    'CASNumber' : '7440-57-5',
    'PubChemCID' : '23985',
    'Category' : ''
},

{
    'StandardName':'Lead',
    'Symbol': 'Pb',
    'Abbreviation' : 'pb',
    'AtomicNumber': 82,
    'AtomicWeight': 207.2,
    'AtomicMass'  : '207.2u',
    'Color' : [112, 128, 144],
    'Block' : 'p',
    'Group' : 14,
    'Period': 6,
    'CrystalStructure': 'Face Centered Cubic',
    'LatticeAngles' : [PI/2, PI/2, PI/2],
    'LatticeConstants': [495.08, 495.08, 495.08],
    'SpaceGroupNumber' : 225,
    'SpaceGroup' : 'Fm_3m',
    'AtomicRadius' : 180,
    'CovalentRadius': 146,
    'VanDerWaalsRadius': 202,
    'ElectronConfiguration' : '[Xe]6s2,4f14,5d10,6p2',
    'ElectronConfigurationString' : '[Xe]6s2,4f14,5d10,6p2',
    'ElectronShellConfiguration': '',
    'IonizationEnergies' : 715.6,
    'QuantumNumbers' : '3P0',
    'CASNumber' : '7439-92-1',
    'PubChemCID' : '5352425',
    'Category' : ''
}

]

ElementData = ElementCollection()

for ele in PeriodicTable:
    ele_obj = Element()
    for key, val in ele.iteritems():
        ele_obj.setData(key, [val])
    ElementData.add(ele_obj)
