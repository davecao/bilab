# -*- coding: utf-8 -*-

from numpy import array

# Literature:
#   Beatriz et al. Covalent radii revisited, Dalton Trans(21):2832-2838.
#      doi:10.1039/b801115j Table 2 1-96
__all__ = ['get_covalent_radius']

READONLY = set()

DTYPE = array(['a']).dtype.char  # 'S' for PY2K and 'U' for PY3K

class ElementProperties(object):
    """Element data field."""

    __slots__ = ['name', 'dtype',  'doc', 'doc_pl', 'meth', 'meth_pl',
                 'ndim', 'none', 'selstr', 'synonym', 'readonly', 'call',
                 'private', 'depr', 'depr_pl', 'desc', 'flags']

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
        #: expected dimension of the data array
        self.ndim = kwargs.get('ndim', 1)
        #: atomic get/set method name
        self.meth = kwargs.get('meth', name.capitalize())
        #: get/set method name in plural form
        self.meth_pl = kwargs.get('meth_pl', self.meth + 's')
        #: :class:`.AtomGroup` attributes to be set None, when ``setMethod``
        #: is called
        self.none = kwargs.get('none')
        #: list of selection string examples
        self.selstr = kwargs.get('selstr')
        #: deprecated method name
        self.depr = kwargs.get('depr')
        #: deprecated method name in plural form
        self.depr_pl = None
        if self.depr is not None:
            self.depr_pl = kwargs.get('depr_pl', self.depr + 's')
        #: synonym used in atom selections
        self.synonym = kwargs.get('synonym')
        #: read-only attribute without a set method
        self.readonly = kwargs.get('readonly', False)
        #: list of :class:`.AtomGroup` methods to call when ``getMethod`` is
        #: called
        self.call = kwargs.get('call', None)
        #: define only _getMethod for :class:`.AtomGroup` to be used by
        #: :class:`.Select` class
        self.private = kwargs.get('private', False)
        #: **True** when there are flags associated with the data field
        self.flags = kwargs.get('flags', False)

        if self.readonly:
            READONLY.add(self.name)
class ElementData(object):
    """ Checmical elements data in periodic table """

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
