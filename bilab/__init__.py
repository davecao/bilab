# -*- coding: utf-8 -*-

__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2014/08/04"
__version__ = '1.0.0'

import os.path
import sys
import warnings
import sysconfig

if sys.version_info[:2] < (2, 7):
    raise Exception('Python version less than 2.7')

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is a required package for ProDy')
else:
    if tuple(map(int, np.__version__.split('.')[:2])) < (1, 8):
        raise ImportError('Numpy v1.8 or later is required.')

from .utilities import PackageLogger, PackageSettings
from .utilities import getPackagePath, joinRepr, tabulate

_PY3K = PY3K = sys.version_info[0] > 2
PY2K = not PY3K

LOGGER = PackageLogger('bilab')
SETTINGS = PackageSettings('bilab', logger=LOGGER)
SETTINGS.load()

# get path of resources, i.e., aaindex1, aadindex2 and aaindex3
#data = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'data'
data = "{}{}{}{}{}".format(
        sysconfig.get_config_var('datarootdir'), 
        os.sep, 'bilab', os.sep, 'data')

__all__ = []

#from . import chemicals
from .chemicals import *
#__all__.extend(chemicals.__all__)
#__all__.append('chemicals')

#from . import utilities
from .utilities import *
#__all__.extend(utilities.__all__)
#__all__.append('utilities')

#from .utilities import PackageLogger
#LOGGER = PackageLogger('.bilab')

#from . import sequence
from .sequence import *
#__all__.extend(sequence.__all__)
#__all__.append('sequence')

#from . import webservice
from .webservice import *
#__all__.extend(webservice.__all__)
#__all__.append('webservice')

#from . import aaprop
from .aaprop import *
#__all__.extend(aaprop.__all__)
#__all__.append('aaprop')

#from . import io
from .io import *
#__all__.extend(io.__all__)
#__all__.append('io')

#from . import ml
from .ml import *
#__all__.extend(io.__all__)
#__all__.append('io')

#from . import optimization
from .optimization import *
#__all__.extend(optimization.__all__)
#__all__.append('optimization')

#from . import optimization
from .geometry import *
#__all__.extend(geometry.__all__)
#__all__.append('geometry')

#LOGGER = PackageLogger('.bilab')
#from . import structure
from .structure import *
#__all__.extend(structure.__all__)
#__all__.append('structure')

#from . import concurrent
from .concurrent import *
#__all__.extend(concurrent.__all__)
#__all__.append('concurrent')


# Try to load bilab
    #scriptDir = os.path.dirname(os.path.realpath(__file__))
    #covalent = import_module("covalent_radius",
    #                         path=scriptDir + os.sep + 'bilab',
    #                         verbose=True)
    #utilities = import_module("utilities",
    #                         path=scriptDir + os.sep + 'bilab',
    #                         verbose=True)
    #bilab = import_module("bilab",
    #                      path = scriptDir + os.sep + 'bilab',
    #                      verbose = opt.verbose)
    # Try to load prody
    #prody = import_module("prody",
    #                      path=opt.prody_path,
    #                      verbose = opt.verbose)
    #print(dir(bilab))
    #bilab.utilities.get_loaded_modules()
# add the following two lines
# so that import_module could add it to global space

import bilab
__all__.append('bilab')

#LOGGER = bilab.utilities.PackageLogger('.bilab')

# prody path
#def import_prody():
#    prody_path = os.path.dirname(os.path.realpath(__file__))
#    print("prody_path: {}".format(prody_path))
#    bilab.utilities.import_module("prody", path=prody_path, verbose = True)
#    try:
#        import prody
#    except ImportError:
#        raise ImportError('Prody is a required package for bilab')
#
#import_prody()
#from . import prody
#from .prody import *


