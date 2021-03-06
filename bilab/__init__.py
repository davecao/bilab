# -*- coding: utf-8 -*-
import os
import sys
import types

import sysconfig
from datetime import datetime

from jinja2 import Environment, FileSystemLoader

if sys.version_info[:2] < (2, 7):
    raise Exception('Python version less than 2.7')

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is a required package')
else:
    if tuple(map(int, np.__version__.split('.')[:2])) < (1, 8):
        raise ImportError('Numpy v1.8 or later is required.')

from .utilities import PackageLogger, PackageSettings
# from .utilities import getPackagePath, joinRepr, tabulate
from jinja2 import Environment, FileSystemLoader

def is_float(n):
    return isinstance(n, float)

TPLPATH = "{}/templates/".format(os.path.dirname(__file__))
jinja2_ENV = Environment(loader=FileSystemLoader(TPLPATH))
jinja2_ENV.add_extension('jinja2.ext.loopcontrols')
jinja2_ENV.tests['Float'] = is_float

__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2014/08/04"
__version__ = '1.0.0'

_TIMESTAMP_FORMAT = '%Y-%m-%d-%H:%M:%S'
_PY3K = PY3K = sys.version_info[0] > 2
PY2K = not PY3K

def is_float(n):
    return isinstance(n, float)

TPLPATH = "{}/templates/".format(os.path.dirname(__file__))
jinja2_ENV = Environment(loader=FileSystemLoader(TPLPATH))
jinja2_ENV.add_extension('jinja2.ext.loopcontrols')
jinja2_ENV.tests['Float'] = is_float
VERSION = __version__

LOGGER = PackageLogger('bilab')
SETTINGS = PackageSettings('bilab', logger=LOGGER)
#SETTINGS.load()

# get path of resources, i.e., aaindex1, aadindex2 and aaindex3
# data = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'data'
# On ubuntu: sysconfig.get_config_var("datarootdir")
#             -- /usr/share
#            sysconfig.get_path("data")
#             -- /usr/local
#     Point to different locations
#  using sysconfig.get_path("data") instead
#             -- /usr/local
# if using egg, have to install data
#  i.e., python setup.py install_data
#    -- /usr/local/lib/python2.7/dist-packages/bilab-1.0.0-py2.7-linux-x86_64.egg/
#      bilab/data/
data = "{}{}{}{}{}".format(
            sysconfig.get_path('data'),
            os.sep, 'share', os.sep, 'bilab', os.sep, 'data')

ff_location = data + os.sep + "data" + os.sep + 'ff' + os.sep

forcefieldsList = {
    'amber91': ff_location + 'amber_parm91',
    'amber94': ff_location + 'amber_parm94',
    'amber99': ff_location + 'amber_parm99',
    'opls':    ff_location + 'opls_parm',
}
__all__ = []

if PY3K:
    from copyreg import __newobj__ as reduce_newobj
    string_types = str,
    integer_types = int,
    class_types = type,
    text_type = str,
    binary_type = bytes
else:
    from copy_reg import __newobj__ as reduce_newobj
    string_types = basestring,
    integer_types = (int, long),
    class_types = (type, types.ClassType),
    text_type = unicode,
    binary_type = str


def get_timestamp():
    return datetime.now().strftime(_TIMESTAMP_FORMAT)

# from . import chemicals
from .chemicals import *
# __all__.extend(chemicals.__all__)
# __all__.append('chemicals')

# from . import utilities
from .utilities import *
# __all__.extend(utilities.__all__)
# __all__.append('utilities')

# from . import sequence
from .sequence import *
# __all__.extend(sequence.__all__)
# __all__.append('sequence')

# from . import webservice
from .webservice import *
# __all__.extend(webservice.__all__)
# __all__.append('webservice')

# from . import aaprop
from .aaprop import *
# __all__.extend(aaprop.__all__)
# __all__.append('aaprop')

# from . import io
from .io import *
# __all__.extend(io.__all__)
# __all__.append('io')

# from . import ml
from .ml import *
# __all__.extend(io.__all__)
# __all__.append('io')

# from . import optimization
from .optimization import *
# __all__.extend(optimization.__all__)
# __all__.append('optimization')

# from . import optimization
from .geometry import *
# __all__.extend(geometry.__all__)
# __all__.append('geometry')

from .graph import *
# __all__.extend(graph.__all__)
# __all__.append('graph')


# from . import structure
from .structure import *
# __all__.extend(structure.__all__)
# __all__.append('structure')

# from . import concurrent
from .concurrent import *
# __all__.extend(concurrent.__all__)
# __all__.append('concurrent')

from .apps import *
