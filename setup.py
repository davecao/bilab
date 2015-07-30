# This script was automatically generated by distutils2
import codecs
import io, os, sys
import glob
from distutils.core import setup,Extension
from distutils.sysconfig import get_config_var

#from distutils.command.install_data import install_data

try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser
try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is a required package')

current_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
#install_data = distutils.command.install_data.install_data
#print("Install data:{}".format(install_data.install_dir))
np_include_dirs = np.__path__[0] + '/core/include'

_UNWANTED_OPTS = frozenset(['-Wstrict-prototypes'])
os.environ['OPT'] = ' '.join(
    _ for _ in get_config_var('OPT').strip().split() if _ not in _UNWANTED_OPTS)
#boost_root = None
#if "BOOST_ROOT" in os.environ:
#    boost_root = os.environ["BOOST_ROOT"]
#else:
#    print("Could not find BOOST_ROOT; terminated!!!")
#    sys.exit(0)

#boost_stage_lib = boost_root + os.path.sep + "stage" + os.path.sep + "lib"
Boost_LIB_PATH = None
print os.environ
if 'BOOST_LIBRARYDIR' in os.environ:
    Boost_LIB_PATH = os.environ['BOOST_LIBRARYDIR']
else:
    print("Could not find BOOST_LIBRARYDIR; terminated!!!")
    sys.exit(0)


distance_wrap = Extension('bilab.geometry._distance_wrap',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = [np_include_dirs],
                    #libraries = [''],
                    #library_dirs = [''],
                    sources = ['bilab/geometry/distance/distance_wrap.c'])
kdtree_lib = Extension('bilab.structure._CKDTree',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = [np_include_dirs],
                    #libraries = [''],
                    #library_dirs = [''],
                    sources = ['bilab/structure/kdtree/KDTree.c', 
                        'bilab/structure/kdtree/KDTreemodule.c'])

CWMatrix = Extension('bilab.linalg._CWMatrix',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    #include_dirs = [boost_root],
                    libraries = ['boost_python'],
                    #library_dirs = [boost_stage_lib],
                    extra_compile_args = ['-ftemplate-backtrace-limit=64'],
                    sources = ['bilab/linalg/matrix/export.cpp']) 
                               #'bilab/linalg/matrix/matrix.cpp'])

bhtsne_wrap = Extension('bilab.ml.NDR.tSNE._bhtsne_wrap',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = [np_include_dirs],
                    library_dirs = [Boost_LIB_PATH],
                    libraries = ['boost_python'],
                    extra_compile_args = ['-ftemplate-backtrace-limit=64', 
                                          '-std=c++11'],
                    sources = ['bilab/ml/NDR/tSNE/bhtsne_export.cpp',
                               'bilab/ml/NDR/tSNE/bhtsne.cpp',
                               'bilab/ml/NDR/tSNE/sptree.cpp'])

def split_multiline(value):
    """Split a multiline string into a list, excluding blank lines."""
    val = value.encode('utf-8')
    return [element for element in
            (line.strip() for line in val.split('\n'))
            if element]


def cfg_to_args(path='setup.cfg'):
    """Compatibility helper to use setup.cfg in setup.py.

    This functions uses an existing setup.cfg to generate a dictionnary of
    keywords that can be used by distutils.core.setup(**kwargs).  It is used
    by generate_setup_py.

    *file* is the path to the setup.cfg file.  If it doesn't exist,
    PackagingFileError is raised.
    """

    # XXX ** == needs testing
    D1_D2_SETUP_ARGS = {"name": ("metadata",),
                        "version": ("metadata",),
                        "author": ("metadata",),
                        "author_email": ("metadata",),
                        "maintainer": ("metadata",),
                        "maintainer_email": ("metadata",),
                        "url": ("metadata", "home_page"),
                        "description": ("metadata", "summary"),
                        "long_description": ("metadata", "description"),
                        "download-url": ("metadata",),
                        "classifiers": ("metadata", "classifier"),
                        "platforms": ("metadata", "platform"),  # **
                        "license": ("metadata",),
                        "requires": ("metadata", "requires_dist"),
                        "provides": ("metadata", "provides_dist"),  # **
                        "obsoletes": ("metadata", "obsoletes_dist"),  # **
                        "package_dir": ("files", 'packages_root'),
                        "packages": ("files",),
                        "scripts": ("files",),
                        "resources": ("files",),
                        "py_modules": ("files", "modules"),  # **
                        }

    MULTI_FIELDS = ("classifiers",
                    "platforms",
                    "requires",
                    "provides",
                    "obsoletes",
                    "packages",
                    "scripts",
                    "py_modules",
                    "extension",
                    "resources"
                    )

    def has_get_option(config, section, option):
        if config.has_option(section, option):
            return config.get(section, option)
        elif config.has_option(section, option.replace('_', '-')):
            return config.get(section, option.replace('_', '-'))
        else:
            return False

    # The real code starts here
    config = RawConfigParser()
    #f = codecs.open(path, encoding='utf-8')
    f = io.open(path,encoding="utf-8")
    try:
        config.readfp(f)
    finally:
        f.close()

    kwargs = {}
    for arg in D1_D2_SETUP_ARGS:
        if len(D1_D2_SETUP_ARGS[arg]) == 2:
            # The distutils field name is different than distutils2's
            section, option = D1_D2_SETUP_ARGS[arg]

        else:
            # The distutils field name is the same thant distutils2's
            section = D1_D2_SETUP_ARGS[arg][0]
            option = arg

        in_cfg_value = has_get_option(config, section, option)
        if not in_cfg_value:
            # There is no such option in the setup.cfg
            if arg == 'long_description':
                filenames = has_get_option(config, section, 'description-file')
                if filenames:
                    filenames = split_multiline(filenames)
                    in_cfg_value = []
                    for filename in filenames:
                        #fp = codecs.open(filename, encoding='utf-8')
                        fp = io.open(filename, encoding="utf-8")
                        try:
                            in_cfg_value.append(fp.read())
                        finally:
                            fp.close()
                    in_cfg_value = '\n\n'.join(in_cfg_value)
            else:
                continue

        if arg == 'package_dir' and in_cfg_value:
            in_cfg_value = {'': in_cfg_value}

        if arg in MULTI_FIELDS:
            # support multiline options
            in_cfg_value = split_multiline(in_cfg_value)
            if arg == 'resources' and in_cfg_value:
                datafiles = []
                for line in in_cfg_value:
                    file_str = line.split('=')[0].strip()
                    datafiles.extend(glob.glob(file_str))
                #kwargs['data_files'] = [('share/bilab',[line.split('=')[0].strip() for line in in_cfg_value])]
                kwargs['data_files'] = [('share/bilab/data',datafiles)]

        kwargs[arg] = in_cfg_value

    return kwargs

# setup
general_settings = cfg_to_args()
# Key: resources has to be removed
general_settings.pop('resources')
#extensions
general_settings['ext_modules'] = [distance_wrap, kdtree_lib, bhtsne_wrap ]
#                                   CWMatrix]
setup(**general_settings)
#setup(**cfg_to_args())
