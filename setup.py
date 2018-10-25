# This script was automatically generated by distutils2
# import codecs
import sys
import os
import io
import re
import platform as plat

import distutils.sysconfig as dsc
# from distutils import core, dir_util
from setuptools import setup, Extension
# from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var

## scipy_distutils Script
# from numpy.distutils.core import setup, Extension
# from numpy.distutils.command import build_src
import Cython
import Cython.Compiler.Main
from Cython.Build import cythonize
from Cython.Distutils import build_ext

# try:
#     from Cython.Distutils import build_ext
#     from Cython.Build import cythonize
# except ImportError:
#     from distutils.command import build_ext
#  setuptools --> pip install_requires
#  from setuptools import setup, Extension

# build_src.Pyrex = Cython
# build_src.have_pyrex = True

# from distutils.ccompiler import show_compilers

from glob import glob
pj = os.path.join

_PY3K = PY3K = sys.version_info[0] > 2

# from distutils.command.install_data import install_data

try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

# include_dirs = []
# library_dirs = []

pyincdir = dsc.get_python_inc(plat_specific=1)
pylibdir = os.path.join('/', *pyincdir.split('/')[:-2] + ['lib'])

try:
    import numpy as np
    np_include_dir = np.get_include()
    # include_dirs.append(np_include_dir)
except ImportError:
    np_include_dir = os.path.join(
                pylibdir.replace('lib/python', 'local/lib/python'),
                'numpy', 'core', 'include')
    print("Unable to import numpy, trying header \n{}".format(
            np_include_dir))
    # raise ImportError('Numpy is a required package')


class FindBoostPackage(object):
    """
    Find the include path and libraries path of Boost library
    """
    def __init__(self, baseDirKey="BOOST_ROOT", minimal_version="1.55.0"):
        self.baseDirKey = baseDirKey
        self.minimal_version = minimal_version
        self.staged_version = False
        self.available = False            # Track availability
        # Need boost python
        self.libext = '.dylib' if sys.platform == 'darwin' else '.so'
        self.found_incs = None
        self.found_lib_paths = None

        # configurable options
        self.baseDir = os.environ[self.baseDirKey]
        self.incDir = None
        self.setupLibrarySettings()

    def setupLibrarySettings(self):
        # Map from library name to header to check for
        self.headerMap = {
            'program_options': 'boost/program_options.hpp',
            'python': 'boost/python.hpp',
            'thread': 'boost/thread.hpp',
            'filesystem': 'boost/filesystem/path.hpp'
        }

        # Map for extra libs needed for config test
        self.extraEnvOptions = {}

        # ---
        #   Build up settings using distutils.sysconfig to
        #   get Python build options
        # ---
        # distutils.sysconfig.get_config_vars()
        # self.python_version = distutils.sysconfig.get_python_version()
        # self.python_version = distutils.sysconfig.get_config_var("VERSION")
        # self.python_inc_dir = distutils.sysconfig.get_python_inc()
        # python_link_share_flags = distutils.sysconfig.get_config_var(
        #                               'LINKFORSHARED')
        # self.python_link_share_flags = "-Wl,-export-dynamic"
        # self.python_lib_path =\
        #   distutils.sysconfig.get_python_lib(standard_lib=True) + "/config"
        # self.python_extra_libs = [
        #           "python"+self.python_version,
        #           "util", "pthread", "dl"]
        # See SHLIBS

        if not self.baseDir:
            print("Boost base dir not set")
            return
        if not os.path.isdir(self.baseDir):    # If we don't have a directory
            print("Boost base dir is not a directory: {}".format(self.baseDir))
            return
        # Try to find inc dir
        if not self.incDir:
            self.incDir = pj(self.baseDir, 'include')
            if os.path.isdir(pj(self.baseDir, 'include')):
                self.incDir = pj(self.baseDir, 'include')
            elif os.path.isdir(pj(self.baseDir)):
                #  stage version?
                self.incDir = pj(self.baseDir)
            if not os.path.isfile(pj(self.incDir, 'boost', 'version.hpp')):
                print("Searching for correct boost include dir...")
                potential_dirs = os.listdir(self.incDir)
                potential_dirs = [d for d in potential_dirs if os.path.isfile(
                                        pj(self.incDir, d, 'boost',
                                            'version.hpp'))]
                potential_dirs.sort()
                if 0 == len(potential_dirs):
                    print("none found.")
                else:
                    self.incDir = pj(self.incDir, potential_dirs[-1])
                    print("found: {}".format(self.incDir))
        if self.incDir and (not os.path.isdir(pj(self.incDir, 'boost'))):
            print(
                "Boost inc dir is not a valid directory: {}".format(
                    pj(self.incDir, 'boost')))
            return
        # Check the version header is there
        version_header = pj(self.incDir, 'boost', 'version.hpp')
        if not os.path.isfile(version_header):
            print("Boost version.hpp header does not exist:{}".format(
                    version_header))
            return
        print("   boost header path: {}".format(self.incDir))
        # --- Check version requirement --- #
        if PY3K:
            ver_file = open(version_header)
        else:
            ver_file = file(version_header)
        ver_match = re.search(
                    "define\s+?BOOST_VERSION\s+?(\d*)", ver_file.read())
        if not ver_match:
            print("could not find BOOST_VERSION in file: {}".format(
                    version_header))
            return
        found_ver_str = int(ver_match.group(1))
        found_ver_str = "{}.{}.{}".format(
            found_ver_str/100000,
            found_ver_str/100 % 1000,
            found_ver_str % 100)

        req_ver = [int(n) for n in self.minimal_version.split('.')]
        found_ver = [int(n) for n in found_ver_str.split('.')]
        print("   boost version:", ".".join([str(x) for x in found_ver]))
        if found_ver < req_ver:
            print("found version is to old: required:{} found:{}" % (
                    self.minimal_version, found_ver_str))
            return
        # Set lists of the options we want
        self.found_incs = self.incDir
        if os.path.isdir(pj(self.baseDir, 'lib')):
            self.found_lib_paths = pj(self.baseDir, 'lib')
        elif os.path.isdir(pj(self.baseDir, 'stage/lib')):
            self.staged_version = True
            self.found_lib_paths = pj(self.baseDir, 'stage/lib')

    def dumpSettings(self):
        # Write out the settings
        print("BoostBaseDir: {}".format(self.baseDir))
        print("BoostIncludeDir:{}".format(self.incDir))
        print("LIBPATH: {}".format(self.found_lib_paths))
        # print("Python settings")
        # print("   inc : {}".format(self.python_inc_dir))
        # print("   link: {}".format(self.python_link_share_flags))
        # print("   lib : {}".format(self.python_extra_libs))
        # print("   lib path: {}".format(self.python_lib_path))

# Remove unwanted options
_UNWANTED_OPTS = frozenset(['-Wstrict-prototypes'])
os.environ['OPT'] = ' '.join(
    _ for _ in get_config_var('OPT').strip().split()
    if _ not in _UNWANTED_OPTS)

fboost = FindBoostPackage()
fboost.dumpSettings()
boost_inc_dir = fboost.found_incs
boost_lib_dir = fboost.found_lib_paths
netcdf_inc_dir = '/opt/local/include'
# if 'Boost_LIBRARY_DIR' in os.environ:
#    Boost_LIBRARY_PATH = os.environ['Boost_LIBRARY_DIR']
# else:
#    print("Could not find Boost_LIBRARY_DIR; terminated!!!")
#    sys.exit(0)

basic_extra_compile_args = ['-ftemplate-backtrace-limit=64','-m64','-std=c++11']
basic_extra_link_args = []
# clang -fopenmp -I <path to omp.h> -L <LLVM OpenMP library path>
lfdfiles = ""
if plat.system() == "Darwin":
    # clang on MacOS
    basic_extra_compile_args += ['-stdlib=libc++']
    basic_extra_link_args += ['-stdlib=libc++']

    lfdfiles = Extension(
        'bilab.io._lfdfiles',
        sources=['bilab/io/_lfdfiles.pyx'],
        include_dirs=[np_include_dir])
elif plat.system() == "Linux":
    # gcc on Linux
    basic_extra_compile_args += ['-static-libstdc++']
    basic_extra_link_args += ['-static-libstdc++']

    lfdfiles = Extension(
        'bilab.io._lfdfiles',
        sources=['bilab/io/_lfdfiles.pyx'],
        include_dirs=[np_include_dir],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'])

# could not compiled by distutils
helanal = Extension(
    name='bilab.geometry._helanal',
    sources=[
        'bilab/geometry/helanal.f95',
        'bilab/geometry/helanal.pyf'])

voroplusplus = Extension(
    name='bilab.geometry.voro.voroplusplus',
    include_dirs=['bilab/geometry/voro/src'],
    sources=[
        'bilab/geometry/voro/voroplusplus.pyx',
        'bilab/geometry/voro/vpp.cpp',
        'bilab/geometry/voro/src/voro++.cc'],
    extra_compile_args=basic_extra_compile_args,
    extra_link_args=basic_extra_link_args,
    language="c++")

crc32 = Extension(
    name='bilab.utilities.crc32',
    include_dirs=[np_include_dir],
    sources=[
        'bilab/utilities/_crc32.pyx'],
    language="c")

distance_wrap = Extension(
    name='bilab.geometry._distance_wrap',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0')],
    include_dirs=[np_include_dir],
    # libraries = [''],
    # library_dirs = [''],
    sources=['bilab/geometry/distance/distance_wrap.c'])

kdtree_lib = Extension(
    name='bilab.structure._CKDTree',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0')],
    include_dirs=[np_include_dir],
    # libraries = [''],
    # library_dirs = [''],
    sources=['bilab/structure/kdtree/KDTree.c',
             'bilab/structure/kdtree/KDTreemodule.c'])

marching_cubes = Extension(
    name="bilab.geometry.isosurface._marching_cubes_cy",
    include_dirs=[np_include_dir],
    sources=[
        'bilab/geometry/isosurface/_marching_cubes_cy.pyx'
    ])

gpmetis_wrap = Extension(
    name="bilab.graph.metis._gpmetis_wrap",
    define_macros=[('MAJOR_VERSION', '5'),
                   ('MINOR_VERSION', '1')],
    include_dirs=['.',
                  np_include_dir,
                  'bilab/graph/metis',
                  'bilab/graph/metis/GKlib',
                  'bilab/graph/metis/libmetis'
                  ],
    sources=[
        #'bilab/graph/metis/GKlib/conf/check_thread_storage.c',
        'bilab/graph/metis/GKlib/b64.c',
        'bilab/graph/metis/GKlib/blas.c',
        'bilab/graph/metis/GKlib/csr.c',
        'bilab/graph/metis/GKlib/error.c',
        'bilab/graph/metis/GKlib/evaluate.c',
        'bilab/graph/metis/GKlib/fkvkselect.c',
        'bilab/graph/metis/GKlib/fs.c',
        'bilab/graph/metis/GKlib/getopt.c',
        'bilab/graph/metis/GKlib/gkregex.c',
        'bilab/graph/metis/GKlib/graph.c',
        'bilab/graph/metis/GKlib/htable.c',
        'bilab/graph/metis/GKlib/io.c',
        'bilab/graph/metis/GKlib/itemsets.c',
        'bilab/graph/metis/GKlib/mcore.c',
        #'bilab/graph/metis/GKlib/memory.c',
        'bilab/graph/metis/GKlib/omp.c',
        'bilab/graph/metis/GKlib/pdb.c',
        'bilab/graph/metis/GKlib/pqueue.c',
        'bilab/graph/metis/GKlib/random.c',
        'bilab/graph/metis/GKlib/rw.c',
        'bilab/graph/metis/GKlib/seq.c',
        'bilab/graph/metis/GKlib/sort.c',
        'bilab/graph/metis/GKlib/string.c',
        'bilab/graph/metis/GKlib/timers.c',
        'bilab/graph/metis/GKlib/tokenizer.c',
        'bilab/graph/metis/GKlib/util.c',
        #'bilab/graph/metis/libmetis/auxapi.c',
        'bilab/graph/metis/libmetis/balance.c',
        'bilab/graph/metis/libmetis/bucketsort.c',
        'bilab/graph/metis/libmetis/checkgraph.c',
        'bilab/graph/metis/libmetis/coarsen.c',
        'bilab/graph/metis/libmetis/compress.c',
        'bilab/graph/metis/libmetis/contig.c',
        'bilab/graph/metis/libmetis/debug.c',
        'bilab/graph/metis/libmetis/fm.c',
        'bilab/graph/metis/libmetis/fortran.c',
        'bilab/graph/metis/libmetis/frename.c',
        'bilab/graph/metis/libmetis/gklib.c',
        'bilab/graph/metis/libmetis/graph.c',
        'bilab/graph/metis/libmetis/initpart.c',
#        'bilab/graph/metis/libmetis/kmetis.c',
        'bilab/graph/metis/libmetis/kwayfm.c',
        'bilab/graph/metis/libmetis/kwayrefine.c',
        'bilab/graph/metis/libmetis/mcutil.c',
        'bilab/graph/metis/libmetis/mesh.c',
        'bilab/graph/metis/libmetis/meshpart.c',
        'bilab/graph/metis/libmetis/minconn.c',
        'bilab/graph/metis/libmetis/mincover.c',
        'bilab/graph/metis/libmetis/mmd.c',
        'bilab/graph/metis/libmetis/ometis.c',
        'bilab/graph/metis/libmetis/options.c',
        'bilab/graph/metis/libmetis/parmetis.c',
#        'bilab/graph/metis/libmetis/pmetis.c',
        'bilab/graph/metis/libmetis/refine.c',
        'bilab/graph/metis/libmetis/separator.c',
        'bilab/graph/metis/libmetis/sfm.c',
        'bilab/graph/metis/libmetis/srefine.c',
        'bilab/graph/metis/libmetis/stat.c',
        'bilab/graph/metis/libmetis/timing.c',
        'bilab/graph/metis/libmetis/util.c',
        'bilab/graph/metis/libmetis/wspace.c',
        'bilab/graph/metis/gpmetis_wrap.pyx'
    ])

# Isosurface = Extension(
#    'bilab.geometry._IsoSurface',
#    define_macros=[('MAJOR_VERSION', '1'),
#                   ('MINOR_VERSION', '0')],
#    include_dirs=[boost_inc_dir],
#    libraries=['boost_python'],
#    library_dirs=[boost_lib_dir],
#    extra_compile_args=['-ftemplate-backtrace-limit=64'],
#    sources=['bilab/geometry/isosurface/export.cpp'])

# CWMatrix = Extension(
#     'bilab.linalg._CWMatrix',
#     define_macros=[('MAJOR_VERSION', '1'),
#                    ('MINOR_VERSION', '0')],
#     # include_dirs = [boost_root],
#     libraries=['boost_python'],
#     # library_dirs = [boost_lib_dir],
#     extra_compile_args=['-ftemplate-backtrace-limit=64'],
#     sources=['bilab/linalg/matrix/export.cpp'])
#     # 'bilab/linalg/matrix/matrix.cpp'])

# using Boost.Python
#bhtsne_wrap = Extension(
#    'bilab.ml.NDR.tSNE._bhtsne_wrap',
#    define_macros=[('MAJOR_VERSION', '1'),
#                   ('MINOR_VERSION', '0')],
#    include_dirs=[np_include_dir, boost_inc_dir],
#    library_dirs=[boost_lib_dir],
#    libraries=['boost_python'],
#    extra_compile_args=['-ftemplate-backtrace-limit=64',
#                        '-std=c++11', '-g'],
#    sources=['bilab/ml/NDR/tSNE/bhtsne_export.cpp',
#             'bilab/ml/NDR/tSNE/bhtsne.cpp',
#             'bilab/ml/NDR/tSNE/sptree.cpp'])

# Cython version
bhtsne_wrap = Extension(
    name='bilab.ml.NDR.tSNE._bhtsne_wrap',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0')],
    include_dirs=[np_include_dir],
    # library_dirs=[boost_lib_dir],
    # libraries=['boost_python'],
    extra_compile_args=basic_extra_compile_args,
    extra_link_args=basic_extra_link_args,
    sources=['bilab/ml/NDR/tSNE/_bhtsne_wrap.pyx',
             'bilab/ml/NDR/tSNE/bhtsne.cpp',
             'bilab/ml/NDR/tSNE/sptree.cpp'],
    language="c++")

netcdf_wrap = Extension(
    name='bilab.structure.minimization._NetCDF',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0'),
                   ('NUMPY', 1)],
    include_dirs=[np_include_dir,
                  #  boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization',
                  'bilab/structure/minimization/MMTK'],
    # library_dirs=[boost_lib_dir],
    # libraries=['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64', '-g'],
    sources=['bilab/structure/minimization/_netcdf.c'])

MMTK_universe = Extension(
    name='bilab.structure.minimization._MMTK_universe',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0'),
                   ('NUMPY', 1)],
    include_dirs=[np_include_dir,
                  #  boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization'],
    # library_dirs=[boost_lib_dir],
    # libraries=['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64', '-g'],
    sources=['bilab/structure/minimization/MMTK_universe.c'])

MMTK_surface = Extension(
    name='bilab.structure.minimization._MMTK_surface',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0'),
                   ('NUMPY', 1)],
    include_dirs=[np_include_dir,
                  #  boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization'],
    # library_dirs=[boost_lib_dir],
    # libraries=['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64', '-g'],
    sources=['bilab/structure/minimization/MMTK_surface.c'])

MMTK_trajectory = Extension(
    name='bilab.structure.minimization._MMTK_trajectory',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0'),
                   ('NUMPY', 1)],
    include_dirs=[np_include_dir,
                  # boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization'],
    # library_dirs = [boost_lib_dir],
    # libraries = ['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64', '-g'],
    sources=['bilab/structure/minimization/MMTK_trajectory.c'])

MMTK_minimization = Extension(
    name='bilab.structure.minimization._MMTK_minimization',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0')],
    include_dirs=[np_include_dir,
                  # boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization'],
    # library_dirs=[boost_lib_dir],
    # libraries=['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64', '-DNUMPY=1', '-g'],
    sources=['bilab/structure/minimization/MMTK_minimization.c'])

MMTK_energy_term = Extension(
    name='bilab.structure.minimization._MMTK_energy_term',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0'),
                   ('NUMPY', 1)],
    include_dirs=[np_include_dir,
                  # boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization'],
    # library_dirs=[boost_lib_dir],
    # libraries=['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64', '-g'],
    sources=['bilab/structure/minimization/_MMTK_energy_term.c'])

MMTK_force_field = Extension(
    name='bilab.structure.minimization._MMTK_forcefield',
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0'),
                   ('SERIAL', None),
                   ('VIRIAL', None),
                   ('MACROSCOPIC', None),
                   ('NUMPY', 1),
                   ('WITH_DPMTA', 1)],
    include_dirs=[np_include_dir,
                  # boost_inc_dir,
                  netcdf_inc_dir,
                  'bilab/structure/minimization',
                  'bilab/structure/minimization/dpmta',
                  'bilab/structure/minimization/dpmta/src',
                  'bilab/structure/minimization/dpmta/mpole'],
    # library_dirs = [boost_lib_dir],
    # libraries = ['boost_python'],
    extra_compile_args=['-ftemplate-backtrace-limit=64',
                        '-O3', '-ffast-math',
                        '-fomit-frame-pointer',
                        '-g'],
    sources=[
        'bilab/structure/minimization/MMTK_forcefield.c',
        'bilab/structure/minimization/bonded.c',
        'bilab/structure/minimization/nonbonded.c',
        'bilab/structure/minimization/ewald.c',
        'bilab/structure/minimization/sparsefc.c',
        'bilab/structure/minimization/dpmta/mpole/mpe_fft.c',
        'bilab/structure/minimization/dpmta/mpole/mpe_misc.c',
        'bilab/structure/minimization/dpmta/mpole/mpe_mpoleC.c',
        'bilab/structure/minimization/dpmta/mpole/mpe_allocC.c',
        'bilab/structure/minimization/dpmta/mpole/mpe_mpoleLJ.c',
        'bilab/structure/minimization/dpmta/mpole/mpe_allocLJ.c',
        'bilab/structure/minimization/dpmta/src/dpmta_serial.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvmkcell.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvmcalc.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvpcalc.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvmkil.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvmkhl.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvcompute.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvmacro.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvscale.c',
        'bilab/structure/minimization/dpmta/src/dpmta_timer.c',
        'bilab/structure/minimization/dpmta/src/dpmta_slvglobals.c',
        'bilab/structure/minimization/dpmta/src/dpmta_distmisc.c'
    ])

mmcif_extra_compile_args = ['-ftemplate-backtrace-limit=64','-std=c++11']
mmcif_extra_link_args = ['-std=c++11']
if plat.system() == "Darwin":
    # clang on MacOS
    mmcif_extra_compile_args += ['-stdlib=libc++']
    mmcif_extra_link_args += ['-stdlib=libc++']
elif plat.system() == "Linux":
    # gcc on Linux
    mmcif_extra_compile_args += ['-static-libstdc++']
    mmcif_extra_link_args += ['-static-libstdc++']

mmcif = Extension(
    name="bilab.io._mmcifio",
    define_macros=[('MAJOR_VERSION', '1'),
                   ('MINOR_VERSION', '0')],
    include_dirs=['.',
                  np_include_dir,
                  'bilab/io/cifparser/mmcif'
                  ],
    sources=[
        'bilab/io/cifparser/mmcif/utils.cpp',
        'bilab/io/cifparser/mmcif/mmcifIO_wrapper.cpp',
        'bilab/io/cifparser/mmcifIO_ext.cpp',
        ],
    extra_compile_args = mmcif_extra_compile_args,
    extra_link_args = mmcif_extra_link_args,
    language="c++")

def split_multiline(value):
    """Split a multiline string into a list, excluding blank lines."""
    val = value.encode('utf-8')
    return [element for element in
            (line.decode('utf-8').strip() for line in val.split(b'\n'))
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
                        "package_data": ("files", "package_data"),
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
                    "resources",
                    "package_data"
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
    # f = codecs.open(path, encoding='utf-8')
    f = io.open(path, encoding="utf-8")
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
                        # fp = codecs.open(filename, encoding='utf-8')
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
            if arg == 'packages' and in_cfg_value:
                if 'package_dir' in kwargs:
                    if kwargs['package_dir']['']:
                        in_cfg_value = [
                            kwargs['package_dir']['']+'.'+pack
                            for pack in in_cfg_value]
            if arg == 'package_data' and in_cfg_value:
                datafiles = {}
                for line in in_cfg_value:
                    split_path = line.split('/')
                    package_strcut_str = ".".join(split_path[:-1]).strip()
                    files = [f.split('/')[-1] for f in glob(line)]
                    datafiles[package_strcut_str] = files
                in_cfg_value = datafiles
                kwargs['install_package_data'] = True

            if arg == 'resources' and in_cfg_value:
                datafiles = []
                for line in in_cfg_value:
                    file_str = line.split('=')[0].strip()
                    files = glob(file_str)
                    dfile = ""
                    if os.path.dirname(files[0]):
                        dfile = os.path.dirname(files[0])
                    # datafiles.extend(glob(file_str))
                    datafiles.append(
                        ('share' + os.sep + 'bilab' + os.sep + dfile,
                         files))
                print(datafiles)
                # kwargs['data_files'] = [
                #   ('share/bilab',[line.split('=')[0].strip()
                #    for line in in_cfg_value])]
                # install into python share lib
                # kwargs['data_files'] = [('share/bilab/data',datafiles)]
                kwargs['data_files'] = datafiles

        kwargs[arg] = in_cfg_value

    return kwargs

# setup
general_settings = cfg_to_args()
# Key: resources has to be removed
general_settings.pop('resources')
# cython
general_settings['cmdclass'] = {'build_ext': build_ext}
# extensions
general_settings['ext_modules'] = [
    distance_wrap,
    kdtree_lib,
    bhtsne_wrap,
    gpmetis_wrap,
    lfdfiles,
    voroplusplus,
    marching_cubes,
#    netcdf_wrap,
#    MMTK_surface,
#    MMTK_minimization,
#    MMTK_trajectory,
#    MMTK_universe,
#    MMTK_energy_term,
#    MMTK_force_field,
#    crc32,
    mmcif
]
#                                   CWMatrix]
general_settings['install_requires'] = [
    'scipy>0.15.0',
    'numpy',
    'ete3',
    'jinja2'
    ]
#    'matplotlib>1.4.0',
#    'numpy>1.9.0']

for k in general_settings:
    if k == "packages" and sys.version_info[:3] < (2, 7, 13):
        general_settings[k] = [ p.encode('ascii', 'ignore') for p in general_settings[k]]
setup(**general_settings)
# setup(**cfg_to_args())