# -*- coding: utf-8 -*-

import sys
import os
import abc
import re

__all__ = ["FindBoostPackage"]

DYNAMIC_LIBS_EXT = {
    "darwin": 'dylib',
    "linux": 'so',
    "win32": 'dll'
}

PY3K = sys.version_info[0] > 2

pj = os.path.join

class FindPackage(abc.ABC):
    """ Base class """
    def __init__(self, env_var):
        self.search_include_paths = [
            "/usr/local/include",
            "/opt/local/include"
        ]
        self.search_library_paths = [
            "/usr/local/lib",
            "/opt/local/lib"
        ]
        self.available = False
        self.environ_var = None

        if env_var in os.environ:
            self.environ_var = os.environ[env_var]
            if os.path.isdir(self.environ_var):
                self.search_include_paths.append(self.environ_var)
                self.search_library_paths.append(self.environ_var)

        self.libext = None
        if sys.platform in DYNAMIC_LIBS_EXT:
           self.libext = DYNAMIC_LIBS_EXT[sys.platform]

    @abc.abstractmethod
    def __call__(self, inc_paths=None, lib_paths=None):
        raise NotImplementedError


class FindBoostPackage(FindPackage):
    """
    Find the include path and libraries path of Boost library

    Search paths:
    Mac, Linux: boost
        BOOST_ROOT
        /opt/include
        /opt/lib
        /usr/local/include
        /usr/local/lib
    """
    def __init__(self,
                 search_inc_paths=None,
                 search_lib_paths=None,
                 eviron_var_name="BOOST_ROOT",
                 minimal_version="1.55.0"):

        super(FindBoostPackage, self).__init__(eviron_var_name)

        self.minimal_version = minimal_version
        self.staged_version = False
        self.available = False  # Track availability
        self.baseDir = None
        self.incDir = None

        # configurable options
        self.headerMap = {
            'program_options': 'boost/program_options.hpp',
            'python': 'boost/python.hpp',
            'thread': 'boost/thread.hpp',
            'filesystem': 'boost/filesystem/path.hpp',
            'version': 'boost/version.hpp'
        }

    def __call__(self, inc_paths=None, lib_paths=None):
        """ search

        Args:
            inc_paths (list):
                a list of paths of include files
            lib_paths (list):
                a list of paths of libraries
        """
        if (inc_paths is not None) and isinstance(inc_paths, list):
            self.search_include_paths.extend(inc_paths)


    def find_boost_version(self, version_header_file):
        """ Check the boost version from version.hpp """
        if not os.path.isfile(version_header_file):
            return
        if PY3K:
            ver_file = open(version_header_file)
        else:
            ver_file = file(version_header_file)
        ver_match = re.search(
                    "define\s+?BOOST_VERSION\s+?(\d*)",
                    ver_file.read())
        if not ver_match:
            return None
        found_ver_str = int(ver_match.group(1))
        found_ver_str = "{}.{}.{}".format(
            found_ver_str/100000,
            found_ver_str/100 % 1000,
            found_ver_str % 100)
        return found_ver_str

    def setupLibrarySettings(self):
        # Map for extra libs needed for config test
        self.extraEnvOptions = {}

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
            print("found version is to old: required:{} found:{}".format(
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
