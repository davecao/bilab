# -*- coding: utf-8 -*-

__all__ = []

import os
import sys
import inspect
import tempfile
from glob import glob
from os.path import abspath, split, join, relpath, splitext
from os.path import sep as dirsep


try:
    import unittest2 as unittest
    from unittest2 import TestCase, skipIf, skipUnless
except ImportError:
    import unittest
    from unittest import TestCase, skipIf, skipUnless

TESTDIR = abspath(split(inspect.getfile(inspect.currentframe()))[0])
TEMPDIR = tempfile.gettempdir()

def test_suite():
    suite = unittest.TestSuite()
    here = os.path.dirname(__file__) or os.curdir
    print("Current dir:{}".format(here))
    for fn in os.listdir(here):
        if fn.startswith("test") and fn.endswith(".py"):
            modname = "build.lib.bilab.tests." + fn[:-3]
            __import__(modname)
            module = sys.modules[modname]
            suite.addTest(module.test_suite())
    return suite
