# distutils: language = c++
#
# isoSurface.pyx : cython interface to CIsoSurface
#
# this file provides a python interface of implementation of marching cube
# using c++.
#
#
# Copyright (c) Wei Cao 2015
# contact: <davecao@bi.a.u-tokyo.ac.jp> or <caotiger@gmail.com>
#

# import dereference and increment operators
from cython.operator cimport dereference as deref, preincrement as inc

cdef extern from "CIsoSurface.h":
    cdef cppclass CIsoSurface[T]:
        CIsoSurface() except +
        void GenerateSurface(T*, T, int, int, int, float, float, float)
        bool IsSurfaceValid()
        void DeleteSurface()
        int GetVolumeLengths(float&, float&, float&)


# Wrapper class
cdef class IsoSurface:
    cdef CIsoSurface *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new CIsoSurface()

    def __dealloc__(self):
        del self.thisptr

    
