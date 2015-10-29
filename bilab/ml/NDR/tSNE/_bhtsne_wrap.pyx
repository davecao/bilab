# distutils: language = c++
#
# isoSurface.pyx : cython interface to bhtsne
#
# this file provides a python interface of implementation of t-SNE
# using c++.
#
#
# Copyright (c) Wei Cao 2015
# contact: <davecao@bi.a.u-tokyo.ac.jp> or <caotiger@gmail.com>
#

import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
cimport cython
from cython.view cimport array as cvarray

DTYPE = np.float64

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# declaration
cdef extern from "bhtsne.h":
    cdef cppclass BHTSNE:
        BHTSNE(double*, int, int, int, double, double, int, bool, bool) except +
        int N, D  # NxD matrix (N = #points, D = dimensionality)
        int no_dims  # m' x n' after mapping
        double* X
        double* Y
        double perplexity
        double theta
        bint verbosity
        bint scaling
        int randseed
        void print_()
        void run()



# Wrapper class
cdef class bhtsne:
    cdef BHTSNE *thisptr  # hold a C++ instance which we're wrapping
#    def __cinit__(self, double* samples, int numOfsamples, int dims, 
#                        int mapped_D, double perplex, double th,
#                        int rseed, bool scale, bool verbose)
    
    @cython.wraparound(False)
    @cython.boundscheck(False) # turn of bounds-checking for entire function
    def __init__(self, np.ndarray samples not None,
                       int N, int D, int no_dims, float perplexity,
                       float theta, int randseed, bint scale, bint verbose):
        assert samples.dtype == DTYPE

        self.thisptr = new BHTSNE(<double*>samples.data, N, D, 
                                  no_dims, perplexity, theta, randseed, 
                                  scale, verbose)


    def __dealloc__(self):
        del self.thisptr

    # disable bounds checks on indexing
    @cython.boundscheck(False)
    # disable negative indexes
    @cython.wraparound(False)
    # disable checking that an object isn't None before calling methods
    @cython.nonecheck(False)
    def run(self):

        # run tsne
        self.thisptr.run()

        # Convert double* Y: self.thisptr.Y to np
        cdef unsigned int rows = self.thisptr.N,
        cdef unsigned int cols = self.thisptr.no_dims
        cdef unsigned int counter = 0

        # Return a numpy array
        #cdef np.ndarray[np.double_t] Y = np.zeros(rows, cols)
        cdef np.ndarray Y = np.empty(shape = (rows, cols), dtype = np.float64, order='c')
        #Y = np.zeros(rows, cols)
        #cdef double[:] Y
        for row in range(rows):
            for col in range(cols):
              Y[row, col] = self.thisptr.Y[counter]
              counter = counter + 1

        return Y


