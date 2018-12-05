# -*- coding: utf-8 -*-
# lfdfiles.pyx

# Copyright (c) 2012-2015, Christoph Gohlke
# Copyright (c) 2012-2015, The Regents of the University of California
# Produced at the Laboratory for Fluorescence Dynamics.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Alternative implementation of some lfdfiles.py functions in Cython.

The lfdfiles.py module for reading many proprietary file formats used to
store experimental data at the Laboratory for Fluorescence Dynamics is
available at http://www.lfd.uci.edu/~gohlke/.

:Author:
  `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`_

:Organization:
  Laboratory for Fluorescence Dynamics, University of California, Irvine

:Version: 2015.02.19

Requirements
------------
* `CPython 2.7 or 3.4 <http://www.python.org>`_
* `Numpy 1.8.2 <http://www.numpy.org>`_
* `Cython 0.22 <http://cython.org/>`_
* A Python distutils compatible C compiler

Install
-------
Use this Cython distutils setup script to build the extension module::

  # setup.py
  # Usage: ``python setup.py build_ext --inplace``
  from distutils.core import setup, Extension
  from Cython.Distutils import build_ext
  setup(name='_lfdfiles',
        cmdclass={'build_ext': build_ext},
        ext_modules=[Extension('_lfdfiles', ['lfdfiles.pyx'],
                               extra_compile_args=['/openmp', '-fopenmp'],
                               extra_link_args=['-fopenmp'])])

"""

from cython.parallel import parallel, prange
cimport cython

import numpy as np
cimport numpy as np

ctypedef np.int8_t int8_t
ctypedef np.int16_t int16_t
ctypedef np.int32_t int32_t
ctypedef np.int64_t int64_t
ctypedef np.uint16_t uint16_t
ctypedef np.uint32_t uint32_t
ctypedef np.uint64_t uint64_t

ctypedef fused uintxx_t:
    uint32_t
    uint64_t

__version__ = "2015.02.19"

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def simfcsfbd_decode(
    uint16_t[::] data,
    int8_t[:, ::1] bins_out,
    uintxx_t[::1] times_out,
    ssize_t[::1] markers_out,
    int windows,
    int pmax,
    int pdiv,
    int harmonics,
    int16_t[:, ::] decoder_table,
    uint16_t tcc_mask,
    uint32_t tcc_shr,
    uint16_t pcc_mask,
    uint32_t pcc_shr,
    uint16_t marker_mask,
    uint32_t marker_shr,
    uint16_t win_mask,
    uint32_t win_shr
    ):
    """Decode flimbox data stream.

    See the lfdfiles.SimfcsFbd documentation for parameter descriptions.

    """
    cdef int maxwindex = <int>decoder_table.shape[1]
    cdef ssize_t maxmarker = markers_out.size
    cdef ssize_t nchannel = bins_out.shape[0]
    cdef ssize_t datasize = bins_out.shape[1]
    cdef ssize_t i
    cdef ssize_t j
    cdef ssize_t c
    cdef uint16_t d
    cdef uintxx_t tcc_max
    cdef uintxx_t t0
    cdef uintxx_t t1
    cdef int m0
    cdef int m1
    cdef int pcc
    cdef int win
    cdef int pmax_win = harmonics * pmax // windows

    if bins_out.shape[0] != decoder_table.shape[0]:
        raise ValueError("shape mismatch between bins and decoder_table")
    if bins_out.shape[1] != times_out.size:
        raise ValueError("shape mismatch between bins and time")

    # calculate cross correlation phase index
    for c in prange(nchannel, nogil=True):
        for i in range(datasize):
            d = data[i]
            pcc = <int>((d & pcc_mask) >> pcc_shr)
            win = <int>((d & win_mask) >> win_shr)
            if win < maxwindex:
                win = <int>(decoder_table[c, win])
                if win >= 0:
                    bins_out[c, i] = <int8_t>(
                        (pmax-1 - (pcc + win*pmax_win) % pmax) // pdiv)
                else:
                    bins_out[c, i] = -1  # no event
            else:
               bins_out[c, i] = -2  # should never happen

    # record up-markers and absolute time
    tcc_max = (tcc_mask >> tcc_shr) + 1
    j = 0
    m0 = <int>(data[0] & marker_mask)
    t0 = (data[0] & tcc_mask) >> tcc_shr
    times_out[0] = 0
    for i in range(1, datasize):
        d = data[i]
        # detect up-markers
        if j < maxmarker:
            m1 = <int>(d & marker_mask)
            if m1 > m0:
                markers_out[j] = i
                j += 1
            m0 = m1
        # cumulative sum of differences of cross correlation time
        t1 = t0
        t0 = (d & tcc_mask) >> tcc_shr
        if t0 > t1:
            times_out[i] = times_out[i-1] + (t0 - t1)
        elif t0 == 0:
            times_out[i] = times_out[i-1] + (tcc_max - t1)
        else:
            # is this supposed to happen?  0 < t0 <= t1
            times_out[i] = times_out[i-1] + (tcc_max - t1) + t0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def simfcsfbd_histogram(
    int8_t[:, ::1] bins,
    uintxx_t[::1] times,
    frame_markers,
    double units_per_sample,
    double scanner_frame_start,
    uint16_t[:, :, :, ::1] hist_out
    ):
    """Calculate histograms from decoded flimbox data and frame markers.

    See the lfdfiles.SimfcsFbd documentation for parameter descriptions.

    """
    cdef ssize_t nframes = hist_out.shape[0]
    cdef ssize_t nchannels = hist_out.shape[1]
    cdef ssize_t framelen = hist_out.shape[2]
    cdef ssize_t nwindows = hist_out.shape[3]
    cdef ssize_t i
    cdef ssize_t j
    cdef ssize_t k
    cdef ssize_t f
    cdef ssize_t c
    cdef ssize_t idx
    cdef uintxx_t t0
    cdef int8_t w

    if bins.shape[0] != hist_out.shape[1]:
        raise ValueError("shape mismatch between bins and hist_out")
    if bins.shape[1] != times.shape[0]:
        raise ValueError("shape mismatch between bins and times")

    units_per_sample = 1.0 / units_per_sample
    for f, (j, k) in enumerate(frame_markers):
        f = f % nframes
        t0 = times[j]
        for c in prange(nchannels, nogil=True):
            for i in range(j, k):
                idx = <ssize_t>(<double>(times[i] - t0) * units_per_sample
                                - scanner_frame_start)
                if idx >= 0 and idx < framelen:
                    w = bins[c, i]
                    if w >= 0:
                        hist_out[f, c, idx, w] += 1