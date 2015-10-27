# borrowed from http://www.lfd.uci.edu/~gohlke/
# -*- coding: utf-8 -*-
# lfdfiles.py

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

"""Read Laboratory for Fluorescence Dynamics (LFD) file formats.

The lfdfiles.py module provides classes for reading many proprietary file
formats used to store experimental data at the `Laboratory for Fluorescence
Dynamics <http://www.lfd.uci.edu/>`_.

:Author:
  `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`_

:Organization:
  Laboratory for Fluorescence Dynamics, University of California, Irvine

:Version: 2015.02.19

Requirements
------------
* `CPython 2.7 or 3.4 <http://www.python.org>`_
* `Numpy 1.8.2 <http://www.numpy.org>`_
* `Lfdfiles.pyx 2015.02.19 <http://www.lfd.uci.edu/~gohlke/>`_
  (recommended for speedup of some functions)

Revisions
---------
2015.02.19
    Initial support for new FBD files containing headers.
2014.12.02
    Read B64, R64, I64 and Z64 files (SimFCS version 4).
2014.10.10
    Read SimFCS FIT files.
2014.04.08
    Read and write CCP4 MAP volume files.
2013.08.10
    Read second harmonics FlimBox data.

Notes
-----
The API is not stable yet and is expected to change between revisions.

Most of the file formats are not well documented and evolve organically.
This implementation is largely based on reverse engineering existing files.
No promise can be made as to the correctness of code and documentation.

Experimental data are often stored in plain binary files with metadata
available in separate, human readable journal files (\*.jrn).

Unless specified otherwise, data are stored in little endian, C contiguous
order.

Software
--------
The following software is referenced in this module:

(1) `SimFCS <http://www.lfd.uci.edu/globals/>`_, a.k.a. Globals for
    Images, is software for fluorescence image acquisition, analysis, and
    simulation, developed by Enrico Gratton.
(2) `Globals <http://www.lfd.uci.edu/globals/>`_, a.k.a. Globals for
    Spectroscopy, is software for the analysis of multiple files from
    fluorescence spectroscopy, developed by Enrico Gratton.
(3) ImObj is software for image analysis, implemented on Win16.
(4) FlimFast is software for frequency-domain, full-field, fluorescence
    lifetime imaging at video rate, developed by Christoph Gohlke.
(5) FLImage is software for frequency-domain, full-field, fluorescence
    lifetime imaging, developed by Christoph Gohlke. Implemented in LabVIEW.
(6) FLIez is software for frequency-domain, full-field, fluorescence
    lifetime imaging, developed by Glen Redford.
(7) Flie is software for frequency-domain, full-field, fluorescence
    lifetime imaging, developed by Peter Schneider. Implemented on a
    Sun UltraSPARC.
(8) FLOP is software for frequency-domain, cuvette, fluorescence lifetime
    measurements, developed by Christoph Gohlke. Implemented in LabVIEW.
(9) `VistaVision <http://www.iss.com/microscopy/software/vistavision.html>`_
    is commercial software for instrument control, data acquisition and data
    processing by ISS Inc.

"""

from __future__ import division, print_function

import os
import sys
import re
import glob
import struct
import zlib
import zipfile
import warnings

import numpy

__version__ = '2015.02.19'
__docformat__ = 'restructuredtext en'
__all__ = (
    'LfdFile', 'LfdFileSequence', 'RawPal', 'SimfcsVpl', 'SimfcsVpp',
    'SimfcsJrn', 'SimfcsBin', 'SimfcsRaw', 'SimfcsInt', 'SimfcsIntPhsMod',
    'SimfcsFit', 'SimfcsCyl', 'SimfcsRef', 'SimfcsBh', 'SimfcsBhz',
    'SimfcsFbf', 'SimfcsFbd', 'SimfcsGpSeq', 'SimfcsMap',
    'SimfcsB64', 'SimfcsI64', 'SimfcsZ64', 'SimfcsR64', 'GlobalsLif',
    'GlobalsAscii', 'FlimfastFlif', 'FlimageBin', 'FlieOut', 'FliezI16',
    'FliezDb2', 'VistaIfli', 'save_map')


class lazyattr(object):
    """Lazy object attribute whose value is computed on first access."""
    __slots__ = ('func', )

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, owner):
        if instance is None:
            return self
        result = self.func(instance)
        if result is NotImplemented:
            return getattr(super(owner, instance), self.func.__name__)
        setattr(instance, self.func.__name__, result)
        return result


class LfdFile(object):
    """Base class for reading single LFD files."""
    _filemode = 'rb'  # file open mode
    _filepattern = r'.*'  # regular expression matching valid file names

    class Error(Exception):
        """Exception to indicate invalid LFD files."""
        def __init__(self, arg='', msg=''):
            if isinstance(arg, LfdFile):
                arg = 'not a valid %s file' % arg.__class__.__name__
            Exception.__init__(self, arg + msg)

    def __init__(self, filename, *args, **kwargs):
        """Open file and read headers and metadata.

        Parameters
        ----------
        validate : bool, optional
            It True, filename must match the _filepattern regular expression.

        """
        self._filepath, self._filename = os.path.split(filename)
        if 'validate' in kwargs:
            if kwargs['validate']:
                self._valid_name()
            del kwargs['validate']
        self._fh = open(filename, self._filemode)
        self._init(*args, **kwargs)
        if self._fh:
            self._pos = self._fh.tell()

    def __str__(self):
        """Return string with information about file."""
        return '\n'.join([
            os.path.join(self._filepath, self._filename).capitalize(),
            (" (%s file)" % self.__class__.__name__),
            self._str()])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self._fh:
            self._fh.close()
            self._fh = None

    def asarray(self, *args, **kwargs):
        """Return data in file as numpy array."""
        if self._fh:
            self._fh.seek(self._pos)
        return self._data(*args, **kwargs)

    def _init(self, *args, **kwargs):
        """Override to validate file and read metadata."""
        pass

    def _data(self, *args, **kwargs):
        """Override to read data from file and return as numpy array."""
        raise NotImplementedError()

    def _str(self):
        """Override to return extra information about file."""
        return ""

    def _valid_name(self):
        """Raise LfdFile.Error if file name does not match _filepattern."""
        name = self._filename
        pattern = self._filepattern
        if not pattern or pattern == r'.*':
            return
        if re.search(pattern, name, re.IGNORECASE) is None:
            raise LfdFile.Error(
                self, ".\n    File name '%s' doesn't match '%s')" % (name,
                                                                     pattern))

    def _decompress_header(self, max_length, max_read=256):
        """Return the first uncompressed bytes of a zlib compressed file."""
        data = self._fh.read(max_read)
        self._fh.seek(0)
        if not data.startswith(b'\x78\x9c'):
            raise LfdFile.Error(self, "not a zlib compressed file")
        decompressor = zlib.decompressobj()
        return decompressor.decompress(data, max_length)

    @lazyattr
    def _fstat(self):
        """Return status of open file."""
        return os.fstat(self._fh.fileno())

    @lazyattr
    def _filesize(self):
        """Return file size in bytes."""
        pos = self._fh.tell()
        self._fh.seek(0, 2)
        size = self._fh.tell()
        self._fh.seek(pos)
        return size


class LfdFileSequence(object):
    """Sequence of LFD files.

    Attributes
    ----------
    files : list of str
        List of file names.
    shape : tuple
        Shape of image sequence detected from file names.

    Examples
    --------
    >>> ims = LfdFileSequence("gpint/v?001.int", readfunc=SimfcsInt)
    >>> ims = ims.asarray()
    >>> ims.shape
    (2, 256, 256)

    """
    _readfunction = None
    _indexpattern = None

    class _ParseError(Exception):
        pass

    def __init__(self, files, readfunc=None, indexpattern=None,
                 asarray='asarray'):
        """Initialize instance from multiple files.

        Parameters
        ----------
        files : str, or sequence of str
            Glob pattern or sequence of file names.
        readfunc : LfdFile class
            File read function or class with ``asarray`` function returning
            numpy array from single file.
        indexpattern : str
            Regular expression pattern that matches sequence indices in the
            file names.
        asarray : str (optional)
            Name of class instance function returning numpy array.
            The default name is 'asarray'.

        """
        if isinstance(files, basestring):
            files = natural_sorted(glob.glob(files))
        if not files or not os.path.isfile(files[0]):
            raise ValueError("no files found")
        self.files = files

        if readfunc is None:
            readfunc = self._readfunction
        if not readfunc:
            raise ValueError("invalid data read function")
        if hasattr(readfunc, asarray):
            _readfunction = readfunc

            def readfunc(fname, *args, **kwargs):
                with _readfunction(fname) as im:
                    return getattr(im, asarray)(*args, **kwargs)

        self._readfunction = readfunc

        if indexpattern is None:
            indexpattern = self._indexpattern
        if indexpattern:
            self._parse(indexpattern)
        else:
            self.shape = (len(files),)
            self._start_index = (0,)
            self._indices = ((i,) for i in range(len(files)))

    def asarray(self, *args, **kwargs):
        """Read image data from all files and return as single numpy array.

        Raise IndexError if image shapes don't match.

        """
        im = self._readfunction(self.files[0], *args, **kwargs)
        result_shape = self.shape + im.shape
        result = numpy.zeros(result_shape, dtype=im.dtype)
        result = result.reshape(-1, *im.shape)
        for index, fname in zip(self._indices, self.files):
            index = [i-j for i, j in zip(index, self._start_index)]
            index = numpy.ravel_multi_index(index, self.shape)
            im = self._readfunction(fname, *args, **kwargs)
            result[index] = im
        result.shape = result_shape
        return result

    def _parse(self, pattern):
        """Get shape from file names."""
        if not pattern:
            raise self._ParseError("invalid pattern")
        pattern = re.compile(pattern, re.IGNORECASE | re.VERBOSE)
        matches = pattern.findall(self.files[0])
        if not matches:
            raise self._ParseError("pattern doesn't match file names")
        indices = []
        for fname in self.files:
            matches = pattern.findall(fname)[-1]
            indices.append([int(m) for m in matches])
        shape = tuple(numpy.max(indices, axis=0))
        start_index = tuple(numpy.min(indices, axis=0))
        shape = tuple(i-j+1 for i, j in zip(shape, start_index))
        if product(shape) != len(self.files):
            warnings.warn("files are missing. Missing data are zeroed")
        self.shape = shape
        self._indices = indices
        self._start_index = start_index

    def __str__(self):
        """Return string with information about image sequence."""
        return "\n".join([
            self.files[0],
            '* files: %i' % len(self.files),
            '* shape: %s' % str(self.shape)])

    def __len__(self):
        return len(self.files)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        pass


class RawPal(LfdFile):
    """Raw color palette.

    PAL files contain a single RGB or RGBA color palette, stored as 256x3 or
    256x4 unsigned bytes in C or Fortran order, without any header.

    Examples
    --------
    >>> with RawPal('rgb.pal') as f:
    ...     print(f.asarray()[100])
    [ 16 255 239]
    >>> with RawPal('rgba.pal') as f:
    ...     print(f.asarray()[100])
    [219 253 187 255]
    >>> with RawPal('rrggbb.pal') as f:
    ...     print(f.asarray()[100])
    [182 114  91]
    >>> with RawPal('rrggbbaa.pal') as f:
    ...     print(f.asarray()[100])
    [182 114  91 170]
    >>> with RawPal('rrggbbaa.pal') as f:
    ...     print(f.asarray('F')[100])
    [182 114  91 170]

    """
    _filepattern = r'.*\.(pal|raw|bin|lut)'

    def _init(self):
        """Verify file size."""
        if self._filesize not in (768, 1024):
            raise LfdFile.Error(self)
        self.dtype = numpy.dtype('u1')

    def _data(self, order=None):
        """Return palette data as (256, 3) or (256, 4) shaped array of uint8.

        Parameters
        ----------
        order : {'C', 'F'}, optional
            Determines whether the data is stored in C (row-major) or
            Fortran (column-major) order. By default the order is
            determined based on size and sum of differences.

        """
        data = numpy.fromfile(self._fh, 'u1').reshape(256, -1)
        if order is None:
            a = data.astype('i4')
            b = a.reshape(-1, 256).T
            if (numpy.sum(numpy.abs(numpy.diff(a, axis=0))) >
                    numpy.sum(numpy.abs(numpy.diff(b, axis=0)))):
                data = data.reshape(-1, 256).T
        elif order == 'F':
            data = data.reshape(-1, 256).T
        elif order != 'C':
            raise ValueError("unknown order", order)
        if data.shape[1] == 4 and numpy.all(data[:, 3] == 0):
            data[:, 3] = 255  # fix transparency
        return data


class SimfcsVpl(LfdFile):
    """SimFCS or ImObj color palette.

    SimFCS VPL files contain a single RGB color palette, stored as 256x3
    unsigned bytes in C or Fortran order, preceded by a 22 or 24 bytes header.

    Attributes
    ----------
    name : str or None
        Name of the palette.

    Examples
    --------
    >>> with SimfcsVpl('simfcs.vpl') as f:
    ...     print(f.asarray()[100])
    [189 210 246]
    >>> with SimfcsVpl('imobj.vpl') as f:
    ...    print(f.asarray()[100])
    [  0 254  27]

    """
    _filepattern = r'.*\.vpl'

    def _init(self):
        """Verify file size and header."""
        if self._filesize == 790:
            if self._fh.read(7) != b'vimage:':
                raise LfdFile.Error(self)
            self.name = bytes2str(stripnull(self._fh.read(15)))
        elif self._filesize == 792:
            self._fh.seek(24)
            self.name = None
        else:
            raise LfdFile.Error(self)
        self.shape = (256, 3)
        self.dtype = numpy.dtype('u1')

    def _data(self):
        """Return palette data as (256, 3) shaped array of uint8."""
        data = numpy.fromfile(self._fh, 'u1', 768)
        if self._filesize == 792:
            return data.reshape(256, 3)
        else:
            return data.reshape(3, 256).T

    def __str__(self):
        """Return name of palette."""
        return self.name


class SimfcsVpp(LfdFile):
    """SimFCS color palettes.

    SimFCS VPP files contain multiple BGRA color palettes, each stored as
    256x4 values of unsigned bytes preceded by a 24 byte Pascal string.

    Attributes
    ----------
    names : list of str
        Names of all palettes in file.

    Examples
    --------
    >>> with SimfcsVpp('simfcs.vpp') as f:
    ...     print(f.asarray('nice.vpl')[100])
    [ 16 255 239 255]

    """
    _filepattern = r'.*\.vpp'

    def _init(self):
        """Read list of palette names from file."""
        self.names = []
        while True:
            # read pascal string
            strlen = ord(self._fh.read(1))
            name = bytes2str(self._fh.read(23)[:strlen].lower())
            if not name.endswith('.vpl'):
                break
            self.names.append(name)
            self._fh.seek(1024, 1)
        if not self.names:
            raise LfdFile.Error(self)
        self.shape = (256, 4)
        self.dtype = numpy.dtype('u1')

    def _data(self, key=0, rgba=True):
        """Return palette data as (256, 4) shaped array of uint8.

        Parameters
        ----------
        key : int or str, optional
            The index or name of the palette to return.
        rgba : bool, optional
            Return RGBA palette if True (default), else BGRA.

        """
        try:
            self.names[key]  # validation
        except TypeError:
            key = self.names.index(key)
        self._fh.seek(key*1048 + 24)
        data = numpy.fromfile(self._fh, 'u1', 1024).reshape(256, 4)
        if rgba:
            data[:, :3] = data[:, 2::-1]
        if numpy.all(data[:, 3] == 0):
            data[:, 3] = 255  # fix transparency
        return data

    def __len__(self):
        return len(self.names)

    def _str(self):
        """Return names of all palettes in file."""
        return ", ".join(self.names)


class SimfcsJrn(LfdFile):
    """SimFCS journal.

    SimFCS JRN files contain metadata for several measurements, stored as
    key, value pairs in an unstructured ASCII format. Records usually start
    with lines of 80 '*' characters. The files do not contain array data.

    The metadata can be accessed as a list of dictionaries.

    Examples
    --------
    >>> with SimfcsJrn('simfcs.jrn', lower=True) as f:
    ...     f[1]['paramters for tracking']['samplimg frequency']
    15625

    """
    _filemode = 'rU'
    _filepattern = r'.*\.jrn'

    # regular expressions of all keys found in journal files
    _keys = r"""
        Image experiment
        Correlation expt
        Card type
        Channel for tracking
        Up down flag
        DATE
        TIME
        Dark \d+
        Movement type
        Mode
        DC threshold
        EWmission
        Emission
        Extension
        Frames integrated
        Frequency domain
        Int scale factor\d+
        Linescan point\s*\d+
        Maximum cycles
        Points per orbit
        Sampli[mn]g frequency
        R harmonics
        R%
        Mod for auto R
        Points per pixel
        Radius
        Dwell time
        Period
        Rperiods
        Sampling freq
        Scan type
        Scanner t-constant
        Scanner time const
        Scanner voltage
        Scanning period
        Scanning radius
        Time/photon mode
        Number or r-periods
        Cycles per particle
        Voltage full scale
        DC threshold
        Z[_-]radius
        Z[_-]Period
        frame_[xy]_range
        wait_end_of_frame
        wait_end_of_line
        [xyz][0o]_frame
        [xyz]_Offset
        [xyz]_Position
        [xyz]_Range\d*
        [PM]\d_\d calibration factor for FLIMbox
        FLIM BOX data
        roi (?:serial|parrallel|parallel) (?:start|width|bin)
        clear cycles
        clearing mode
        amplifier gain
        on-chip multiplier gain
        pmode
        shutter open mode
        exposure time
        frame time
        total time
        frames written
        """
    _keys = '|'.join(i.strip().replace(' ', '[ ]')
                     for i in _keys.splitlines() if i.strip())
    _keys = re.compile('(%s)' % _keys)
    _skip = re.compile(r'\s*\(.*\)')  # ignore parenthesis in values

    def _init(self, lower=False):
        """Read journal file and parse into list of dictionaries.

        Parameters
        ----------
        lower : bool
            Convert keys to lower case.

        """
        firstline = self._fh.readline()
        if not (firstline.startswith('*'*80) or
                firstline.startswith('roi')):
            raise LfdFile.Error("not a SimFCS journal file")
        content = self._fh.read()
        self._filesize = len(firstline) + len(content)
        if not firstline.startswith('***'):
            content = firstline + '\n' + content
        self._lower = lower = (lambda x: x.lower()) if lower else (lambda x: x)
        self._records = []
        for record in content.split('*'*80):
            recdict = {}
            record = record.split('COMMENTS', 1)
            if len(record) > 1:
                record, comments = record
            else:
                record, comments = record[0], ''
            recdict[lower('COMMENTS')] = comments.strip('* :;=\n\r\t')
            record = re.split("[*]{5}([\w\s]+)[*]{5}", record)
            for key, value in zip(record[1::2], record[2::2]):
                newdict = {}
                key = lower(key.strip('* :;=\n\r\t'))
                value = self._parse_journal(
                    value, self._keys, result=newdict, lower=lower)
                recdict[key] = newdict
            self._parse_journal(
                record[0], self._keys, result=recdict, lower=lower)
            self._records.append(recdict)
        self.close()

    def _data(self):
        """Raise ValueError."""
        raise ValueError("file does not contain array data")

    @staticmethod
    def _parse_journal(journal, repattern, result=None, lower=None):
        """Return dictionary of keys and values in journal string."""
        if result is None:
            result = {}
        if lower is None:
            def lower(x):
                return x
        keyval = re.split(repattern, journal, maxsplit=0)
        keyval = list(s.strip('* :;=\n\r\t') for s in keyval)
        val = [as_type(re.sub(SimfcsJrn._skip, '', v)) for v in keyval[2:-1:2]]
        key = [lower(k) for k in keyval[1:-1:2]]
        result.update(zip(key, val))
        return result

    def _str(self):
        """Return string with information about file."""
        lower = self._lower
        comments = lower('COMMENTS')
        result = []
        for i, record in enumerate(self._records):
            result.extend([
                "Record %i" % i,
                format_dict(record, excludes=[comments, '_']),
                lower("* COMMENTS:"),
                record[comments], '',
            ])
        return '\n'.join(result)

    def __getitem__(self, key):
        """Return selected record."""
        return self._records[key]

    def __len__(self):
        """Return number of records."""
        return len(self._records)

    def __iter__(self):
        """Return iterator over records."""
        return iter(self._records)


class SimfcsBin(LfdFile):
    """SimFCS raw binary data.

    SimFCS BIN and RAW files contain homogeneous data of any type and shape,
    stored C-contiguously in little endian order.
    A common format is: shape=(256, 256), dtype='uint16'.

    Attributes
    ----------
    shape : tuple of int
        Shape of data array contained in file.
    dtype : numpy dtype
        Type of data contained in file.

    Examples
    --------
    >>> with SimfcsBin('NPR_PMMA1001.bin', (256, 256), 'uint16') as f:
    ...     print(f.asarray()[751, 127, 127])
    1

    """
    _filepattern = r'.*\.(bin|raw)'

    def _init(self, shape, dtype):
        """Validate file size is a multiple of provided shape and type.

        Parameters
        ----------
        shape : int or sequence of int
            Shape of array to read from file.
        dtype : numpy dtype
            Desired output data-type.

        """
        self.shape = tuple(numpy.array(shape, int))
        self.dtype = numpy.dtype(dtype)
        size = product(self.shape) * self.dtype.itemsize
        if self._filesize % size:
            raise LfdFile.Error(self, "file size mismatch.")
        if self._filesize // size > 1:
            self.shape = (self._filesize // size, ) + self.shape

    def _data(self):
        """Return data as array of specified shape and type."""
        data = numpy.fromfile(self._fh, '<' + self.dtype.char)
        return data.reshape(*self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* shape: %s\n* dtype: %s" % (self.shape, self.dtype)

SimfcsRaw = SimfcsBin


class SimfcsInt(LfdFile):
    """SimFCS intensity image.

    SimFCS INT files contain a single intensity image, stored as 256x256
    float32 or uint16 (older format). The measurement extension is usually
    encoded in the file name.

    Attributes
    ----------
    shape : tuple of int
        Shape of data array contained in file. Always (256, 256).
    dtype : numpy dtype
        Type of data contained in file (float32 or uint16).


    Examples
    --------
    >>> with SimfcsInt('simfcs2036.int') as f:
    ...     print(f.asarray()[255, 255])
    3.0
    >>> with SimfcsInt('simfcs1006.int') as f:
    ...     print(f.asarray()[255, 255])
    9

    """
    _filepattern = r'.*\.int'

    def _init(self):
        """Validate file size is 256 KB"""
        if self._filesize == 262144:
            self.dtype = numpy.dtype('<f4')
        elif self._filesize == 131072:
            self.dtype = numpy.dtype('<u2')
        else:
            raise LfdFile.Error(self, "file size mismatch")
        self.shape = (256, 256)

    def _data(self):
        """Return data as 256x256 array of float32 or uint16."""
        return numpy.fromfile(self._fh, self.dtype).reshape(256, 256)

    def _str(self):
        """Return additional information about file."""
        return "* height: 256\n* width: 256\n* dtype: %s" % self.dtype


class SimfcsIntPhsMod(LfdFile):
    """SimFCS lifetime component images.

    SimFCS INT, PHS and MOD files contain fluorescence lifetime image data
    from frequency-domain measurements.
    Three 256x256 float32 images are stored in separate files:
    intensity (``*.int``), phase (``*.phs``) and modulation (``*.mod``).
    Phase values are in degrees, modulation in percent.
    The measurement extension and channel is often encoded in the file name.

    Examples
    --------
    >>> with SimfcsIntPhsMod('simfcs_1000.phs') as f:
    ...     print(f.asarray().mean(-1).mean(-1))
    [ 5.72  0.    0.05]

    """
    _filepattern = r'.*\.(int|phs|mod)'

    def _init(self):
        """Validate file size is 256 KB and other components exist."""
        if self._filesize != 262144:
            raise LfdFile.Error(self, "file size mismatch")
        for ext in ('int', 'phs', 'mod'):
            fname = os.path.join(self._filepath, self._filename[:-3] + ext)
            if not os.path.exists(fname):
                raise LfdFile.Error(self, "%s not found" % fname)
        self.dtype = numpy.dtype('<f4')
        self.shape = (256, 256)

    def _data(self):
        """Return image data as (3, 256, 256) shaped array of float32."""
        data = []
        for ext in ('int', 'phs', 'mod'):
            fname = os.path.join(self._filepath, self._filename[:-3] + ext)
            data.append(numpy.fromfile(fname, self.dtype))
        return numpy.array(data, '=f4').reshape(3, *self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* height: 256\n* width: 256\n* dtype: %s" % self.dtype


class SimfcsFit(LfdFile):
    """SimFCS fit data.

    SimFCS FIT files contain results from image scan analysis.
    The fit parameters are stored as a 1024x16 float64 array, followed by
    a 8 bytes buffer and the intensity image used for the fit, stored as
    a 256x256 float32 array.

    The 16 fit parameters are::

        W0
        Background
        Pixel size
        Triplet amplitude
        G1
        D1 (um2/s)
        G2
        D2 (um2/s)
        Exp aplitude
        Exp time/Ch1 int
        Triplet rate/Ch2 int
        Fraction vesicle
        Radius vesicle
        Velocity modulus
        Velocity x
        Velocity y

    Examples
    --------
    >>> with SimfcsFit('simfcs.fit') as f:
    ...     p_fit, dc_ref = f.asarray(size=7)
    ...     print(p_fit[6, 1, 1], dc_ref[128, 128])
    0.936775481992 20.23

    """
    _filepattern = r'.*\.fit'
    # type of data in file
    _record_t = numpy.dtype([('p_fit', '<f8', (1024, 16)), ('_', '<f8'),
                             ('dc_ref', '<f4', (256, 256))])
    # parameter labels
    labels = ('W0', 'Background', 'Pixel size', 'Triplet amplitude',
              'G1', 'D1 (um2/s)', 'G2', 'D2 (um2/s)', 'Exp aplitude',
              'Exp time/Ch1 int', 'Triplet rate/Ch2 int',
              'Fraction vesicle', 'Radius vesicle', 'Velocity modulus',
              'Velocity x', 'Velocity y')

    def _init(self):
        """Validate file size is 384 KB."""
        if self._filesize != 393224:
            raise LfdFile.Error(self, "size mismatch")

    def _data(self, size=32):
        """Return fit parameters and intensity image as numpy arrays.

        The array shape of the fit parameters will be (size, size, 16).

        Parameters
        ----------
        size : int
            Number of rows and columns of fit parameters array (default: 32).

        """
        if not 0 < size < 33:
            raise ValueError("size out of range [1..32]")
        p_fit = numpy.fromfile(self._fh, '<f8', 1024*16).reshape(1024, 16)
        p_fit = p_fit[:size*size].reshape(size, size, 16)
        self._fh.seek(8, 1)
        dc_ref = numpy.fromfile(self._fh, '<f4', 65536).reshape(256, 256)
        return p_fit, dc_ref


class SimfcsCyl(LfdFile):
    """SimFCS orbital tracking data.

    SimFCS CYL files contain intensity data from orbital tracking measurements,
    stored as a uint16 array of shape (2 channels, number of orbits,
    256 points per orbit).

    The number of channels and points per orbit can be read from the
    associated jounal file.

    Attributes
    ----------
    shape : tuple of 3 int
        Shape of data array contained in file:
        0: number or channels (1 or 2).
        1: number of orbits.
        2: number of samples per orbit (1 to 256).

    Examples
    --------
    >>> with SimfcsCyl('simfcs.cyl') as f:
    ...     print(f.shape, f.asarray()[0, -1, :4])
    (2, 1829, 256) [58 61 64 65]

    """
    _filepattern = r'.*\.cyl'

    def _init(self, points_per_orbit=256, channels=2):
        """Verify file size is multiple of 1024."""
        if channels > 2 or channels < 1:
            raise ValueError("channels out of range [1..2]")
        if points_per_orbit > 256 or points_per_orbit < 1:
            raise ValueError("points_per_orbit out of range [1..256]")
        if self._filesize % 1024:
            raise LfdFile.Error(self)
        self.shape = (channels, self._filesize // 1024, points_per_orbit)
        self.dtype = numpy.dtype('<u2')

    def _data(self):
        """Return data as (channels, -1, points_per_orbit) array of uint16."""
        a = numpy.fromfile(self._fh, self.dtype)
        return a.reshape(2, -1, 256)[:self.shape[0], :, :self.shape[2]]

    def _str(self):
        """Return additional information about file."""
        return "* channels: %i\n* orbits: %i\n* points_per_orbit: %i" % (
            self.shape)


class SimfcsRef(LfdFile):
    """SimFCS referenced fluorescence lifetime images.

    SimFCS REF files contain referenced fluorescence lifetime image data.
    Five 256x256 float32 images are stored consecutively:
    0) dc - intensity
    1) ph1 - phase of 1st harmonic
    2) md1 - modulation of 1st harmonic
    3) ph2 - phase of 2nd harmonic
    4) md2 - modulation of 2nd harmonic
    Phase values are in degrees, the modulation values are normalized.
    Phase and modulation values may be NaN.

    Examples
    --------
    >>> with SimfcsRef('simfcs.ref') as f:
    ...     print(f.asarray()[:, 255, 255])
    [ 301.33   44.71    0.62   68.13    0.32]

    """
    _filepattern = r'.*\.ref'

    def _init(self):
        """Verify file size is 1280 KB."""
        if not self._filesize == 1310720:
            raise LfdFile.Error(self)
        self.shape = (5, 256, 256)
        self.dtype = numpy.dtype('<f4')

    def _data(self):
        """Return images as (5, 256, 256) shaped array of float32."""
        return numpy.fromfile(self._fh, self.dtype).reshape(self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* images: %i\n* height: %i\n* width: %i" % self.shape


class SimfcsBh(LfdFile):
    """SimFCS Becker and Hickl fluorescence lifetime histogram.

    SimFCS B&H files contain time-domain fluorescence lifetime histogram data,
    acquired from Becker and Hickl(r) TCSPC cards, or converted from other
    data sources.
    The data are stored consecutively as 256 bins of 256x256 float32 images.
    SimFCS BHZ files are zipped B&H files: a Zip archive containing a single
    B&H file.
    BHZ files are occasionally used to store consecutive 256x256 float32
    images, e.g. volume data.

    Attributes
    ----------
    iszip : bool
        True, if file is a ZIP/BHZ file.
    shape : tuple of 3 int
        Shape of data array contained in file. Usually (256, 256, 256).

    Examples
    --------
    >>> with SimfcsBh('simfcs.b&h') as f:
    ...     print(f.asarray()[59, 1, 84])
    12.0
    >>> with SimfcsBh('simfcs.bhz') as f:
    ...     print(f.asarray()[59, 1, 84])
    12.0

    """
    _filepattern = r'.*\.(b&h|bhz)'

    def _init(self):
        """Verify file size or Zip file contains B&H file."""
        self.iszip = self._filename.lower().endswith('bhz')
        if self.iszip:
            with zipfile.ZipFile(self._fh) as zf:  # requires Python >= 2.7
                try:
                    filesize = zf.filelist[0].file_size
                except (zipfile.BadZipfile, IndexError, AttributeError):
                    raise LfdFile.Error(self)
        else:
            filesize = self._filesize
        if filesize % 262144:
            raise LfdFile.Error(self)
        self.shape = (filesize // 262144, 256, 256)
        self.dtype = numpy.dtype('<f4')

    def _data(self):
        """Return image data as (256, 256, 256) shaped array of float32."""
        if self.iszip:
            with zipfile.ZipFile(self._fh) as zf:
                data = zf.read(zf.filelist[0])
            data = numpy.fromstring(data, self.dtype)
        else:
            data = numpy.fromfile(self._fh, self.dtype)
        return data.reshape(self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* bins: %i\n* height: %i\n* width: %i" % self.shape

SimfcsBhz = SimfcsBh


class SimfcsB64(LfdFile):
    """SimFCS integer intensity image(s).

    SimFCS B64 files contain one or more square intensity images.
    The image data are stored as int16 contiguously after one int32 defining
    the image size in x and y dimensions.
    The measurement extension is usually encoded in the file name.

    Attributes
    ----------
    shape : tuple of int
        Shape of data array contained in file.
    dtype : numpy dtype
        Type of data contained in file. Always int16.

    Examples
    --------
    >>> with SimfcsB64('simfcs.b64') as f:
    ...     print(f.shape, f.asarray()[101, 255, 255])
    (102, 256, 256) 0

    """
    _filepattern = r'.*\.b64'

    def _init(self, dtype='<i2', maxsize=1024):
        """Read file header."""
        size = struct.unpack('<i', self._fh.read(4))[0]
        if maxsize < size < 2:
            raise LfdFile.Error(self, "image size out of range")
        self.shape = (size, size)
        self.dtype = numpy.dtype(dtype)
        # determine number of images in file
        size = product(self.shape) * self.dtype.itemsize
        fsize = self._filesize - 4
        if fsize % size:
            raise LfdFile.Error(self, "file size mismatch.")
        if fsize // size > 1:
            self.shape = (fsize // size, ) + self.shape

    def _data(self):
        """Return image data as 2D or 3D array of int16."""
        data = numpy.fromfile(self._fh, '<' + self.dtype.char)
        return data.reshape(*self.shape)

    def _str(self):
        """Return additional information about file."""
        return "shape: %s\n* dtype: %s" % (self.shape, self.dtype)


class SimfcsI64(LfdFile):
    """SimFCS compressed intensity image.

    SimFCS I64 files contain a single square intensity image, stored as a
    zlib deflate compressed stream of one int32 (defining the image size in
    x and y dimensions) and the float32 image data.
    The measurement extension is usually encoded in the file name.

    Attributes
    ----------
    shape : tuple of int
        Shape of data array contained in file.
    dtype : numpy dtype
        Type of data contained in file. Always float32.

    Examples
    --------
    >>> with SimfcsI64('simfcs1000.i64') as f:
    ...     print(f.shape, f.asarray()[511, 511])
    (512, 512) 0.0

    """
    _filepattern = r'.*\.i64'

    def _init(self, dtype='<f4', maxsize=1024):
        """Read file header."""
        if 67108864 < self._filesize < 32:  # limit to 64 MB
            raise LfdFile.Error(self, "file size out of range")
        size = struct.unpack('<i', self._decompress_header(4))[0]
        if maxsize < size < 2:
            raise LfdFile.Error(self, "image size out of range")
        self.shape = (size, size)
        self.dtype = numpy.dtype(dtype)

    def _data(self):
        """Return data as 2D array of float32."""
        bufsize = product(self.shape) * self.dtype.itemsize + 4
        data = zlib.decompress(self._fh.read(), 15, bufsize)
        return numpy.fromstring(data[4:], self.dtype).reshape(*self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* height: %i\n* width: %i\n* dtype: %s" % (
            self.shape[0], self.shape[1], self.dtype)


class SimfcsZ64(LfdFile):
    """SimFCS compressed image stack.

    SimFCS Z64 files contain stacks of square images such as intensity volumes
    or time-domain fluorescence lifetime histograms acquired from
    Becker and Hickl(r) TCSPC cards.
    The data are stored as zlib deflate compressed stream of two int32
    (defining the image size in x and y dimensions and the number of images)
    and a maximum of 256 square float32 images.

    Attributes
    ----------
    dtype : numpy dtype
        Type of data contained in file. Always float32.
    shape : tuple of 3 int
        Shape of data array contained in file.

    Examples
    --------
    >>> with SimfcsZ64('simfcs.z64') as f:
    ...     print(f.shape, f.asarray()[128, 128, 128])
    (256, 1024, 1024) 0.0

    """
    _filepattern = r'.*\.z64'

    def _init(self, dtype='<f4', maxsize=(256, 1024)):
        """Read file header."""
        if self._filesize < 32:
            raise LfdFile.Error(self, "file size out of range")
        size, inum = struct.unpack('<ii', self._decompress_header(8))[:2]
        if maxsize[-1] < size < 2 or maxsize[0] < inum < 2:
            raise LfdFile.Error(self, "image size out of range")
        self.shape = (inum, size, size)
        self.dtype = numpy.dtype(dtype)

    def _data(self):
        """Return data as 3D array of float32."""
        bufsize = product(self.shape) * self.dtype.itemsize + 8
        data = zlib.decompress(self._fh.read(), 15, bufsize)
        data = numpy.fromstring(data[8:], '<' + self.dtype.char)
        return data.reshape(*self.shape)

    def _str(self):
        """Return additional information about file."""
        return "*depth: %i\n * height: %i\n* width: %i\n* dtype: %s" % (
            self.shape[0], self.shape[1], self.shape[2], self.dtype)


class SimfcsR64(LfdFile):
    """SimFCS compressed referenced fluorescence lifetime images.

    SimFCS R64 files contain referenced fluorescence lifetime images.
    The data are stored as a zlib deflate compressed stream of one int32
    (defining the image size in x and y dimensions) and five square float32
    images:
    0) dc - intensity
    1) ph1 - phase of 1st harmonic
    2) md1 - modulation of 1st harmonic
    3) ph2 - phase of 2nd harmonic
    4) md2 - modulation of 2nd harmonic
    Phase values are in degrees, the modulation values are normalized.
    Phase and modulation values may be NaN.

    Attributes
    ----------
    dtype : numpy dtype
        Type of data contained in file. Always float32.
    shape : tuple of 3 int
        Shape of data array contained in file.

    Examples
    --------
    >>> with SimfcsR64('simfcs.r64') as f:
    ...     print(f.shape, f.asarray()[:, 100, 200])
    (5, 256, 256) [   0.25   23.22    0.64  104.33    2.12]

    """
    _filepattern = r'.*\.r64'

    def _init(self, dtype='<f4', maxsize=1024):
        """Read file header."""
        if self._filesize < 32:
            raise LfdFile.Error(self, "file size out of range")
        size = struct.unpack('<i', self._decompress_header(4))[0]
        if maxsize < size < 2:
            raise LfdFile.Error(self, "image size out of range")
        self.shape = (5, size, size)
        self.dtype = numpy.dtype(dtype)

    def _data(self):
        """Return data as 3D array of float32."""
        bufsize = product(self.shape) * self.dtype.itemsize + 4
        data = zlib.decompress(self._fh.read(), 15, bufsize)
        data = numpy.fromstring(data[4:], '<' + self.dtype.char)
        return data.reshape(*self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* height: %i\n* width: %i\n* dtype: %s" % (
            self.shape[1], self.shape[2], self.dtype)


class SimfcsMap(LfdFile):
    """SimFCS volume data.

    SimFCS MAP files contain 3D volume data stored in CCP4 map format
    used by the Electron Microscopy Data Bank.
    <http://emdatabank.org/mapformat.html>
    <http://www.ccp4.ac.uk/html/maplib.html>
    <ftp://ftp.wwpdb.org/pub/emdb/doc/map_format/EMDB_mapFormat_v1.0.pdf>

    Attributes
    ----------
    shape : tuple of 3 int
        Shape of data array contained in file.
    start : tuple of 3 int
        Position of first section, row, and column (voxel grid units).
    cell_interval : tuple of 3 int
        Intervals per unit cell repeat along Z, Y, X.
    cell_length : tuple of 3 float
        Unit Cell repeats along Z, Y, X in Angstroms.
    cell_angle : tuple of 3 float
        Unit Cell angles (alpha, beta, gamma, in degrees).
    map_src : tuple of 3 int
        Relationship of Z, Y, X axes to sections, rows, columns
    density_min, density_max, density_mean: float
        Minimum, maximum, average density.
    density_rms : float
        Rms deviation of map from mean density.
    spacegoup : int
        IUCr space group number (1-230).
    skew_matrix, skew_translation : ndarray or None
        Skew matrix and translation vector (if any, else None).
    symboltable : list of byte strings
       Symmetry records as defined in International Tables.

    Examples
    --------
    >>> with SimfcsMap('simfcs.ccp4') as f:
    ...     print(f.asarray()[100, 100, 100])
    1.0

    """
    _filepattern = r'.*\.(map|ccp4)'
    _dtypes = {0: 'i1', 1: 'i2', 2: 'f4', 4: 'q8', 5: 'i1'}

    def _init(self):
        """Read CCP4 file header and symboltable."""
        header = self._fh.read(1024)
        if (header[208:212] != b'MAP '):
            raise LfdFile.Error(self, " %s" % header[:32])
        try:
            (nc, nr, ns,
             mode,  # data type
             ncstart, nrstart, nsstart,
             nx, ny, nz,
             x_length, y_length, z_length,
             alpha, beta, gamma,
             mapc, mapr, maps,
             self.density_min, self.density_max, self.density_mean,
             self.spacegoup,
             nsymbt,  # number of bytes used for storing symmetry operators
             skew_matrix_flag,
             S11, S12, S13, S21, S22, S23, S31, S32, S33,
             T1, T2, T3,
             # extra,
             map_,  # b'MAP '
             machst,  # machine stamp,
             self.density_rms,
             nlabl,  # number of labels used
             L0, L1, L2, L3, L4, L5, L6, L7, L8, L9
             ) = struct.unpack('3ii3i3i3f3f3i3fiii9f3f60x4s4sfi'
                               '80s80s80s80s80s80s80s80s80s80s', header)
        except struct.error as e:
            raise LfdFile.Error(self, e)
        try:
            # machst = header[212:216]
            byteorder = {b'DA\x00\x00': '<', b'\x11\x11\x00\x00': '>'}[machst]
        except KeyError:
            byteorder = '='
            warnings.warn("SimfcsMap: unknown machine stamp: %s" % machst)
        try:
            self.dtype = byteorder + SimfcsMap._dtypes[mode]
        except KeyError:
            raise LfdFile.Error(self, "unknown mode: %s" % mode)
        self.shape = ns, nr, nc
        self.start = nsstart, nrstart, ncstart
        self.cell_interval = nz, ny, nx
        self.cell_length = z_length, y_length, x_length
        self.cell_angle = alpha, beta, gamma
        self.map_src = maps, mapr, mapc
        if skew_matrix_flag != 0:
            self.skew_translation = numpy.array([T1, T2, T3], 'float64')
            self.skew_matrix = numpy.array([[S11, S12, S13],
                                            [S21, S22, S23],
                                            [S31, S32, S33]], 'float64')  # .T?
        else:
            self.skew_translation = None
            self.skew_matrix = None
        if 0 <= nlabl <= 10:
            self.labels = [stripnull(lbl) for lbl in
                           (L0, L1, L2, L3, L4, L5, L6, L7, L8, L9)[:nlabl]]
        else:
            self.labels = []
        if nsymbt < 0 or nsymbt % 80:
            raise LfdFile.Error(self, "invalid symbol table size: %i" % nsymbt)
        self.symboltable = [stripnull(self._fh.read(80))
                            for i in range(nsymbt // 80)]

    def _data(self, memmap=False):
        """Return volume data as numpy array.

        Parameters
        ----------
        memmap : bool
            If True, use numpy.memmap to read array.

        """
        if memmap:
            return numpy.memmap(self._fh, dtype=self.dtype, mode='r',
                                offset=self._pos, shape=self.shape)
        data = numpy.fromfile(self._fh, self.dtype, product(self.shape))
        return data.reshape(self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* shape: %s\n* length: %s\n* dtype: %s" % (
            self.shape, self.cell_length, self.dtype)


def save_map(filename, data, start=(0, 0, 0), cell_interval=None,
             cell_length=None, cell_angle=(90, 90, 90), map_src=(3, 2, 1),
             density=None, density_rms=0, spacegroup=1,
             skew_matrix=None, skew_translation=None, symboltable=b'',
             labels=[b'Created by lfdfiles.py']):
    """Save 3D volume data to CCP4 MAP formatted file.

    Parameters
    ----------
    filename : str
        Name of file to write.
    data : 3D array_like
        Input volume.

    Refer to the SimfcsMap attributes for other parameters.

    Examples
    --------
    >>> data = numpy.arange(1000000).reshape(100, 100, 100).astype('f4')
    >>> save_map('test.ccp4', data)

    """
    data = numpy.asarray(data)
    if data.ndim != 3:
        raise ValueError("data must be 3 dimensional")
    try:
        mode = {'i1': 0, 'i2': 1, 'f4': 2, 'q8': 4}[data.dtype.str[-2:]]
    except KeyError:
        raise ValueError("dtype not supported by MAP format")
    if cell_interval is None:
        cell_interval = data.shape
    if cell_length is None:
        cell_length = data.shape  # ?
    if density is None:
        density = numpy.min(data), numpy.max(data), numpy.mean(data)
    if skew_matrix is None or skew_translation is None:
        skew_matrix_flag = 0
        S = numpy.zeros((3, 3))
        T = numpy.zeros(3)
    else:
        skew_matrix_flag = 1
        S = numpy.array(skew_matrix)
        T = numpy.array(skew_translation)
    assert(S.shape == (3, 3))
    assert(T.shape == (3,))
    labels = list(labels)[:10] if labels else []
    nlabl = len(labels)
    labels = [l[:79] for l in labels]
    labels.extend(b'' for _ in range(10 - nlabl))

    header = struct.pack(
        '3ii3i3i3f3f3i3fiii9f3f60x4s4sfi80s80s80s80s80s80s80s80s80s80s',
        data.shape[2], data.shape[1], data.shape[0],
        mode,
        start[2], start[1], start[0],
        cell_interval[2], cell_interval[1], cell_interval[0],
        cell_length[2], cell_length[1], cell_length[0],
        cell_angle[0], cell_angle[1], cell_angle[2],
        map_src[2], map_src[1], map_src[0],
        density[0], density[1], density[2],
        spacegroup,
        len(symboltable),
        skew_matrix_flag,
        S[0, 0], S[0, 1], S[0, 2],
        S[1, 0], S[1, 1], S[1, 2],
        S[2, 0], S[2, 1], S[2, 2],
        T[0], T[1], T[2],
        # extra,
        b'MAP ',
        {'little': b'DA\x00\x00', 'big': b'\x11\x11\x00\x00'}[sys.byteorder],
        density_rms,
        len(labels),
        labels[0], labels[1], labels[2], labels[3], labels[4],
        labels[5], labels[6], labels[7], labels[8], labels[9],
    )
    with open(filename, 'wb') as fh:
        fh.write(header)
        fh.write(symboltable)
        data.tofile(fh)


class SimfcsFbf(LfdFile):
    """FlimBox firmware.

    SimFCS FBF files contain FlimBox device firmwares, stored in binary form
    following a NULL terminated ASCII string containing properties and
    description.

    The properties (lower cased) can be accessed via dictionary interface.

    Examples
    --------
    >>> with SimfcsFbf('FlimExt10_Ver1_2_2.fbf') as f:
    ...     f['windows'], f['channels'], f['secondharmonic'], 'extclk' in f
    (16, 2, 0, True)

    """
    _filepattern = r'.*\.(fbf)'

    def _init(self):
        """Read and parse NULL terminated header string."""
        try:
            # the first 1024 bytes contain the header
            header, _ = self._fh.read(1024).split(b'\x00', 1)
        except ValueError:
            raise LfdFile.Error(self)
        header = bytes2str(header)
        self._fh.seek(len(header) + 1)
        try:
            header, comment = header.rsplit('/', 1)
            comment = [comment]
        except ValueError:
            comment = []
        self._settings = {}
        for (name, value, unit) in re.findall(
                r'([a-zA-Z\s]*)([.\d]*)([a-zA-Z\d]*)/', header):
            name = name.strip().lower()
            if not name:
                name = {'w': 'windows', 'ch': 'channels'}.get(unit, None)
                unit = ''
            if not name:
                comment.append(name + value + unit)
                continue
            if unit == 'MHz':
                unit = 1000000
            try:
                if unit:
                    value = int(value) * int(unit)
                else:
                    value = int(value)
                unit = 0
            except ValueError:
                pass
            self._settings[name] = (value + unit) if value != '' else True
        self._settings['comment'] = '/'.join(reversed(comment))

    def _data(self):
        """Return firmware as binary string."""
        return self._fh.read()

    def _str(self):
        """Return string with header settings."""
        return format_dict(self._settings)

    def __getitem__(self, key):
        return self._settings[key]

    def __contains__(self, key):
        return key in self._settings

    def __len__(self):
        return len(self._settings)

    def __iter__(self):
        return iter(self._settings)


class SimfcsFbd(LfdFile):
    """FlimBox data.

    SimFCS FDB files contain raw data from the FlimBox device, storing a
    stream of 16 bit integers (data words) that can be decoded to photon
    arrival windows, channels, and times.

    The measurement's frame size, pixel dwell time, number of sampling
    windows, and scanner type are encoded in the last four letters of the
    file name.
    Newer FBD files, where the 3rd character in the file name tag is '0',
    start with the first 1kB of the firmware file used for the measurement,
    followed by 31kB containing a binary record with measurement settings.

    It depends on the application and its setting how to interpret the
    decoded data, e.g. as time series, line scans, or image frames of FCS
    or digital frequency domain fluorescence lifetime measurements.

    The data word format depends on the device's firmware.
    A common layout is::

        |F|E|D|C|B|A|9|8|7|6|5|4|3|2|1|0|  data word bits
                            |-----------|  pcc (cross correlation phase)
                        |---------------|  tcc (cross correlation time)
                      |-|                  marker (indicates start of frame)
        |-------------|                    index into decoder table

    The data word can be decoded into a cross correlation phase histogram
    index (shown for the 1st harmonics)::

        bin = (pmax-1 - (pcc + win * (pmax//windows)) % pmax) // pdiv

    ``bin``
        Cross correlation phase index (phase histogram bin number).
    ``pcc``
        Cross correlation phase (counter).
    ``pmax``
        Number of entries in cross correlation phase histogram.
    ``pdiv``
        Divisor to reduce number of entries in phase histogram.
    ``win``
        Arrival window.
    ``windows``
        Number of sampling windows.

    Attributes
    ----------
    frame_size : int
        Number of pixels/samples in one line scan, excluding retrace.
    pixel_dwell_time : float
        Number of microseconds the scanner remains at each pixel/sample.
    windows : int
        Number of sampling windows used by flimbox.
    channels : int
        Number of channels used by flimbox.
    harmonics : int
        First or second harmonics.
    pmax : int
        Number of entries in cross correlation phase histogram.
    pdiv : int
        Divisor to reduce number of entries in phase histogram.
    laser_frequency : float
        Laser frequency in Hz. Default is 20000000 Hz, the internal
        flimbox frequency.
    correction_factor : float
        Factor to correct dwell_time/laser_frequency when the scanner
        clock is not known exactly. Default is 1.0.
    scanner : str
        Acquisition software or hardware.
    scanner_line_add : int
        Number of pixels/samples added to each line (for retrace).
    scanner_line_start : int
        Index of first valid pixel/sample in scan line.
    scanner_frame_start : int
        Index of first valid pixel/sample after marker.

    Examples
    --------
    >>> with SimfcsFbd('simfcs$CBCO.fbd') as f:
    ...     bins, times, markers = f.asarray(word_count=500000,
    ...                                      skip_words=1900000)
    >>> print(bins[0, :2], times[:2], markers)
    [53 51] [ 0 42] [ 44097 124815]
    >>> hist = [numpy.bincount(b[b>=0]) for b in bins]
    >>> numpy.argmax(hist[0])
    53

    """
    _filepattern = r'.*\.fbd'

    _attributes = (
        'frame_size', 'windows', 'channels', 'harmonics', 'pmax', 'pdiv',
        'pixel_dwell_time', 'laser_frequency', 'correction_factor',
        'scanner', 'scanner_line_add', 'scanner_line_start',
        'scanner_frame_start')

    _frame_size = {
        # Map 1st character in file name tag to image frame size
        'A': 64, 'B': 128, 'C': 256, 'D': 320, 'E': 512,
        'F': 640, 'G': 800, 'H': 1024}

    _flimbox_settings = {
        # Map 3rd character in file name tag to (windows, channels, harmonics)
        # '0': file contains header
        'A': (2, 2, 1),
        'B': (4, 2, 1),
        'C': (8, 2, 1),
        'F': (8, 4, 1),
        'D': (16, 2, 1),
        'E': (32, 2, 1),
        'H': (64, 1, 1),
        # second harmonics. P and Q might be switched in some files?
        'N': (2, 2, 2),
        'O': (4, 2, 2),  # frequency = 40000000 ?
        'P': (8, 2, 2),
        'G': (8, 4, 2),
        'Q': (16, 2, 2),
        'R': (32, 2, 2),
        'S': (64, 1, 2),
    }

    _histogram_settings = {
        # Map (windows, channels, harmonics) to (pmax, pdiv)
        (2, 2, 1): (256, 4),
        (4, 2, 1): (256, 4),
        (8, 2, 1): (64, 1),
        (8, 4, 1): (64, 1),
        (16, 2, 1): (64, 1),
        # (32, 2, 1): (256, 4), ?
        # (64, 1, 1): (64, 1), ?
        # second harmonics
        (2, 2, 2): (256, 4),
        (4, 2, 2): (128, 2),
        (8, 2, 2): (32, 1),
        (8, 4, 2): (32, 1),
        (16, 2, 2): (32, 1),
        # (32, 2, 2): (128, 2), ?
        # (64, 1, 2): (32, 1), ?
    }

    _scanner_settings = {
        # Map 4th and 2nd character in file name tag to pixel_dwell_time,
        # scanner_line_add, scanner_line_start, and scanner_frame_start.
        # As of SimFCS ~2011.
        # These values may not apply any longer and need to be overridden.
        'S': {
            'name': 'Native SimFCS, 3-axis card',
            'A': (4, 198, 99, 0), 'J': (100, 20, 10, 0),
            'B': (5, 408, 204, 0), 'K': (128, 16, 8, 0),
            'C': (8, 256, 128, 0), 'L': (200, 10, 5, 0),
            'D': (10, 204, 102, 0), 'M': (256, 8, 4, 0),
            'E': (16, 128, 64, 0), 'N': (500, 4, 2, 0),
            'F': (20, 102, 51, 0), 'O': (512, 4, 2, 0),
            'G': (32, 64, 32, 0), 'P': (1000, 2, 1, 0),
            'H': (50, 40, 20, 0), 'Q': (1024, 2, 1, 0),
            'I': (64, 32, 16, 0), 'R': (2000, 1, 0, 0)},
        'O': {
            'name': 'Olympus FV 1000, NI USB',
            'A': (2, 10, 8, 0), 'F': (20, 56, 55, 0),
            'B': (4, 10, 8, 0), 'G': (40, 28, 20, 0),
            'C': (8, 10, 8, 0), 'H': (50, 12, 10, 0),
            'D': (10, 112, 114, 0), 'I': (100, 12, 9, 0),
            'E': (12.5, 90, 80, 0), 'J': (200, 10, 89, 0)},
        'Y': {
            'name': 'Zeiss LSM510',
            'A': (6.39, 344, 22, 0), 'D': (51.21, 176, 22, 414),
            'B': (12.79, 344, 22, 0), 'E': (102.39, 88, 0, 242),
            'C': (25.61, 344, 22, 0), 'F': (204.79, 12, 10, 0)},
        'Z': {
            'name': 'Zeiss LSM710',
            'A': (6.39, 344, 2, 0), 'D': (51.21, 176, 8, 414),
            'B': (12.79, 344, 2, 0), 'E': (102.39, 88, 10, 242),
            'C': (25.61, 344, 2, 0), 'F': (204.79, 12, 10, 0)},
        'I': {
            'name': 'ISS Vista slow scanner',
            'A': (4, 112, 73, 0), 'H': (32, 37, 28, 0),
            'B': (6, 112, 73, 0), 'I': (40, 28, 19, 0),
            'C': (8, 112, 73, 0), 'J': (64, 18, 14, 0),
            'D': (10, 112, 73, 0), 'K': (200, 6, 4, 0),
            'E': (12.5, 90, 59, 0), 'L': (500, 6, 4, 0),
            'F': (16, 73, 48, 0), 'M': (1000, 6, 4, 0),
            'G': (20, 56, 37, 0)},
        'V': {
            'name': 'ISS Vista fast scanner',
            'A': (4, 112, 73, 0), 'H': (32, 37, 28, 0),
            'B': (6, 112, 73, 0), 'I': (40, 28, 19, 0),
            'C': (8, 112, 73, 0), 'J': (64, 21, 14, 0),
            'D': (10, 112, 73, 0), 'K': (100, 12, 8, 0),
            'E': (12.5, 90, 59, 0), 'L': (200, 6, 4, 0),
            'F': (16, 73, 48, 0), 'M': (500, 6, 4, 0),
            'G': (20, 56, 37, 0), 'N': (1000, 6, 4, 0)},
        'T': {
            # used with new file format only?
            'name': 'IOTech scanner card'}
    }

    _header_t = [
        # Binary header starting at offset 1024 in files with $xx0x names.
        # This is written as a memory dump of a Delphi record, hence the pads
        ('owner', '<i4'),  # must be 0
        ('pixel_dwell_time_index', '<i4'),
        ('frame_size_index', '<i4'),
        ('line_length', '<i4'),
        ('points_end_of_frame', '<i4'),
        ('x_starting_pixel', '<i4'),
        ('line_integrate', '<i4'),
        ('scanner_index', '<i4'),
        ('synthesizer_index', '<i4'),
        ('windows_index', '<i4'),
        ('channels_index', '<i4'),
        ('_pad1', '<i4'),
        ('line_time', '<f8'),
        ('frame_time', '<f8'),
        ('scanner', '<i4'),
        ('_pad2', '<f4'),
        ('laser_frequency', '<f8'),
        ('laser_factor', '<f8'),
        ('frames_to_average', '<i4'),
        ('enabled_before_start', '<i4'),
        ('mypcalib', '<16f4'),
        ('mymcalib', '<16f4'),
        ('mypcalib1', '<16f4'),
        ('mymcalib1', '<16f4'),
        ('mypcalib2', '<16f4'),
        ('mymcalib2', '<16f4'),
        ('mypcalib3', '<16f4'),
        ('mymcalib3', '<16f4'),
        ('h1', '<i4'),
        ('h2', '<i4'),
        ('process_enable', 'b'),
        ('integrate', 'b'),
        ('detect_maitai', 'b'),
        ('trigger_on_up', 'b'),
        ('line_scan', 'b'),
        ('circular_scan', 'b'),
        ('average_frames_on_reading', 'b'),
        ('show_each_frame', 'b'),
        ('write_each_frame', 'b'),
        ('write_one_big_file', 'b'),
        ('subtract_background', 'b'),
        ('normalize_to_frames', 'b'),
        ('second_harmonic', 'b'),
        ('_pad3', '3b'),
        ('bin_pixel_by_index', '<i4'),
        ('jitter_level', '<i4'),
        ('phaseshift1', '<i4'),
        ('phaseshift2', '<i4'),
        ('acquire_item_index', '<i4'),
        ('acquire_number_of_frames', '<i4'),
        ('show_channel_index', '<i4')]

    _header_pixel_dwell_time = (4, 5, 8, 10, 16, 20, 32, 50, 64, 100,
                                128, 200, 256, 500, 512, 1000, 1, 2)
    _header_frame_size = (64, 128, 256, 320, 512, 640, 800, 1024)
    _header_bin_pixel_by = (1, 2, 4, 8)
    _header_windows = (4, 8, 16, 32, 64)
    _header_channels = (1, 2, 4)  # '2 to 4 ch', '2 ch 8w Spartan6'
    _header_scanner_types = (
        'None',
        'Native SimFCS, 3-axis card',
        'Olympus FV 1000, NI USB',
        'Zeiss LSM510',
        'Zeiss LSM710',
        'ISS Vista slow scanner',
        'ISS Vista fast scanner',
        'M2 laptop only',
        'IOTech scanner card')
    _header_synthesizer_names = (
        'Internal FLIMbox frequency',
        'Fianium',
        'Spectra Physics MaiTai',
        'Spectra Physics Tsunami',
        'Coherent Chameleon')

    @lazyattr
    def decoder_settings(self):
        """Return dictionary of parameters to decode flimbox data stream.

        Returns
        -------
        decoder_table : ndarray of int16 and shape (channels, window indices)
            Decoder table, mapping channel and window indices to actual
            arrival windows.
        tcc_mask, tcc_shr : int
            Binary mask and number of bits to right shift in order to extract
            cross correlation time from data word.
        pcc_mask, pcc_shr : int
            Binary mask and number of bits to right shift in order to extract
            cross correlation phase from data word.
        marker_mask, marker_shr: int
            Binary mask and number of bits to right shift in order to extract
            markers from data word.
        win_mask, win_shr : int
            Binary mask and number of bits to right shift in order to extract
            index into lookup table from data word.

        """
        return dict(zip(
            ('decoder_table', 'tcc_mask', 'tcc_shr', 'pcc_mask', 'pcc_shr',
             'marker_mask', 'marker_shr', 'win_mask', 'win_shr'),
            getattr(self, '_w%ic%i' % (self.windows, self.channels))()))

    @staticmethod
    def _w4c2():
        # Return parameters to decode 4 windows, 2 channels flimbox data.
        a = [[-1, 0, -1, 1, -1, 2, -1, 3, -1, 0, -1, -1, 1, 0, -1, 2,
              1, 0, 3, 2, 1, 0, 3, 2, 1, -1, 3, 2, -1, -1, 3, -1],
             [-1, -1, 0, -1, 1, -1, 2, -1, 3, 0, -1, -1, 0, 3, -1, 0,
              1, 2, 2, 1, 2, 1, 1, 2, 3, -1, 0, 3, -1, -1, 3, -1]]
        a = numpy.array(a, 'int16')
        return a, 0x3ff, 0, 0xff, 0, 0x400, 10, 0xf800, 11

    @staticmethod
    def _w8c2():
        # Return parameters to decode 8 windows, 2 channels flimbox data.
        a = numpy.zeros((2, 81), 'int16') - 1
        a[0, 1:9] = range(8)
        a[1, 9:17] = range(8)
        a[:, 17:] = numpy.mgrid[0:8, 0:8].reshape(2, -1)[::-1, :]
        return a, 0xff, 0, 0x3f, 0, 0x100, 8, 0xffff, 9

    @staticmethod
    def _w8c4():
        # Return parameters to decode 8 windows, 4 channels flimbox data.
        # Not tested.
        a = numpy.zeros((4, 128), 'int16') - 1
        for i in range(128):
            win = (i & 0b0000111)
            ch0 = (i & 0b0001000) >> 3
            ch1 = (i & 0b0010000) >> 4
            ch2 = (i & 0b0100000) >> 5
            ch3 = (i & 0b1000000) >> 6
            if ch0 + ch1 + ch2 + ch3 != 1:
                continue
            if ch0:
                a[0, i] = win
            elif ch1:
                a[1, i] = win
            elif ch2:
                a[2, i] = win
            elif ch3:
                a[3, i] = win
        return a, 0xff, 0, 0x3f, 0, 0x100, 8, 0xffff, 9

    @staticmethod
    def _w16c1():
        # Return parameters to decode 16 windows, 1 channel flimbox data.
        # Not tested.
        a = numpy.zeros((1, 32), 'int16') - 1
        for i in range(32):
            win = (i & 0b11110) >> 1
            ch0 = (i & 0b00001)
            if ch0:
                a[0, i] = win
        return a, 0xff, 0, 0x3f, 0, 0x100, 8, 0xffff, 11

    @staticmethod
    def _w16c2():
        # Return parameters to decode 16 windows, 2 channels flimbox data.
        # Not tested.
        a = numpy.zeros((2, 64), 'int16') - 1
        for i in range(64):
            win = (i & 0b111100) >> 2
            ch0 = (i & 0b000010) >> 1
            ch1 = (i & 0b000001)
            if ch0 + ch1 != 1:
                continue
            if ch0:
                a[0, i] = win
            elif ch1:
                a[1, i] = win
        return a, 0xff, 0, 0x3f, 0, 0x100, 8, 0xffff, 10

    @staticmethod
    def _w32c2():
        # Return parameters to decode 32 windows, 2 channels flimbox data.
        # TODO
        raise NotImplementedError()

    @staticmethod
    def _w64c1():
        # Return parameters to decode 64 windows, 1 channel flimbox data.
        # TODO
        raise NotImplementedError()

    def _init(self, code=None, **kwargs):
        """Initialize attributes from file name code and additional parameters.

        Parameters
        ----------
        code : str (optional)
            Four character string, encoding frame size (1st char),
            pixel dwell time (2nd char), number of sampling windows (3rd char),
            and scanner type (4th char).
            By default this will be extracted from the file name.
        kwargs : dict
            Optional named parameters to override instance attributes.

        """
        if not code:
            try:
                code = re.search(r'.*\$([A-Z0]{4,4})\.fbd',
                                 self._filename, re.IGNORECASE).group(1)
            except AttributeError:
                pass
        if not code:
            code = 'CFCS'  # old flimbox file ?
            warnings.warn("failed to detect settings from file name;"
                          "assuming a '%s' file" % code)

        if code[2] == '0':
            # New FBD file format with header
            # the first 1024 bytes contain the start of a flimbox firmware file
            with SimfcsFbf(self._fh.name) as fbf:
                self._fbf = fbf._settings
            # the next 31kB contain the binary file header
            self._fh.seek(1024, 0)
            self._header = numpy.fromfile(self._fh, self._header_t, 1)[0]
            hdr = self._header
            if hdr['owner'] != 0:
                raise ValueError("unknown header format `%i`" % hdr['owner'])
            self.frame_size = self._header_frame_size[hdr['frame_size_index']]
            self.scanner = self._header_scanner_types[hdr['scanner_index']]
            self.pixel_dwell_time = self._header_pixel_dwell_time[
                hdr['pixel_dwell_time_index']]
            self.scanner_line_add = int(hdr['line_length']) - self.frame_size
            self.scanner_line_start = int(hdr['x_starting_pixel'])
            self.scanner_frame_start = 0
            self.windows = self._header_windows[hdr['windows_index']]
            self.channels = self._header_channels[hdr['channels_index']]
            self.harmonics = (1, 2)[hdr['second_harmonic']]
            self.laser_frequency = float(hdr['laser_frequency'])
            self.correction_factor = float(hdr['laser_factor'])
            assert self.windows == self._fbf['windows']
            assert self.channels == self._fbf['channels']
            assert hdr['second_harmonic'] == self._fbf['secondharmonic']
            # encoded data start at offset 32768
            self._data_offset = 32768
        else:
            self.frame_size = self._frame_size[code[0]]
            self.scanner = self._scanner_settings[code[3]]['name']
            (self.pixel_dwell_time,
             self.scanner_line_add,
             self.scanner_line_start,
             self.scanner_frame_start
             ) = self._scanner_settings[code[3]][code[1]]
            (self.windows,
             self.channels,
             self.harmonics
             ) = self._flimbox_settings[code[2]]
            self.laser_frequency = 20000000
            self.correction_factor = 1.0
            self._data_offset = 0

        i = self.windows, self.channels, self.harmonics
        try:
            self.pmax, self.pdiv = self._histogram_settings[i]
        except IndexError:
            raise NotImplementedError(
                "can't handle %i windows, %i channels, %i harmonics" % i)

        # override attributes with those specified by user
        for k, v in kwargs.items():
            assert k in self._attributes
            setattr(self, k, v)

        assert self.harmonics in (1, 2)
        self.laser_frequency *= self.harmonics

    def units_per_sample(self):
        """Return number of flimbox units per scanner sample."""
        return ((self.pixel_dwell_time / 1000000) *
                (self.pmax / (self.pmax - 1) *
                 self.laser_frequency * self.correction_factor))

    def decode(self, data=None, word_count=-1, skip_words=0, max_markers=65536,
               **kwargs):
        """Read flimbox data stream from file and return decoded data.

        Parameters
        ----------
        data : ndarray or None
            Flimbox data stream. If None (default), read data from file.
        word_count : int (optional)
            Number of data words to process (default: -1, all words).
        skip_words : int (optional)
            Number of data words to skip at beginning of stream (default: 0).
        max_markers : int (optional)
            Maximum number of markers expected in data stream (default: 65536).

        Returns
        -------
        bins : ndarray of int8 and shape (channels, size)
            The cross correlation phase index for all channels and data points.
            A value of -1 means no photon was counted.
        times : ndarray of uint64 or uint32
            The times in flimbox counter units at each data point.
            A potentially huge array.
        markers : ndarray of ssize_t
            The indices of up markers in the data stream, usually indicating
            frame starts.

        """
        if data is None:
            self._fh.seek(self._data_offset + skip_words*2, 0)
            data = numpy.fromfile(self._fh, dtype='<u2', count=word_count)
        elif skip_words or word_count != -1:
            if word_count < 0:
                data = data[skip_words: word_count]
            else:
                data = data[skip_words: skip_words+word_count]

        bins_out = numpy.empty((self.channels, data.size), dtype=numpy.int8)
        times_out = numpy.empty(data.size, dtype=numpy.uint64)
        markers_out = numpy.zeros(max_markers, dtype=numpy.intp)

        simfcsfbd_decode(
            data, bins_out, times_out, markers_out,
            self.windows, self.pmax, self.pdiv, self.harmonics,
            **self.decoder_settings)

        return bins_out, times_out, markers_out[markers_out > 0]

    def frames(self, decoded, select_frames=None, aspect_range=(0.8, 1.2),
               frame_cluster=0, **kwargs):
        """Return shape and start/stop indices of scanner frames.

        If unable to detect any frames using the default settings, try to
        determine a correction factor from clusters of frame durations.

        Parameters
        ----------
        decoded : tuple of 3 ndarrays, or None
            Times and markers as returned by the decode function.
            If None, the decode function is called.
        select_frames : slice (optional)
            Specifies which image frames to return (default: all frames).
        aspect_range : tuple (optional)
            Minimum and maximum aspect ratios of valid frames. The default
            lets 1:1 aspects pass.
        frame_cluster : int
            Index of the frame duration cluster to use when calculating
            the correction factor.

        Returns
        -------
        shape : tuple
            Dimensions of scanner frame.
        frame_markers : list of tuples of 2 int
            Start and stop indices of detected image frames.

        """
        if decoded is None:
            decoded = self.decode(**kwargs)
        times, markers = decoded[-2:]

        line_len = self.frame_size + self.scanner_line_add
        line_time = line_len * self.units_per_sample()
        frame_durations = numpy.ediff1d(times[markers])

        frame_markers = []
        if aspect_range:
            # detect frame markers assuming correct settings
            line_num = sys.maxsize
            for i, duration in enumerate(frame_durations):
                lines = duration / line_time
                aspect = self.frame_size / lines
                if aspect_range[0] < aspect < aspect_range[1]:
                    frame_markers.append((int(markers[i]),
                                          int(markers[i+1]) - 1))
                else:
                    continue
                if line_num > lines:
                    line_num = lines
            line_num = int(round(line_num))

        if not frame_markers:
            # calculate frame duration clusters, assuming few clusters that
            # are narrower and more separated than cluster_size.
            cluster_size = 1024
            clusters = []
            cluster_indices = []
            for d in frame_durations:
                d = int(d)
                for i, c in enumerate(clusters):
                    if abs(d - c[0]) < cluster_size:
                        cluster_indices.append(i)
                        c[0] = min(c[0], d)
                        c[1] += 1
                        break
                else:
                    cluster_indices.append(len(clusters))
                    clusters.append([d, 1])
            clusters = list(sorted(clusters, key=lambda x: x[1], reverse=True))
            # possible correction factors, assuming square frame shape
            line_num = self.frame_size
            correction_factors = [c[0]/(line_time*line_num) for c in clusters]
            # select specified frame cluster
            frame_cluster = min(frame_cluster, len(correction_factors)-1)
            self.correction_factor = correction_factors[frame_cluster]
            frame_markers = [(int(markers[i]), int(markers[i+1])-1)
                             for i, c in enumerate(cluster_indices)
                             if c == frame_cluster]
            msg = ["no frames detected with default settings."
                   "Using square shape and correction factor %.5f." % (
                       self.correction_factor)]
            if len(correction_factors) > 1:
                msg.append("The most probable correction factors are: %s" % (
                    ", ".join("%.5f" % i for i in correction_factors[:4])))
            warnings.warn('\n'.join(msg))

        if not isinstance(select_frames, slice):
            select_frames = slice(select_frames)
        frame_markers = frame_markers[select_frames]
        if not frame_markers:
            raise ValueError("no frames selected")
        return (line_num, line_len), frame_markers

    def asimage(self, decoded, frames, integrate_frames=1, square_frame=True,
                **kwargs):
        """Return image histograms from decoded data and detected frames.

        This function may fail to produce expected results when settings
        were recorded incorrectly, scanner and flimbox frequencies were out
        of sync, or scanner settings were changed during acquisition.

        Parameters
        ----------
        decoded : tuple of 3 ndarrays, or None
            Bins, times, and markers as returned by the decode function.
            If None, the decode function is called.
        frames : tuple or None
            Scanner_shape and frame_markers as returned by the frames function.
            If None, the frames function is called.
        integrate_frames : int
            Specifies which frames to sum. By default (1), all frames are
            summed into one. If 0, no frames are summed.
        square_frame : bool
            If True (default), return square image (frame_size x frame_size),
            else return full scanner frame.

        Returns
        -------
        result : 5D ndarray of uint16
            Image histogram of shape (number of frames, channels in bins
            array, detected line numbers, frame_size, histogram bins).

        """
        if decoded is None:
            decoded = self.decode(**kwargs)
        bins, times, markers = decoded
        if frames is None:
            frames = self.frames(decoded, **kwargs)
        scanner_shape, frame_markers = frames
        # an extra line to scanner frame to allow clipping indices
        scanner_shape = scanner_shape[0] + 1, scanner_shape[1]
        # allocate output array of scanner frame shape
        shape = (integrate_frames if integrate_frames else len(frame_markers),
                 bins.shape[0],  # channels
                 scanner_shape[0] * scanner_shape[1],
                 self.pmax // self.pdiv)
        result = numpy.zeros(shape, dtype='uint16')
        # calculate frame data histogram
        simfcsfbd_histogram(
            bins, times, frame_markers, self.units_per_sample(),
            self.scanner_frame_start, result)
        # reshape frames and slice valid region
        result = result.reshape(shape[:2] + scanner_shape + shape[-1:])
        if square_frame:
            result = result[..., :self.frame_size, self.scanner_line_start:
                            self.scanner_line_start+self.frame_size, :]
        return result

    def _data(self, **kwargs):
        """Read flimbox data stream from file and return decoded data."""
        return self.decode(data=None, **kwargs)

    def _str(self):
        """Return additional information about file."""
        return "\n".join("* %s: %s" % (a, getattr(self, a))
                         for a in self._attributes)


def simfcsfbd_decode(data, bins_out, times_out, markers_out,
                     windows, pmax, pdiv, harmonics, decoder_table,
                     tcc_mask, tcc_shr, pcc_mask, pcc_shr,
                     marker_mask, marker_shr, win_mask, win_shr):
    """Decode flimbox data stream.

    See the documentation of the SimfcsFbd class for parameter descriptions
    and the lfdfiles.pyx file for a faster implementation.

    """
    tcc = data & tcc_mask  # cross correlation time
    if tcc_shr:
        tcc >>= tcc_shr
    times_out[:] = tcc
    times_out[times_out == 0] = (tcc_mask >> tcc_shr) + 1
    times_out[0] = 0
    times_out[1:] -= tcc[:-1]
    del tcc
    numpy.cumsum(times_out, out=times_out)

    markers = data & marker_mask
    markers = numpy.diff(markers.view('int16'))
    markers = numpy.where(markers > 0)[0]  # trigger up
    markers += 1
    size = min(len(markers), len(markers_out))
    markers_out[:size] = markers[:size]
    del markers

    if win_mask != 0xffff:  # window index
        win = data & win_mask
        win >>= win_shr
    else:
        win = data >> win_shr
    win = decoder_table.take(win, axis=1)
    nophoton = win == -1
    win *= pmax // windows * harmonics
    pcc = data & pcc_mask  # cross correlation phase
    if pcc_shr:
        pcc >>= pcc_shr
    win += pcc
    del pcc
    win %= pmax
    win += 1 - pmax
    numpy.negative(win, win)
    if pdiv > 1:
        win //= pdiv
    win = win.astype('int8')
    win[nophoton] = -1
    bins_out[:] = win


def simfcsfbd_histogram(bins, times, frame_markers, units_per_sample,
                        scanner_frame_start, out):
    """Calculate histograms from decoded flimbox data and frame markers.

    See the documentation of the SimfcsFbd class for parameter descriptions
    and the lfdfiles.pyx file for a much faster implementation.

    """
    warnings.warn("patience, please.")
    nframes, nchannels, frame_length, nwindows = out.shape
    for f, (j, k) in enumerate(frame_markers):
        f = f % nframes
        t = times[j:k] - times[j]
        t /= units_per_sample  # index into flattened array
        if scanner_frame_start:
            t -= scanner_frame_start
        t = t.astype('uint32')
        numpy.clip(t, 0, frame_length-1, out=t)
        for c in range(nchannels):
            d = bins[c, j:k]
            for w in range(nwindows):
                x = numpy.where(d == w)[0]
                x = t.take(x)
                x = numpy.bincount(x, minlength=frame_length)
                out[f, c, :, w] += x

try:
    from bilab.io._lfdfiles import simfcsfbd_histogram, simfcsfbd_decode
except ImportError:
    warnings.warn("failed to import _lfdfiles module; using slow alternative.")


class SimfcsGpSeq(LfdFileSequence):
    """SimFCS generalized polarization image sequence.

    SimFCS GP sequences contain intensity images from two channels, stored in
    separate SimfcsInt files with consecutive names.

    Examples
    --------
    >>> ims = SimfcsGpSeq("gpint/v*.int").asarray()
    >>> ims.shape
    (2, 135, 256, 256)

    """
    _readfunction = SimfcsInt
    _indexpattern = r'(\d)(\d+)\.int'


class GlobalsLif(LfdFile):
    """Globals binary lifetime data.

    Globals LIF files contain array and meta data of multiple frequency-domain
    cuvette lifetime measurement, stored as consecutive 472 byte records.
    The number of frequencies per record is limited to 25. The format was
    also used by ISS software.

    Examples
    --------
    >>> with GlobalsLif("globals.lif") as lif:
    ...     print(len(lif), lif[42]['date'], lif[42].asarray().shape)
    43 1987.8.8 (5, 11)

    """
    _filepattern = r'.*\.lif'

    _record_t = numpy.dtype([
        ('_title_len', 'u1'),
        ('title', 'a80'),
        ('number', 'i2'),
        ('frequency', [('_len', 'u1'), ('str', 'a6')], 25),
        ('phase', 'i2', 25),
        ('modulation', 'i2', 25),
        ('deltap', 'i2', 25),
        ('deltam', 'i2', 25),
        ('nanal', 'i2'),
        ('date', 'i2', 3),
        ('time', 'i2', 3)])

    class Record(dict):
        def asarray(self):
            return numpy.array((
                self['frequency'],
                self['phase'],
                self['modulation'],
                self['deltap'],
                self['deltam']))

        def __str__(self):
            return format_dict(self)

    def _init(self):
        """Verify file size and read all records."""
        if (self._filesize % 472
                or self._filesize == 0
                or self._filesize // 472 > 1024):
            raise LfdFile.Error(self)
        self._records = records = []
        for rec in numpy.rec.fromfile(self._fh, self._record_t):
            number = int(rec['number'])
            if number == 0:
                continue
            elif number > 25:
                warnings.warn('corrupted record')
                continue
            record = self.Record()
            record['number'] = number
            record['title'] = bytes2str(rec['title'][:rec['_title_len']])
            record['nanal'] = int(rec['nanal'])
            record['date'] = "%i.%i.%i" % tuple(rec['date'])
            record['time'] = "%i:%i:%i" % tuple(rec['time'])
            record['frequency'] = numpy.array(
                [float(f[:i].strip()) for i, f in rec['frequency'][:number]],
                dtype='float64')
            record['phase'] = rec['phase'][:number] / 100.0
            record['modulation'] = rec['modulation'][:number] / 100.0
            record['deltap'] = rec['deltap'][:number] / 100.0
            record['deltam'] = rec['deltam'][:number] / 100.0
            records.append(record)

    def _data(self, key):
        """Return freq, phi, mod, dp, dm of selected record as numpy array."""
        return self._records[key].asarray()

    def _str(self):
        return '\n'.join("Record %i\n%s" % (i, format_dict(r))
                         for i, r in enumerate(self._records))

    def __getitem__(self, key):
        """Return selected record."""
        return self._records[key]

    def __len__(self):
        """Return number of records."""
        return len(self._records)

    def __iter__(self):
        """Return iterator over records."""
        return iter(self._records)


class GlobalsAscii(LfdFile):
    """Globals ASCII lifetime data.

    Globals ASCII files contain array and meta data of a single frequency
    domain lifetime measurement, stored as human readable ASCII string.
    Consecutive measurements are stored in separate files with increasing file
    extension numbers. The format is also used by ISS and FLOP software.

    The metadata can be accessed via dictionary getitem interface.

    Examples
    --------
    >>> with GlobalsAscii('FLOP.001') as f:
    ...    print(f['experiment'], f.asarray().shape)
    LIFETIME (5, 20)

    """
    _filemode = 'rU'
    _filepattern = r'.*\.(\d){3}'

    def _init(self, lower=True):
        """Read file and parse into dictionary and data array.

        Parameters
        ----------
        lower : bool (True)
            Convert keys to lower case and replace spaces by underscores.

        """
        if not self._fh.read(5) == 'TITLE':
            raise LfdFile.Error("not a Globals ASCII data file")

        self._fh.seek(0)
        content = self._fh.read()
        self._filesize = len(content)
        self._lower = lower = (lambda x: x.lower().replace(' ', '_')
                               if lower else (lambda x: x))
        # parse keys and values
        self._record = record = {}
        matches = re.findall(r"(.*?)(?:\((.*)\))?:\s(.*)",
                             content, re.IGNORECASE)
        for key, unit, value in matches:
            key = lower(key.strip())
            unit = unit.strip()
            value = as_type(value.strip())
            record[key] = value
            if unit:
                record[key + '_unit'] = unit
        # extract data array
        data = re.search(r"DATA:(.*)[\r\n]([-+\.\d\s\r\n]*)ENDDATA",
                         content, re.IGNORECASE)
        headers = [d.strip().title() for d in data.group(1).split(', ')]
        record[lower('DATA')] = headers
        data = numpy.fromstring(data.group(2), sep=' ')
        self.__data = data.reshape((len(headers), -1), order='F')
        record[lower('DATA_SHAPE')] = self.__data.shape
        self.close()

    def _data(self):
        """Return freq, phase, modulation, deltap, deltam as numpy array."""
        return self.__data.copy()

    def _str(self):
        """Return string with information about file."""
        return format_dict(self._record)

    def __getitem__(self, key):
        """Return value of key in record."""
        return self._record[key]


class VistaIfli(LfdFile):
    """VistaVision fluorescence lifetime image.

    VistaVision IFLI files contain phasor and lifetime images for several
    time points, channels, slices, and frequencies from analog or digital
    frequency domain fluorescence lifetime measurements.
    After a header of 1024 bytes, the images are stored as two 7 dimensional
    arrays of float32 and shape (time, channel, Z, Y, X, frequency, sample).
    The phasor array has three samples: the average intensity (DC) and the
    real (g) and imaginary (s) parts of the phasor.
    The lifetime array has two samples, the lifetimes calculated from phase
    and modulation (apparent single lifetimes).

    Attributes
    ----------
    header : numpy.recarray
        File header.
    shape : tuple of 6 int
        The number of times, channels, Z, Y, X coordinates, and frequencies.

    Examples
    --------
    >>> ifli = VistaIfli('vista.ifli')
    >>> print(ifli.header.frequencies)
    [ 48000000.  96000000.]
    >>> ifli.asarray().shape
    (1, 2, 1, 128, 128, 2, 3)
    >>> ifli.close()

    """
    _filepattern = r'.*\.ifli'

    def _header_t(self, number_frequencies):
        return numpy.dtype([
            ('signature', 'a12'),  # 'VistaFLImage'
            ('version', 'u1'),
            ('channel_bits', 'u1'),
            ('compression', 'u1'),
            ('dimensions', 'u2', 5),  # TCZYX
            ('boundaries', 'f4', 6),  # x0, x1, y0, y1, z0, z1
            ('coordinate', 'u1'),  # unknown, px, um, mm
            ('pixel_sampling_time', 'f4'),  # milliseconds
            ('pixel_interval_time', 'f4'),
            ('line_interval_time', 'f4'),
            ('frame_interval_time', 'f4'),
            ('number_frequencies', 'i4'),
            ('frequencies', 'f4', number_frequencies),
            ('cross_frequency', 'f4'),
            ('reference_lifetime', 'f4'),  # ns
            ('reference_phasor', 'f4', (number_frequencies, 3))])

    def _init(self):
        """Read header and metadata from file."""
        if self._fh.read(12) != b'VistaFLImage':
            raise LfdFile.Error(self)
        self._fh.seek(66)
        numfreq = struct.unpack('<i', self._fh.read(4))[0]
        self._fh.seek(0)
        self.header = h = numpy.rec.fromfile(
            self._fh, self._header_t(numfreq), 1, byteorder='<')[0]
        assert h['version'] == 1
        assert h['compression'] == 0
        self.shape = tuple(
            int(i) for i in reversed(h['dimensions'])) + (numfreq, )
        self.axes = 'TCZYXF'

    def __getitem__(self, key):
        """Return header attribute."""
        return self.header[key]

    @lazyattr
    def channel_indices(self):
        """Return indices of valid channels."""
        bits = self.header['channel_bits']
        dims = self.header['dimensions'][-2]
        return tuple(i for i in range(dims) if bits & 2**i)

    def phasor(self):
        """Return average intensity and phasor coordinates from file.

        The returned array is of type float32 and shape (time, channel, Z, Y,
        X, frequency, sample). The three samples are the average intensity (DC)
        and the real (g) and imaginary (s) parts of the phasor.

        """
        self._fh.seek(1024)
        shape = self.shape + (3, )
        phasor = numpy.fromfile(self._fh, '<f4', product(shape))
        phasor.shape = shape
        return phasor

    def lifetime(self):
        """Return apparent single lifetime data from file.

        The returned array is of type float32 and shape (time, channel, Z, Y,
        X, frequency, sample). It has two samples, the lifetimes calculated
        from phase and modulation.

        """
        self._fh.seek(1024 + product(self.shape) * 3)
        shape = self.shape + (2, )
        lifetime = numpy.fromfile(self._fh, '<f4', product(shape))
        lifetime.shape = shape
        return lifetime

    def _data(self):
        """Return phasor data from file."""
        return self.phasor()


class FlimfastFlif(LfdFile):
    """FlimFast fluorescence lifetime image.

    FlimFast FLIF files contain camera images and metadata of frequency-domain
    fluorescence lifetime measurements.
    A 640 bytes header is followed by a variable number of uint16 images,
    each preceded by a 64 bytes record.

    Attributes
    ----------
    header : numpy.recarray
        File header.
    records : numpy.recarray
        Image headers.
    shape : tuple of 3 int
        Shape of data array in file.

    Examples
    --------
    >>> flif = FlimfastFlif('flimfast.flif')
    >>> flif.header.frequency
    80.652
    >>> flif.records['phase'][31]
    348.75
    >>> flif.asarray()[31, 219, 299]
    366
    >>> flif.close()

    """
    _filepattern = r'.*\.flif'

    _header_t = numpy.dtype([
        ('magic', 'a8'),  # '\211FLF\r\n0\n'
        ('creator', 'a120'),
        ('date', 'a32'),
        ('comments', 'a351'),
        ('_', 'u1'),
        ('fileprec', '<u2'),
        ('records', '<u2'),
        ('phases', '<i4'),
        ('width', '<i4'),
        ('height', '<i4'),
        ('dataprec', '<i4'),
        ('background', '<i4'),
        ('camframes', '<i4'),
        ('cambin', '<i4'),
        ('roileft', '<i4'),
        ('roitop', '<i4'),
        ('frequency', '<f4'),
        ('ref_tauphase', '<f4'),
        ('measured_phase', '<f4'),
        ('measured_mod', '<f4'),
        ('start', '<u4'),
        ('duration', '<u4'),
        ('phaseoffset', '<f4'),
        ('ref_taumod', '<f4')])  # ('padding', 'i4', 14)

    _record_t = numpy.dtype([
        ('index', '<i4'),
        ('order', '<i4'),
        ('phase', '<f4'),
        ('integrated', '<i4'),
        ('time', '<u4')])  # ('padding', 'u4', 11)

    def _init(self):
        """Read header and record metadata from file."""
        h = self.header = numpy.rec.fromfile(self._fh, self._header_t, 1,
                                             byteorder='<')[0]
        h["creator"] = stripnull(h.creator)
        h["date"] = stripnull(h.date)
        h["comments"] = stripnull(h.comments)
        if ((h.magic != b'\211FLF\r\n0\n')
                or (4097 < h.width) or (h.width < 1)
                or (4097 < h.height) or (h.height < 1)
                or (1025 < h.phases) or (h.phases < 3)
                or (1025 < h.records) or (h.records < 2)):
            raise LfdFile.Error(self)
        self.records = numpy.recarray((h.phases, ), self._record_t)
        stride = 11*4 + 2*h.width*h.height  # padding + image
        self._fh.seek(14*4, 1)  # header padding
        for i in range(h.phases):
            self.records[i] = numpy.rec.fromfile(self._fh, self._record_t, 1,
                                                 byteorder='<')[0]
            self._fh.seek(stride, 1)
        self.shape = (h.records, h.height, h.width)
        self.dtype = numpy.dtype('<u2')

    def _data(self):
        """Return images as (records, height, width) shaped uint16 array."""
        p, h, w = self.shape
        data = numpy.empty((p, h * w), 'uint16')
        self._fh.seek(640 + 64)
        for i in range(p):
            data[i] = numpy.fromfile(self._fh, self.dtype, h * w)
            self._fh.seek(64, 1)
        data.shape = self.shape
        return data

    def _str(self):
        """Return file header as string."""
        return "\n".join(("* %s: %s" % (name, getattr(self.header, name)))[:79]
                         for name in self._header_t.names[1:]
                         if not name.startswith('_'))


class FlimageBin(LfdFile):
    """FLImage fluorescence lifetime image.

    FLImage BIN files contain referenced fluorescence lifetime image data
    from frequency-domain measurements.
    Three 300x220 big endian float32 images are stored in separate files:
    intensity (``*.int.bin``), phase (``*.phi.bin``) and modulation
    (``*.mod.bin``), respectively single apparent lifetimes from phase
    (``*.tph.bin``) and modulation (``\*.tmd.bin``).
    Phase values are in degrees, modulation in percent.

    Examples
    --------
    >>> with FlimageBin('flimage.int.bin') as f:
    ...     print(f.asarray(True)[:, 219, 299])
    [   1.23  111.8    36.93]

    """
    _filepattern = r'.*\.(int|mod|phi|tmd|tph)\.bin'

    def _init(self):
        """Verify file size."""
        if not self._filesize == 264000:
            raise LfdFile.Error(self)
        self.shape = (220, 300)
        self.dtype = numpy.dtype('>f4')

    def _data(self, component=False):
        """Return images as (220, 300) shaped array of float32.

        Parameters
        ----------
        component : bool, optional
           If True, return image data from all component files.

        """
        if component:
            data = []
            for ext in ('int.bin', 'phi.bin', 'mod.bin', 'tph.bin', 'tmd.bin'):
                fname = os.path.join(self._filepath, self._filename[:-7] + ext)
                try:
                    data.append(numpy.fromfile(fname, self.dtype))
                except IOError:
                    pass
        else:
            data = numpy.fromfile(self._fh, self.dtype)
        return numpy.array(data, '=f4').reshape(-1, 220, 300).squeeze()

    def _str(self):
        return "* height: %i\n* width: %i" % self.shape


class FlieOut(LfdFile):
    """Flie fluorescence lifetime image.

    Flie OUT files contain referenced fluorescence lifetime image data
    from frequency-domain measurements.
    Three 300x220 big endian float32 images are stored in separate files:
    intensity (``off_*.out``), phase (``phi_*.out``), and modulation
    (``mod_*.out``). Phase values are in degrees, modulation in percent.
    No metadata are available.

    Examples
    --------
    >>> with FlieOut('off_flie.out') as f:
    ...     print(f.asarray(True)[:, 219, 299])
    [ 91.85  28.24  69.03]

    """
    _filepattern = r'(off|phi|mod)_.*\.out'

    def _init(self):
        """Verify file size."""
        if not self._filesize == 264000:
            raise LfdFile.Error(self)
        self.shape = (220, 300)
        self.dtype = numpy.dtype('>f4')

    def _data(self, component=False):
        """Return image data as (220, 300) shaped array of float32.

        Parameters
        ----------
        component : bool, optional
           If True, return image data from all off, phi, mod component files.

        """
        if component:
            data = []
            for pre in ('Off_', 'Phi_', 'Mod_'):
                fname = os.path.join(self._filepath, pre + self._filename[4:])
                try:
                    data.append(numpy.fromfile(fname, self.dtype))
                except IOError:
                    pass
        else:
            data = numpy.fromfile(self._fh, self.dtype)
        return numpy.array(data, '=f4').reshape(-1, 220, 300).squeeze()

    def _str(self):
        return "* height: %i\n* width: %i" % self.shape


class FliezI16(LfdFile):
    """FLIez integer image.

    FLIez I16 files contain camera images, usually for one phase cycle of
    frequency-domain fluorescence lifetime measurements.
    A number of 256x256 uint16 intensity images is stored consecutively.
    No metadata are available.

    Attributes
    ----------
    shape : tuple of 3 int
        Shape of data array contained in file.

    Examples
    --------
    >>> with FliezI16('fliez.i16') as f:
    ...     print(f.shape, f.asarray()[::8, 108, 104])
    (32, 256, 256) [401 538 220 297]

    """
    _filepattern = r'.*\.i16'

    def _init(self):
        """Verify file size is 128 KB."""
        if self._filesize % 131072:
            raise LfdFile.Error(self)
        self.shape = (self._filesize // 131072, 256, 256)
        self.dtype = numpy.dtype('<u2')

    def _data(self):
        """Return images as (-1, 256, 256) shaped array of uint16."""
        return numpy.fromfile(self._fh, self.dtype).reshape(self.shape)

    def _str(self):
        return "* images: %i\n* height: %i\n* width: %i" % self.shape


class FliezDb2(LfdFile):
    """FLIez double image.

    FLIez DB2 files contain a sequence of images from fluorescence lifetime
    measurements. The modality of the data stored in different files varies:
    phase intensities, average intensities, phase or modulation.
    After a header specifying the 3D data shape, images are stored
    consecutively as float64. No other metadata are available.

    Attributes
    ----------
    shape : tuple of int
        Shape of data array contained in file.

    Examples
    --------
    >>> with FliezDb2('fliez.db2') as f:
    ...     print(f.shape, f.asarray()[8, 108, 104])
    (32, 256, 256) 234.0

    """
    _filepattern = r'.*\.db2'

    def _init(self):
        """Read data shape and verify file size."""
        shape = struct.unpack('<iii', self._fh.read(12))
        if self._filesize - 12 != product(shape) * 8:
            raise LfdFile.Error(self)
        self.shape = shape[::-1]
        self.dtype = numpy.dtype('<f8')

    def _data(self):
        """Return images as 3D array of float64."""
        return numpy.fromfile(self._fh, self.dtype).reshape(self.shape)

    def _str(self):
        """Return additional information about file."""
        return "* images: %i\n* height: %i\n* width: %i" % self.shape


def stripnull(s):
    """Return string truncated at the first NULL character."""
    try:
        return s[:s.index(b'\x00')]
    except ValueError:
        return s


def product(iterable):
    """Return product of sequence of numbers.

    Equivalent of functools.reduce(operator.mul, iterable, 1).

    >>> product([2**8, 2**30])
    274877906944
    >>> product([])
    1

    """
    prod = 1
    for i in iterable:
        prod *= i
    return prod


def as_type(value, types=None):
    """Return argument as one of types if possible."""
    if types is None:
        types = int, float, bytes2str
    for typ in types:
        try:
            return typ(value)
        except (ValueError, TypeError, UnicodeEncodeError):
            pass
    return value


def format_dict(adict, prefix='', indent='  ', bullets=('* ', '* '),
                excludes=('_', ), linelen=79):
    """Return pretty-print of nested dictionary."""
    result = []
    for k, v in sorted(adict.items(), key=lambda x: x[0].lower()):
        if any(k.startswith(e) for e in excludes):
            continue
        if isinstance(v, dict):
            v = '\n' + format_dict(v, prefix=prefix+indent, excludes=excludes)
            result.append(prefix + bullets[1] + '%s: %s' % (k, v))
        else:
            result.append(
                (prefix + bullets[0] + '%s: %s' % (k, v))[:linelen].rstrip())
    return '\n'.join(result)


def natural_sorted(iterable):
    """Return human sorted list of strings.

    >>> natural_sorted(['f1', 'f2', 'f10'])
    ['f1', 'f2', 'f10']

    """
    def sortkey(x):
        return [(int(c) if c.isdigit() else c) for c in re.split(numbers, x)]
    numbers = re.compile(r'(\d+)')
    return sorted(iterable, key=sortkey)


if sys.version_info[0] > 2:
    def bytes2str(x):
        return str(x, 'cp1252')
    basestring = str
else:
    bytes2str = str


if __name__ == "__main__":
    import doctest
    os.chdir('data')
    numpy.set_printoptions(suppress=True, precision=2)
    doctest.testmod()
