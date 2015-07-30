from __future__ import division, print_function, absolute_import

import warnings
import numpy as np

from bilab.ml.NDR.tSNE import _bhtsne_wrap
import collections
import operator
import sys
import types

PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str
    integer_types = int,
    class_types = type,
    text_type = str
    binary_type = bytes

    MAXSIZE = sys.maxsize

else:
    string_types = basestring,
    integer_types = (int, long)
    class_types = (type, types.ClassType)
    text_type = unicode
    binary_type = str

    if sys.platform.startswith("java"):
        # Jython always uses 32 bits.
        MAXSIZE = int((1 << 31) - 1)
    else:
        # It's possible to have sizeof(long) != sizeof(Py_ssize_t).
        class X(object):
            def __len__(self):
                return 1 << 31
        try:
            len(X())
        except OverflowError:
            # 32-bit
            MAXSIZE = int((1 << 31) - 1)
        else:
            # 64-bit
            MAXSIZE = int((1 << 63) - 1)
            del X

# Replacement for lazy loading stuff in upstream six.  See gh-2764
if PY3:
    import builtins
    import functools
    reduce = functools.reduce
    zip = builtins.zip
    xrange = builtins.range
else:
    import __builtin__
    import itertools
    builtins = __builtin__
    reduce = __builtin__.reduce
    zip = itertools.izip
    xrange = __builtin__.xrange


if PY3:
    def get_unbound_function(unbound):
        return unbound

    Iterator = object

    def callable(obj):
        return any("__call__" in klass.__dict__ for klass in type(obj).__mro__)
else:
    def get_unbound_function(unbound):
        return unbound.im_func

    class Iterator(object):

        def next(self):
            return type(self).__next__(self)

    callable = callable

def _copy_array_if_base_present(a):
    """
    Copies the array if its base points to a parent array.
    """
    if a.base is not None:
        return a.copy()
    elif np.issubsctype(a, np.float32):
        return np.array(a, dtype=np.double)
    else:
        return a


def _copy_arrays_if_base_present(T):
    """
    Accepts a tuple of arrays T. Copies the array T[i] if its base array
    points to an actual array. Otherwise, the reference is just copied.
    This is useful if the arrays are being passed to a C function that
    does not do proper striding.
    """
    l = [_copy_array_if_base_present(a) for a in T]
    return l


def _convert_to_bool(X):
    if X.dtype != np.bool:
        X = X.astype(np.bool)
    if not X.flags.contiguous:
        X = X.copy()
    return X


def _convert_to_double(X):
    if X.dtype != np.double:
        X = X.astype(np.double)
    if not X.flags.contiguous:
        X = X.copy()
    return X


def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u

def bhtsne(samples, no_dims=2, perplexity=30.0, theta=0.5, randseed=-1, verbose=False):
    N = samples.shape[0]
    Y = np.zeros((N, no_dims), dtype=np.double)
    _bhtsne_wrap.BHTSNE(_convert_to_double(samples), 
                        N, no_dims, perplexity, theta, randseed, verbose)
    return Y
