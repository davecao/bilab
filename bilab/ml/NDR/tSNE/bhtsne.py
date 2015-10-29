from __future__ import division, print_function, absolute_import

import warnings
import numpy as np
import ctypes
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


def _convert_to_double_p(X):
    X = _convert_to_double(X)
    return X.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u


def _mean_and_std(X, axis=0, with_mean=True, with_std=True):
    """Compute mean and std deviation for centering, scaling.
    Zero valued std components are reset to 1.0 to avoid NaNs when scaling.
    """
    X = np.asarray(X, dtype=float)
    Xr = np.rollaxis(X, axis)

    if with_mean:
        mean_ = Xr.mean(axis=0)
    else:
        mean_ = None

    if with_std:
        std_ = Xr.std(axis=0)
        std_ = _handle_zeros_in_scale(std_)
    else:
        std_ = None

    return mean_, std_


def _handle_zeros_in_scale(scale):
    ''' Makes sure that whenever scale is zero, we handle it correctly.
    This happens in most scalers when we have constant features.'''

    # if we are fitting on 1D arrays, scale might be a scalar
    if np.isscalar(scale):
        if scale == 0:
            scale = 1.
    elif isinstance(scale, np.ndarray):
        scale[scale == 0.0] = 1.0
        scale[~np.isfinite(scale)] = 1.0
    return scale


def _scale(X, axis=0, with_mean=True, with_std=True, copy=True):
    X = np.asarray(X)
    mean_, std_ = _mean_and_std(X, axis, with_mean=with_mean,
                                with_std=with_std)
    if copy:
        X = X.copy()
    # Xr is a view on the original array that enables easy use of
    # broadcasting on the axis in which we are interested in
    Xr = np.rollaxis(X, axis)
    if with_mean:
        Xr -= mean_
        mean_1 = Xr.mean(axis=0)
        # Verify that mean_1 is 'close to zero'. If X contains very
        # large values, mean_1 can also be very large, due to a lack of
        # precision of mean_. In this case, a pre-scaling of the
        # concerned feature is efficient, for instance by its mean or
        # maximum.
        if not np.allclose(mean_1, 0):
            warnings.warn("Numerical issues were encountered "
                          "when centering the data "
                          "and might not be solved. Dataset may "
                          "contain too large values. You may need "
                          "to prescale your features.")
            Xr -= mean_1
    if with_std:
        Xr /= std_
        if with_mean:
            mean_2 = Xr.mean(axis=0)
            # If mean_2 is not 'close to zero', it comes from the fact that
            # std_ is very small so that mean_2 = mean_1/std_ > 0, even if
            # mean_1 was close to zero. The problem is thus essentially due
            # to the lack of precision of mean_. A solution is then to
            # substract the mean again:
            if not np.allclose(mean_2, 0):
                warnings.warn("Numerical issues were encountered "
                              "when scaling the data "
                              "and might not be solved. The standard "
                              "deviation of the data is probably "
                              "very close to 0. ")
                Xr -= mean_2
    return X


def _bhtsne_call_wrap(samples, no_dims=2, perplexity=30.0, theta=0.5,
                      randseed=-1, sampling=100, scale=False, verbose=False):
    N = samples.shape[0]
    D = samples.shape[1]

    # Y = np.zeros((N, no_dims), dtype=np.double)
    # tsne = _bhtsne_wrap.bhtsne(samples, Y, N, D, no_dims, perplexity,
    #                           theta, randseed, scale, verbose)
    tsne = _bhtsne_wrap.bhtsne(samples, N, D, no_dims, perplexity,
                               theta, randseed, scale, verbose)
    # print(tsne)
    Y = tsne.run()
    Y_scale = _scale(Y)
    Y = np.cosh(Y_scale)
    return Y


def bhtsne(samples, no_dims=1, perplexity=30.0, theta=0.5, randseed=-1,
           sampling=100, scale=False, verbose=False):
    Y = _bhtsne_call_wrap(samples, no_dims=no_dims, perplexity=perplexity,
                          theta=theta, randseed=randseed,
                          sampling=sampling, verbose=False)

    for i in xrange(sampling):
        if verbose:
            sys.stdout.write("Sampling at iter {}/{}\r".format(i+1, sampling))
            sys.stdout.flush()

        res = _bhtsne_call_wrap(samples, no_dims=no_dims,
                                perplexity=perplexity,
                                theta=theta, randseed=randseed,
                                sampling=sampling, verbose=False)
        Y = np.column_stack((Y, res))
    # find the maximum and minimum in rows
    col_min = np.amin(Y, axis=1)
    col_max = np.amax(Y, axis=1)
    # get their summation
    sum_col_del = col_max + col_min
    # Remove the max and min
    sum_cols = np.sum(Y, axis=1) - sum_col_del

    # Compute the average on row
    NumofSampplingDimensions = Y.shape[1] - 2
    final_mapping = sum_cols/NumofSampplingDimensions

    return final_mapping
