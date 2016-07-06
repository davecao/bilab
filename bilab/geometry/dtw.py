# -*- coding: utf-8 -*-

"""
   Implementation of Dynamic Time Wraping
"""
import numpy as np
#i mport bilab


def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u


def _traceback(D):
    """
        traceback for DTW matrix
    """
    D = np.asarray(D, dtype=None, order='c').squeeze()
    D = np.atleast_2d(D)
    if D.ndim > 2:
        raise ValueError("Input vector should be 2-D.")
    wrap_path = [(0, 0)]
    i, j = np.array(D.shape) - 1
    while (i > 0 and j > 0):
        ind = np.argmin((D[i-1, j-1], D[i-1, j], D[i, j-1]))
        if ind == 0:
            i = i - 1
            j = j - 1
        elif ind == 1:
            i = i - 1
        elif ind == 2:
            j = j - 1
        wrap_path.append((i, j))
    return wrap_path


def stdDTW(u, v, dist=lambda x, y: abs(x-y)):
    """

    Standard DTW algorithm with the complexity O(N^2)

    Args:
        u (np array):
        v (np array):

    Return:
        distance and warp path between u and v
    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    M = len(u) + 1
    N = len(v) + 1
    s = np.zeros([M, N])
    s[:] = np.inf
    s[0, 0] = 0
    wrap_path = []
    for i in range(1, M):
        for j in range(1, N):
            cost = abs(u[i-1] - v[j-1])
            s[i, j] = cost + min(s[i-1, j],    # insertion
                                 s[i, j-1],    # deletion
                                 s[i-1, j-1])  # match
    wrap_path = _traceback(s)  # tuple
    return (s[M-1, N-1]/sum(s.shape), wrap_path)


def ExpandedResWindow(lowResPath, shrinkU, shrinkV, X, Y, radius):
    pass


def constrainedDTW(u, v, window):
    pass


def fastDTW(u, v, radius=1):
    """
        Implementation of Fast DTW algorithm

    Reference:

    Stan Salvador & Philip Chan, FastDTW: Toward Accurate Dynamic Time Warping
    in Linear Time and Space. KDD Workshop on Mining Temporal and
    Sequential Data, pp. 70-80, 2004

    Args:
        u (np array):
        v (np array):

    Return:
        similarity distance and wrap path between u and v
    """

    u = _validate_vector(u)
    v = _validate_vector(v)
    M = len(u)
    N = len(v)
    minSize = radius + 2
    resolution_factor = 2  # for shrinking data
    if M < minSize or N < minSize:
        # call Standard DTW
        return stdDTW(u, v)
    else:
        # Recursive case
        # Project the warp path from a coarser resolution onto
        # the current current resolution. Run DTW only along the
        # projected path (and also ‘radius’ cells from the projected path).
        shrinkU = u[:M/resolution_factor]
        shrinkV = v[:N/resolution_factor]
        lowResPath = fastDTW(shrinkU, shrinkV, radius)
        window = ExpandedResWindow(
                        lowResPath, shrinkU, shrinkV, u, v, radius)

        return constrainedDTW(u, v, window)
