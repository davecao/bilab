# -*- coding: utf-8 -*-

import numpy as np
import sys

__all__ = ['SVDaffineTransformation',
           'UmeyamaTransformation',
           'ransac']


def SVDaffineTransformation(src, dst, ax=0):
    """ Kabasch algorithm """
    assert src.shape == dst.shape
    # get mean
    src_mean = src.mean(ax)
    dst_mean = dst.mean(ax)
    # centerize
    src_demean = src - src_mean
    dst_demean = dst - dst_mean

    matrix = np.dot(src_demean.T, dst_demean)

    U, s, Vh = np.linalg.svd(matrix)
    Id = np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, np.sign(np.linalg.det(matrix))]])
    rotation = np.dot(Vh.T, np.dot(Id, U.T))
    translation = dst_mean - np.dot(src_mean, rotation.T)
    return rotation, translation


def UmeyamaTransformation(src, dst, with_scaling=True):
    # n x d: n samples with d-dimensions
    # Reference: Eigen/Core/src/Umeyama.h
    if src.shape != dst.shape:
        print("The size of matrice of src and dst should be same.")
        sys.exit(1)
    ax = 0
    dx = 1
    dim = src.shape[dx]
    n = src.shape[ax]
    one_over_n = 1.0 / n
    # get mean
    src_mean = src.mean(ax)
    dst_mean = dst.mean(ax)
    # centerize
    src_centerized = src - src_mean
    dst_centerized = dst - dst_mean

    # variance Equation (36)-(37)
    src_variance = np.sum(src_centerized**2, axis=1)
    src_variance = src_variance.mean(ax)

    # matrix Equation (38)
    sigma = one_over_n * np.dot(dst_centerized.T, src_centerized)

    # SVD
    U, d, Vt = np.linalg.svd(sigma)
    # Prepare result: dim is 2,3 generally
    Qmatrix = np.eye(dim+1)

    # reorgnize the diag S in Equation (39)
    S = np.ones(dim)
    det_sigma = np.linalg.det(sigma)
    if det_sigma < 0:
        S[dim-1] = -1
    # Determine the optimum transfromation parameters
    # when rank(simga) >= d-1
    threshold = 1e-12 * d[0]  # the maximum of singular value
    r = sum([1 for i in d if i > threshold])

    det_UV = np.linalg.det(U) * np.linalg.det(Vt.T)

    if r == dim-1:  # Condition in Equation (43)
        if det_UV > 0:
            Qmatrix[0:dim, 0:dim] = np.dot(U, Vt)
        else:
            s = S[dim-1]
            S[dim-1] = -1
            Qmatrix[0:dim, 0:dim] = np.dot(np.dot(U, np.diag(S)), Vt)
            S[dim-1] = s
    else:
        # Optimum transfromation parameters
        # Equation (40), (41)
        Qmatrix[0:dim, 0:dim] = np.dot(np.dot(U, np.diag(S)), Vt)

    if with_scaling:
        # Equation (42): scaling
        c = 1/src_variance * np.dot(d, S)
        # Equation (41): translation
        Qmatrix[0:dim, dim] = dst_mean
        Qmatrix[0:dim, dim] -= c*np.dot(Qmatrix[0:dim, 0:dim], src_mean)
        Qmatrix[0:dim, 0:dim] *= c
    else:
        # last column in Q
        Qmatrix[0:dim, dim] = dst_mean
        Qmatrix[0:dim, dim] -= np.dot(Qmatrix[0:dim, 0:dim], src_mean)

    return Qmatrix


def random_partition(n, n_data):
    """return n random rows of data (and also the other len(data)-n rows)"""
    all_idxs = np.arange(n_data)
    np.random.shuffle(all_idxs)
    idxs1 = all_idxs[:n]
    idxs2 = all_idxs[n:]
    return idxs1, idxs2


def ransac(src, dst, n, num_iter, threshold, min_inliers):
    """
    Random sample consensus (RANSAC) is an iterative method to estimate
    parameters of a mathematical model from a set of observed data which
    contains outliers.
    See https://en.wikipedia.org/wiki/RANSAC

    Args:
        n - the minimum number of data values required to fit the model
        num_iter - the number of iterations
        max_iter - maximum number of iterations
        min_iter - minimum number of iterations
        threshold:
    """

    max_inliers = 0
    inliers = None

    for i in range(num_iter):
        maybe_idxs, test_idxs = random_partition(n, src.shape[0])
        # Compute transformation using these n points
        Q = UmeyamaTransformation(src[maybe_idxs, :], dst[maybe_idxs, :])
        # Apply transformation to src
        src_trans = np.dot(Q, src)
        # find the inliers
        dist = np.sum((dst - src_trans)**2, axis=0)
        num_inliers = np.count_nonzero(dist < threshold)
        if max_inliers < num_inliers:
            max_inliers = num_inliers
            inliers = np.nonzero(dist < threshold)[0]
    if max_inliers < min_inliers:
        return None, []
    # reestimation
    Q = UmeyamaTransformation(src[maybe_idxs, inliers],
                              dst[maybe_idxs, inliers])
    print("Finished ransac with {} inliers".format(max_inliers))
    return Q, inliers
