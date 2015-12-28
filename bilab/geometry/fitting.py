# -*- coding: utf-8 -*-

import numpy as np
import logging

from bilab.concurrent import ThreadPool, ThreadTask, ThreadWorker,\
                             CallBackResults
logging.basicConfig(
    level=logging.DEBUG,
    format='[%(levelname)s] (%(threadName)-10s) %(funcName)s %(message)s')


def global_eval(x, W):
    # P = I -W*W^T projection matrix
    # P: 3x3, W: 3x1
    # W: unit-length direction
    W = np.asmatrix(W)
    P = np.identity(3) - W.T*W

    # pprint.pprint(W)
    # S: 3x3
    S = np.matrix([[0, -1 * W[0, 2], W[0, 1]],
                   [W[0, 2], 0, -1 * W[0, 0]],
                   [-1*W[0, 1], W[0, 0], 0]])
    # 3x3
    A = np.zeros((3, 3))
    # 3x1
    B = np.zeros(3).T
    # number of samples: N X 3
    n = x.shape[0]
    # Y: 3xN = P 3x3 x.T 3xN
    x = np.asmatrix(x)
    if x.all() == 0:
        logging.debug("x is zero")

    Y = P * x.T

    # 1xN
    sqrlength = np.sum(np.square(Y), 0)
    # print("There are %d atoms" % n)
    for i in range(n):
        A = A + np.outer(Y[:, i], Y[:, i].T)
        B = B + sqrlength[:, i] * Y[:, i].T
    # get mean value
    ave_sqrlength = np.sum(sqrlength) / n
    # 3x3, 3x3, 3x3 = 3x3
    Ahat = -1 * S * A * S

    # 3x3, 3x1 = 3x1
    denominator = np.trace(Ahat * A)

    PC = (Ahat * B.T) / denominator
    # print 'value: {}'.format(denominator)
    # sqrlength: 1xn, ave_sqrlength: 1x1, PC:1x3 Y:3xn
    term = sqrlength - ave_sqrlength - 2 * PC.T * Y
    error = np.sum(np.square(term)) / n
    # PC->3xN
    diff = np.repeat(PC, n, axis=1) - Y
    r_sqr = np.sum(np.sum(np.square(diff), 0)) / n
    return (error, PC, r_sqr, W)


def fit_multi(data, imax, jmax, num_threads=12,
              CheckInterval=0.5, verbose=False):
    """
        Multiple threaded fitting
    """
    minError = np.inf
    w_direct = np.zeros(3)
    c_center = np.zeros(3)
    r_sqr = 0
    half_pi = np.pi / 2
    two_pi = 2 * np.pi

    cbk = CallBackResults(10)
    cbk.curr_w = w_direct
    cbk.minError = minError
    cbk.w_direct = w_direct
    cbk.c_center = c_center
    cbk.r_sqr = r_sqr

    def th_callback(data, res=cbk):
        try:
            res.lock.acquire()
            error, curr_c, curr_rsqr, w_direct = data
            if error < res.minError:
                res.minError = error
                res.w_direct = w_direct
                res.c_center = curr_c
                res.r_sqr = curr_rsqr
            c_d = np.asarray(res.curr_w).flatten().tolist()
            w_d = np.asarray(res.w_direct).flatten().tolist()
        finally:
            res.lock.release()

    thread_pool = ThreadPool(num_threads,
                             CheckInterval=CheckInterval,
                             verbose=verbose)

    for j in range(jmax):
        phi = half_pi * j / jmax
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        for i in range(imax):
            theta = two_pi * i / imax
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)
            curr_w = np.matrix([cos_theta * sin_phi,
                                sin_theta * sin_phi,
                                cos_phi])

            taskid = "task_{}_{}".format(j, i)
            cbk.curr_w = curr_w
            thread_pool.add_task(taskid, global_eval, data, curr_w,
                                 callback=th_callback)
    thread_pool.wait_completion()
    thread_pool.shutdown()

    return (cbk.w_direct,
            cbk.c_center,
            float(cbk.r_sqr),
            float(cbk.minError))


def fit_single(data, imax, jmax, verbose=False):
    """
        Single threaded fitting
    """
    minError = np.inf
    w_direct = np.zeros(3)
    c_center = np.zeros(3)
    r_sqr = 0
    half_pi = np.pi / 2
    two_pi = 2 * np.pi

    for j in range(jmax):
        phi = half_pi * j / jmax
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        for i in range(imax):
            theta = two_pi * i / imax
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)
            curr_w = np.matrix([cos_theta*sin_phi,
                                sin_theta*sin_phi,
                                cos_phi])
            error, curr_c, curr_rsqr, w = global_eval(data, curr_w)

            if error < minError:
                minError = error
                w_direct = curr_w
                c_center = curr_c
                r_sqr = curr_rsqr
            if verbose:
                c_d = np.asarray(curr_w).flatten().tolist()
                w_d = np.asarray(w_direct).flatten().tolist()

    return (w_direct, c_center, float(r_sqr), float(minError))


def cylinder_fitting(points, imax=64, jmax=64,
                     description=None,
                     verbose=False,
                     num_threads=1,
                     CheckInterval=0.5):
    """
    Fitting a cylinder to a set of points:

    The cylinder radius is r, point C lies at the axis and has
    unit-length direction W.

    .. note::
        x = X - mean(X) --> sample mean is zero

    Error function:

    .. math::
        E(r^{2}, C, W) = min \\sum_{i=1}^{n}(r_{i}^{2} - r^{2})^{2}

    contains 6 parameters:
    **one** for the squared radius :math:`r^{2}`,
    :math:`r^{2} = \\frac{1}{n}\\sum_{i=1}^{n}r_{i}^{2}`,
    **three** for the point C and **two** for the unit-length direction W

    Args:
        points (array) : nx3 matrix

        imax (integer) : grid value for x axis

        jmax (integer) : grid value for y axis

        description (str) : option value

        verbose (bool) : show info in detail

    Returns:
        r_sqr (float) : the square of r

        C (array)  : center 3x1

        W (array) : direction 3x1

    """
    if not type(num_threads) == int:
        raise ValueError("Error: option num_threads should be an integer")

    minError = np.inf
    data_mat = np.asmatrix(points)
    # zero mean
    data_mean = data_mat.mean(0)
    sample = data_mat - data_mean
    w_direct = np.zeros(3)
    c_center = np.zeros(3)
    r_sqr = 0
    half_pi = np.pi / 2
    two_pi = 2 * np.pi

#    if verbose:
#        print "phi  theta  r2 error {}x{}".format(imax,jmax)

    if num_threads == 1:
        w_direct, c_center, r_sqr, minError = \
            fit_single(sample, imax, jmax, verbose=verbose)
    else:
        w_direct, c_center, r_sqr, minError = fit_multi(
            sample, imax, jmax,
            num_threads=num_threads,
            CheckInterval=CheckInterval,
            verbose=verbose)

    c_data = c_center + data_mean.T

    return (np.array(w_direct).flatten().tolist(),
            np.array(c_data).flatten().tolist(),
            float(r_sqr), float(minError))
