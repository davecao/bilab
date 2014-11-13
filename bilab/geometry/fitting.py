# -*- coding: utf-8 -*-

import numpy as np

__all__ = ['cylinder_fitting']

def global_eval(x, W):
    # P = I -W*W^T projection matrix
    # P: 3x3, W: 3x1
    # W: unit-length direction
    W = np.asmatrix(W)
    P = np.identity(3) - W.T*W

    #pprint.pprint(W)
    # S: 3x3
    S = np.matrix([ [         0, -1 * W[0,2],     W[0,1] ],
                       [     W[0,2],          0, -1 * W[0,0]],
                       [-1 * W[0,1],     W[0,0],           0]])
    # 3x3
    A = np.zeros((3,3))
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
    #print("There are %d atoms" % n)
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
    #print 'value: {}'.format(denominator)
    # sqrlength: 1xn, ave_sqrlength: 1x1, PC:1x3 Y:3xn
    term = sqrlength - ave_sqrlength - 2 * PC.T * Y
    error = np.sum(np.square(term)) / n
    # PC->3xN
    diff = np.repeat(PC, n, axis = 1) - Y
    r_sqr = np.sum(np.sum(np.square(diff), 0)) / n
    return (error, PC, r_sqr)

def cylinder_fitting(points, imax=64, jmax=64, description=None,verbose=False):
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
    
    minError = np.inf
    data_mat = np.asmatrix(points)
    # zero mean
    data_mean = data_mat.mean(0)
    sample = data_mat - data_mean
    w_direct = np.zeros(3)
    c_center = np.zeros(3)
    r_sqr    = 0
    half_pi = np.pi / 2
    two_pi  = 2 * np.pi
    if verbose:
        print "phi  theta  r2 error {}x{}".format(imax,jmax)
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
            error, curr_c, curr_rsqr = global_eval(sample, curr_w)
            #if description is None:
            #    print("err: %s" % error)
            #else:
            #    print("%s err:%s" % (description,error))

            if error < minError:
                minError = error
                w_direct = curr_w
                c_center = curr_c
                r_sqr = curr_rsqr
                if verbose:
                    print "{0:8.3f} {1:8.3f} {2:8.3f} {3:8.3f}".format(phi,
                              theta, r_sqr, error)
    c_data = c_center + data_mean.T
    return (np.array(w_direct).flatten().tolist(),
            np.array(c_data).flatten().tolist(),
            float(r_sqr), float(minError))
