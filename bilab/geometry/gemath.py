# -*- coding: utf-8 -*-
import math
import numpy as np

__all__ = [ 'unit', 
            'degree', 
            'mean', 
            'normal_vector_from_matrix', 
            'projection_vector' ]

def unit(vec):
    """ unitify a given vector 

    .. ipython:: python

        import bilab
        vec = np.array([1, 2, 3])
        bilab.geometry.unit(vec)

    Args:
        **vec** (array): a n-dimensional vector

    Returns:
        an array of normalized n-dimensional vector

    """
    return np.asarray(vec) / np.linalg.norm(vec)

def degree(vec1, vec2):
    """
    Calculate angle between two vectors

    .. ipython:: python

        vec1 = np.array([1, 2, 3])
        vec2 = np.array([4, 5, 6])
        bilab.geometry.degree(vec1, vec2)

    Args:

        **vec1** (array or list) : a n-dimensional vector  

        **vec2** (array or list) : a n-dimensional vector  

    Returns:
        return an acute angle between two given vectors

    """
    v1 = np.array(vec1, dtype=np.float64, copy=True).squeeze()
    v2 = np.array(vec2, dtype=np.float64, copy=True).squeeze()

    if v1.shape != v2.shape:
        raise ValueError("Two input vectors are not in the same shape")

    unit_vec1 = unit(v1)
    unit_vec2 = unit(v2)

    # The first version:  following is not work as the two parallel vectors
    #angle_rad = numpy.arccos(numpy.dot(vec1, vec2) /
    #    (numpy.linalg.norm(vec1) * numpy.linalg.norm(vec2)))
    #angle_deg = angle_rad * 180 / numpy.pi

    # The revised version:
    #angle_rad = numpy.arccos(numpy.dot(unit_vec1, unit_vec2))
    #angle_deg = angle_rad * 180 / numpy.pi
    #if numpy.isnan(angle):
    #    if (unit_vec1 == unit_vec2).all():
    #        return 0.0
    #    else
    #        return 180.0
    # The second revised with clip
    angle_rad = np.arccos(np.clip(np.dot(unit_vec1, unit_vec2), -1, 1))
    angle_deg = angle_rad * 180 / np.pi
    return  180 - angle_deg if angle_deg > 90 else angle_deg

def mean(data, raw=True):
    """ Calculate the mean of a vector/ matrix

    .. Note:
        Return an raw array if the dimension is one

    .. ipython:: python
        vec1 = np.array([[1, 2, 3],[4, 5, 6]])
        bilab.geometry.mean(vec1)

    Args:
        data (array or list): a numpy array or a list

    Returns:
        a numpy array

    """
    d = np.array(data, dtype=np.float64, copy=True)
    if d.ndim == 1 and raw:
        return d
    elif d.ndim == 0:
        return d.mean()
    return d.mean(0)

def normal_vector_from_matrix(matrix, eps=1e-8):
    """
    Calculate the normal vector from point matrix through svd decomposition

    Args:
        matrix(numpy array): nx3

    Returns:
        normal (1x3 list): the normal vector defined by points

    """
    if not isinstance(matrix, (np.ndarray, np.generic)) or matrix.ndim != 2:
        raise ValueError("{}:Inappropriate argument value, it should be a 2d numpy array".format('normal_vector_from_matrix'))

    M = np.matrix(matrix, dtype=np.float64, copy=True)

    # Centralized by column
    centeroid = np.mean(M, axis=0)
    M -= centeroid

    # covariance matrix (3xn x nx3)= 3x3
    cov_matrix = M.T * M

    # SVD decomposition
    U, s, V = np.linalg.svd(cov_matrix)
    i = np.where(abs(np.real(s)) < eps)[0]

    # normal: unit eigenvector corresponding to lowest eigenvalue
    normal = np.squeeze(np.asarray(np.real(V[:, i[0]])))

    return normal

def projection_vector(vec_u, vec_v):
    """
    Project the vector 'u' to the vector 'v'
                        dot(U, V)       V
         Prj(u to v) = ----------- * --------
                         norm(V)      norm(V)
    .. math::
        Prj(u to v) = \\frac{dot(u,v)}{norm(v)}\\times\\frac{v}{norm(v)}
    
    Args:
        vec_u (array or list, 1x3): vector U
        vec_v (array or list, 1x3): vector V

    Returns
        a coordinate of the projected point

    """
    u = np.array(vec_u, dtype=np.float64, copy=True)
    v = np.array(vec_v, dtype=np.float64, copy=True)

    magnitude_v = np.linalg.norm(v)
    norm_v = v / magnitude_v
    dot_uv = np.dot(u, v)
    
    return (dot_uv / magnitude_v * norm_v).tolist()
