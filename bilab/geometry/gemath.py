# -*- coding: utf-8 -*-
import math
import numpy as np

__all__ = [ 'unit', 'degree', 'mean', 'normal_vector_from_matrix' ]

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

        **vec1** (array) : a n-dimensional vector  

        **vec2** (array) : a n-dimensional vector  

    Returns:
        return an acute angle between two given vectors

    """
    unit_vec1 = unit(vec1)
    unit_vec2 = unit(vec2)

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

def mean(data):
    """ Calculate the mean of a vector/ matrix 
    Args:
        data (array or list): a numpy array or a list

    Returns:
        a numpy array
    """
    d = np.array(data, dtype=np.float64, copy=True)
    if d.ndim < 2:
        return d
    return d.mean(0)

def normal_vector_from_matrix(matrix, eps=1e-8):
    """
    Calculate the normal vector from point matrix through svd decomposition

    Args:
        matrix(numpy array): nx3

    Returns:
        normal (1x3 array): the normal vector defined by points

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
    normal = np.real(V[:, i[0]]).squeeze()

    return normal

def projection_vector(vec_u, vec_v):
    """
        Project the vector 'u' to the vector 'v'

        Args:
           vec_u (array, 1x3): vector U
           vec_v (array, 1x3): vector V
        return a coordinate of the projected point

                        dot(U, V)       V
         Prj(u to v) = ----------- * --------
                         norm(V)      norm(V)
    """
    magnitude_v = np.linalg.norm(vec_v)
    norm_v = vec_v / magnitude_v
    dot_uv = np.dot(vec_u, vec_v)
    
    return (dot_uv / magnitude_v * norm_v).tolist()
