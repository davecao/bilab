# -*- coding: utf-8 -*-

import numpy as np

__all__ = [ 'unit', 'degree' ]

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

