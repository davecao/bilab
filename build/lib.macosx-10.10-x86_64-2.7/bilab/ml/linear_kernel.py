# -*- coding: utf-8 -*-

__all__ = []

import sys
from bilab.ml import Kernel

__all__ = ['Linear_kernel']

class Linear_kernel(Kernel):
    """ docstring for class

    """

    def __init__(self):
        """Initialization.

        Args:
           param (type):

        Kwargs:

        """
        super(Linear_kernel, self).__init__()

    def __str__(self):
        """ Serialize """
        return "Kernel.name: {}".format(self.__class__.__name__)

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        """kernel function

            If the input two vector of same size it returns the dot product of
            them.
            If a matrix is given, it will return a kernel matrix.

        Args:
            u (numpy array): mxn matrix or a n-dim vector,
                num of samples in rows and features in columns
            v (numpy array, optional): mxn matrix or a n-dim vector
                num of samples in rows and features in columns

        Kwargs:

        Returns:
            a float value or a matrix.

        Raises:

        """
        argc = len(args)
        K = np.inf

        if argc == 2:
            # Two variables are given
            u = args[0]
            v = args[1]
            if isinstance(u, (np.ndarray, np.generic)) and \
                isinstance(v, (np.ndarray, np.generic)):
                #check shapes of u,v
                if u.shape == v.shape:
                    K = np.dot(u, v)
                else:
                    raise ValueError("The sizes of the two matrice "
                        "are not same")
            else:
                raise TypeError("Inappropriate argument type for {}"
                    .format(self.__class__.__name__))
        elif argc == 1:
            u = args[0]
            if isinstance(u, (np.ndarray, np.generic)):
                K = u * u.T
            else:
                raise TypeError("Inappropriate argument type for {}"
                    .format(self.__class__.__name__))
        else:
            print("Incorrect number of arguments")
            sys.exit(1)
        return K
