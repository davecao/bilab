# -*- coding: utf-8 -*-

import sys
from bilab.ml import Kernel

__all__ = ['Polynomial_kernel']

class Polynomial_kernel(Kernel):
    """ docstring for class

    """

    def __init__(self, gamma=None, coef=1, degree=3):
        """Initialization.

        Args:

        Kwargs:

        """
        super(Polynomial_kernel, self).__init__()
        self.gamma = gamma
        self.coef = coef
        self.degree = degree

    def __str__(self):
        """ Serialize """
        des  = "Kernel.name:{}\n".format(self.__class__.__name__)
        des += "Kernel.gamma:{}\n".format(self.gamma)
        des += "Kernel.coef:{}\n".format(self.coef)
        des += "Kernel.degree:{}".format(self.degree)
        return des

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        """kernel function

            If the input two vector of same size it returns the dot product of
            them.
            If a matrix is given, it will return a kernel matrix.
            (gamma*u'*v + coef)^degree
        Args:
            u (numpy array): mxn matrix or a n-dim vector,
                num of samples in rows and features in columns
            v (numpy array, optional): mxn matrix or a n-dim vector
                num of samples in rows and features in columns

        Kwargs:
            gamma (float): the reciprocal of the number of features of a sample.
            coef  (float): default is 1
            degree (int) : default is 3

        Returns:
            a float value or a matrix.

        Raises:

        """
        argc = len(args)
        K = np.inf
        gamma = kwargs.pop('gamma', None)
        coef = kwargs.pop('coef', 1)
        degree = kwargs.pop('degree', 3)

        if argc == 2:
            # Two variables are given
            u = args[0]
            v = args[1]
            if isinstance(u, (np.ndarray, np.generic)) and \
                isinstance(v, (np.ndarray, np.generic)) and \
                u.shape[0] == 1 and v.shape[0] == 1:
                #check shapes of u,v
                if u.shape == v.shape:
                    num_of_features = u.size if u.ndim == 1 else u.shape[1]
                    if gamma is None:
                        gamma = 1.0 / num_of_features

                    K = np.power(gamma * np.dot(u, v) + coef, degree)
                else:
                    raise ValueError("The sizes of the two matrice "
                        "are not same")
                # save parameters actually used
                self.gamma = gamma
                self.coef = coef
                self.degree = degree
            else:
                raise TypeError("Inappropriate argument type for {}"
                    .format(self.__class__.__name__))
        elif argc == 1:
            u = args[0]
            if isinstance(u, (np.ndarray, np.generic)):
                num_of_features = u.size if u.ndim == 1 else u.shape[1]
                if gamma is None:
                    gamma = 1.0 / num_of_features

                K = np.power(gamma * u * u.T + coef, degree)
                self.gamma = gamma
                self.coef = coef
                self.degree = degree
            else:
                raise TypeError("Inappropriate argument type for {}"
                    .format(self.__class__.__name__))
        else:
            print("Incorrect number of arguments")
            sys.exit(1)
        return K
