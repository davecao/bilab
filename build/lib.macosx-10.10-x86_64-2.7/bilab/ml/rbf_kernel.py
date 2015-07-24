# -*- coding: utf-8 -*-

import sys
import numpy as np

from bilab.ml import Kernel

__all__ = ['RBF_kernel']

class RBF_kernel(Kernel):
    """ docstring for class

    """

    def __init__(self, gamma=None):
        """Initialization.

        Args:
           gamma (float):

        Kwargs:

        """
        super(RBF_kernel, self).__init__()
        self.gamma = gamma

    def __str__(self):
        """ Serialize """
        return "Kernel.name:{} \nKernel.gamma:{}\n".format(self.__class__.__name__, self.gamma)

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        """kernel function

            If the input is two vectors of same size, it returns a float number
                exp(-gamma* dot((u-v), (u-v).T)).

            If the input is only a matrix, it returns a kernel matrix
                exp(-gamma*|u - v|^2)
        Args:
            u (numpy array): mxn matrix or a n-dim vector,
                num of samples in rows and features in columns
            v (numpy array, optional): mxn matrix or a n-dim vector
                num of samples in rows and features in columns

        Kwargs:
            gamma (float): the parameter for rbf basis function.
                exp(-gamma*|sample1 - sample2|^2). If it is not set,
                the reciprocal of the number of features will be used as
                gamma.

        Returns:
            a float value or a matrix.

        Raises:

        """
        argc = len(args)
        K = np.inf
        gamma = kwargs.get('gamma')
        if argc == 2:
            # Two variables are given
            u = args[0]
            v = args[1]

            if isinstance(u, (np.ndarray, np.generic)) and \
                isinstance(v, (np.ndarray, np.generic)):

                # count num of feaures in column
                num_of_features = u.size if u.ndim == 1 else u.shape[1]

                #check gamma
                if self.gamma is None and gamma is None:
                    # set gamma to 1/num_of_features
                    gamma = 1 / num_of_features
                elif gamma is None and self.gamma is not None:
                    gamma = self.gamma

                # check whether the sizes of the twos are same or not
                if u.shape == v.shape:
                    diff_uv = u - v
                    coefficient = -1 * gamma
                    K = np.exp(coefficient * np.dot(k, k.T))
                else:
                    raise ValueError("The sizes of the two matrice "
                        "are not same")
            else:
                raise TypeError("Inappropriate argument type for {}"
                    .format(self.__class__.__name__))
        elif argc == 1:
            # the input must be a matrix
            # see kernel-methods.net
            #   K
            u = args[0]
            if isinstance(u, (np.ndarray, np.generic)):
                # count num of feaures in column
                num_of_features = u.size if u.ndim == 1 else u.shape[1]
                num_of_samples = u.shape[0] if u.ndim == 2 else u.ndim
                # check gamma
                if self.gamma is None and gamma is None:
                    # set gamma to 1/num_of_features
                    gamma = 1.0 / num_of_features
                elif gamma is None and self.gamma is not None:
                    gamma = self.gamma
                #print("num_of_samples:{}x{}".format(u.shape[0], u.shape[1]))
                K = u * u.T
                # must convert diag to matrix
                d = np.asmatrix(np.diag(K)).T
                K = K - np.ones((num_of_samples, 1)) * d.T
                K = K - d * np.ones((1, num_of_samples))
                K = np.exp(gamma * K)
            else:
                raise TypeError("Inappropriate argument type for {}"
                    .format(self.__class__.__name__))
        #else:
        #    print("RBF kernel:Incorrect number of arguments")
        #    sys.exit(1)
        # save parameter
        self.gamma = gamma
        return K

    def _print_np_matrix(self, d):
        """ Pretty print numpy matrix """
        np.set_printoptions(precision=4, suppress=True)
        print("{}".format(d))
        # reset to default
        np.set_printoptions(precision=8, suppress=False)
