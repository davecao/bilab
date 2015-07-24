# -*- coding: utf-8 -*-

import numpy as np


from bilab.tests import unittest
from bilab.ml import Kernel, RBF_kernel

def test_Kernel():
    u = np.asmatrix([np.zeros(10), np.ones(10) * 10])
    ans = np.matrix('1 0;0 0')
    kernel = Kernel(RBF_kernel)
    x = kernel(u)
    if u == x:
        print("Test Kernel: RBF strategy is OK")
    else:
        print("Test Kernel: RBF strategy is Failed")

if __name__ == '__main__':
    unittest.main(defaultTest="test_Kernel")
