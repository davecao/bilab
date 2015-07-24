# -*- coding: utf-8 -*-

from bilab import LOGGER

try:
    from bilab.linalg._CWMatrix import CWMatrix as CWMatrix
except ImportError:
    raise ImportError('CWMatrix module could not be imported. '
                      'Maybe _CWMatrix.so no found.')

__all__ = ['CWMatrix']
