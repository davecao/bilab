# -*- coding: utf-8 -*-

from functools import  wraps

from Scientific.IO.FortranFormat import FortranFormat, FortranLine
from Scientific.IO.TextFile import TextFile
from Scientific.DictWithDefault import DictWithDefault

# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    line = method.func_code.co_firstlineno
    filename = method.func_code.co_filename
    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError('Abstract method %s at File "%s", line %s'
            'should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented


class ForceFields(object):

    """ Store the support vectors
        Base class

    """
    def __init__(self):
        """ Initialization """
        super(ForceFields, self).__init__()

    def __str__(self):
        return ""

    def __repr__(self):
        return ""

    def __getitem__(self, item):
        return self.__dict__[item]

    @abstractmethod
    def load(self, params):
        """
            Load force field parameters from a file

        Args:
            params (str): a file name of a force field

        Returns:
            status (bool): True for successfully load a parameters from a file
        """

