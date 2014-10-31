# -*- coding: utf-8 -*-

from functools import  wraps

__all__ = ['AbstractModel','Model']

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

class AbstractModel(object):
    """ Store the support vectors
        Base class

    """
    def __init__(self):
        """ Initialization """
        super(AbstractModel, self).__init__()

    def __str__(self):
        return ""

    def __repr__(self):
        return ""

    def __getitem__(self, item):
        return self.__dict__[item]

    @abstractmethod
    def load(self, model):
        """
            Load precomputed model from a file

        Args:
            model (str): a file name of the svm model

        Returns:
            status (bool): True for successfully load a model from a file
        """

    @abstractmethod
    def classify(self, sample):
        """
            For a given sample, do the classification

        Args:
            sample (array): input sample in np array

        Returns:
            decision (float): decision values for each input sample
        """


class Model(AbstractModel):
    """ A strategy class
    """
    def __init__(self, implementation=None):
        super(Model, self).__init__()
        self.__implementation = implementation
        self.__imp_name = self.__implementation.__class__.__name__

    def __str__(self):
        return self.__implementation.__str__()

    def load(self,model):
        status = self.__implementation.load(model)
        return status

    def classify(self, sample):
        return self.__implementation.classify(sample)

    def __getattr__(self, attr_name):
        try:
            return self.__implementation.__getattribute__(attr_name)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr_name))
