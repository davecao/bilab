# -*- coding: utf-8 -*-
# from __future__ import absolute_import

from functools import wraps
# from abc import ABCMeta, abstractmethod
from bilab import PY3K
from bilab import forcefieldsList
from bilab.structure.forcefields import AmberParamParser
# from bilab.io.FortranFormat import (FortranFormat, FortranLine)

__all__ = ['ForceField', 'AmberForceField']


def abstractmethod(method):
    # decorator borrowed from Mozilla mxr
    if PY3K:
        line = method.__code__.co_firstlineno
        filename = method.__code__.co_filename
    else:
        line = method.func_code.co_firstlineno
        filename = method.func_code.co_filename

    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError(
            'Abstract method %s at File "%s", line %s'
            ' should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented


class Register(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(Register, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, '_registry'):
            cls._registry = set()
        cls._registry.add(cls)
        # Remove base classes
        cls._registry -= set(bases)
        # print("{} in {} for {}".format("FFRegister",
        #    "__init__", cls.__name__))

    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls._registry)

    def __str__(cls):
        if cls in cls._registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])


class ForceField(object):
    """
    interface of Force Field
        Base class

        Note: multiple inherited subclasses may cause metaclass conflict
    """
    __metaclass__ = Register
    # __metaclass__ = ABCMeta

    def __init__(self, *args, **kwargs):
        """ Initialization """
        super(ForceField, self).__init__()

#    @abstractmethod
#    def __call__(cls, clsname, *args, **kwargs):
#        """ ForceField callable
#            Generate an object of specified subclass name
#        """
        # return cls.__metaclass__._registry(clsname)(*args, **kwargs)
#        print("{} in {}".format("ForceField", "__call__"))
        # print("{}".format(self))
        # print("{}".format(self.__class__))
        # return self.__class__.__Factory(clsname)(*args, **kwargs)
#        try:
#            for c in cls.__registry:
#                if c.__name__ == clsname:
#                    return c(*args, **kwargs)
#                #cls.__metaclass__.registry[name]
#        except ValueError:
#            raise ValueError( "{} is not registered yet or wrong name".format(
#                               clsname)#)

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass

    @classmethod
    def Factory(cls, name):
        """
        usage:
            ff = bilab.structure.forcefields.ForceField()
            # create an object of class AmberForceField
            amberFF = ff.Factory("AmberForceField")(*args, **kwargs)
        """
        # print("classmethod Factory is called by {}".format(cls.__name__))
        try:
            for c in cls._registry:
                if c.__name__ == name:
                    return c()
                # cls.__metaclass__.registry[name]
        except ValueError:
            raise ValueError(
                "{} is not registered yet or wrong name".format(name))


class AmberForceField(ForceField):
    """
        Amber Force Field
    """
    def __init__(self, *args, **kwargs):
        try:
            self.ff_name = kwargs.get("ff_name", "amber94")
            self.data_parser = kwargs.get("parser", AmberParamParser())
            self.ff_dat_loc = forcefieldsList[self.ff_name]
        except KeyError:
            raise KeyError(
                "Wrong name of the force field: {}".format(self.ff_name))
        # ff_params is a dict
        #
        # print("{}".format(type(self.data_parser)))
        self.ff_params = self.data_parser.parse(self.ff_dat_loc)

    def __str__(self):
        return "{} from {}".format(self.ff_name, self.ff_dat_loc)

    def __repr__(self):
        return self.ff_name




#ForceField.register(AmberForceField)
