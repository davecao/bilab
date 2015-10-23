# -*- coding: utf-8 -*-
#from __future__ import absolute_import

from functools import wraps
#from abc import ABCMeta, abstractmethod
from bilab.io.FortranFormat import FortranFormat, FortranLine
from bilab import PY3K
from bilab import forcefieldsList

__all__ = ['FFParserInterface', 'AmberParamParser']

# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    if PY3K:
        line = method.__code__.co_firstlineno
        filename = method.__code__.co_filename
    else:
        line = method.func_code.co_firstlineno
        filename = method.func_code.co_filename
    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError('Abstract method %s at File "%s", line %s'
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
        cls._registry -= set(bases) #Remove base classes
        #print("{} in {} for {}".format("FFRegister", "__init__", cls.__name__))
    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls._registry)

    def __str__(cls):
        if cls in cls._registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])

class FFParserInterface(object):
    """ 
        Interface for parsing force field data file 
    """
    def __init__(self):
        super(FFParserInterface, self).__init__()

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def parse(self, file):
        pass

#class ParamDataParser(object):
#    """
#        A parser's register class for parsing force field parameter file
#    """
#    __metaclass__ = Register
#
#    def __init__(self, *args, **kwargs):
#        super(ParamDataParser, self).__init__()
##        parser_name = kwargs.get("parser_name", "AmberParamParser")
##        self.parser = self.__get_parser(parser_name)
#        #print("{}".format(type(self.parser)))
#
#    def __str__(self):
#        return "{}".format(self.parser.__class__.__name__)
#
#    def __repr__(self):
#        return self.__str__()
#
#    def __getattr__(self, attr):
#        #delegated to the object
#        return getattr(self.parser, attr)
#
##    def __get_parser(self, name):
##        try:
##            for c in self._registry:
##                if c.__name__ == name:
##                    return c()
##        except ValueError:
##            raise ValueError("{} is not registered yet or wrong name".format(name))
#    
##    def parse(self, filename):
##        return self.__getattr__(self.parser, "parse")(filename)

class AmberParamParser(FFParserInterface):
    """
        Parser for Amber Force Field parameter file
    """
    def __init__(self, *args, **kwargs):
        super(AmberParamParser, self).__init__()

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__str__()

    def parse(self, filename):
        print("parse file {}".format(filename))
