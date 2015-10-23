# -*- coding: utf-8 -*-
#from __future__ import absolute_import

from functools import wraps
#from abc import ABCMeta, abstractmethod
from bilab.io.FortranFormat import FortranFormat, FortranLine
from bilab import PY3K
from bilab import forcefieldsList

__all__ = ['ForceField', 'AmberForceField']

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

class ParamDataParser(object):
    """
        A parser's register class for parsing force field parameter file
    """
    __metaclass__ = Register

    def __init__(self, *args, **kwargs):
        super(ParamDataParser, self).__init__()
#        parser_name = kwargs.get("parser_name", "AmberParamParser")
#        self.parser = self.__get_parser(parser_name)
        #print("{}".format(type(self.parser)))

    def __str__(self):
        return "{}".format(self.parser.__class__.__name__)

    def __repr__(self):
        return self.__str__()

    def __getattr__(self, attr):
        #delegated to the object
        return getattr(self.parser, attr)

#    def __get_parser(self, name):
#        try:
#            for c in self._registry:
#                if c.__name__ == name:
#                    return c()
#        except ValueError:
#            raise ValueError("{} is not registered yet or wrong name".format(name))
    
#    def parse(self, filename):
#        return self.__getattr__(self.parser, "parse")(filename)

class AmberParamParser(FFParserInterface, ParamDataParser):
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


class ForceField(object):
    """ 
        interface of Force Field 
        Base class

        Note: multiple inherited subclasses may cause metaclass conflict
    """
    __metaclass__ = Register
    #__metaclass__ = ABCMeta

    def __init__(self, *args, **kwargs):
        """ Initialization """
        super(ForceField, self).__init__()

#    @abstractmethod
#    def __call__(cls, clsname, *args, **kwargs):
#        """ ForceField callable 
#            Generate an object of specified subclass name
#        """
        #return cls.__metaclass__._registry(clsname)(*args, **kwargs)
#        print("{} in {}".format("ForceField", "__call__"))
        #print("{}".format(self))
        #print("{}".format(self.__class__))
        #return self.__class__.__Factory(clsname)(*args, **kwargs)
#        try:
#            for c in cls.__registry:
#                if c.__name__ == clsname:
#                    return c(*args, **kwargs)
#                #cls.__metaclass__.registry[name]
#        except ValueError:
#            raise ValueError("{} is not registered yet or wrong name". format(clsname)#)

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
        #print("classmethod Factory is called by {}".format(cls.__name__))
        try:
            for c in cls._registry:
                if c.__name__ == name:
                    return c()
                #cls.__metaclass__.registry[name]
        except ValueError:
            raise ValueError("{} is not registered yet or wrong name".format(name))

class AmberForceField(ForceField):
    """
        Amber Force Field
    """
    def __init__(self, *args, **kwargs):
        try:
            self.ff_name = kwargs.get("ff_name", "amber94")
            self.data_parser = kwargs.get("parser", AmberParamParser())

                                #ParamDataParser("AmberParamParser"))
            self.ff_dat_loc = forcefieldsList[self.ff_name]
        except KeyError:
            raise KeyError("Wrong name of the force field: {}".format(self.ff_name))
        # ff_params is a dict
        # 
        #print("{}".format(type(self.data_parser)))
        self.ff_params = self.data_parser.parse(self.ff_dat_loc)

    def __str__(self):
        return "{} from {}".format(self.ff_name, self.ff_dat_loc)
    
    def __repr__(self):
        return self.ff_name




#ForceField.register(AmberForceField)
