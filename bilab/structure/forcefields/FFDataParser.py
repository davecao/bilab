# -*- coding: utf-8 -*-
#from __future__ import absolute_import

from functools import wraps
#from abc import ABCMeta, abstractmethod
from bilab.io.FortranFormat import FortranFormat, FortranLine
from bilab import PY3K
from bilab import forcefieldsList
from bilab.structure.forcefields import AmberAtomType, AmberBondParameters, \
            AmberBondAngleParameters, AmberDihedralParameters, \
            AmberImproperParameters,  AmberHbondParameters,\
            AmberLJParameterSet

# Replace enumrate with i, item in zip(range(len(foo)), foo)
if PY3K:
    from itertools import zip 
else:
    from itertools import izip as zip


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

    Amber Params structure: in Order 
        Title:  one line
        1. AtomType
            FortranFormat('A2,1X,F10.2')
            FortranLine(file.readline()[:-1], format)
            -- blank line --
        2. Hydrophylic atom types
            FortranFormat('20(A2,2X)') <-- Hydrophylic atom 
            -- blank line --
        3. BondParameters
            FortranFormat('A2,1X,A2,2F10.2')
            -- blank line --
        4. AngleParameters (bond angles)
            FortranFormat('A2,1X,A2,1X,A2,2F10.2')
            -- blank line --
        5. DihedralParameters
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2')
            -- blank line --
        6. ImproperParamters
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2')
            -- blank line --
        7. HydrongenBondParameters
            FortranFormat('2X,A2,2X,A2,2X,2F10.2')
            -- blank line --
        8. Lennard-Jones Parameters
            FortranFormat('2X,A2,6X,3F10.6')
            -- blank line --
        END
    """
    def __init__(self, *args, **kwargs):
        super(AmberParamParser, self).__init__()

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__str__()

    def parse(self, filename):

        AtomType, Hydrophilic, Bond, BondAngle, Dihedral, ImproperParams, HydrongenBond, LJParams = range(8)
        state_ = -1
        title = ""
        atomType_dict = {}
        bond_dict = {}
        bondangle_dict = {}
        dihedral_dict = {}
        impropParams_arr = {}
        hbond_arr = {}
        ljParams_arr = {}
        formats = [
            FortranFormat('A2,1X,F10.2'),
            FortranFormat('20(A2,2X)'),
            FortranFormat('A2,1X,A2,2F10.2'),
            FortranFormat('A2,1X,A2,1X,A2,2F10.2'),
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2'),
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2'),
            FortranFormat('2X,A2,2X,A2,2X,2F10.2'),
            FortranFormat('2X,A2,6X,3F10.6')
        ]
        with open(filename, 'r') as fhandle:
            # line by line
            line_format = formats[0]
            #for lineno, line in zip(range(len(fhandle)), fhandle):
            for line in fhandle:
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    continue

                if state_ == -1 :
                    # read title 
                    title = line.strip()
                    state_ = AtomType

                if state_ == AtomType:
                    con = FortranLine(line.strip(), line_format)
                    #atomType_arr.append(AmberAtomType(con[0].strip(), con[1]))
                    atomType_dict[con[0].strip()] = \
                            AmberAtomType(con[0].strip(), con[1])

                if state_ == Hydrophylic:
                    con = FortranLine(line.strip(), line_format)
                    for at in con:
                        atomType_dict[at.strip()].hydrophylic = True

                if state_ == Bond:
                    con = FortranLine(line.strip(), line_format)
                    name1, name2 = sorted(con[0].strip(), con[1].strip())
                    bond_dict[(name1, name2)] = \
                        AmberBondParameters(atomType_dict[name1],
                                            atomType_dict[name2],
                                            con[2], con[3])
                if state_ == BondAngle:
                    con = FortranLine(line.strip(), line_format)
                    name1, name2, name3 = sorted([con[0].strip(), 
                                                  con[1].strip(), 
                                                  con[2].strip()])
                    bondangle_dict[(name1, name2, name3)] = \
                        AmberBondAngleParameters(atomType_dict[name1],
                                                 atomType_dict[name2],
                                                 atomType_dict[name3],
                                                 con[3].strip(),
                                                 con[4].strip())
                if state_ == Dihedral:
                    con = FortranLine(line.strip(), line_format)
                    name1, name2, name3, name4 = sorted([con[0].strip(), 
                                                  con[1].strip(), 
                                                  con[2].strip(),
                                                  con[3].strip()])

                    p = AmberDihedralParameters(self.atom_types[name1],
                                            self.atom_types[name2],
                                            self.atom_types[name3],
                                            self.atom_types[name4],
                                            con[4], con[5], con[6], con[7])
