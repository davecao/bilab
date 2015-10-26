# -*- coding: utf-8 -*-
# from __future__ import absolute_import

from functools import wraps
# from abc import ABCMeta, abstractmethod
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
                ' should be implemented by a concrete class' % (
                    repr(method), filename, line))
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
        #                                 "__init__", cls.__name__))

    def __iter__(cls):
        # Meta methods, called on class objects:
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

# class ParamDataParser(object):
#     """
#         A parser's register class for parsing force field parameter file
#     """
#     __metaclass__ = Register
#
#     def __init__(self, *args, **kwargs):
#         super(ParamDataParser, self).__init__()
#         parser_name = kwargs.get("parser_name", "AmberParamParser")
#         self.parser = self.__get_parser(parser_name)
#         #print("{}".format(type(self.parser)))
#
#     def __str__(self):
#         return "{}".format(self.parser.__class__.__name__)
#
#     def __repr__(self):
#         return self.__str__()
#
#     def __getattr__(self, attr):
#         #delegated to the object
#         return getattr(self.parser, attr)
#
#     def __get_parser(self, name):
#         try:
#             for c in self._registry:
#                 if c.__name__ == name:
#                     return c()
#         except ValueError:
#             raise ValueError("{} is not registered yet or wrong name".format(
#                               name))
#
#     def parse(self, filename):
#        return self.__getattr__(self.parser, "parse")(filename)


class AmberParamParser(FFParserInterface):
    """
        Parser for Amber Force Field parameter file

    Amber Params structure: in Order
        1. Title:
            FortranFormat('20A4')
        2. AtomType
            FortranFormat('A2,1X,F10.2')
            FortranLine(file.readline()[:-1], format)
            -- blank line --
        3. Hydrophilic atom types
            FortranFormat('20(A2,2X)') <-- Hydrophyilic atom
            -- blank line --
        4. BondParameters (Bond length)
            FortranFormat('A2,1X,A2,2F10.2')
            -- blank line --
        5. AngleParameters (Bond angle)
            FortranFormat('A2,1X,A2,1X,A2,2F10.2')
            -- blank line --
        6. DihedralParameters
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2')
            -- blank line --
        7. ImproperParamters (Improper Dihedral Parameters)
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2')
            -- blank line --
        8. HydrongenBondParameters (H-bond 10-12 Potential Parameters)
            FortranFormat('2X,A2,2X,A2,2X,2F10.2')
            -- blank line --
        9. Lennard-Jones AtomType for Non-bonded 6-12 Potential
            FortranFormat(20(A2,2X))
            -- blank line --
        10. Lennard-Jones Parameters (Non-bonded 6-12 Potential)
            FortranFormat('A4,6X,A2')
            -- blank line --
        END
        Caution: FortranFormat('2X,A2,6X,3F10.6')
            'SK' : Slater-Kirkwood parameters are input.
            'RE' : van der Waals radius and the potential well depth
                    parameters are read.
            'AC' : The 6-12 potential coefficients are read.

        10A.  SK,  Slater-Kirkwood if SK in 10
        10B.  RE,  vdW
        10C.  AC
    """
    def __init__(self, *args, **kwargs):
        super(AmberParamParser, self).__init__()

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__str__()

    def parse(self, filename, modifications=[]):
        """
            Parse force field parameters

        Args:
            filename (str) : file name of the parameter file
                            internally, all data files are stored in
                            the {sysconfig.get_config_var('datarootdir')}
                            /share/bilab/data/ff
        Kwargs:
            modifications (tuple of 2 elements) :
                    first is file name
                    second is the name of the modifcation
                modifications of force field parameters
                see  http://ambermd.org/formats.html#topo.cntrl
        """
        def _readAtomTypes(fhandle):
            while True:
                line = fhandle.next()
                if line.isspace():
                    break
                con = FortranLine(line.strip(), formats[0])
                atomType_dict[con[0].strip()] = AmberAtomType(
                                con[0].strip(), con[1])
            return _readAtomTypes

        def _readBondParameters(fhandle):
            while True:
                line = fhandle.next()
                if line.isspace():
                    break
                con = FortranLine(line.strip(), formats[2])
                name1, name2 = sorted(con[0].strip(), con[1].strip())
                bond_dict[(name1, name2)] = AmberBondParameters(
                                            atomType_dict[name1],
                                            atomType_dict[name2],
                                            con[2], con[3])
            return _readBondParameters

        def _readAngleParameters(fhandle):
            while True:
                line = fhandle.next()
                if line.isspace():
                    break
                con = FortranLine(line.strip(), formats[3])
                name1, name2, name3 = sorted(
                            [con[0].strip(), con[1].strip(), con[2].strip()])
                bondangle_dict[(name1, name2, name3)] = \
                    AmberBondAngleParameters(
                        atomType_dict[name1], atomType_dict[name2],
                        atomType_dict[name3], con[3].strip(), con[4].strip())
            return _readAngleParameters

        def _readDihedralParameters(fhandle):
            app = None
            while True:
                line = fhandle.next()
                if line.isspace():
                    break
                con = FortranLine(line.strip(), fromats[4])
                name1, name2, name3, name4 = sorted(
                        [con[0].strip(), con[1].strip(),
                         con[2].strip(), con[3].strip()])
                if app is not None:
                    app.addTerm(con[4], con[5], con[6], con[7])
                    if con[7] >= 0:
                        app = None
                else:
                    p = AmberDihedralParameters(
                        self.atom_types[name1], self.atom_types[name2],
                        self.atom_types[name3], self.atom_types[name4],
                        con[4], con[5], con[6], con[7])
                    if name1 == 'X' and name4 == 'X':
                        dihedral2_dict[(name2, name3)] = p
                    else:
                        dihedral_dict[(name1, name2, name3, name4)] = p
                    if con[7] < 0:
                        app = p
            return _readDihedralParameters

        def _readImproperParameters(fhandle):
            while True:
                line = fhandle.next()
                if line.isspace():
                    break
                con = FortranLine(line.strip(), formats[5])
                name1 = con[0].strip()
                name2 = con[1].strip()
                name3 = con[2].strip()
                name4 = con[3].strip()
                if name1 == 'X':
                    if name2 == 'X':
                        name1 = name3
                        name2 = name4
                        name3 = 'X'
                        name4 = 'X'
                    else:
                        name1 = name3
                        name2, name3 = sorted([name2, name4])
                        name4 = 'X'
                else:
                    name1, name2, name3, name4 = \
                       (name3, ) + sorted[(name1, name2, name4)]
                p = AmberImproperParameters(atomType_dict[name1],
                                            atomType_dict[name2],
                                            atomType_dict[name3],
                                            atomType_dict[name4],
                                            con[4], con[5], con[6], con[7])
                if name4 == 'X':
                    if name3 == 'X':
                        impropers2_dict[(name1, name2)] = p
                    else:
                        impropers1_dict[(name1, name2, name3)] = p
                else:
                    impropers_dict[(name1, name2, name3, name4)] = p

            return _readImproperParameters

        def _readHbondParameters(fhandle):
            while True:
                line = fhandle.next()
                if line.isspace():
                    break
                con = FortranLine(line.strip(), formats[6])
                name1 = con[0]
                name2 = con[1]
                name1, name2 = sorted([name1, name2])
                hbonds_dict[(name1, name2)] = \
                    AmberHbondParameters(atomType_dict[name1],
                                         atomType_dict[name2],
                                         con[2], con[3])
            return _readHbondParameters

        def _readLJParameters(fhandle):
            while True:
                line = fhandle.next()
                con = FortranLine(line.strip(), line_format)
                if con[0] == 'END':
                    break
                set_name = con[0]
                ljpar_ = AmberLJParameterSet(set_name, con[1])
                ljparams_dict[set_name] = ljpar_
                # 10A
                while True:
                    line = lineIter.next()
                    if line.isspace():
                        break
                    con = FortranLine(line.strip(), FortranFormat(
                                    '2X,A2,6X,3F10.6'))
                    name = con[0].strip()
                    ljpar_.addEntry(name, con[1], con[2], con[3])
            return _readLJParameters

        state_ = 0
        title = ""
        atomType_dict = {}
        bond_dict = {}
        bondangle_dict = {}
        dihedral_dict = {}
        dihedral2_dict = {}
        impropers_dict = {}
        impropers1_dict = {}
        impropers2_dict = {}
        hbonds_dict = {}
        lj_equivalent_dict = {}
        ljparams_dict = {}
        formats = [
            FortranFormat('A2,1X,F10.2'),  # AtomTypes
            FortranFormat('20(A2,2X)'),    # Hydrophylic types
            FortranFormat('A2,1X,A2,2F10.2'),  # Bond length
            FortranFormat('A2,1X,A2,1X,A2,2F10.2'),  # Bond angle
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2'),  # dihedral
            FortranFormat('A2,1X,A2,1X,A2,1X,A2,I4,3F15.2'),  # improper
            FortranFormat('2X,A2,2X,A2,2X,2F10.2'),  # H-bond
            FortranFormat('20(A2,2X)'),  # LJ AtomTypes
            FortranFormat('A4,6X,A2')  # LJ Parameters
        ]
        with open(filename, 'r') as fhandle:
            # line by line
            line_format = formats[0]
            # using iterator
            lineIter = iter(fhandle)
            # 1. title
            title = lineIter.next()
            line_format = formats[state_]
            # 2. AtomType
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                atomType_dict[con[0].strip()] = AmberAtomType(
                            con[0].strip(), con[1])
            # 3. hydrophilic
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                for at in con:
                    at = at.strip()
                    if not at:
                        continue
                    try:
                        atomType_dict[at].hydrophilic = True
                    except KeyError:
                        print("Hydrophilic atom {} is not"
                              " existed in Atom types".format(at))
            # 4. Bond length
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                name1, name2 = sorted(con[0].strip(), con[1].strip())
                bond_dict[(name1, name2)] = AmberBondParameters(
                                            atomType_dict[name1],
                                            atomType_dict[name2],
                                            con[2], con[3])
            # 5. Bond Angle
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                name1, name2, name3 = sorted(
                            [con[0].strip(), con[1].strip(), con[2].strip()])
                bondangle_dict[(name1, name2, name3)] = \
                    AmberBondAngleParameters(
                        atomType_dict[name1], atomType_dict[name2],
                        atomType_dict[name3], con[3].strip(), con[4])
            # 6. Dihedral angle
            app = None
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                name1, name2, name3, name4 = sorted(
                        [con[0].strip(), con[1].strip(),
                         con[2].strip(), con[3].strip()])
                if app is not None:
                    app.addTerm(con[4], con[5], con[6], con[7])
                    if con[7] >= 0:
                        app = None
                else:
                    p = AmberDihedralParameters(
                        self.atom_types[name1], self.atom_types[name2],
                        self.atom_types[name3], self.atom_types[name4],
                        con[4], con[5], con[6], con[7])
                    if name1 == 'X' and name4 == 'X':
                        dihedral2_dict[(name2, name3)] = p
                    else:
                        dihedral_dict[(name1, name2, name3, name4)] = p
                    if con[7] < 0:
                        app = p
            # 7. Improper Dihedral Parameters
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                name1 = con[0].strip()
                name2 = con[1].strip()
                name3 = con[2].strip()
                name4 = con[3].strip()
                if name1 == 'X':
                    if name2 == 'X':
                        name1 = name3
                        name2 = name4
                        name3 = 'X'
                        name4 = 'X'
                    else:
                        name1 = name3
                        name2, name3 = sorted([name2, name4])
                        name4 = 'X'
                else:
                    name1, name2, name3, name4 = \
                       (name3, ) + sorted[(name1, name2, name4)]
                p = AmberImproperParameters(atomType_dict[name1],
                                            atomType_dict[name2],
                                            atomType_dict[name3],
                                            atomType_dict[name4],
                                            con[4], con[5], con[6], con[7])
                if name4 == 'X':
                    if name3 == 'X':
                        impropers2_dict[(name1, name2)] = p
                    else:
                        impropers1_dict[(name1, name2, name3)] = p
                else:
                    impropers_dict[(name1, name2, name3, name4)] = p
            # 8. H-bond 10-12 Potential
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                name1 = con[0]
                name2 = con[1]
                name1, name2 = sorted([name1, name2])
                hbonds_dict[(name1, name2)] = \
                    AmberHbondParameters(atomType_dict[name1],
                                         atomType_dict[name2],
                                         con[2], con[3])
            # 9. Atom Types for Non-bonded 6-12 Potential
            while True:
                line = lineIter.next()
                if line.isspace():
                    state_ += 1
                    line_format = formats[state_]
                    break
                con = FortranLine(line.strip(), line_format)
                name1 = con[0]
                for s in con[1:]:
                    name2 = s.strip()
                    lj_equivalent_dict[name2] = name1
            # 10. LJ parameters (END)
            while True:
                line = lineIter.next()
                con = FortranLine(line.strip(), line_format)
                if con[0] == 'END':
                    break
                set_name = con[0]
                ljpar_ = AmberLJParameterSet(set_name, con[1])
                ljparams_dict[set_name] = ljpar_
                # 10A
                while True:
                    line = lineIter.next()
                    if line.isspace():
                        break
                    con = FortranLine(line.strip(), FortranFormat(
                                    '2X,A2,6X,3F10.6'))
                    name = con[0].strip()
                    ljpar_.addEntry(name, con[1], con[2], con[3])
        mod_title = ""
        for mod, ljname in modifications:
            with open(mod, 'r') as fhandle:
                # lineIter = iter(fhandle)
                mod_title = fhandle.next()
                lineIter.next()  # blank line
                while True:
                    line = lineIter.next()
                    kw = line.strip()[:4]
                    if kw == 'MASS':
                        _readAtomTypes(fhandle)
                    elif kw == 'BOND':
                        _readBondParameters(fhandle)
                    elif kw == 'ANGL':
                        _readAngleParameters(fhandle)
                    elif kw == 'DIHE':
                        _readDihedralParameters(fhandle)
                    elif kw == 'IMPR':
                        _readImproperParameters(fhandle)
                    elif kw == 'HBON':
                        _readHbondParameters(fhandle)
                    elif kw == 'NONB':
                        _readLJParameters(fhandle)
