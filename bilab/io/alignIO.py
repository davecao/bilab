# -*- coding: utf-8 -*-
from __future__ import print_function

from functools import  wraps
from types import StringTypes, FileType
import os, sys
import re
import pprint
import numpy as np

from bilab.sequence import Sequence
from bilab import PY3K

#PY2 = sys.version_info[0] == 2
#PY3 = sys.version_info[0] == 3

if PY3K:
    #On Python 3, this will be a unicode StringIO
    from io import StringIO
else:
    # On Python 2 this will be a (bytes) string based handle.
    # Note this doesn't work as it is unicode based:
    # from io import StringIO
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO

__all__=["AlignIO", "MultiFastaIO"]

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

class ClassRegistry(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(ClassRegistry, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, 'registry'):
            cls.registry = set()
        cls.registry.add(cls)
        cls.registry -= set(bases) #Remove base classes
    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls.registry)

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])

class AlignIO(object):
    """
    .. note::
        This class uses two patterns, composite and registry
    """
    __metaclass__ = ClassRegistry

    def __init__(self, ConcreteIO=None, *args, **kwargs):
        """Initialization.

        Args:
           IO (class):

        Kwargs:
        """
        super(AlignIO, self).__init__()
        self.__ConcreteIO = None

        if ConcreteIO:
            # instance, developped by a user or subclasses of Kernel
            self.__ConcreteIO = ConcreteIO()
        elif self.__isstr(ConcreteIO):
            # is string
            self.__ConcreteIO = self.__create(ConcreteIO, *args, **kwargs)

    def __create(self, clsname, *args, **kwargs):
        """ Create an instance for a given class in str or name
        """
        obj = None
        for cls in self.registry:
            if clsname == cls.__name__:
                obj = cls(args)
            elif clsname == cls:
                obj = cls(args)
        if obj:
            # Initialize object
            obj.__init__(*args, **kwargs)
        else:
            print("Unknown class {}".format())
        return obj

    def __str__(self):
        #OK:    return self.__IO.__str__()
        #Failed return self.__getattr__(self.__IO, '__str__')
        IO_func_str_ = self.__getattr__(self.__ConcreteIO, '__str__')
        return IO_func_str_()

    def __repr__(self):
        #return self.__IO.__repr__()
        IO_func_str_ = self.__getattr__(self.__ConcreteIO, '__repr__')
        return IO_func_str_()

    def __call__(self, *args, **kwargs):
        return self.__ConcreteIO(*args, **kwargs)

    def __getattr__(self, obj, attr):
        """ get delegation to the object """
        #return getattr(obj, attr)
        try:
            return self.__ConcreteIO.__getattribute__(attr)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr))

    def __isstr(self, s):
        try:
            return isinstance(s, basestring)
        except NameError:
            return isinstance(s, str)

class IOBase(object):
    """ Base class for Conrete IO classes """
    def __init__(self, handle, *args, **kwargs):
        """Initialization.
        Args:
           *args:
        Kwargs:
            **kwargs:
        """
        super(IOBase, self).__init__()
        if isinstance(handle, FileType):
            #file handle
            self.__handle = handle
        elif isinstance(handle, StringTypes):
            # is string
            try:
                self.__handle = open(handle, 'r')
            except IOError as (errno, strerr):
                print ("I/O error({0}): {1}".formt(errno, strerr))
            except:
                print("Unexpected error:{0}".format(sys.exec_info()[0]))
        else:
            raise IOError("Unknown argument: a file handle or a string")

    @abstractmethod
    def parse(self, handle, alphabet=None):
        """
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """

class MultiFastaIO(IOBase):
    """ Derivated class 
    To read multiple alignments in Fasta format
    e.g.
        >id1:...
        ACDEF....
        DDD--...
        >id2:...
        -CD-F...
        DD-D-...
    """
    def __init__(self, handle):
        super(MultiFastaIO, self).__init__(handle)
        self.__handle = handle

    def parse(self, alphabet=None, isAligned=False):
        """ Impletementation of the abstract method defined in IOBase
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """
        handle = self.__handle
        def build_seq(seq, alphabet, header, header_lineno, 
                      comments, isAligned=False):
            try:
                # Create a bilab.sequence object
                title_fields = re.split("\||:", header)
                seq_id = title_fields[0]
                if comments:
                    header += '\n' + '\n'.join(comments)
                s = Sequence(seq, alphabet, name=seq_id, description=header)
            except ValueError:
                raise ValueError(
                    "Failed to parse the file at the line %d: "
                    "Character not in alphabet:%s" %(header_lineno, alphabet))
            return s

        # loop file handle
        if not isinstance(handle, file):
            print("IOError: argument is not a file handle")
            sys.exit(0)
        seqs = []
        header = []
        seq_str = ""
        for lineno, line in enumerate(handle):
            line = line.strip("\n")
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    # Create Sequence object
                    s = build_seq(seq_str, alphabet, header,
                                  header_lineno, comments, isAligned=isAligned)
                    seqs.append(s)
                    seq_str = ""
                    header =None
                header = line[1:]
                header_lineno = lineno
            else:
                seqs_str += line
        return seqs
