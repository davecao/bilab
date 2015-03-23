# -*- coding: utf-8 -*-
from __future__ import print_function

from functools import  wraps
from types import StringTypes, FileType, ClassType, InstanceType
import os, sys
import re
import pprint
import numpy as np

from bilab.sequence import Sequence, Alphabet
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

#__all__=["AlignIO", "MultiFastaIO", "ClustalWIO"]
__all__=["AlignIO"]

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

    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           IO (class):

        Kwargs:
        """
        super(AlignIO, self).__init__()
        ConcreteIO_cls = kwargs.pop("ConcreteIO", None)
        if ConcreteIO_cls is None:
            return
#        if isinstance(ConcreteIO, ClassType):
#            # developped by a user or subclasses of Kernel
#            if ConcreteIO in self.registry:
#                self.__ConcreteIO = ConcreteIO(*args, **kwargs)
#        elif isinstance(ConcreteIO, InstanceType):

#            self.__ConcreteIO = ConcreteIO
#        elif self.__isstr(ConcreteIO):
#            # is string
        #self.__ConcreteIO = 
        self.__ConcreteIO = self.__create(ConcreteIO_cls, *args, **kwargs)

    def __create(self, clsname, *args, **kwargs):
        """ Create an instance for a given class in str or name
        """
        obj = None
        for cls in self.registry:
            if clsname == cls.__name__ or clsname == cls:
                obj = cls(*args)
        if obj:
            # Initialize object
            obj.__init__(*args, **kwargs)
        else:
            print("Unknown class {}".format(clsname))
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

    def __getattr__(self, attr):
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
            self.handle = handle
        elif isinstance(handle, StringTypes):
            # is string
            try:
                self.handle = open(handle, 'r')
            except IOError as (errno, strerr):
                print ("I/O error({0}): {1}".formt(errno, strerr))
            except:
                print("Unexpected error:{0}".format(sys.exec_info()[0]))
        else:
            raise IOError("Unknown argument: a file handle or a string")
        return self

    @abstractmethod
    def parse(self, alphabet=None, isAligned=False):
        """
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """

class MultiFastaIO(IOBase, AlignIO):
    """ Derived class 
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

    def parse(self, alphabet=None, isAligned=False):
        """ Impletementation of the abstract method defined in IOBase
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """
        handle = self.handle
        alphabet = Alphabet(alphabet) 
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
        seq_str = ""
        comments = []
        header = None
        header_lineno = 0
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
                    comments = []
                header = line[1:]
                header_lineno = lineno
            elif line.startswith(';'):
                comments.append(line[1:])
            else:
                seq_str += line
        # Store last one
        s = build_seq(seq_str, alphabet, header,
                                  header_lineno, comments, isAligned=isAligned)
        seqs.append(s)
        return seqs

class ClustalWIO(IOBase, AlignIO):
    """
    Derived class
    To read multiple alignments in clustalw format
    e.g.
CLUSTAL W (1.81) multiple sequence alignment


CXCR3_MOUSE       --------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
BLR_HUMAN         --------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
CXCR1_HUMAN       --------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
CXCR4_MURINE      -----------------------------------YTSDN---------YSGSGDYDSNK
                                                     :  :          :..     ..
 
CXCR3_MOUSE       -SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
BLR_HUMAN         --LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
CXCR1_HUMAN       -SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
CXCR4_MURINE      -EPC-RDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVAD
                             :.  .:   * ::** .::**  *  ::   :   * *: : ::*::**

CXCR3_MOUSE       VLLVLTLPLWAVDAA-VQWVFGPGLCKVAGALFNINFYAGAFLLACISFDRYLSIVHATQ
BLR_HUMAN         LLLVFILPFAVAEGS-VGWVLGTFLCKTVIALHKVNFYCSSLLLACIAVDRYLAIVHAVH
CXCR1_HUMAN       LLFALTLPIWAASKV-NGWIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATR
CXCR4_MURINE      LLFVITLPFWAVDAM-ADWYFGKFLCKAVHIIYTVNLYSSVLILAFISLDRYLAIVHATN
                  :*:.: **: ...     * :*  ***..  :  :*:*.. ::** *:.****:****..
    

    """
    class Token(object):
        """Represents the items returned by a file scanner, normally processed
        by a parser.
        
        Attributes :
        o typeof    -- a string describing the kind of token
        o data      -- the value of the token
        o lineno    -- the line of the file on which the data was found (if known)
        o offset    -- the offset of the data within the line (if known)
        """
        __slots__ = [ 'typeof', 'data', 'lineno', 'offset']
        def __init__(self, typeof, data=None, lineno=-1, offset=-1) :
            self.typeof = typeof
            self.data = data
            self.lineno = lineno
            self.offset = offset
     
        def __repr__(self) :
            return stdrepr( self) 

        def __str__(self):
            coord = str(self.lineno)
            if self.offset != -1 : coord += ':'+str(self.offset)
            coord = coord.ljust(7)
            return (coord+ '  '+ self.typeof +' : ').ljust(32)+ str(self.data or '')

    def __init__(self, handle):
        super(ClustalWIO, self).__init__(handle)


    def parse(self, alphabet=None, isAligned=False):
        """ Impletementation of the abstract method defined in IOBase
            referred to clustal_io.py in weblogo
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """

        def scan(handle):
            """Scan a clustal format MSA file and yield tokens.
            The basic file structure is
            
            begin_document
                header?     
               (begin_block
                   (seq_id seq seq_index?)+
                   match_line?
               end_block)*
            end_document     
    
            Usage:
                for token in scan(clustal_file):
                    do_something(token)
            """
            header, body, block = range(3)
            yield self.Token("begin")
            leader_width = -1
            state = header
            for L, line in enumerate(handle):
                if state==header :
                    if line.isspace() : continue
                    m = header_line.match(line)
                    state = body
                    if m is not None :
                        yield self.Token("header", m.group() )
                        continue
                    # Just keep going and hope for the best.
                    #else :
                        #raise ValueError("Cannot find required header")

                if state == body :
                    if line.isspace() : continue
                    yield self.Token("begin_block")
                    state = block
                    # fall through to block
                
                if state ==  block:
                    if line.isspace() :
                        yield self.Token("end_block")
                        state = body
                        continue
                    
                    m = match_line.match(line)
                    if m is not None :
                        yield self.Token("match_line", line[leader_width:-1])
                        continue
             
                    m = seq_line.match(line) 
                    if m is None: 
                        raise ValueError("Parse error on line: %d (%s)" % (L,line))
                    leader_width = len(m.group(1))
                    yield self.Token("seq_id", m.group(1).strip() )
                    yield self.Token("seq", m.group(2).strip() )
                    if m.group(3)  :
                        yield self.Token("seq_num", m.group(3)) 
                    continue

                # END state blocks. If I ever get here something has gone terrible wrong
                raise RuntimeError()
            
            if state==block:
                 yield self.Token("end_block")
            yield self.Token("end")     
            return

        handle = self.handle
        header_line = re.compile(r'(CLUSTAL.*)$')
        # (sequence_id) (Sequence) (Optional sequence number)
        seq_line   = re.compile(r'(\s*\S+\s+)(\S+)\s*(\d*)\s*$')
        # Saved group includes variable length leading space.
        # Must consult a seq_line to figure out how long the leading space is since
        # the maximum CLUSTAL ids length (normally 10 characters) can be changed.
        match_line = re.compile(r'([\s:\.\*]*)$')
        alphabet = Alphabet(alphabet) 
        seq_ids = []
        seqs = []
        block_counts = 0
        data_len = 0
        for token in scan(handle):
            if token.typeof== "begin_block":
                block_count = 0
            elif token.typeof == "seq_id":
                if len(seqs) <= block_count :
                    seq_ids.append(token.data)
                    seqs.append([])
            elif token.typeof == "seq":
                if not alphabet.alphabetic(token.data) :
                    raise ValueError(
                        "Character on line: %d not in alphabet: %s : %s" % (
                        token.lineno, alphabet, token.data) )
                seqs[block_count].append(token.data)
                if block_count==0 :
                    data_len = len(token.data) 
                elif data_len != len(token.data) :
                    raise ValueError("Inconsistent line lengths")
                    
                block_count +=1
        seqs = [Sequence("".join(s), alphabet, name=i) for s, i in zip(seqs, seq_ids)]
        return seqs
