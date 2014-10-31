# -*- coding: utf-8 -*-
"""This module defines the sequence class
"""
import os
import sys
import re
import math
import pprint
import numpy as np

from weakref import WeakValueDictionary
from types import *
from bilab.aaprop import AAindex, AAindex1Parser

__all__ = ['Sequence']

class Sequence(object):
    """ 
    Sequence class for single protein sequence

    usage:

    >>> seq = bilab.sequence.Sequence()
    >>> seq.id = 'test'
    >>> seq.sequence = 'AAA'

    """

#        Note: __slots__ works with instance variables


    _instances = WeakValueDictionary()
    #__slots__ = ['id', 'sequence', '__dict__']

    @property
    def Count(self):
        return len(self._instances)

    def __init__(self, *args, **kwargs):
        
        super(Sequence, self).__init__()
        # weak count
        #print("{}".format(id(self)))
        self._obj_id = id(self)
        # initial class
        arg_len = len(args)
        if arg_len == 2:
            self._id = args[0]
            try:
                sequence = ''.join(args[1].split())
                _ = sequence.isalpha()
            except AttributeError:
                raise ValueError('sequence must be a string.')
            else:
                if not _:
                    raise ValueError('not a valid sequence')
            self._sequence = sequence
            self._length = len(sequence)
        else:
            raise ValueError("Wrong number of arguments.")
        self._instances[self._obj_id] = self

    def __str__(self):
        return ">{0}\n{1}".format(self._id, self._sequence)

    def __repr__(self):
        return self.__str__();

    def __getitem__(self, pos):
        return self._sequence[pos]

#    def __del__(self):
#        """ Count references """
#        if self.Count == 0:
#            print 'Last Counter object deleted'
#        else:
#            print self.Count, ' Counter objects remaining'

    def get_sequence(self):
        return self._sequence

    def get_id(self):
        return self._id

    def get_length(self):
        return self._length

    def get_aa_freq(self):
        """ aa frequency """
        if isinstance(self._sequence, bytes):
            self._sequence = self._sequence.decode("utf-8")
        enc = 'utf-16' + ('le' if sys.byteorder == 'little' else 'be')
        # convert into int
        a = np.frombuffer(self._sequence.encode(enc), dtype = np.uint16)
        # count occurrences
        counts = np.bincount(a)
        counts = [(unichr(i), v) for i, v in enumerate(counts) if v]
        return dict(counts)

    def to_numeric(self, scale, smooth=False, func=None):
        """Generate a numeric representation by an AAindex scale
        
            Replace each amino acids in the sequence with the corresponding  
            value of an amino acid scale (e.g., Hydrophobic scale)

        Args:
            scale (string or dict): find the scale name in AAindex db  
                                    if a string is given. Or user defined  
                                    scale in dictionary object.

            smooth (bool): if true, it will apply the :param func to the list  
                           obtained after converting the sequence to numeric  
                           form. Otherwise, the raw numeric form will be  
                           returned.

            func (function): a smoothing function if smooth option is True.

        Kwargs:

        Returns: 
            A list of values

        """
        aa1prop = AAindex(parser = AAindex1Parser)
        numeric_rep = []
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        sequence = std_aa_only(self._sequence)
                # check parameters:
        # if smooth is true and func should be given too
        if smooth and func is None:
            print("if smooth is True, func must be given")
            sys.exit(1)
        elif smooth and type(func) is not FunctionType:
            raise TypeError("The given parameter func must be a function")

        if type(scale) is StringType:
            scale_val = aa1prop.get_scale(scale)
        elif type(scale) is DictType:
            scale_val = scale

        if scale_val is not None:
            for aa in sequence:
                if aa.upper() in scale_val:
                    numeric_rep.append(scale_val[aa.upper()])
                else:
                    # set to zero if amino acid is not found 
                    numeric_rep.append(0.0)
        else:
            return None
        if smooth:
            return func(numeric_rep)
        return numeric_rep

