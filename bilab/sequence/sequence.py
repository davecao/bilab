# -*- coding: utf-8 -*-
import os
import sys
import re
import math
import pprint
import numpy as np

from string import maketrans
from weakref import WeakValueDictionary
from types import *
from bilab.aaprop import AAindex, AAindex1Parser
from bilab.sequence import Alphabet, generic_alphabet,\
                           protein_alphabet, nucleic_alphabet

"""
    This module defines the sequence class,
    revised from weblogo library
"""

__all__ = ['Sequence']


class Sequence(str):
    """
    Sequence class for single protein sequence inherited from str class

    usage:
    >>> import bilab
    >>> from bilab.sequence import Alphabet, generic_alphabet, protein_alphabet, nucleic_alphabet
    >>> pr=bilab.sequence.Alphabet('ACSDGF')
    >>> pr='ACSDGF'
    >>> s1=bilab.sequence.Sequence(pr,name='test',alphabet=protein_alphabet)
    >>> print(repr(s1))
    >test
    ACSDGF
    """
    def __new__(cls, obj,
                alphabet=generic_alphabet,
                name=None,
                description=None,
                isAligned=False):
        self = str.__new__(cls, obj)
        if alphabet is None:
            alphabet = generic_alphabet
        if not isinstance(alphabet, Alphabet):
            alphabet = Alphabet(alphabet)
        if not alphabet.alphabetic(self):
            raise ValueError("Sequence not alphabetic %s, '%s'" %
                             (alphabet, self))
        self._alphabet = alphabet
        self.name = name
        self.description = description
        self.isAligned = isAligned
        return self

    # BEGIN PROPERTIES
    # Make alphabet constant
    @property
    def alphabet(self):
        return self._alphabet
    # END PROPERTIES

    def ords(self):
        """ Convert sequence into an array of integers
        in the range [0, len(alphabet)]
        """
        return self._alphabet.ords(self)

    def tally(self, alphabet=None):
        """ Counts the occurrences of alphabetic characters

        Args:
            alphabet -- an optional alternative alphabet

        returns:
            A list of charcter counts in alphabetic order
        """
        if not alphabet:
            alphabet = self.alphabet
        L = len(alphabet)
        counts = [0, ] * L
        ords = alphabet.ords(self)
        for n in ords:
            if n < L:
                counts[n] += 1
        return counts

    def __getslice__(self, i, j):
        cls = self.__class__
        return cls(str.__getslice__(self, i, j), self.alphabet)

    def __getitem__(self, key):
        # cls = self.__class__
        return (str.__getitem__(self, key), self.alphabet)

    def __add__(self, other):
        # self + other
        cls = self.__class__
        return cls(str.__add__(self, other), self.alphabet)

    def __radd__(self, other):
        # other + self, where other is a superclass of self
        cls = self.__class__
        return cls(str.__add__(self, other), self.alphabet)

    def __eq__(self, other):
        if not hasattr(other, "alphabet"):
            return False
        if self.alphabet != other.alphabet:
            return False
        return str.__eq__(self, other)

    def __ne__(self, other):
        return not self.__eq__(other)


#    def __str__(self):
#        """ self.name and alphabet will be not found when creating
#            an object. This is caused by:
#       __new__  -  if not alphabet.alphabetic(self) :
#            where self is a string
#       """
#        return ">{0}\n{1}".format(self.name, self.alphabet.letters())

    def __repr__(self):
        # return self.__str__();
        return ">{0}\n{1}".format(self.name, self)

#------- Not finished ------------------------

    def join(self, string_list):
        cls = self.__class__
        return cls(super(Sequence, self).join(string_list), self.alphabet)

    def reverse(self):
        """ reverse an order of the sequence
            Return a new string
        """
        cls = self.__class__
        return cls(self[::-1], self.alphabet)

    def ungap(self):
        return self.remove('-.~')

    def remove(self, delchars):
        """ Return a new sequence without specififed chars """
        cls = self.__class__
        cleanseq = ''.join(char for char in str(self)
                           if char not in set(delchars))
        return cls(cleanseq.translate(maketrans('', '')), self.alphabet)

    def translator(self, frm='', to='', delete='', keep=None):
        # Python Cookbook Recipe 1.9
        # Chris Perkins, Raymond Hettinger
        # usage:
        #    >>std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        #    >>std_aa_only("ADDGSFGHXXX")
        #    ADDGSFGH  <--- XXX will be delete
        #    >>
        if len(to) == 1:
            to = to * len(frm)
        trans = maketrans(frm, to)
        if keep is not None:
            allchars = maketrans('', '')
            # delete is expanded to delete everything except
            # what is mentioned in set(keep)-set(delete)
            delete = allchars.translate(
                        allchars,
                        keep.translate(allchars, delete))

        def translate(s):
            return s.translate(trans, delete)
        return translate

    def lower(self):
        """ Return a lowercase string """
        cls = self.__class__
        trans = maketrans('ABCDEFGHIJKLMNOPQRSTUVWXYZ',
                          'abcdefghijklmnopqrstuvwxyz')
        return cls(str(self).translate(trans), self.alphabet)

    def upper(self):
        """ Return a upppercase string """
        cls = self.__class__
        trans = maketrans('abcdefghijklmnopqrstuvwxyz',
                          'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        return cls(str(self).translate(trans), self.alphabet)

    def mask(self, letters='abcdefghijklmnopqrstuvwxyz', mask='X'):
        """
        Replace all occurrences of letters with the mask character.
        The default is to replace all lower case letters with 'X'.
        """
        LL = len(letters)
        if len(mask) != 1:
            raise ValueError("Mask should be single character")
        to = mask * LL
        trans = maketrans(letters, to)
        cls = self.__class__
        return cls(str(self).translate(trans), self.alphabet)

    def words(self, k, alphabet=None):
        """
        Return an iteration over all subwords of length k in the sequence.
        If an optional alphabet is provided, only words from that alphabet
        are returned.

        >>> list(Sequence("abcabc").words(3))
        ['abc', 'bca', 'cab', 'abc']
        """
        if len(self) < k:
            return

        # An optimization. Chopping up strings is faster.
        seq = self.alphabet.normalize(self).tostring()
        # seq = self.tostring()

        for i in range(0, len(seq)-k+1):
            word = seq[i:i+k]
            if alphabet is None or alphabet.alphabetic(word):
                yield word

    def word_count(self, k, alphabet=None):
        """ k-string continue string count

            Sequence("abcabc").word_count(3)
            [('abc', 2), ('bca', 1), ('cab', 1)]
        """
        from itertools import groupby
        words = sorted(self.words(k, alphabet))
        return [(item, sum(1 for n in group))
                for item, group in groupby(words)]


#    def __del__(self):
#        """ Count references """
#        if self.Count == 0:
#            print 'Last Counter object deleted'
#        else:
#            print self.Count, ' Counter objects remaining'

#    def get_aa_freq(self):
#        """ aa frequency """
#        if isinstance(self._sequence, bytes):
#            self._sequence = self._sequence.decode("utf-8")
#        enc = 'utf-16' + ('le' if sys.byteorder == 'little' else 'be')
#        # convert into int
#        a = np.frombuffer(self._sequence.encode(enc), dtype = np.uint16)
#        # count occurrences
#        counts = np.bincount(a)
#        counts = [(unichr(i), v) for i, v in enumerate(counts) if v]
#        return dict(counts)

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
        aa1prop = AAindex(parser=AAindex1Parser)
        numeric_rep = []
        std_aa_only = self.translator(keep="ACDEFGHILKMNPQRSTVWY")
        sequence = std_aa_only(self)
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


