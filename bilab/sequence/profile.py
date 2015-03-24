# -*- coding: utf-8 -*-
"""This module defines a profile class
"""
from __future__ import print_function
from bilab.sequence import Sequence, Alphabet

__all__ = ['Profile']

class Profile(object):
    """Generating a profile from multiple sequence alignment

    """
    def __init__(self, sequenceList):
        super(Profile, self).__init__()
        self.sequenceList = sequenceList # Array of bilab.sequence.Sequence
        # Check obj in sequenceList ?

    def __ords(self, alphabet):
        """ Convert sequences into 2d array of ordinals
        """
        if not alphabet:
            raise ValueError("Unknown alphabet")
        s_ords = []
        for s in self.sequenceList:
            s_ords.append(s.ords())
        return s_ords

    def generate(self, alphabet=None):
        """ Generate profile """
        if (alphabet is None) or not isinstance(alphabet, Alphabet):
            raise ValueError("alphabet is not set correctly")
        N = len(alphabet)
        ords = self.__ords(alphabet)
        L = len(ords[0])
        counts = [ [0,]*N for l in range(0,L)]

        for o in ords :
            if len(o)!=L : 
                raise ValueError("Sequences are of incommensurate lengths.")
            for j,n in enumerate(o) :
                if n<N : 
                    counts[ j][n] +=1
        return counts
