# -*- coding: utf-8 -*-
"""This module defines a profile class
"""
from __future__ import print_function

import numpy as np
from bilab.sequence import Sequence, Alphabet

__all__ = ['Profile']

class Profile(object):
    """Generating a profile from multiple sequence alignment

    """
    def __init__(self, sequenceList, alphabet=None):
        super(Profile, self).__init__()
        self.sequenceList = sequenceList # Array of bilab.sequence.Sequence
        # Check obj in sequenceList ?
        self.alphabet = alphabet
        self.gap = '-'

    def __gap_unknonw_index(self, alphabet=None):
        """ Add ordinal numbers for gap and letters which are not included in
            alphabet
        """
        if not alphabet:
            alphabet = self.alphabet
        if (alphabet is None) or not isinstance(alphabet, Alphabet):
            raise ValueError("alphabet is not set correctly")
        alphabet_length = len(alphabet)
        if '-' not in alphabet:
            gapIndex = alphabet_length + 1
            unknown_Index = gapIndex + 1
        else:
            unknown_Index = alphabet_length + 1
        return (gapIndex, unknown_Index)

    def __ords(self, alphabet, gapIndex=None, unknownIndex=None):
        """ Convert sequences into 2d array of ordinals
        """
        if not alphabet:
            raise ValueError("Unknown alphabet")
        s_ords = []
        for seq in self.sequenceList:
            temp_ords = []
            for aa in seq:
                residue = aa[0]
                if residue in alphabet:
                    temp_ords.append(alphabet.ord(residue))
                elif residue == self.gap:
                    # if gap is not in alphabet
                    temp_ords.append(gapIndex)
                else:
                    # unknonw residue
                    temp_ords.append(unknownIndex)
            s_ords.append(temp_ords)
        return s_ords

    def generate(self, alphabet=None):
        """ Generate frequency profile 
        Args:
            alphabet - Amino acids count for each position
        Kwargs:

        Returns:
            a numpy 2D array of frequencies of specified amino acids 
            row means positions of amino acid in the sequence
            column means frequencies of specified amino acid occurred in the 
            order of the given alphabet.
         note: if the value divided by zero, it will return zero, not nan 
        """
        if not alphabet:
            alphabet = self.alphabet
        if (alphabet is None) or not isinstance(alphabet, Alphabet):
            raise ValueError("alphabet is not set correctly")
        (gapIndex, unknown_Index) = self.__gap_unknonw_index(alphabet)
        N = len(alphabet)
        ords = self.__ords(alphabet, 
                            gapIndex=gapIndex, 
                            unknownIndex=unknown_Index)
        L = len(ords[0])
        #counts = [ [0,]*N for l in range(0,L)]
        # L positions and N kinds of amino acids
        counts = np.zeros((L, N), dtype=float)
        for o in ords :
            if len(o)!=L : 
                raise ValueError("Sequences are of incommensurate lengths.")
            for j,n in enumerate(o) :
                if n<N : 
                    counts[ j][n] +=1
        # sum at each position
        aa_sum = np.sum(counts, axis=1)
        # 20 aa at each position
        # if divisor is zero, 
        # disable RuntimeWarining: invalid value encountered in divide
        old_seterr = np.geterr()
        np.seterr(divide='ignore',invalid='ignore')
        fre_aa = np.nan_to_num(counts/np.asmatrix(aa_sum).T).flatten().tolist()
        np.seterr(**old_seterr)
        return fre_aa

    def profileHMM(self, alphabet=None):
        """
        Generate a profile HMM
        """
        # get frequencies of alphabet occurred in positions
        fre = self.generate(alphabet)
