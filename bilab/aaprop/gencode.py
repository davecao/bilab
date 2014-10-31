# -*- coding: utf-8 -*-

__all__ = ['base2aa', 'str2nuc', 'nuc2aa']

nuc2aa_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
"""
.. autofunction:: base2aa

"""
def base2aa(codon):
    """ Get 1-letter amino acid from nucleic acid
    """
    return nuc2aa_table.get(codon.upper(), 'x')

def str2nuc(dna, frame):
    codons = []
    for i in range(frame - 1, len(dna)-2, 3):
        codon = dna[i:i+3]
        codons.append(codon)
    return codons

def nuc2aa(dna, frame=1):
    """ Convert a base to an amino acid

    Args:
        dna (str): the name of base in three letters

    Returns:
        amino acids (str): one letter of amino acids

    Example:

    >>> from bilab.aaprop import nuc2aa
    >>> print nuc2aa('ATA')
    I
    >>> print nuc2aa(''ATAATATTC')
    IIF

    """

    codons = str2nuc(dna, frame)
    amino_acids = ''
    for codon in codons:
        amino_acids = amino_acids + base2aa(codon)
    return amino_acids

