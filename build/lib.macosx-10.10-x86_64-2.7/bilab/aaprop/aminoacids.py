# -*- coding: utf-8 -*-

import math

from types import *

__all__ = ['properties', 'get_aaprop_distance']

"""
    properties:
       key: amino acids in three letter
       value: tuple
             0: 1-letter code,
             1: monoisotopic mass
             2: average mass
             3: volume

    M3L: N-trimethyl-lysine in 154 pdbs or so. 2014/08/29
    MSE: Seleno-methionine  in 7761 pdbs or so. 2014/08/29
    CAS: S-(dimethylarsanyl)-L-cysteine in 119 pdbs or so. 2014/08/29

    Volume:
        Zamyatnin, A.A., Protein volume in solution,
        Prog. Biophys. Mol. Biol., 24, 107-123 (1972) PMID: 4566650.

    Distance:
        order:ACDEFGHIKLMNPQRSTVWY

        Schneider and Wrede gave 20x20 these distance values.

        Reference:
        G. Schneider and P. Wrede
        The rational design of amino acid sequences by artificial
        neural networks and simulated molecular evolution: de novo design
        of an idealized leader peptidase cleavage site.
        Biophys J. Feb 1994; 66(2 Pt 1): 335–344.
"""

#properties = {
#    "GLY" : ("G",  57.02146,  57.0519,  60.1),
#    "ALA" : ("A",  71.03711,  71.0788,  88.6),
#    "LEU" : ("L", 113.08406, 113.1594, 166.7),
#    "ILE" : ("I", 113.08406, 113.1594, 166.7),
#    "ARG" : ("R", 156.10111, 156.1875, 173.4),
#    "LYS" : ("K", 128.09496, 128.1741, 168.6),
#    "MET" : ("M", 131.04049, 131.1926, 162.9),
#    "CYS" : ("C", 103.00919, 103.1388, 108.5),
#    "TYR" : ("Y", 163.06333, 163.1760, 193.6),
#    "THR" : ("T", 101.04768, 101.1051, 116.1),
#    "PRO" : ("P",  97.05276,  97.1167, 112.7),
#    "SER" : ("S",  87.03203,  87.0782,  89.0),
#    "TRP" : ("W", 186.07931, 186.2132, 227.8),
#    "ASP" : ("D", 115.02694, 115.0886, 111.1),
#    "GLU" : ("E", 129.04259, 129.1155, 138.4),
#    "ASN" : ("N", 114.04293, 114.1038, 114.1),
#    "GLN" : ("Q", 128.05858, 128.1307, 143.8),
#    "PHE" : ("F", 147.06841, 147.1766, 189.9),
#    "HIS" : ("H", 137.05891, 137.1411, 153.2),
#    "VAL" : ("V",  99.06841,  99.1326, 140.0),
#    "M3L" : ("K",float("NaN"), 189.2800, float("NaN")),
#    "MSE" : ("M",float("NaN"), 196.1100, float("NaN")),
#    "CAS" : ("C",float("NaN"), 225.1400, float("NaN"))
#}
properties = {
    "Code": {
            "GLY" : "G",
            "ALA" : "A",
            "LEU" : "L",
            "ILE" : "I",
            "ARG" : "R",
            "LYS" : "K",
            "MET" : "M",
            "CYS" : "C",
            "TYR" : "Y",
            "THR" : "T",
            "PRO" : "P",
            "SER" : "S",
            "TRP" : "W",
            "ASP" : "D",
            "GLU" : "E",
            "ASN" : "N",
            "GLN" : "Q",
            "PHE" : "F",
            "HIS" : "H",
            "VAL" : "V"},

    "Mass" : {"G":  57.0519,
              "A":  71.0788,
              "L": 113.1594,
              "I": 113.1594,
              "R": 156.1875,
              "K": 128.1741,
              "M": 131.1926,
              "C": 103.1388,
              "Y": 163.1760,
              "T": 101.1051,
              "P":  97.1167,
              "S":  87.0782,
              "W": 186.2132,
              "D": 115.0886,
              "E": 129.1155,
              "N": 114.1038,
              "Q": 128.1307,
              "F": 147.1766,
              "H": 137.1411,
              "V":  99.1326},
    "Volume": {
              "G": 60.1,
              "A": 88.6,
              "L":166.7,
              "I":166.7,
              "R":173.4,
              "K":168.6,
              "M":162.9,
              "C":108.5,
              "Y":193.6,
              "T":116.1,
              "P":112.7,
              "S": 89.0,
              "W":227.8,
              "D":111.1,
              "E":138.4,
              "N":114.1,
              "Q":143.8,
              "F":189.9,
              "H":153.2,
              "V":140.0},
    "distance":[
        [0    , 0.112, 0.819, 0.827, 0.540, 0.208, 0.696, 0.407, 0.891, 0.406,
         0.379, 0.318, 0.191, 0.372, 1    , 0.094, 0.220, 0.273, 0.739, 0.552],
        [0.114, 0    , 0.847, 0.838, 0.437, 0.32 , 0.66 , 0.304, 0.887, 0.301,
         0.277, 0.324, 0.157, 0.341, 1    , 0.176, 0.233, 0.167, 0.639, 0.457],
        [0.729, 0.742, 0    , 0.124, 0.924, 0.697, 0.435, 0.847, 0.249, 0.841,
         0.819, 0.56 , 0.657, 0.584, 0.295, 0.667, 0.649, 0.797, 1    , 0.836],
        [0.79 , 0.788, 0.133, 0    , 0.932, 0.779, 0.406, 0.86 , 0.143, 0.854,
         0.83 , 0.599, 0.688, 0.598, 0.234, 0.726, 0.682, 0.824, 1    , 0.837],
        [0.508, 0.405, 0.977, 0.918, 0    , 0.69 , 0.663, 0.128, 0.903, 0.131,
         0.169, 0.541, 0.42 , 0.459, 1    , 0.548, 0.499, 0.252, 0.207, 0.179],
        [0.206, 0.312, 0.776, 0.807, 0.727, 0    , 0.769, 0.592, 0.894, 0.591,
         0.557, 0.381, 0.323, 0.467, 1    , 0.158, 0.272, 0.464, 0.923, 0.728],
        [0.896, 0.836, 0.629, 0.547, 0.907, 1    , 0    , 0.848, 0.566, 0.842,
         0.825, 0.754, 0.777, 0.716, 0.697, 0.865, 0.834, 0.831, 0.981, 0.821],
        [0.403, 0.296, 0.942, 0.891, 0.134, 0.592, 0.652, 0    , 0.892, 0.013,
         0.057, 0.457, 0.311, 0.383, 1    , 0.443, 0.396, 0.133, 0.339, 0.213],
        [0.889, 0.871, 0.279, 0.149, 0.957, 0.9  , 0.438, 0.899, 0    , 0.892,
         0.871, 0.667, 0.757, 0.639, 0.154, 0.825, 0.759, 0.882, 1    , 0.848],
        [0.405, 0.296, 0.944, 0.892, 0.139, 0.596, 0.653, 0.013, 0.893, 0    ,
         0.062, 0.452, 0.309, 0.376, 1    , 0.443, 0.397, 0.133, 0.341, 0.205],
        [0.383, 0.276, 0.932, 0.879, 0.182, 0.569, 0.648, 0.058, 0.884, 0.062,
         0    , 0.447, 0.285, 0.372, 1    , 0.417, 0.358, 0.12 , 0.391, 0.255],
        [0.424, 0.425, 0.838, 0.835, 0.766, 0.512, 0.78 , 0.615, 0.891, 0.603,
         0.588, 0    , 0.266, 0.175, 1    , 0.361, 0.368, 0.503, 0.945, 0.641],
        [0.22 , 0.179, 0.852, 0.831, 0.515, 0.376, 0.696, 0.363, 0.875, 0.357,
         0.326, 0.231, 0    , 0.228, 1    , 0.196, 0.161, 0.244, 0.72 , 0.481],
        [0.512, 0.462, 0.903, 0.861, 0.671, 0.648, 0.765, 0.532, 0.881, 0.518,
         0.505, 0.181, 0.272, 0    , 1    , 0.461, 0.389, 0.464, 0.831, 0.522],
        [0.919, 0.905, 0.305, 0.225, 0.977, 0.928, 0.498, 0.929, 0.141, 0.92 ,
         0.908, 0.69 , 0.796, 0.668, 0    , 0.86 , 0.808, 0.914, 1    , 0.859],
        [0.1  , 0.185, 0.801, 0.812, 0.622, 0.17 , 0.718, 0.478, 0.883, 0.474,
         0.44 , 0.289, 0.181, 0.358, 1    , 0    , 0.174, 0.342, 0.827, 0.615],
        [0.251, 0.261, 0.83 , 0.812, 0.604, 0.312, 0.737, 0.455, 0.866, 0.453,
         0.403, 0.315, 0.159, 0.322, 1    , 0.185, 0    , 0.345, 0.816, 0.596],
        [0.275, 0.165, 0.9  , 0.867, 0.269, 0.471, 0.649, 0.135, 0.889, 0.134,
         0.12 , 0.38 , 0.212, 0.339, 1    , 0.322, 0.305, 0    , 0.472, 0.31],
        [0.658, 0.56 , 1    , 0.931, 0.196, 0.829, 0.678, 0.305, 0.892, 0.304,
         0.344, 0.631, 0.555, 0.538, 0.968, 0.689, 0.638, 0.418, 0    , 0.204],
        [0.587, 0.478, 1    , 0.932, 0.202, 0.782, 0.678, 0.23 , 0.904, 0.219,
         0.268, 0.512, 0.444, 0.404, 0.995, 0.612, 0.557, 0.328, 0.244, 0]
    ]
}

def isstandard_aminoacid(aa, order="ACDEFGHIKLMNPQRSTVWY"):
    """ Check the given string representing one of 20 standard amino acids

    Args:
        aa (str) : The name of an amino acid in 3-letter or 1-letter

    Returns
        found (bool).
            True:  if the given string is a standard amino acid.
            False: otherwise
        None. if the given amino acid has a wrong name
    """
    aa_upper = aa.upper()
    found = False
    coding = properties['Code']

    if type(aa1) is not StringType:
        print("The first argument should be a string. Return None")
    else:
        if len(aa_upper) == 3 or len(aa_upper) == 1:
            found =  aa_upper in coding or aa_upper in coding.values()
        else:
            print("Improper argument {}".format(aa))

    return found

def get_aaprop_distance(aa1, aa2, not_found = "NaN"):
    """ Get the distance parameter for a pair of amino acids

    Args:
       aa1 (str):  The name of an amino acid.
       aa2 (str):  The name of an amino acid.

    Kwargs:

    Returns:
       float.  The distance value between the two amino acids in the
               Schneider and Wrede's distance matrix.

        None.  If the given amino acids could not be found.

    Example:

    >>> print get_aaprop_distance('ALA', 'LYS')
    0.891
    >>> print get_aaprop_distance('A', 'K')
    0.891

    """
    order = "ACDEFGHIKLMNPQRSTVWY"
    if isstandard_aminoacid(aa1) and isstandard_aminoacid(aa2):
        try:
            aa1inx = order.index[properties['Code'][aa1]]
            aa2inx = order.index[properties['Code'][aa2]]
        except:
            pass
        else:
            aa1inx = order.index[aa1]
            aa2inx = order.index[aa2]
        return properties['distance'][aa1inx][aa2inx]
    else:
        return None
