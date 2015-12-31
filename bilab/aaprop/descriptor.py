# -*- coding: utf-8 -*-
import sys
import numpy as np

from types import *
from bilab.utilities import translator
from bilab.aaprop.aminoacids import properties

__author__ = "Wei Cao"
__copyright__ = "Copyright 2014, Bioinformation Engineering Laboratory"
__version__ = "1.0.0-dev"
__maintainer__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2014/08/28"
__status__ = "Dev"


__all__ = [
           'AutoCorrelatedDescriptor',
           'DipeptideDescriptor',
           'TripeptideDescriptor',
           'PseudoAACompoDescriptor',
           'QuasiOrderDescriptor',
           'UnifiedHilbertDescriptor',
           'ZpCurveDescriptor'
           ]

"""
    A vistor pattern
"""


class prDescriptor(object):
    """ base class for descriptors """
    def __init__(self):
        """
        prDescriptor is a base class, implemented in a visitor pattern

        Note:

        Args:
          visitee (class) : The class is visited.

        """
        super(prDescriptor, self).__init__()

    def generate(self, visitee, *args, **kwargs):
        """
        """
        meth = None
        for cls in visitee.__class__.__mro__:
            method_name = 'visitee_' + cls.__name__
            meth = getattr(self, method_name, None)
            if meth:
                break
        if not meth:
            meth = self.generic_visit

    def generic_visit(self, visitee, *args, **kwargs):
        print('generic_visit ' + visitee.__class__.__name__)


class AutoCorrelatedDescriptor(prDescriptor):
    """
    ***Auto-correlation function***

    .. math::
        Ri = \\frac{1}{L-j}\\sum_{i=1}^{L-j} h_{i}h_{i+j}

        h_{i} = \\frac{h_{i}-average(h)}{\\sigma}

    where, j=1,2, ..., m and :math:`h_{i}` is an aaindex value of
    ith residue of the protein.

    .. note::
        Zhang, C. T., Lin, Z. S., Zhang, Z., and Yan, M. (1998).
        Prediction of the helix/strand content of globular proteins
        based on their primary sequences.
        **Protein Eng.** 11, 971-979.

    """

    def generate(self, visitee, *args, **kwargs):
        """Generate feature of the sequence


        Args:
            visitee: the instance of class Sequence

        Kwargs:
            scale (string): AAindex name (string), default is 'KYTJ820101'
            m (int): m terms of autocorrelation function. default is 12
            verbose: show verbose info (bool)

        Returns:
            A list of values as a feature vector

        """
        scale_name = kwargs.pop('scale', 'KYTJ820101')
        m = kwargs.pop('m', 12)
        verbose = kwargs.pop('verbose', False)
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        seq = std_aa_only(visitee.get_sequence())
        length = len(seq)
        if length < m:
            print("The length of the given sequence is \
                    shorter than {}".format(m))
            sys.exit(0)
        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        # generate a numeric representation of the sequence
        val = np.asarray(visitee.to_numeric(scale_name))
        ave_h = val.mean()
        sigma = val.std()
        h = (val - val.mean()) / val.std()
        f_val = []

        for j in range(m):
            hi = h[0:h.size - j]
            hj = h[j:]
            fj = np.correlate(hi, hj)/(h.size - j)
            f_val.extend(fj.flatten())

        return f_val


class DipeptideDescriptor(prDescriptor):
    """

    **Description:**
        A protein sequence is represented by a 400-D feature vector.
        Each element of 20x20 matrix is the occurrence of AA pair given
        a protein sequence.

    Args:
        0. numeric values of a protein sequence (in list)

    Kwargs:
        order = **ACDEFGHILKMNPQRSTVWY**

    Returns:
        A list of values as a feature vector

    """
    def generate(self, visitee, *args, **kwargs):
        """ Generate feature of the sequence

        Note:

        Args:
            visitee (class): the instance of class Sequence

        Kwargs:
            order (string): the order of amino acids
            verbose (bool): show verbose info

        Returns:
            A list of values as a feature vector

        """
        order = kwargs.pop('order', "ACDEFGHILKMNPQRSTVWY")
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        verbose = kwargs.pop('verbose', False)
        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        seq = std_aa_only(visitee.get_sequence())
        length = len(seq)
        rows = seq[:length - 1]
        cols = seq[1:]
        feature_vec = np.zeros((len(order), len(order)))
        for aai, aaj in zip(rows, cols):
            inx1 = order.index(aai)
            inx2 = order.index(aaj)
            feature_vec[inx1, inx2] += 1
        occurrences = np.reshape(feature_vec, feature_vec.size)
        frequencies = occurrences / np.sum(occurrences)
        return frequencies


class TripeptideDescriptor(prDescriptor):
    """
    **Description:**
        A protein sequence is represented by a 8000-D feature vector.
        Each element of 20x20x20 matrix is the occurrence of AA pair given
        a protein sequence.

      args:
        0: numeric values of a protein sequence (in list)

      kwargs:
        order = **ACDEFGHILKMNPQRSTVWY**

    """

    def generate(self, visitee, *args, **kwargs):
        """Generate feature of the sequence

        Note:

        Args:
            visitee (class): the instance of class Sequence

        Kwargs:
            order (string): the order of amino acids
            verbose (bool): show verbose info

        Returns:
            A list of values as a feature vector

        """
        order = kwargs.pop('order', "ACDEFGHILKMNPQRSTVWY")
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        verbose = kwargs.pop('verbose', False)
        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        seq = std_aa_only(visitee.get_sequence())
        length = len(seq)
        rows = seq[:length - 2]
        cols = seq[1:length - 1]
        hight = seq[2:]

        feature_vec = np.zeros((len(order), len(order), len(order)))
        for aai, aaj, aah in zip(rows, cols, hight):
            inx1 = order.index(aai)
            inx2 = order.index(aaj)
            inx3 = order.index(aah)
            feature_vec[inx1, inx2, inx3] += 1
        occurrences = np.reshape(feature_vec, feature_vec.size)
        frequencies = occurrences / np.sum(occurrences)
        return frequencies


class PseudoAACompoDescriptor(prDescriptor):
    """
    **Pseudo AA Composition**

    .. math::
        \\theta_{j} = \\frac{1}{L-j}\\sum_{i=1}^{L-j}C_{i, i+j}

        C(R_{i}, R_{j})=\\frac{1}{3}
        {[H1(R_{j})-H1(R_{i})]^{2} +
        [H2(R_{j})-H2(R_{i})]^{2} + [H3(R_{j})-H3(R_{i})]^{2}}

    where j=1, 2, ..., L-1 (L: the length of protein X)

    :math:`H1(R_{i})`:the hydrophobicity value of :math:`R_{i}`
    from Hopp and Woods.
    :math:`H2(R_{i})`:the hydrophilicity value of :math:`R_{i}`
    from  Tanford
    :math:`H3(R_{i})`:the average mass of the amino acid :math:`R_{i}`

    .. note::
        Chou, K. C. (2001).
        Prediction of protein cellular attributes
        using pseudo-amino acid composition.
        **Proteins** 43, 246-255

    """

    def generate(self, visitee, *args, **kwargs):
        """Generate feature of the sequence

        Note:

        Args:
            visitee (class): the instance of class Sequence

        Kwargs:
            nterms(int) : the number of features will be generated.
                          default is 20 (used in Chou's paper)

            h1 (string or dict) : the name of a hydrophobicity scale.
                          default is Hopp and Woods.

            h2 (string or dict) : the name of a hydrophobicity scale.
                          default is Tanford

            h3 (string or dict) : the name of side-chain mass scale

            verbose (bool): show verbose info

        Returns:
            A list of values as a feature vector

        """
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        nterms = kwargs.pop('nterms', 20)
        h1_scale = kwargs.pop('h1', 'HOPT810101')
        h2_scale = kwargs.pop('h2', 'NOZY710101')
        h3_scale = kwargs.pop('h3', properties['Mass'])
        verbose = kwargs.pop('verbose', False)

        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        seq = std_aa_only(visitee.get_sequence())

        h1_val = np.asarray(visitee.to_numeric(h1_scale))
        h2_val = np.asarray(visitee.to_numeric(h2_scale))
        h3_val = np.asarray(visitee.to_numeric(h3_scale))

        h1 = (h1_val - h1_val.mean()) / h1_val.std()
        h2 = (h2_val - h2_val.mean()) / h2_val.std()
        h3 = (h3_val - h3_val.mean()) / h3_val.std()

        length = len(h1)
        if length < nterms:
            print ("The sequence is too shorter to generate {} features"
                   .format(nterms))
            sys.exit(0)

        theta = []

        for j in range(nterms):

            h1_Ri = h1[0: h1.size - j]
            h1_Rj = h1[j:]
            h1_val = (h1_Rj - h1_Ri)**2

            h2_Ri = h2[0: h2.size - j]
            h2_Rj = h2[j:]
            h2_val = (h2_Rj - h2_Ri)**2

            h3_Ri = h3[0: h3.size - j]
            h3_Rj = h3[j:]
            h3_val = (h3_Rj - h3_Ri)**2

            f_val = 1/3 * (h1_val + h2_val + h3_val)/(length - j)

            theta.extend(f_val)
        return theta


class QuasiOrderDescriptor(prDescriptor):
    """
    QuasiOrderDescriptor

    **Method**:

    :math:`\\tau(j)` is called jth rank sequence-order-coupling number

    .. math::
        \\tau(j)=\\frac{1}{L-j}J_{i, i+j}

        J_{i, k} = D(R_{i},R_{k}) * D(R_{i},R_{k})

    where j=1,2,...,L-1 (L: length of protein X) and the coupling factor
    :math:`J_{i,k}` is a function of amino acids
    :math:`R_{i}` and :math:`R_{k}`
    :math:`D(R_{i},R_{k})` is the distance between :math:`R_{i}` and
    :math:`R_{k}` in 20x20 table shown in Schneider and Wrede's work.

    .. note::
        Chou, K. C. (2000b).
        Prediction of protein subcellular locations by
        incorporating quasi-sequence-order effect.
        **Biochem. Biophys. Res. Commun.** 19,477-83.
    """
    def generate(self, visitee, *args, **kwargs):
        """Generate feature of the sequence

        Note:

        Args:
            visitee (class): the instance of class Sequence

        Kwargs:
            nterms : the number of features will be generated.
                          default is 12 (used in Chou's paper)
        Returns:
            A list of values as a feature vector
        """
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        nterms = kwargs.pop('nterms', 12)
        verbose = kwargs.pop('verbose', False)

        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        seq = std_aa_only(visitee.get_sequence())
        length = len(seq)
        Tau = []

        for j in range(nterms):
            J_Ri = seq[0: length - j]
            J_Rj = seq[j: length]
            J_val = 0
            for (aa1, aa2) in zip(J_Ri, J_Rj):
                dist_aa1_aa2 = bilab.aaprop.get_aaprop_distance(aa1, aa2)
                dist_aa2_aa1 = bilab.aaprop.get_aaprop_distance(aa2, aa1)
                J_val += dist_aa1_aa2 * dist_aa2_aa1
            Tau.append(J_val / length - j)
        return Tau


class UnifiedHilbertDescriptor(prDescriptor):
    """
    UnifiedHilbertDescriptor

    **Method**:

    unified attribute vector in Hilbert Space :math:`X=(a1, a2, . . . , a20)'`;
    a(i) is square root of the occurrence frequency of one of 20 native
    amino acids.

    .. note ::
      Feng, Z. P. and Zhang, C. T. (2001).
      Prediction of the subcellular location of
      prokaryotic proteins based on
      the hydrophobic index of the amino acids.
      **Int. J. Biol. Macromol.** 14, 255-261.

    """
    def generate(self, visitee, *args, **kwargs):
        """Generate feature of the sequence

        Args:
            visitee (class): the instance of class Sequence

        Kwargs:
            order (str): a string of amino acids in 1-letter code
                    default is "ACDEFGHILKMNPQRSTVWY".
                    e.g., order = "ACDEFGHILKMNPQRSTVWY"

        Returns:
            A list of values as a feature vector

        """
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        order = kwargs.pop('order', "ACDEFGHILKMNPQRSTVWY")
        verbose = kwargs.pop('verbose', False)
        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        seq = std_aa_only(visitee.get_sequence())

        # initialize a dict
        # occurrence = {k:0 for k in order}
        occurrence = np.zeros(len(order))
        for aa in seq:
            try:
                occurrence[order.index(aa)] += 1
            except:
                pass

        return np.sqrt(occurrence).tolist()


class ZpCurveDescriptor(prDescriptor):
    """
    **Description**:

    20 kinds of amino acids are divided into 4 groups.

        Apolar-->A

        Polar -->P

        Positively charged --> Cp

        Negatively charged --> Cn

    .. math::
        X_{n}=( A_{n}+P_{n} )  -( Cp_{n}+Cn_{n} );

        Y_{n}=( A_{n}+Cp_{n}) -( P_{n}+Cn_{n} );

        Z_{n}=( A_{n}+Cn_{n}) -( P_{n}+Cp_{n} );

    where n:1...L (L is the number of residues of the sequence)

    .. note::
        Feng, Z. P. and Zhang, C. T. (2002).
        A graphic representation of protein primary structure
        and its application in predicting subcellular
        locations of prokaryotic proteins.
        **Int. J. Biochem. Cell Biol.** 34, 298-307.

    """
    def generate(self, visitee, *args, **kwargs):
        """Generate feature of the sequence

        Args:

            visitee (class): the instance of class Sequence

        Kwargs:

            order (str): amino acids in 1-letter code will be counted.

            category (dict):  grouping 20 standard amino acids

                'key' is the name of a category and must be the four strings,
                        'Nonpolar', 'Polar', 'Positively_charged',
                        'Negatively_charged'.

                'value' is a string containing 1-letter amino acids

                default grouping::

                    category =  {

                        'Nonpolar': "AVLIPMFWCG",

                        'Polar' : 'STYNQ',

                        'Positively_charged':'DE',

                        'Negatively_charged':'KRH'
                    }

            axis (int) : default is 0

                0: return a 1d list from a 3d array.
                   i.e., [x1, y1, z1, x2, y2, z2, ... , xn, yn, zn]

                1: return a 1d list of x only

                2: return a 1d list of y only

                3: return a 1d list of z only

        Returns:

            A list of values as a feature vector.
            The returns is controlled by axis option
            see :param: axis

        """
        std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
        order = kwargs.pop('order', "ACDEFGHILKMNPQRSTVWY")
        verbose = kwargs.pop('verbose', False)
        default_keys = ['Nonpolar', 'Polar',
                        'Positively_charged',
                        'Negatively_charged']
        category = kwargs.pop(
            'category', {'Nonpolar':"AVLIPMFWCG", 'Polar' : 'STYNQ',
                    'Positively_charged':'DE', 'Negatively_charged':'KRH'}
                )
        axis = kwargs.pop('axis', 0)
        if type(category) is not DictType:
            raise ValueError("Inappropriate argument for category")
        elif category.keys() != default_keys:
            raise ValueError(
                "Keys of category is not correct." +
                "It must contain 'Nonpolar', 'Polar'," +
                "'Positively_charged' and 'Negatively_charged'.")
        if verbose:
            print("{} Visited by {}".format(
                    visitee.__class__.__name__,
                    self.__class__.__name__))
        seq = std_aa_only(visitee.get_sequence())
        A = 0
        P = 0
        C_p = 0
        C_n = 0
        f_val = []
        for aa in seq:
            if category['Nonpolar'].find(aa):
                A += 1
            elif category['Polar'].find(aa):
                P += 1
            elif category['Positively_charged'].find(aa):
                C_p += 1
            elif category['Negatively_charged'].find(aa):
                C_n += 1
            else:
                if verbose:
                    print(
                        "{} is not grouped in the given category".format(aa))
            if axis == 0:
                x = (A + P) - (C_p + C_n)
                y = (A + C_p) - (P + C_n)
                z = (A + C_n) - (P + C_p)
                f_val.extend((x, y, z))
            elif axis == 1:
                x = (A + P) - (C_p + C_n)
                f_val.append(x)
            elif axis == 2:
                y = (A + C_p) - (P + C_n)
                f_val.append(y)
            elif axis == 3:
                z = (A + C_n) - (P + C_p)
                f_val.append(z)
        return f_val
