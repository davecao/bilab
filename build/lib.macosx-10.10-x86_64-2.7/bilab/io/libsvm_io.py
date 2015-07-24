# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pprint

from types import *

__all__ = ['AbstractDataParser',
           'LIBSVM_datafile_parser']

class AbstractDataParser(object):
    """ base class for parsing file
    """

    def __init__(self):
        """Initialization.

        """
        super(AbstractDataParser, self).__init__()


    def fromFile(self, filename, *args, **kwargs):
        """class method: fromFile

           Read data from a file

        Args:
           filename (str):

        Kwargs:

        Returns:
            a numpy matrix: samples in row, features in column; the first column
            is labels (+1 or -1).
            i.e.,

             [
               [label, feature1, ... , featureN],

               [label, feature1, ... , featureN],
               ...

             ]

        Raises:
            IOError, ValueError

        Usage:

        >>> print fromFile(filename)

        """
        # force subclass to implement
        raise NotImplementedError("{} method hasn't been implemented yet.")

    def fromString(self, str, *args, **kwargs):
        """class method: fromString

        Read data from a String stream

        Args:
            str (str):  a string stream

        Kwargs:

        Returns:

            a numpy matrix: samples in row, features in column; the first column
            is labels (+1 or -1).

            i.e.,
             [

               [label, feature1, ... , featureN],

               [label, feature1, ... , featureN],

               ...

             ]

        Raise:
            ValueError

        """
        # Do not force subclass to implement
        pass

class LIBSVM_datafile_parser(AbstractDataParser):
    """ File parser for libsvm formatted data file

    """
    def __init__(self):
        """ Initialization

        """
        super(LIBSVM_datafile_parser, self).__init__()

    def __call__(self, sample):
        if os.path.isfile(sample):
            labels, features = self.fromFile(sample)
        elif type(sample) is StringType:
            labels, features = self.fromString(sample)
        return (labels, features)

    def _dictRecord(self, filename):
        """Read data from a file

        Args:
            filename (str): a string of a file name

        Returns:
            a tuple of (nrows, ncols, labels, features):
                   labels and features are in same size.

            nrows: i.e., the number of samples.

            ncols: i.e., the number of featurs in a sample.

            label: an array of +1 or -1
                    i.e.,
                    [1, 1, ... , -1, -1, ... , -1]

            features: an array of dictionary of features:
                      the index is taken as the key
                    e.g.,
                      [

                        {1:1.0, ... , n:10.0},

                        {1:2.0, ... , n:11.0},

                        ...

                      ]

        """
        if os.path.isfile(filename):
            # load into a dictionary
            labels = np.zeros(0, dtype=np.int8)
            features = []
            length = 0
            max_index_of_feature = 0
            with open(filename) as fhandle:
                for line in fhandle:
                    length += 1
                    line = line.split(None, 1)

                    if len(line) == 1:
                        print("Ignore line:{} in {} for no features found"
                            .format(length, filename))
                        continue

                    label, data = line
                    feature_inx = [ x.split(":") for x in data.split()]
                    indexes = [ int(inx[0]) for inx in feature_inx]
                    if max_index_of_feature == 0:
                        max_index_of_feature = max(indexes)
                    elif max_index_of_feature != max(indexes):
                        #The lengths of feature vectors are different
                        raise ValueError("Different length of the " +
                            "feature vector found in line {}".format(length))

                    # in one line
                    #repeated = True if not \
                    #        [ x for x in indexes if indexes.count(x) > 1] \
                    #        else False

                    if not [ x for x in indexes if indexes.count(x) > 1]:
                        # No repeated indexes in the feature line
                        f_val = {}
                        for index, value in feature_inx:
                            f_val[int(index)] = float(value)
                        labels = np.append(labels, int(label))
                        features.append(f_val)
                    else:
                        # The feature have duplicated indexes
                        raise ValueError(
                            'Duplicated indexes found in {}: line {}'
                            .format(filename, length))

        else:
            raise ValueError('{} is not a filename or does not look'
            'like a valid string'.format(filename))
        nrows = len(labels)
        ncols = max_index_of_feature
        return (nrows, ncols, labels, features)

    def fromFile(self, filename, *args, **kwargs):
        """Parse libsvm formatted data file

        Args:
           filename (str): a file or a string

        Kwargs:
            verbose (bool) :  show detail info

        Returns:
            a tuple of (labels, featurs):

            labels:  a numpy array of (+1 or -1).

            features: a numpy dense matrix of samples in row, features in column;
            
            format::

                [
                    [feature1, ... , featureN],

                    [feature1, ... , featureN],
                    
                    ...
                ]

        Raises:

           IOError, ValueError

        Usage:

        >>> print parse_libsvm_formatted_data(filename)

        """
        verbose = kwargs.pop('verbose', False)
        nsamples, nfeatures, tags, features = self._dictRecord(filename)
        labels = tags
        data = np.zeros((nsamples, nfeatures))

        if verbose:
            positives = len(labels[labels>0])
            negatives = len(labels[labels<0])
            print("{} positive samples and "
                  " {} negative samples with {} features in each"
                .format(positives, negatives, nfeatures))

        for row, feature in enumerate(features):
            for col in feature:
                #The numpy array is zero-based array
                data[row, col-1] = feature[col]

        return (labels, data)
