# -*- coding: utf-8 -*-

import numpy as np
import os.path

from io import FileIO
from types import *

from bilab.io import AbstractDataParser
from bilab.ml import Kernel, Model
from bilab.optimization import AbstractSolver

__all__ = ['SVM']

class SVM(object):
    """ SVM class

    implement a facade pattern

    """

    def __init__(self, parser=None,
                       solver=None,
                       outfile=None):
        """Initialization.

        Args:
            parser (AbstractDataParser):
                    load foramtted data
            kernel (object):
                    kernel methods for mapping samples to the higher dimensional
                    space.
            solver (object):
                    solve QP problem
        Kwargs:

        """
        super(SVM, self).__init__()
#        self.parser = self.__initial2obj(parser)
#        self.kernel = self.__initial2obj(kernel)
#        self.solver = self.__initial2obj(solver)
        self.model = None
        self.ready_for_training = True # Ready
        self.ready_for_classifying = False #
        self.train_file_location = None

        if outfile is None:
            self.outfile = "svm.model"
        else:
            self.outfile = outfile

        if isinstance(parser, AbstractDataParser):
            self.parser = parser
        else:
            self.parser = parser()

#        if isinstance(kernel, Kernel):
#            self.kernel = kernel
#        else:
#            self.kernel = kernel()

        if isinstance(solver, AbstractSolver):
            self.solver = solver
        else:
            self.solver = solver()

        if not (self.parser or self.solver):
            self.ready_for_training = False


#    def __initial2obj(self, inType):
#        #print("inType is {}".format(type(inType)))
#        clsname = inType.__class__.__name__
#        print("Create an instance of {}".format(clsname))
#        if type(inType) is ClassType:
#        if isinstance(inType, object):
#            print("Create instance of {} from {}"
#                 .format(inType.__class__.__name__, type(inType)))
#            return inType()
#        if isinstance(inType, object):
#            print("{} is an instance already".format(inType.__class__.__name__))
#            return inType

    def train(self, sample, verbose=True):
        """ Class method train

        training svm
        
        Args:
            sample (str): a file name of training samples
        
        """
        if self.ready_for_training:
            #print("parser is {}".format(type(self.parser)))
            # loading data

            self.train_file_location = os.path.abspath(sample)
            if verbose:
                print("Reading training dataset:{}".format(self.train_file_location))
            labels, features = self.parser(sample)
            if self.outfile is None:
                self.outfile = sample + 'model'
            if verbose:
                print("Create the kernel matrix")
            # create Qmatrix
            #Qmatrix = self.kernel(np.asmatrix(features))
            if verbose:
                print("Solve the QP problem")
            # training
            self.solver.solve(labels, features)
            #self.solver.show()
            if verbose:
                print("Training Finished")
        else:
            if self.parser is None:
                print("SVM: parser is not set yet.")
            if self.solver is None:
                print("SVM: solver is not set yet.")

        self.model = self.solver.getModel()
        #print("Model:{}".format(self.model.__str__()))

    def classify(self, sample):
        """ Class method classify

        Returns:

            descision_vals(nx1 array):

        """
        if self.model:
            descision_vals = self.model.classify(sample)
        return descision_vals

    def load(self, model):
        """ Load a mode from a file """
        status = False
        self.model = Model()
        if os.path.isfile(model):
            fhandle = FileIO(model, 'r')
            for line in fhandle:
                line = line.strip()
                if line.startswith("#"):
                    continue
                #Model
                if line.startswith("Model"):
                    model_line = line.split(":")
                    print("model name: {}".format(model_line[1]))

                #kernel
                if line.startswith("Kernel"):
                    kernel_line = line.split(":")
                    prop_name = kernel_line[0].split('.')[1]
                    prop_val = kernel_line[1]
                    print("kernel {}: {}".format(prop_name, prop_val))

        return status

    def save(self, filename=None):
        """ Save model to file
        """
        if filename is not None:
            self.outfile = filename

        with open(self.outfile,'w+') as fhandle:
            # save kernel parameters
            fhandle.write("# train dataset: {}\n".format(self.train_file_location))
            #fhandle.write(self.solver.__str__())
            fhandle.write(self.model.__str__())
