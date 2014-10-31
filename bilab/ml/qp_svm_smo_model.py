# -*- coding: utf-8 -*-

import os
import numpy as np

from bilab.ml import Kernel
from bilab.ml import AbstractModel

__all__ = ['QP_svm_smo_model']

class QP_svm_smo_model(AbstractModel):
    """ Store a svm model
    """
    def __init__(self, kernel=None,
                       threshold = np.inf,
                       categories=None,
                       sv = [],
                       sv_coef = [],
                       total = None):
        super(QP_svm_smo_model, self).__init__()
        self.name = self.__class__.__name__
        self.kernel = kernel
        self.categories = categories
        self.sv = sv
        self.sv_coef = sv_coef
        self.b = threshold
        self.total = total
        self.kernel = kernel
        self.ready =False

    def __str__(self):
        txt = "Model.name:" + self.__class__.__name__ + '\n'
        txt += self.kernel.__str__()
        txt += "SV.total:{}\n".format(self.total)
        txt += "SV.category:"
        for key in self.categories:
            txt += " {}:{}".format(key, len(self.categories[key]))
        txt += "\n"
        txt += "SV.b:{}\n".format(self.b)
        txt += "SV.sv:\n"

        for i in range(len(self.sv_coef)):
            #get coefficient
            coeff = self.sv_coef[i]
            features = self.sv[i]
            # loop over sv
            feature = ''
            for inx, x in np.ndenumerate(features):
                feature += "{}:{:.6f} ".format(inx[1] + 1, x)
            txt += "{:.3f} {}\n".format(coeff, feature)
        return txt

    def __repr__(self):
        return self.__str__()


    def check_status(self):
        getReady = self.ready
        if not getReady:
            if self.sv and self.sv_coef and (not np.isinf(self.b)) and self.kernel:
                self.ready = True
                getReady = True
        return getReady

    def load(self, model):
        """ Load a model from a saved file

        Args:
            model (str): a file name

        Returns:
            status (bool): True if successfully load a model
        """
        status = False
        if os.path.isfile(model):
            with open(model, 'r') as fhandle:
                for line in fhandle:
                    if line.startwith("#"):
                        continue
                    #kernel
                    if line.startwith("Kernel"):
                        kernel_line = line.split(":")

            status = True
        return status

    def classify(self, sample, decision_func=None, verbose=True):
        """ Prediction

        Args:
            sample(nxm np array): an input samples for prdiction

        Returns:
            descision value(float)

        """
        if self.check_status():
            if verbose:
                print("{} is ready for prediction".format(self.name))

