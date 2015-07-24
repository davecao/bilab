# -*- coding: utf-8 -*-
""" This module contains solvers for the quadratic programing(QP) problem.

    Sequential minimal optimization (smo) :
       it is used in training support vector machines, invented by
       John platt in 1998 at Microsoft Research.

"""
import sys
import time
import numpy as np

from bilab.ml import Kernel, Model, QP_svm_smo_model
from bilab.optimization import AbstractSolver

from subprocess import Popen, PIPE, check_call
from threading import Thread

__all__ = ['QP_svm_smo_solver']

class QP_svm_smo_solver(AbstractSolver):
    """ A QP solver implements the SMO algorithm
        It is an iterative algorithm. Although it is guaranteed to converge,
        the selection of the pair of multipliers affects the speed of
        convergence.

    Reference:
        Platt, John (1998), Sequential Minimal Optimization:
        A Fast Algorithm for Training Support Vector Machines

        Chang, Chih-Chung; Lin, Chih-Jen (2011). "LIBSVM: A library for
        support vector machines".
        ACM Transactions on Intelligent Systems and Technology 2 (3).
    """
    def __init__(self,
                kernel,
                qm = None,
                weights = None,
                epsilon = 0.001,
                tolerance=0.001,
                cost = 1.0,
                probability_estimates = False,
                parallel=False):
        """Initialization.

        Args:

        Kwargs:
            epsilon (float) : evaluate agrange multipliers default is 0.001
            tolerance (float) : default is 0.001 (for KKT)
            cost (float): parameter C for classification SVM (C-SVM),
                     default is 1.0.

            probability_estimates (bool) : estimate probability if True.
                     default is False.
            parallel (bool) : using multiple processors if True.
                     default is False.

        """
        super(QP_svm_smo_solver, self).__init__()
        self.kernel = kernel
        self.epsilon = float(epsilon)
        self.tolerance = float(tolerance)
        self.cost = float(cost)
        self.probability_estimates = bool(probability_estimates)
        self.parallel = bool(parallel)
        #self.sample_size = Qmatrix.shape[0]
        self.sv_index = [] #  indexes of support vectors
        self.convergence = False # check the convergence
        self.b = 0.0 #  threshold
        self.labels = None
        self.samples = None
        self.sv = [] # array of support vectors
                    # each of SVs constains three elements
                    # (alpha, label, features)
                    # where features is a np 1d array
        self.sv_status = {} # the number of support vectors in two categories
                            # e.g., {'1':70, '-1':80}
        # check kernel
        if isinstance(kernel, Kernel):
            self.kernel = kernel
        else:
            self.kernel = kernel()

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        txt = self.kernel.__str__()
        txt += "Opt.epsilon:{}\n".format(self.epsilon)
        txt += "Opt.tolerance:{}\n".format(self.tolerance)
        txt += "Opt.cost:{}\n".format(self.cost)
        txt += "SV.total:{}\n".format(len(self.sv_index))
        txt += "SV.category:"
        for key in self.sv_status:
            txt += " {}:{}".format(key, len(self.sv_status[key]))
        txt += "\n"
        txt += "SV.b:{}\n".format(self.b)
        txt += "SV.sv:\n"
        for alpha,label, features in self.sv:
            # loop over sv
            feature = ''
            for inx, x in np.ndenumerate(features):
                feature += "{}:{:.6f} ".format(inx[1] + 1, x)
            txt += "{:.3f} {}\n".format(alpha*label, feature)
        return txt

    def getModel(self):
        """ Save to model """
        total_sv = len(self.sv_index)
        svm_model = QP_svm_smo_model(kernel=self.kernel,
                                threshold = self.b,
                                categories=self.sv_status,
                                total = total_sv)
        for alpha,label, features in self.sv:
            svm_model.sv_coef.append(alpha * label)
            svm_model.sv.append(features)
        return Model(svm_model)

    def solve(self, labels, samples, *args,  **kwargs):
        """Class method

            Implementation of QP solver

        Args:
            labels (n-dim vector): labels of n samples (1 or -1)

            samples (n x m matrix): n samples with m attributes.

        Kwargs:
            weights (nx1 matrix): weights for each sample used for the
                    unbalanced dataset. default is 1.
            qm (n x n matrix): cached kernel matrix, where
                                n is the number of training samples.

        """
        self.weights = kwargs.pop('weights', None)
        self.Qmatrix = kwargs.pop('qm', None)
        if self.Qmatrix is None:
            self.Qmatrix = self.kernel(np.asmatrix(samples))

        self.__check_input(labels, samples)

        if self.sample_size > sys.maxsize/100:
            maxiter = sys.maxsize
        else:
            maxiter = 100 * self.sample_size
        self.maxiter = np.max([10000000, maxiter])
#        p = ProgressBar(1,"progress", interval=5)

        if self.parallel:
            self.__qp_smo_parallel_solver()
        else:
#            p.start()
            self.__qp_smo_platt()
#            p.join()
        # create a model


    def __check_input(self, labels, samples):

        """ Check the parameters
        """

        if isinstance(labels, (np.ndarray, np.generic)) and labels.ndim == 1:
            self.labels = labels
        else:
            raise TypeError("labels: Inappropriate argument type for {}"
                    .format(self.__class__.__name__))

        if isinstance(self.Qmatrix, (np.ndarray, np.generic)):
            self.Qmatrix = np.asmatrix(self.Qmatrix)
        else:
            raise TypeError("Qmatrix: Inappropriate argument type for {}"
                    .format(self.__class__.__name__))

        if isinstance(samples, (np.ndarray, np.generic)):
            self.samples = np.asmatrix(samples)
            self.sample_size = self.samples.shape[0]
        else:
            raise TypeError("samples: Inappropriate argument type for {}"
                    .format(self.__class__.__name__))

        #labels size
        if labels.size != self.sample_size:
            raise ValueError("The size of labels and samples are not same.")

        if self.weights is None:
            self.weights = np.ones(self.Qmatrix.shape[0])

        elif isinstance(weights, list):
            self.weights = np.array(weights)

        elif isinstance(weights, (np.ndarray, np.generic)):
            if weights.shape[0] == self.Qmatrix.shape[0]:
                self.weights = weights
            else:
                raise ValueError("The size of weights "
                    "is not same as the Qmatrix")
        # initialize alphas
        self.alpha = np.zeros(self.sample_size)
        # initialize error cache
        self.error_cache = np.zeros(self.sample_size) # error
        uniq_labels = np.unique(self.labels)
        for l in uniq_labels:
            self.sv_status[l] = []

    def __qp_smo_platt(self):
        """
            Implementation based on the John Platt's pseudo code

            The two inner functions: takeStep and examineAll

        """
        target = self.labels # 1 or -1 of samples' label
        # inner function as closure, so take b and error_cache as class members
        #b = self.b # threshold
        #error_cache = np.zeros(self.sample_size)

        def check_lagrange_multiplier(alpha):
            """ Check condition for alpha
                returns:
                True if alpha > eps and alpha < C-eps, otherwise return False.
            """
            state = False
            bound = self.cost - self.epsilon
            if alpha > self.epsilon and alpha < bound:
                state = True
            return state

        def svm_output(i):
            """ calculate svm output for sample i
                u_{i} = \sum_{j=1}^{N}(y_{j}* \alpha_{j} * K(x_{j},x_{i}) - b
                    equ(10)
            """
            val = 0.0
            for j in range(self.sample_size):
                if self.alpha[j] == 0.0:
                    continue
                val += target[j] * self.alpha[j] * self.Qmatrix[j, i]
            val -= self.b
            return val

        def takeStep(i1, i2):
            """ closure method __takeStep

            Args:
                i1: the index of a sample
                i2: the index of a sample

            """
            if i1 == i2:
                return 0
            b = self.b

            alpha1 = self.alpha[i1]
            alpha2 = self.alpha[i2]

            y1 = target[i1]
            y2 = target[i2]
            # E1 = SVMoutput on point[i1] - y1 (check in error cache)
            E1 = svm_output(i1) - y1
            E2 = svm_output(i2) - y2
            s = y1 * y2
            # compute bound: L and H
            # if target y1 != y2 equ(13)
            #   L = max(0, alpha2-alpha1), H = min (C, C + alpha2 - aplha1)
            # else               equ(14)
            #   L = max(0,alpha2 + aplha1 -C), H = min(C, alpha2 + aplha1)
            if y1 != y2:
                L = np.max([0, alpha2 - alpha1])
                H = np.min([self.cost, self.cost + alpha2 - alpha1])
            else:
                L = np.max([0, alpha1 + alpha2 - self.cost])
                H = np.min([self.cost, alpha1 + alpha2])
            if L == H:
                return 0

            k11 = self.Qmatrix[i1, i1]
            k12 = self.Qmatrix[i1, i2]
            k22 = self.Qmatrix[i2, i2]
            eta = k11 + k22 - 2*k12
            if eta > 0:
                # update alpha
                # \alpha_{2}^{new} = \alpha_{2} + \frac{y_{2}*(E1-E2)}{\eta}
                # where \E_{i} = u_{i} -y_{i}
                a2 = alpha2 + y2 * (E1-E2)/eta
                # constrained minimum is found by clipping the unconstrained
                # minimum to the ends of the line segment:
                if a2 < L:
                    a2 = L
                elif a2 > H:
                    a2 = H
                #print("takeStep: {}, {} eta={:.3} E1={:.3} E2={:.3} L={:.3} H={:.3} a2={:.3}"
                #    .format(i1, i2, eta, E1, E2, L, H, a2))
            else:
                # under unusual circumstances, eta will not be positive.
                # A negative eta will occur if the kernel K does not obey
                #    Mercer's condition, which can cause the objective function
                #    to become indefinite.
                # A zero eta: caused by more than one training example has
                #    the same input vector x.
                # Lobj = objective function at a2 = L
                # Hobj = objective function at a2 = H
                # objective function:
                #   W(a1,a2) = a1 + a2
                #              - 0.5 * k11 * a1 * a1
                #              - 0.5 * k22 * a2 * a2
                #              - s   * k12 * a1 * a2
                #              -y1*a1*v1 -y2*a2*v2 + W_constant,
                #   where Kij = k(xi, xj)
                #    vi = \sum_{j=3}^{l} y_{j}*a_{j}^{old}*K_{ij}
                #       = f^{old}(xi) + b^{old}
                #         - y1*a1^{old}*k_{1i}
                #         - y2*a2^{old}*k_{2i}
                # a2 = L
                #a1 = L
                #a2 = alpha2 + y1*y2*(alpha1 - a1)
                a2 = L
                a1 = alpha1 + s * (alpha2 - a2)

                #f1 = y1*(E1 + b) - a1*k11 - s*a2*k12
                #f2 = y2*(E2 + b) - s*a1*k12 - a2*k22

                #L1 = a1 + s*(a2 - L)
                #Lobj = L1*f1 + L*f2 + 0.5*L1*L1*k11 + 0.5*L*L*k22 + s*L*L1*k12
                v1 = svm_output(i1) + b - y1*alpha1*k11 - y2*alpha2*k12
                v2 = svm_output(i2) + b - y1*alpha1*k12 - y2*alpha2*k22
                Lobj = a1 + a2 - 0.5*k11*a1*a1 - 0.5*k22*a2*a2 \
                       - s*k12*a1*a2 - y1*a1*v1 - y2*a2*v2
                #a1 = H
                #a2 = alpha2 + y1*y2*(alpha1 - a1)
                a2 = H
                a1 = alpha1 + s * (alpha2 - a2)
                #H1 = a1 + s*(a2 - H)
                #Hobj = H1*f1 + H*f2 + 0.5*H1*H1*k11 + 0.5*H*H*k22 + s*H*H1*k12
                v1 = svm_output(i1) + b - y1*alpha1*k11 - y2*alpha2*k12
                v2 = svm_output(i2) + b - y1*alpha1*k12 - y2*alpha2*k22
                Hobj = a1 + a2 - 0.5*k11*a1*a1 - 0.5*k22*a2*a2 \
                       - s*k12*a1*a2 - y1*a1*v1 - y2*a2*v2

                if Lobj > Hobj + self.epsilon:
                    a2 = L
                elif Lobj < Hobj - self.epsilon:
                    a2 = H
                else:
                    a2 = alpha2
                #print("takeStep: {}, {} eta={:.3} E1={:.3} E2={:.3} L={:.3} H={:.3} a2={:.3}"
                #    .format(i1, i2, eta, E1, E2, L, H, a2))
            #if a2 < self.tolerance:
            #    a2 = 0
            #elif a2 > (self.cost - self.tolerance):
            #    a2 = self.cost
            if abs(a2-alpha2) < self.epsilon * (a2 + alpha2 + self.epsilon):
                return 0
            a1 = alpha1 + s * (alpha2 - a2)
            old_b = b
            # Update threshold to reflect change in Lagrange multipliers
            #if a1 > self.epsilon and a1 < (self.cost - self.epsilon):
            if check_lagrange_multiplier(a1):
                # a1 is not at the bound
                b += E1 + y1*(a1 - alpha1)*k11 + y2*(a2-alpha2)*k12
                self.error_cache[i1] = E1
            #elif a2 > self.epsilon and a2 < (self.cost - self.epsilon):
            elif check_lagrange_multiplier(a2):
                b += E2 + y1*(a1 - alpha1)*k12 + y2*(a2-alpha2)*k22
                self.error_cache[i2] = E2
            else:
                b += 0.5*\
                    (E1 + y1*(a1 - alpha1)*k11 + y2*(a2-alpha2)*k12 + \
                     E2 + y1*(a1 - alpha1)*k12 + y2*(a2-alpha2)*k22)

            # update weight vectors to reflect change in a1 and a2, if SVM is linear
            # update error cache using new Largrange multipliers
            for m in range(self.sample_size):
                if m == i1 or m == i2:
                    continue
                #if self.alpha[m] > self.epsilon and \
                #   self.alpha[m] < (self.cost - self.epsilon):
                if check_lagrange_multiplier(self.alpha[m]):
                   self.error_cache[m] += y2*(a2-alpha2)*self.Qmatrix[i2, m] + \
                                     y1*(a1-alpha1)*self.Qmatrix[i1, m] + \
                                     old_b - b
            # store a1 in the alpha array
            self.alpha[i1] = a1
            # store a2 in the alpha array
            self.alpha[i2] = a2
            self.b = b

        def examineExample(i2):
            """ closure method examineExample

            Karush-Kuhn-Tucker conditions
                alpha_{i} = 0 , y_{i} * u_{i} >= 1
                0< alpha_{i} < C , y_{i} * u_{i} = 1
                alpha_{i} = C , y_{i} * u_{i} <= 1

            where u_{i} is the output of the SVM for the ith training example.

            Args:
                i2: the index of a sample

            """
            y2 = target[i2]
            alpha2 = self.alpha[i2]
            if check_lagrange_multiplier(alpha2):
                E2 = self.error_cache[i2]
            else:
                E2 = svm_output(i2) - y2
            r2 = E2 * y2

            max_Ej = 0.0
            max_j = -1

            # Check KKT condition
            #if ((r2 < -1 * self.tolerance) and alpha2 < self.cost ) or \
            #    (r2 > self.tolerance and alpha2 > 0):
            if ((r2 < -1 * self.tolerance) and alpha2 < self.cost -self.epsilon)\
                 or \
                (r2 > self.tolerance and alpha2 > self.epsilon):
                # i2 violates the KKT conditions, need to be optimized.
                #1. number of non-zero & non-C alpha > 1
                num = np.count_nonzero((self.alpha > 0) & \
                                       (self.alpha != self.cost))
                if num > 1:
                    #print("second choice heuristic")
                    offset = np.random.random_integers(0, self.sample_size - 1)
                    for j in range(self.sample_size):
                        pos = (j + offset) % self.sample_size
                        if check_lagrange_multiplier(self.alpha[pos]):
                            Ej = self.error_cache[pos]
                            E_diff = Ej - E2
                            if E_diff > max_j:
                                max_Ej = E_diff
                                max_j = pos
                    if max_j >= 0:
                        if takeStep(max_j, i2) == 1:
                            return 1
                #2. loop over all non-zero and non-C alpha,
                # starting at a random point
                #alpha_selection = (self.alpha > 0) & (self.alpha != self.cost)
                #working_set = [ idx for idx, val in enumerate(alpha_selection) \
                #                if val]
                offset = np.random.random_integers(0, self.sample_size - 1)
                for j in range(self.sample_size):
                    pos = (j + offset) % self.sample_size
                    if check_lagrange_multiplier(self.alpha[pos]):
                        if takeStep(pos, i2) == 1:
                            return 1
                #3.loop over all possible i1, starting at a random point
                offset = np.random.random_integers(0, self.sample_size - 1)
                for j in range(self.sample_size):
                        pos = (j + offset) % self.sample_size
                        if not check_lagrange_multiplier(self.alpha[pos]):
                            if takeStep(pos, i2) == 1:
                                return 1
            return 0

        """ Main process """
        numChanged = 0
        examineAll = 1
        loop = 0
        while numChanged > 0 or examineAll:
            numChanged = 0
            if loop > self.maxiter:
                print("Reach to max iteration {}".format(self.maxiter))
                break
            else:
                loop += 1
            if examineAll:
                # loop over all training examples
                for i in range(self.sample_size):
                    numChanged += examineExample(i)
            else:
                # loop over examples where alpha is not 0 and not C
                alpha_selection = (self.alpha > 0) & (self.alpha != self.cost)
                working_set = [ idx for idx, val in enumerate(alpha_selection) \
                                if val]
                print("There are {} alphas not 0 and not C".format(len(working_set)))
                for j in working_set:
                    numChanged += examineExample(j)
            if examineAll == 1:
                examineAll = 0
            elif numChanged == 0:
                examineAll = 1
            print("Iter {}: numChanged:{} examineAll:{}"
                .format(loop, numChanged, examineAll))
        #end of while
        if loop <= self.maxiter:
            self.convergence = True
        # store sv
        for i in range(self.sample_size):
            if self.alpha[i] != 0.0:
                self.sv_index.append(i)
                self.sv.append((self.alpha[i], self.labels[i], self.samples[i]))
                self.sv_status[self.labels[i]].append(i)
        # recalculate weights


    def __qp_smo_parallel_solver(self):
        """ A parallel QP solver implements the SMO algorithm

        Reference:
            Luca Zanni (2006). Parallel Software for Training Large Scale
            Support Vector Machines on Multiprocessor Systems
        """
        pass

    def is_symmetry(self, x):
        """ Check symmetry of a matrix """
        return np.all((x - x.T) == 0)

    def is_positive_definite(self, x):
        """ Check whether a matrix is positive definited or not
            by Cholesky decomposition
            (Much faster than using eigen vals)
        """
        try:
            np.linalg.cholesky(x)
        except  np.linalg.LinAlgError:
            raise np.linalg.LinAlgError("Matrix does not appear to be"
                " positive definited")

    def show(self):
        if self.convergence:
            print("#Optimization status: convergence.")
        print("NumofSVs:{}", len(self.sv_index))

    def show_globals(self):
        tmp = globals().copy()
        print("Globals varaibles:")
        for k,v in tmp.items():
            if not k.startswith('__') and k!='tmp' and k!='In' and k!='Out' \
                and not hasattr(v, '__call__'):
                print("{} : {}".format(k, v))

        # in a single line
        #[print(k,'  :  ',v,'\n') for k,v in tmp.items() \
        #    if not k.startswith('_') and k!='tmp' and k!='In' and k!='Out' \
        #    and not hasattr(v, '__call__')]

    def show_locals(self):
        tmp = locals().copy()
        print("locals varaibles:")
        for k,v in tmp.items():
            if not k.startswith('__') and k!='tmp' and k!='In' and k!='Out' \
                and not hasattr(v, '__call__'):
                print("{} : {}".format(k, v))

        # in a single line
        #[print(k,'  :  ',v,'\n') for k,v in tmp.items() \
        #    if not k.startswith('_') and k!='tmp' and k!='In' and k!='Out' \
        #    and not hasattr(v, '__call__')]
