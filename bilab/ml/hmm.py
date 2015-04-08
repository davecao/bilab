# -*- coding: utf-8 -*-

import numpy as np

__all__ = ['HMM']

# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    line = method.func_code.co_firstlineno
    filename = method.func_code.co_filename
    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError('Abstract method %s at File "%s", line %s'
            'should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented

class ClassRegistry(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(ClassRegistry, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, 'registry'):
            cls.registry = set()
        cls.registry.add(cls)
        cls.registry -= set(bases) #Remove base classes
    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls.registry)

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])


class HMM(object):
    """Base class for various HMM

    """
    __metaclass__ = ClassRegistry

    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           ConcreteHMM (class):

        Kwargs:
        """
        super(HMM, self).__init__()
        Concrete_cls = kwargs.pop("Concrete", None)
        if Concrete_cls is None:
            return
        self.__Concrete = self.__create(Concrete_cls, *args, **kwargs)

    def __create(self, clsname, *args, **kwargs):
        """ Create an instance for a given class in str or name
        """
        obj = None
        for cls in self.registry:
            if clsname == cls.__name__ or clsname == cls:
                obj = cls(*args)
        if obj:
            # Initialize object
            obj.__init__(*args, **kwargs)
        else:
            print("Unknown class {}".format(clsname))
        return obj

    def __str__(self):
        func_str_ = self.__getattr__(self.__Concrete, '__str__')
        return func_str_()

    def __repr__(self):
        #return self.__IO.__repr__()
        func_str_ = self.__getattr__(self.__Concrete, '__repr__')
        return func_str_()

    def __call__(self, *args, **kwargs):
        return self.__Concrete(*args, **kwargs)

    def __getattr__(self, attr):
        """ get delegation to the object """
        #return getattr(obj, attr)
        try:
            return self.__Concrete.__getattribute__(attr)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr))

    def __isstr(self, s):
        try:
            return isinstance(s, basestring)
        except NameError:
            return isinstance(s, str)

class DiscreteHMM(HMM):
    """
    A Discrete HMM - The most basic implementation of a Hidden Markov Model,
    where each hidden state uses a discrete probability distribution for
    the physical observations.
    
    Args:
        hiddenStates -  the alphabets of hidden states (N), list or array
                         e.g. hiddenStates = ('Rainy', 'Sunny')
        observedStates - the alphabets of observable symbols (M)
                         e.g. observedStates = ('walk', 'shop', 'clean')

    Kwargs:
        transitionProbMatrix - hidden states transition probability matrix 
                               ([NxN] numpy array)
                                   Rainy  Sunny
                            Rainy  0.7     0.3
                            Sunny  0.4     0.6
        emissionProbMatrix -  denoting each state's distribution ([NxM] numpy array)
                                  walk     shop     clean
                        Rainy      0.1      0.4      0.5
                        Sunny      0.6      0.3      0.1
        initStatesProb- initial state's Probability mass function 
                        ([N] numpy array).
                       (Rainy, Sunny) = np.array([0.6, 0.4])

    Additional attributes:
      precision - a numpy element size denoting the precision
      verbose - a flag for printing progress information, mainly when learning
    """

    def __init__(self, hiddenStates, observerStates, 
                transitionProbMatrix=None, 
                emissionProbMatrix=None, 
                initStatesProb=None, 
                initPMFtype=None,
                precision=np.double, 
                verbose=False):
        super(DiscreteHMM, self).__init__():

        # Number of hidden states
        self.numOfStates = len(hiddenStates)
        # Number of observations
        self.numOfObservations = len(observedStates)

        self.transitionProbMatrix = transitionProbMatrix
        self.emissionProbMatrix = emissionProbMatrix
        self.initStatesProb = initStatesProb
        self.precision = precision
        self.verbose = verbose
        self.setInitialProb()

    def setInitialProb(self):
        if self.initPMFtype == 'uniform':
            self.stateDistribution = 
                np.ones( 
                    (self.numOfHiddenStates), dtype=self.precision) * 
                    (1.0/self.numOfHiddenStates)
            self.hiddenTransProbMatrix = np.ones( 
                (self.numOfHiddenStates, self.numOfHiddenStates), 
                    dtype=self.precision)*(1.0/self.numOfHiddenStates)
            self.stateDistribution = np.ones( 
                (self.numOfHiddenStates,self.numOfObsStates), 
                    dtype=self.precision)*(1.0/self.numOfObsStates)


class HMM_Interface(object):
    """ 
    Interface of all subclasses of HMM
    """
    def __init__(self, n, m, precision=np.double, verbose=False):
        self.n = n
        self.m = m
        self.precision = precision
        self.verbose = verbose

    def _eta(self, t, T):
        """
        Governs how each sample in the time series should be weighed.
        This is the default case where each sample has the same weigh, 
        i.e: this is a 'normal' HMM.
        """
        return 1

    def _alpha(self, observations):
        """
        Calculates 'alpha' the forward variable.
    
        The alpha variable is a numpy array indexed by time, then state (TxN).
        alpha[t][i] = the probability of being in state 'i' after observing the 
        first t symbols.
        """
        alpha = numpy.zeros((len(observations),self.n),dtype=self.precision)

    def _beta(self, observations):
        """
        Calculates 'beta' the backward variable.
        
        The beta variable is a numpy array indexed by time, then state (TxN).
        beta[t][i] = the probability of being in state 'i' and then observing the
        symbols from t+1 to the end (T).
        """
    def _xi(self, observations, alpha=None, beta=None):
        """
        Calculates 'xi', a joint probability from the 'alpha' and 'beta' variables.
        
        The xi variable is a numpy array indexed by time, state, and state (TxNxN).
        xi[t][i][j] = the probability of being in state 'i' at time 't', and 'j' at
        time 't+1' given the entire observation sequence.
        """
    def _gamma(self, xi, seqlen):
        """
        Calculates 'gamma' from xi.
        
        Gamma is a (TxN) numpy array, where gamma[t][i] = the probability of being
        in state 'i' at time 't' given the full observation sequence.
        """
    def _updatemodel(self,new_model):
        '''
        Replaces the current model parameters with the new ones.
        '''
        self.pi = new_model['pi']
        self.A = new_model['A']

    def _calcstats(self,observations):
        '''
        Calculates required statistics of the current model, as part
        of the Baum-Welch 'E' step.
        
        Deriving classes should override (extend) this method to include
        any additional computations their model requires.
        
        Returns 'stat's, a dictionary containing required statistics.
        '''
        stats = {}
        
        stats['alpha'] = self._calcalpha(observations)
        stats['beta'] = self._calcbeta(observations)
        stats['xi'] = self._calcxi(observations,stats['alpha'],stats['beta'])
        stats['gamma'] = self._calcgamma(stats['xi'],len(observations))
        
        return stats

    def _reestimate(self,stats,observations):
        '''
        Performs the 'M' step of the Baum-Welch algorithm.
        
        Deriving classes should override (extend) this method to include
        any additional computations their model requires.
        
        Returns 'new_model', a dictionary containing the new maximized
        model's parameters.
        '''        
        new_model = {}
        
        # new init vector is set to the frequency of being in each step at t=0 
        new_model['pi'] = stats['gamma'][0]
        new_model['A'] = self._reestimateA(observations,stats['xi'],stats['gamma'])
        
        return new_model

    def _baumwelch(self,observations):
        """
        An EM(expectation-modification) algorithm devised by Baum-Welch. Finds a local maximum
        that outputs the model that produces the highest probability, given a set of observations.
        
        Returns the new maximized model parameters
        """
        # inner functions

        # E step - calculate statistics
        stats = self._calcstats(observations)

        # M step
        return self._reestimate(stats,observations)

    def _forwardbackward(self, observations, cache=False):
        """
        Forward-Backward procedure is used to efficiently calculate the probability of the observation, given the model - P(O|model)
        alpha_t(x) = P(O1...Ot,qt=Sx|model) - The probability of state x and the observation up to time t, given the model.
        
        The returned value is the log of the probability, i.e: the log likehood model, give the observation - logL(model|O).
        
        In the discrete case, the value returned should be negative, since we are taking the log of actual (discrete)
        probabilities. In the continuous case, we are using PDFs which aren't normalized into actual probabilities,
        so the value could be positive.
        """

    def _viterbi(self, observations):
        """
         Find the best state sequence (path) using viterbi algorithm - a method of dynamic programming,
        very similar to the forward-backward algorithm, with the added step of maximization and eventual
        backtracing.
        
        delta[t][i] = max(P[q1..qt=i,O1...Ot|model] - the path ending in Si and until time t,
        that generates the highest probability.
        
        psi[t][i] = argmax(delta[t-1][i]*aij) - the index of the maximizing state in time (t-1), 
        i.e: the previous state.

        """
    def decode(self, observations):
        '''
        Find the best state sequence (path), given the model and an observation. i.e: max(P(Q|O,model)).
        
        This method is usually used to predict the next state after training. 
        '''        
        # use Viterbi's algorithm. It is possible to add additional algorithms in the future.
        return self._viterbi(observations)

    def train(self, observations, iterations=1,epsilon=0.0001,thres=-0.001):
        """ 
        Updates the HMMs parameters given a new set of observed sequences.
        
        observations can either be a single (1D) array of observed symbols, or when using
        a continuous HMM, a 2D array (matrix), where each row denotes a multivariate
        time sample (multiple features).
        
        Training is repeated 'iterations' times, or until log likelihood of the model
        increases by less than 'epsilon'.
        
        'thres' denotes the algorithms sensitivity to the log likelihood decreasing
        from one iteration to the other.
        """
        # Mapping
        for i in xrange(iterations):
            # EM algorithm model
            model = self._baumwelch(observations)
            # log likelihood of the model
            prob_ = self._forwardbackward(observations, cache=True)
            # update the model with the new estimation
            self.pi = new_model['pi']
            self.A  = new_model['A']
            # log liklihood of the new model
            prob_new = self._forwardbackward(observations, cache=False)
            # verbose
            if (self.verbose):
                print("iter: {}, L(model|O) = {}, L(model_new|O) = {}, "
                    " converging = {}".format ( prob_new-prob_old > thres )
            if ( abs(prob_new-prob_old) < epsilon ):
                # converged
                break
