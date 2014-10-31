# -*- coding: utf-8 -*-

__all__ = [ 'DivergenceError', 'AbstractSolver']

class DivergenceError(Exception):
    """Exception raised if the entropy dual has no finite minimum.
    """
    def __init__(self, message):
        self.message = message
        Exception.__init__(self)

    def __str__(self):
        return repr(self.message)

class AbstractSolver(object):
    """ Abstract class for QP solvers """
    def __init__(self):
        """ Initialization """
        super(AbstractSolver, self).__init__()

    def solve(self, labels, samples, **kwargs):
        """ Class method

        Args: at least two arguments
            labels (n-dim vector): labels of n samples (1 or -1)
            samples (n x m matrix): n samples with m attributes.

        """
        # force subclass to implement
        raise NotImplementedError("{} method hasn't been implemented yet."
            .format('solve'))

    def getModel(self):
        # force subclass to implement
        raise NotImplementedError("{} method hasn't been implemented yet."
            .format('getModel'))
