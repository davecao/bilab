#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2015-12-16 16:45:54
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-02 19:29:46
# -------------------------------------------------------
#  Refer to the implementation of basinhopping in scipy
# Reference:
#
# Stefan Goedecker
# Minima hopping: An efficient search method for the global minimum of the
# potential energy surface of complex molecular systems.
# J. Chem. Phys., Vol. 120, 9911 (2004)

import numpy as np
import scipy.optimize as opt
import collections
from scipy.optimize import OptimizeResult


class Storage(object):
    """
    Class used to store the lowest energy structure
    """
    def __init__(self, x, f):
        self._add(x, f)

    def _add(self, x, f):
        self.x = np.copy(x)
        self.f = f

    def update(self, x, f):
        if f < self.f:
            self._add(x, f)
            return True
        else:
            return False

    def get_lowest(self):
        return self.x, self.f


class MinimaHoppingRunner(object):
    """This class implements the core of the minimahopping algorithm.
    x0 : ndarray
        The starting coordinates.
    minimizer : callable
        The local minimizer, with signature ``result = minimizer(x)``.
        The return value is an `optimize.OptimizeResult` object.
    step_taking : callable
        This function displaces the coordinates randomly.  Signature should
        be ``x_new = step_taking(x)``.  Note that `x` may be modified in-place.
    accept_tests : list of callables
        To each test is passed the kwargs `f_new`, `x_new`, `f_old` and
        `x_old`.  These tests will be used to judge whether or not to accept
        the step.  The acceptable return values are True, False, or ``"force
        accept"``.  If the latter, then this will override any other tests in
        order to accept the step.  This can be used, for example, to forcefully
        escape from a local minimum that ``minimahopping`` is trapped in.
    disp : bool, optional
        Display status messages.
    """
    def __init__(self, x0, minimizer, step_taking, accept_tests,
                 alpha1=0.9803921569, alpha2=1.02,
                 beta1=1.05, beta2=1.05, beta3=0.9523809524, disp=False):
        self.x = np.copy(x0)
        self.minimizer = minimizer
        self.step_taking = step_taking
        self.accept_tests = accept_tests
        self.disp = disp

        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.beta1 = beta1
        self.beta2 = beta2
        self.beta3 = beta3

        self.nstep = 0

        # initialize return object
        self.res = OptimizeResult()
        self.res.minimization_failures = 0

        # do initial minimization
        minres = minimizer(self.x)
        if not minres.success:
            self.res.minimization_failures += 1
            if self.disp:
                print("warning: minmahopping: local minimization failure")
        self.x = np.copy(minres.x)
        self.energy = minres.fun
        if self.disp:
            print("minimahopping step %d: f %g" % (self.nstep, self.energy))

        # initialize storage class
        self.storage = Storage(self.x, self.energy)

        if hasattr(minres, "nfev"):
            self.res.nfev = minres.nfev
        if hasattr(minres, "njev"):
            self.res.njev = minres.njev
        if hasattr(minres, "nhev"):
            self.res.nhev = minres.nhev

    def _monte_carlo_step(self):
        """Do one monte carlo iteration
        Randomly displace the coordinates, minimize, and decide whether
        or not to accept the new coordinates.
        """
        # Take a random step.  Make a copy of x because the step_taking
        # algorithm might change x in place
        x_after_step = np.copy(self.x)
        x_after_step = self.step_taking(x_after_step)

        # do a local minimization
        minres = self.minimizer(x_after_step)
        x_after_quench = minres.x
        energy_after_quench = minres.fun
        if not minres.success:
            self.res.minimization_failures += 1
            if self.disp:
                print("warning: minimahopping: local minimization failure")
        if hasattr(minres, "nfev"):
            self.res.nfev += minres.nfev
        if hasattr(minres, "njev"):
            self.res.njev += minres.njev
        if hasattr(minres, "nhev"):
            self.res.nhev += minres.nhev

        # accept the move based on self.accept_tests. If any test is false,
        # than reject the step.  If any test returns the special value, the
        # string 'force accept', accept the step regardless.  This can be used
        # to forcefully escape from a local minimum if normal minima hopping
        # steps are not sufficient.
        accept = True
        for test in self.accept_tests:
            testres = test(f_new=energy_after_quench, x_new=x_after_quench,
                           f_old=self.energy, x_old=self.x)
            if isinstance(testres, bool):
                if not testres:
                    accept = False
            elif isinstance(testres, str):
                if testres == "force accept":
                    accept = True
                    break
                else:
                    raise ValueError("accept test must return bool or string "
                                     "'force accept'. Type is", type(testres))
            else:
                raise ValueError("accept test must return bool or string "
                                 "'force accept'. Type is", type(testres))

        # Report the result of the acceptance test to the take step class.
        # This is for adaptive step taking
        if hasattr(self.step_taking, "report"):
            self.step_taking.report(accept, f_new=energy_after_quench,
                                    x_new=x_after_quench, f_old=self.energy,
                                    x_old=self.x)

        return x_after_quench, energy_after_quench, accept

    def one_cycle(self):
        """Do one cycle of the minimahopping algorithm
        """
        self.nstep += 1
        new_global_min = False

        xtrial, energy_trial, accept = self._monte_carlo_step()

        if accept:
            self.energy = energy_trial
            self.x = np.copy(xtrial)
            new_global_min = self.storage.update(self.x, self.energy)

        # print some information
        if self.disp:
            self.print_report(energy_trial, accept)
            if new_global_min:
                print("found new global minimum on step %d with function"
                      " value %g" % (self.nstep, self.energy))

        # save some variables as MinimaHoppingRunner attributes
        self.xtrial = xtrial
        self.energy_trial = energy_trial
        self.accept = accept

        return new_global_min

    def print_report(self, energy_trial, accept):
        """print a status update"""
        xlowest, energy_lowest = self.storage.get_lowest()
        print("minimahopping step %d: f %g trial_f %g accepted %d "
              " lowest_f %g" % (self.nstep, self.energy, energy_trial,
                                accept, energy_lowest))


class AdaptiveStepsize(object):
    """
    Class to implement adaptive stepsize.
    This class wraps the step taking class and modifies the stepsize to
    ensure the true acceptance rate is as close as possible to the target.
    Parameters
    ----------
    takestep : callable
        The step taking routine.  Must contain modifiable attribute
        takestep.stepsize
    accept_rate : float, optional
        The target step acceptance rate
    interval : int, optional
        Interval for how often to update the stepsize
    factor : float, optional
        The step size is multiplied or divided by this factor upon each
        update.
    verbose : bool, optional
        Print information about each update
    """
    def __init__(self, takestep, accept_rate=0.5, interval=50, factor=0.9,
                 verbose=True):
        self.takestep = takestep
        self.target_accept_rate = accept_rate
        self.interval = interval
        self.factor = factor
        self.verbose = verbose

        self.nstep = 0
        self.nstep_tot = 0
        self.naccept = 0

    def __call__(self, x):
        return self.take_step(x)

    def _adjust_step_size(self):
        old_stepsize = self.takestep.stepsize
        accept_rate = float(self.naccept) / self.nstep
        if accept_rate > self.target_accept_rate:
            # We're accepting too many steps.  This generally means we're
            # trapped in a basin.  Take bigger steps
            self.takestep.stepsize /= self.factor
        else:
            # We're not accepting enough steps.  Take smaller steps
            self.takestep.stepsize *= self.factor
        if self.verbose:
            print("adaptive stepsize: acceptance rate %f target %f new "
                  "stepsize %g old stepsize %g" % (
                    accept_rate,
                    self.target_accept_rate, self.takestep.stepsize,
                    old_stepsize))

    def take_step(self, x):
        self.nstep += 1
        self.nstep_tot += 1
        if self.nstep % self.interval == 0:
            self._adjust_step_size()
        return self.takestep(x)

    def report(self, accept, **kwargs):
        "called by minimahopping to report the result of the step"
        if accept:
            self.naccept += 1


class RandomDisplacement(object):
    """
    Add a random displacement of maximum size, stepsize, to the coordinates
    update x inplace
    """
    def __init__(self, stepsize=0.5):
        self.stepsize = stepsize

    def __call__(self, x):
        x += np.random.uniform(-self.stepsize, self.stepsize, np.shape(x))
        return x


class MinimizerWrapper(object):
    """
    wrap a minimizer function as a minimizer class
    """
    def __init__(self, minimizer, func=None, **kwargs):
        self.minimizer = minimizer
        self.func = func
        self.kwargs = kwargs

    def __call__(self, x0):
        if self.func is None:
            return self.minimizer(x0, **self.kwargs)
        else:
            return self.minimizer(self.func, x0, **self.kwargs)


class Metropolis(object):
    """
    Metropolis acceptance criterion

    http://bigdft.org/images/b/b4/2011-10_SG_Minhop.pdf
    Minima hopping can be transformed into basin hopping by 3 steps:
    1. Thresholding with E_diff is replaced by Metropolis step with
        exp(−E/E_diff ): no significant loss of performance

    2. A constant variable value of E_diff is replaced by a constant value:
        factor of 2 loss in performance

    3. MD based escapes with variable kinetic energy are replaced by random
        moves: orders of magnitude loss in performance
        – Either a very large temperature E_diff has to be allowed to escape
          soon from current minimum because random displacement go over high
          barriers and therefore into high energy local minima (BEP principle)
        – Or most accepted configurations are identical to current minimum
          and escapes happen rarely.
    """
    def __init__(self, T, Tmax=1000, alpha1=0.9803921569, alpha2=1.02,
                 beta1=1.05, beta2=1.05, beta3=0.9523809524):

        # self.beta = 1.0 / T
        self.E_diff = 0.5
        self.Tmax = Tmax
        self.temperature = T

        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.beta1 = beta1
        self.beta2 = beta2
        self.beta3 = beta3

    def accept_reject(self, energy_new, energy_old):
        # exp(deltaE/kbT)
        # w = min(1.0, np.exp(-(energy_new - energy_old) * self.beta))
        rand = np.random.rand()

        # ? exp(-E/E_diff) ?
        # Warniing: exp might be overflow if E/E_diff became very small
        # E = energy_new
        # E_diff = energy_new - energy_old
        # boltzmann_f = 1.0 / (self.E_diff*self.temperature)
        boltzmann_f = 1.0 / self.E_diff
        # boltzmann_f = 1.0 / self.temperature
        # w = min(1.0, np.exp(-(energy_new - energy_old)*boltzmann_f))
        w = np.exp(-(energy_new - energy_old)*boltzmann_f)
        # accepted = (energy_new - energy_old) < self.E_diff
        accepted = w >= rand
        if accepted:
            print("Accepted: Enew - {:5.3f} deltaE - {:5.3f} w - {:5.3f} rand - {:5.3f} Ediff - {:5.3f}".format(
                  energy_new, energy_new - energy_old, w, rand, self.E_diff))
            # Accept new minimum
            # alpha1 < 1 decrease thresholding
            # beta1 > 1  increase T
            # beta2 > 1  increase T
            self.E_diff *= self.alpha1
            # self.temperature *= self.beta1
            # if energy_new == energy_old:
            #    self.temperature *= self.beta2
        else:
            # Reject new minimum
            # alpha2 > 1 increase thresholding
            # beta3 < 1 decrease T
            print("Rejected: Enew - {:5.3f} deltaE - {:5.3f} w - {:5.3f} rand - {:5.3f} Ediff - {:5.3f}".format(
                  energy_new, energy_new - energy_old, w, rand, self.E_diff))
            self.E_diff *= self.alpha2
            # self.temperature *= self.beta3
            if self.E_diff <= 0.0:
                self.E_diff = 0.5
            elif self.E_diff > 1:
                self.E_diff = 1
        return accepted

    def __call__(self, **kwargs):
        """
        f_new and f_old are mandatory in kwargs
        """
        return bool(self.accept_reject(kwargs["f_new"],
                    kwargs["f_old"]))


def minimahopping(func, x0, niter=100,
                  T=1.0,
                  stepsize=0.5,
                  alpha1=0.9803921569,
                  alpha2=1.02,
                  beta1=1.05,
                  beta2=1.05,
                  beta3=0.9523809524,
                  minimizer_kwargs=None,
                  take_step=None,
                  accept_test=None,
                  callback=None,
                  interval=50,
                  disp=False,
                  niter_success=None):
    """
    Pseudocode in the original papaer
    ------------------------------------------------------------------------
    initialize a current minimum ‘Mcurrent’
    MDstart
       ESCAPE TRIAL PART
        start a MD trajectory with kinetic energy Ekinetic from
        current minimum ‘Mcurrent’. Once the potential reaches the
        mdmin-th minimum along the trajectory stop MD and optimize geometry
        to find the closest local minimum ‘M’
        if (‘M’ equals ‘Mcurrent’) then
            Ekinetic = Ekinetic * betal (betal>1)
            goto MDstart
        else if (‘M’ equals a minimum visited previously) then
            Ekinetic = Ekinetic * beta2 (beta2>1)
            goto MDstart
        else if (‘M’ equals a new minimum) then
            Ekinetic = Ekinetic * beta3 (beta3<1)
        endif
       DOWNWARD PREFERENCE PART
        if (energy(‘M’)–energy(‘Mcurrent’)<Ediff) then
            accept new minimum: ‘Mcurrent’=‘M’
            add ‘Mcurrent’ to history list
            Ediff = Ediff * alpha1 (alpha1<1)
        else if rejected
            Ediff = Ediff * alpha2 (alpha2>1)
        endif
        goto MDstart
    ----------------------------------------------------------------------
    Args:  Standard parameter set presented in the original paper
        alpha1 (float): =1/alpha2, which decreases E_diff when a new minimum is
                        accepted.
        alpha2 (float): = 1.02, which increases the threshold on rejection.
        beta1 (float): = 1.05, increases the kinetic energy when a MD escape
                        trial fails.
        beta2 (float): = 1.05,
                    increases E_kin when the new minimum is already known.
        beta3 (float): = 1/1.05, decreases the kinetic energy if the minimum
                        is unknown.

    The following parameters are from minimahopping
    func : callable ``f(x, *args)``
        Function to be optimized.  ``args`` can be passed as an optional item
        in the dict ``minimizer_kwargs``
    x0 : ndarray
        Initial guess.
    niter : integer, optional
        The number of minima hopping iterations
    T : float, optional
        The "temperature" parameter for the accept or reject criterion.  Higher
        "temperatures" mean that larger jumps in function value will be
        accepted.  For best results ``T`` should be comparable to the
        separation
        (in function value) between local minima.
    stepsize : float, optional
        initial step size for use in the random displacement.
    minimizer_kwargs : dict, optional
        Extra keyword arguments to be passed to the minimizer
        ``scipy.optimize.minimize()`` Some important options could be:
            method : str
                The minimization method (e.g. ``"L-BFGS-B"``)
            args : tuple
                Extra arguments passed to the objective function (``func``) and
                its derivatives (Jacobian, Hessian).
    take_step : callable ``take_step(x)``, optional
        Replace the default step taking routine with this routine.  The default
        step taking routine is a random displacement of the coordinates, but
        other step taking algorithms may be better for some systems.
        ``take_step`` can optionally have the attribute ``take_step.stepsize``.
        If this attribute exists, then ``minimahopping`` will adjust
        ``take_step.stepsize`` in order to try to optimize the global minimum
        search.
    accept_test : callable, ``accept_test(f_new=f_new, x_new=x_new, f_old=fold,
        x_old=x_old)``, optional
        Define a test which will be used to judge whether or not to accept the
        step.  This will be used in addition to the Metropolis test based on
        "temperature" ``T``.  The acceptable return values are True,
        False, or ``"force accept"``.  If the latter, then this will
        override any other tests in order to accept the step.  This can be
        used, for example, to forcefully escape from a local minimum that
        ``minimahopping`` is trapped in.
    callback : callable, ``callback(x, f, accept)``, optional
        A callback function which will be called for all minima found.  ``x``
        and ``f`` are the coordinates and function value of the trial minimum,
        and ``accept`` is whether or not that minimum was accepted. This can be
        used, for example, to save the lowest N minima found.  Also,
        ``callback`` can be used to specify a user defined stop criterion by
        optionally returning True to stop the ``minimahopping`` routine.
    interval : integer, optional
        interval for how often to update the ``stepsize``
    disp : bool, optional
        Set to True to print status messages
    niter_success : integer, optional
        Stop the run if the global minimum candidate remains the same for this
        number of iterations.
    """
    x0 = np.array(x0)

    # set up minimizer
    if minimizer_kwargs is None:
        minimizer_kwargs = dict()
    wrapped_minimizer = MinimizerWrapper(opt.minimize, func,
                                         **minimizer_kwargs)

    # set up step taking algorithm
    if take_step is not None:
        if not isinstance(take_step, collections.Callable):
            raise TypeError("take_step must be callable")
        # if take_step.stepsize exists then use AdaptiveStepsize to control
        # take_step.stepsize
        if hasattr(take_step, "stepsize"):
            take_step_wrapped = AdaptiveStepsize(take_step, interval=interval,
                                                 verbose=disp)
        else:
            take_step_wrapped = take_step
    else:
        # use default
        displace = RandomDisplacement(stepsize=stepsize)
        take_step_wrapped = AdaptiveStepsize(displace, interval=interval,
                                             verbose=disp)

    # set up accept tests
    if accept_test is not None:
        if not isinstance(accept_test, collections.Callable):
            raise TypeError("accept_test must be callable")
        accept_tests = [accept_test]
    else:
        accept_tests = []
    # use default
    metropolis = Metropolis(T, alpha1=alpha1, alpha2=alpha2)
    accept_tests.append(metropolis)

    if niter_success is None:
        niter_success = niter + 2

    mh = MinimaHoppingRunner(x0, wrapped_minimizer, take_step_wrapped,
                             accept_tests,
                             alpha1=alpha1,
                             alpha2=alpha2,
                             beta1=beta1,
                             beta2=beta2,
                             beta3=beta3, disp=disp)

    # start main iteration loop
    count = 0
    message = ["requested number of minimahopping iterations completed"
               " successfully"]
    for i in range(niter):
        new_global_min = mh.one_cycle()

        if isinstance(callback, collections.Callable):
            # should we pass a copy of x?
            val = callback(mh.xtrial, mh.energy_trial, mh.accept)
            if val is not None:
                if val:
                    message = ["callback function requested stop early by"
                               "returning True"]
                    break

        count += 1
        if new_global_min:
            count = 0
        elif count > niter_success:
            message = ["success condition satisfied"]
            break

    # prepare return object
    lowest = mh.storage.get_lowest()
    res = mh.res
    res.x = np.copy(lowest[0])
    res.fun = lowest[1]
    res.message = message
    res.nit = i + 1
    return res
