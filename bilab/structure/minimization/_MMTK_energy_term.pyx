# Energy term interface for Python code
#
# Written by Konrad Hinsen
# last revision: 2007-1-11
#

include 'MMTK/forcefield.pxi'

import bilab.structure.minimization
#import MMTK
#from Scientific import N
#from numpy.oldnumeric import *
#from numpy import ndarray as array_type
import numpy as N

cdef class PyEnergyTerm(EnergyTerm):

    cdef public object universe
    cdef object configuration

    def __init__(self, name, universe):
        EnergyTerm.__init__(self, universe, name, (name,))
        self.eval_func = <void *>PyEnergyTerm.c_evaluate
        self.universe = universe
        self.configuration = None

    def __getattr__(self, name):
        if name == "name":
            return self.evaluator_name
        elif name == "term_names":
            names = []
            for i in range(self.nterms):
                names.append(self.term_names[i])
            return tuple(names)
        else:
            raise AttributeError, name

    cdef void c_evaluate(self, EnergyEvaluator eval,
                         energy_spec *input, energy_data *energy) except *:
        cdef int do_gradients, do_fc
        do_gradients = energy.gradients != NULL
        do_fc = energy.force_constants != NULL
        if self.configuration is None:
            self.configuration = bilab.structure.minimization.Configuration(self.universe,
                                                    <object>input.coordinates)
        else:
            self.configuration.array = <object>input.coordinates
        try:
            energy_dict = self.evaluate(self.configuration,
                                        do_gradients, do_fc)
        except:
            energy.error = 1
            raise
        energy.energy_terms[self.index] = energy_dict['energy']
        try:
            energy.energy_terms[self.virial_index] = \
                energy.energy_terms[self.virial_index] + energy_dict['virial']
        except KeyError:
            energy.virial_available = 0
        if do_gradients:
            try:
                N.add(<object>energy.gradients,
                      energy_dict['gradients'].array,
                      <object>energy.gradients)
            except (KeyError, AttributeError):
                energy.error = 1
                raise
        if do_fc:
            try:
                N.add(<object>energy.force_constants,
                      energy_dict['force_constants'].array,
                      <object>energy.force_constants)
            except (KeyError, AttributeError):
                energy.error = 1
                raise
