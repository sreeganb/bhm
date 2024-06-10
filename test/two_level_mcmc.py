#!/usr/bin/env python
"""MCMC protocol for the Bayesian hierarchical modeling
   This involves two stages of parameter estimation.
   The first stage is a simple MCMC step, whereas the 
   second stage involves the joint posterior estimation of the 
   local and global parameters. Use the core IMP library for 
   most of the MCMC steps."""
import IMP
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.restraints


class TwoLevelMCMC:
    def __init__(self):
        self.m = IMP.Model()

    def _create_global_particle(self):
        # Create a nuisance particle to denote the end to end distance of the two strings
        # Give it a non-informative prior (Jeffrey's prior maybe)
        self.chi_is_sampled = True
        chiminnuis = 1.0
        chimaxnuis = 30.0
        chiinit = 15.0
        chimin = 0.01
        chimax = 100.0
        #chitrans = 0.5
        chi = IMP.pmi.tools.SetupNuisance(self.model, chiinit, chiminnuis, chimaxnuis, self.chi_is_sampled).get_particle()
        self.rs.add_restraint(IMP.isd.UniformPrior(self.model, chi, 10000.0, chimax, chimin))
        
        self.global_particle = IMP.Particle(self.m)
        self.global_particle.set_name('global')
        self.global_particle.add_attribute(IMP.FloatKey('mean'), 0)
        self.global_particle.add_attribute(IMP.FloatKey('precision'), 1)
        