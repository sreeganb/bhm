#!/usr/bin/env python
"""MCMC protocol for the Bayesian hierarchical modeling
   This involves two stages of parameter estimation.
   The first stage is a simple MCMC step, whereas the 
   second stage involves the joint posterior estimation of the 
   local and global parameters. Use the core IMP library for 
   most of the MCMC steps."""
import IMP
import pymc
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.restraints
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm

class TwoLevelMCMC:
    def __init__(self, data):
        self.data = data
        self.model = pymc.Model()
        self.mean = pymc.Normal('mean', mu=0, tau=1)
        self.precision = pymc.Gamma('precision', alpha=1, beta=1)
        self.observed = pymc.Normal('observed', mu=self.mean, tau=self.precision, observed=self.data)
        self.mcmc = pymc.MCMC(self.model)
    
    def run_mcmc(self, iterations, burn, thin):
        self.mcmc.sample(iter=iterations, burn=burn, thin=thin)
    
    def get_samples(self):
        samples_mean = self.mcmc.trace('mean')[:]
        samples_precision = self.mcmc.trace('precision')[:]
        return samples_mean, samples_precision

# Create an instance of the TwoLevelMCMC class
data = [1, 2, 3, 4, 5]
mcmc = TwoLevelMCMC(data)
mcmc.run_mcmc(iterations=1000, burn=100, thin=10)
samples_mean, samples_precision = mcmc.get_samples()
print(samples_mean)