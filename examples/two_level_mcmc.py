"""MCMC protocol for the Bayesian hierarchical modeling
   This involves two stages of parameter estimation.
   The first stage is a simple MCMC step, whereas the 
   second stage involves the joint posterior estimation of the 
   local and global parameters. Use the core IMP library for 
   most of the MCMC steps."""

from __future__ import print_function, division
import IMP
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.analysis
import IMP.pmi.io
import IMP.pmi.alphabets
import IMP.rmf
import IMP.isd
import IMP.pmi.dof
import numpy as np
import numpyro
from numpyro.infer import MCMC, NUTS
import numpyro.distributions as dist
import matplotlib.pyplot as plt
import jax.random as random
import pymc as pm

class TwoLevelMCMC:
    def __init__(self, data):
        self.data = data

    
    def create_mu(self):
        """Create a nuisance on the length of the helix."""
        lengthinit = 20.0
        self.lengthissampled = True
        lengthminnuis = 4.0
        lengthmaxnuis = 40.0
        length = IMP.pmi.tools.SetupNuisance(self.model, lengthinit,
                                             lengthminnuis, lengthmaxnuis,
                                             self.lengthissampled
                                             ).get_particle()
        self.rs.add_restraint(
            IMP.isd.LognormalRestraint(length, 15.0, 6.0))
        
    def model(self):
        with pm.Model() as model:
            # Define priors for global parameters
            global_mean = pm.Normal("global_mean", mu=0, sd=10)
            global_variance = pm.HalfNormal("global_variance", sd=10)
            # Define priors for local parameters
            local_mean = pm.Normal("local_mean", mu=0, sd=10)
            local_variance = pm.HalfNormal("local_variance", sd=10)
            # Define the global Gaussian process
            kernel = pm.Exponential("kernel", lam=1)
            global_gp = pm.GaussianProcess("global_gp", kernel=kernel, observed=self.data["x"], sigma=global_variance)
            # Define the local Gaussian process
            local_gp = pm.GaussianProcess("local_gp", kernel=kernel, observed=self.data["x"], sigma=local_variance)
            # Define the likelihood
            obs = pm.Normal("obs", mu=global_mean + local_mean, sd=global_gp + local_gp, observed=self.data["y"])
        return model
    
    def run_mcmc(self, num_warmup, num_samples, num_chains):
        model = self.model()
        with model:
            trace = pm.sample(num_samples, tune=num_warmup, chains=num_chains)
        # Plot the resultsts
        pm.plot_posterior(trace, var_names=["global_mean", "local_mean"])
        plt.show()

# Define an instance of the class
data = {"x": np.array([1, 2, 3, 4, 5]), "y": np.array([2, 4, 6, 8, 10])}
mcmc_instance = TwoLevelMCMC(data)
mcmc_instance.run_mcmc(1000, 1000, 1)