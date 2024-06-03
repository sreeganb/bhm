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
from operator import itemgetter
from collections import defaultdict
import numpy as np
import random
import numpyro
from numpyro.infer import MCMC, NUTS
import numpyro.distributions as dist
import matplotlib.pyplot as plt
import jax.random as random
class TwoLevelMCMC:
    def __init__(self, data):
        self.data = data

    def model(self):
        # Define priors for global parameters
        global_mean = numpyro.sample("global_mean", dist.Normal(0, 10))
        global_variance = numpyro.sample("global_variance", dist.HalfNormal(10))

        # Define priors for local parameters
        local_mean = numpyro.sample("local_mean", dist.Normal(0, 10))
        local_variance = numpyro.sample("local_variance", dist.HalfNormal(10))

        # Define the global Gaussian process
        kernel = numpyro.sample("kernel", dist.Exponential(1))
        global_gp = numpyro.sample("global_gp", dist.GaussianProcess(kernel, self.data["x"], global_variance))

        # Define the local Gaussian process
        local_gp = numpyro.sample("local_gp", dist.GaussianProcess(kernel, self.data["x"], local_variance))

        # Define the likelihood
        numpyro.sample("obs", dist.Normal(global_mean + local_mean, global_gp + local_gp), obs=self.data["y"])

    def run_mcmc(self, num_warmup, num_samples, num_chains):
        kernel = NUTS(self.model)
        mcmc = MCMC(kernel, num_warmup=num_warmup, num_samples=num_samples, num_chains=num_chains)
        mcmc.run(random.PRNGKey(0))
        mcmc.print_summary()

        # Get the posterior samples
        samples = mcmc.get_samples()

        # Plot the results
        plt.figure()
        plt.plot(self.data["x"], self.data["y"], "ro", label="Data")
        plt.plot(self.data["x"], samples["global_mean"].mean(axis=0), label="Global Mean")
        plt.plot(self.data["x"], samples["local_mean"].mean(axis=0), label="Local Mean")
        plt.fill_between(self.data["x"], samples["global_mean"].mean(axis=0) - 2 * samples["global_mean"].std(axis=0),
                         samples["global_mean"].mean(axis=0) + 2 * samples["global_mean"].std(axis=0), alpha=0.3)
        plt.fill_between(self.data["x"], samples["local_mean"].mean(axis=0) - 2 * samples["local_mean"].std(axis=0),
                         samples["local_mean"].mean(axis=0) + 2 * samples["local_mean"].std(axis=0), alpha=0.3)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.show()

# Define an instance of the class
data = {"x": np.array([1, 2, 3, 4, 5]), "y": np.array([2, 4, 6, 8, 10])}
mcmc_instance = TwoLevelMCMC(data)
mcmc_instance.run_mcmc(1000, 1000, 1)
