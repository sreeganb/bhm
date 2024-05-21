#!/usr/bin/env python

import torch
import pyro
import pyro.distributions as dist
import matplotlib.pyplot as plt

class HMC:
    def __init__(self):
        self.mean = torch.tensor([0., 0.])
        self.cov = torch.tensor([[1., 0.8], [0.8, 1.]])
    
    def model(self):
        x = pyro.sample("x", dist.MultivariateNormal(self.mean, self.cov))
    
    def run(self):
        num_samples = 1000
        num_burnin_steps = 100
        
        hmc_kernel = pyro.infer.mcmc.HMC(self.model)
        mcmc_run = pyro.infer.mcmc.MCMC(hmc_kernel, num_samples=num_samples, warmup_steps=num_burnin_steps)
        mcmc_run.run()
        
        samples = mcmc_run.get_samples()
        return samples

if __name__ == "__main__":
    hmc = HMC()
    samples = hmc.run()
    plt.scatter(samples['x'][:, 0], samples['x'][:, 1])
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.title('Sampled Distribution')
    plt.show()