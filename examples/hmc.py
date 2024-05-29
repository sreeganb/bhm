#!/usr/bin/env python

import torch
import pyro
import numpy as np
import pyro.distributions as dist
import matplotlib.pyplot as plt

# Create a numpy array
data = np.array([[1, 2], [3, 4], [5, 6]])
print(data)
# Convert numpy array to torch tensor
tensor = torch.from_numpy(data)

# Convert torch tensor to pyro-compatible tensor
pyro_tensor = tensor.float()
print(pyro_tensor)

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

    class TwoLevelHMC(HMC):
        def __init__(self):
            super().__init__()
            self.mean_mean = torch.tensor([0., 0.])
            self.mean_cov = torch.tensor([[1., 0.5], [0.5, 1.]])
        
        def model(self):
            mean = pyro.sample("mean", dist.MultivariateNormal(self.mean_mean, self.mean_cov))
            x = pyro.sample("x", dist.MultivariateNormal(mean, self.cov))
        
    if __name__ == "__main__":
        hmc = TwoLevelHMC()
        samples = hmc.run()
        plt.scatter(samples['x'][:, 0], samples['x'][:, 1])
        plt.xlabel('x1')
        plt.ylabel('x2')
        plt.title('Sampled Distribution')
        plt.show()