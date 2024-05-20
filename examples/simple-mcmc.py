#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Define the target distribution
def target_distribution(x):
    # Define your target distribution here
    # Example: Gaussian distribution with mean 0 and standard deviation 1
    return np.exp(-0.5 * x**2) / np.sqrt(2 * np.pi)

# Define the MCMC function
def mcmc_sampling(initial_state, num_iterations):
    # Initialize the current state
    current_state = initial_state
    # Initialize the accepted samples list
    accepted_samples = []
    
    # Run MCMC iterations
    for iteration in range(num_iterations):
        # Propose a new state
        proposed_state = np.random.normal(current_state, 1)
        
        # Calculate the acceptance probability
        acceptance_prob = min(1, target_distribution(proposed_state) / target_distribution(current_state))
        
        # Accept or reject the proposed state
        if np.random.uniform() < acceptance_prob:
            current_state = proposed_state
            accepted_samples.append(current_state)
    
    return accepted_samples

# Set the initial state and number of iterations
initial_state = 0
num_iterations = 100

# Run the MCMC sampling
samples = mcmc_sampling(initial_state, num_iterations)

# Print the final samples
print("Final samples:", samples)
import matplotlib.pyplot as plt

# Plot the final sampled distribution
plt.hist(samples, bins=30, density=True, alpha=0.5)
plt.xlabel('Samples')
plt.ylabel('Probability Density')
plt.title('Final Sampled Distribution')
plt.show()