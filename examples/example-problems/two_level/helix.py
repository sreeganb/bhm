import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pymc as pm
import arviz as az

# Set random seed for reproducibility
np.random.seed(10101)

# True parameters for simulation
pi_true = 0.6  # Probability that a monomer is in the helical state
mu_helix_true = 15
sigma_helix_true = 2
mu_glob_true = 6
sigma_glob_true = 3

# Simulate monomer data
N_monomer = 800
monomer_states = np.random.binomial(1, pi_true, size=N_monomer)
monomer_obs = np.where(
    monomer_states == 1,
    np.random.normal(mu_helix_true, sigma_helix_true, size=N_monomer),
    np.random.normal(mu_glob_true, sigma_glob_true, size=N_monomer)
)

# Simulate dimer data
# Dimer states:
# 0: globular-globular
# 1: globular-helical
# 2: helical-globular
# 3: helical-helical

# Interaction terms (delta) for each dimer state
delta_true = np.array([0.0, 0.0, 0.0, 0.0])  # Adjust if interactions exist

# Probabilities for each dimer state
pi_dimer_true = np.array([
    (1 - pi_true) ** 2,   # globular-globular
    (1 - pi_true) * pi_true,  # globular-helical
    pi_true * (1 - pi_true),  # helical-globular
    pi_true ** 2           # helical-helical
])

# Dimer means and standard deviations
mu_dimer_true = np.array([
    mu_glob_true + delta_true[0],
    (mu_glob_true + mu_helix_true) / 2 + delta_true[1],
    (mu_glob_true + mu_helix_true) / 2 + delta_true[2],
    mu_helix_true + delta_true[3]
])

sigma_dimer_true = np.array([
    sigma_glob_true,
    (sigma_glob_true + sigma_helix_true) / 2,
    (sigma_glob_true + sigma_helix_true) / 2,
    sigma_helix_true
])

# Simulate dimer data
N_dimer = 200
dimer_states = np.random.choice(4, size=N_dimer, p=pi_dimer_true)
dimer_obs = np.random.normal(mu_dimer_true[dimer_states], sigma_dimer_true[dimer_states])

# Plot the simulated data
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
sns.histplot(monomer_obs, bins=15, kde=True)
plt.title('Simulated Monomer End-to-End Distances')
plt.xlabel('Distance')
plt.subplot(1, 2, 2)
sns.histplot(dimer_obs, bins=15, kde=True)
plt.title('Simulated Dimer End-to-End Distances')
plt.xlabel('Distance')
plt.tight_layout()
plt.show()

# Hierarchical Bayesian Model
with pm.Model() as hierarchical_model:
    # Hyperpriors for the monomer means
    mu_mu = pm.Normal('mu_mu', mu=10, sigma=5)
    sigma_mu = pm.HalfNormal('sigma_mu', sigma=5)
    
    # Monomer state probabilities
    pi = pm.Beta('pi', alpha=2, beta=2)
    
    # Monomer parameters
    mu_glob = pm.Normal('mu_glob', mu=mu_mu, sigma=sigma_mu)
    sigma_glob = pm.HalfNormal('sigma_glob', sigma=5)
    
    mu_helix = pm.Normal('mu_helix', mu=mu_mu, sigma=sigma_mu)
    sigma_helix = pm.HalfNormal('sigma_helix', sigma=5)
    
    # Monomer mixture components
    monomer_components = pm.Normal.dist(
        mu=[mu_glob, mu_helix],
        sigma=[sigma_glob, sigma_helix]
    )
    
    # Monomer mixture model
    monomer_obs_ = pm.Mixture(
        'monomer_obs',
        w=[1 - pi, pi],
        comp_dists=monomer_components,
        observed=monomer_obs
    )
    
    # Dimer state probabilities based on monomer probabilities
    pi_dimer = pm.Deterministic('pi_dimer', pm.math.stack([
        (1 - pi) ** 2,        # globular-globular
        (1 - pi) * pi,        # globular-helical
        pi * (1 - pi),        # helical-globular
        pi ** 2               # helical-helical
    ]))
    
    # Interaction terms for dimer states
    delta = pm.Normal('delta', mu=0, sigma=1, shape=4)
    
    # Dimer means and standard deviations
    mu_dimer = pm.Deterministic('mu_dimer', pm.math.stack([
        mu_glob + delta[0],
        (mu_glob + mu_helix) / 2 + delta[1],
        (mu_glob + mu_helix) / 2 + delta[2],
        mu_helix + delta[3]
    ]))
    
    sigma_dimer = pm.Deterministic('sigma_dimer', pm.math.stack([
        sigma_glob,
        (sigma_glob + sigma_helix) / 2,
        (sigma_glob + sigma_helix) / 2,
        sigma_helix
    ]))
    
    # Dimer mixture components
    dimer_components = pm.Normal.dist(mu=mu_dimer, sigma=sigma_dimer)
    
    # Dimer mixture model
    dimer_obs_ = pm.Mixture(
        'dimer_obs',
        w=pi_dimer,
        comp_dists=dimer_components,
        observed=dimer_obs
    )
    
    # Sampling
    trace = pm.sample(
        draws=3000, 
        tune=2000, 
        target_accept=0.9, 
        cores=1, 
        random_seed=42
    )

# Summarize the results
summary = az.summary(trace, var_names=[
    'pi', 'mu_glob', 'sigma_glob', 'mu_helix', 'sigma_helix', 'delta'
])
print(summary)

# Plot posterior distributions
az.plot_posterior(trace, var_names=[
    'pi', 'mu_glob', 'sigma_glob', 'mu_helix', 'sigma_helix', 'delta'
])
plt.show()