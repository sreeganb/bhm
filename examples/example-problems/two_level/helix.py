import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pymc as pm
import pymc.math as pm_math  # For PyMC 4.x

# Set random seed for reproducibility
np.random.seed(42)

# True parameters for simulation
pi_true = 0.6  # Probability that a monomer is in the helical state
mu_helix_true = 15
sigma_helix_true = 2
mu_glob_true = 7
sigma_glob_true = 3

# Simulate monomer data
N_monomer = 500
monomer_states = np.random.binomial(1, pi_true, size=N_monomer)
monomer_obs = np.where(
    monomer_states == 1,
    np.random.normal(mu_helix_true, sigma_helix_true, size=N_monomer),
    np.random.normal(mu_glob_true, sigma_glob_true, size=N_monomer)
)

# Simulate dimer data
# Possible dimer states:
# 0: globular-globular
# 1: globular-helical
# 2: helical-globular
# 3: helical-helical

# Interaction terms for each dimer state
delta_true = np.array([1.5, 1.0, 1.0, 0.5])

# Create all possible combinations of monomer states
monomer_pairs = np.array(np.meshgrid([0, 1], [0, 1])).T.reshape(-1, 2)
dimer_states_possible = monomer_pairs[:, 0] + monomer_pairs[:, 1]*2  # Encode states as 0-3

# Probabilities for each dimer state
pi_dimer_true = np.array([
    (1 - pi_true)**2,     # globular-globular
    (1 - pi_true)*pi_true,  # globular-helical
    pi_true*(1 - pi_true),  # helical-globular
    pi_true**2            # helical-helical
])
pi_dimer_true /= pi_dimer_true.sum()  # Normalize

# Means and standard deviations for dimer states
mu_dimer_true = np.array([
    2*mu_glob_true + delta_true[0],
    mu_glob_true + mu_helix_true + delta_true[1],
    mu_helix_true + mu_glob_true + delta_true[2],
    2*mu_helix_true + delta_true[3]
])
sigma_dimer_true = np.array([
    np.sqrt(2*sigma_glob_true**2),
    np.sqrt(sigma_glob_true**2 + sigma_helix_true**2),
    np.sqrt(sigma_helix_true**2 + sigma_glob_true**2),
    np.sqrt(2*sigma_helix_true**2)
])

# Simulate dimer data
N_dimer = 400
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

with pm.Model() as hierarchical_model:
    # Hyperpriors for the monomer means
    mu_mu = pm.Normal('mu_mu', mu=10, sigma=5)  # Adjusted mean and sigma

    # Hyperpriors for the monomer standard deviations
    mu_sigma = pm.HalfNormal('mu_sigma', sigma=2)

    # Prior for monomer state probability π
    pi = pm.Beta('pi', alpha=2.0, beta=2.0)

    # Monomer helical state parameters
    mu_helix = pm.Normal('mu_helix', mu=mu_mu, sigma=2)        # Narrower sigma
    sigma_helix = pm.HalfNormal('sigma_helix', sigma=2)

    # Monomer globular state parameters
    mu_glob = pm.Normal('mu_glob', mu=mu_mu, sigma=2)          # Narrower sigma
    sigma_glob = pm.HalfNormal('sigma_glob', sigma=2)
    
    # Monomer mixture components
    monomer_components = [
        pm.Normal.dist(mu=mu_glob, sigma=sigma_glob),
        pm.Normal.dist(mu=mu_helix, sigma=sigma_helix)
    ]
    
    # Monomer mixture model
    monomer_obs_ = pm.Mixture(
        'monomer_obs',
        w=[1 - pi, pi],
        comp_dists=monomer_components,
        observed=monomer_obs
    )
    
    # Prior for dimer state probabilities π_k
    alpha_dirichlet = np.ones(4)
    pi_dimer = pm.Dirichlet('pi_dimer', a=alpha_dirichlet)
    
    # Priors for interaction terms δ_k
    delta = pm.Normal('delta', mu=1.0, sigma=0.5, shape=4)
    
    # Dimer state means and standard deviations
    mu_dimer = pm.Deterministic('mu_dimer', pm.math.stack([
        2 * mu_glob + delta[0],          # globular-globular
        mu_glob + mu_helix + delta[1],   # globular-helical
        mu_helix + mu_glob + delta[2],   # helical-globular
        2 * mu_helix + delta[3]          # helical-helical
    ]))

    sigma_dimer = pm.Deterministic('sigma_dimer', pm.math.stack([
        pm.math.sqrt(2 * sigma_glob**2),
        pm.math.sqrt(sigma_glob**2 + sigma_helix**2),
        pm.math.sqrt(sigma_helix**2 + sigma_glob**2),
        pm.math.sqrt(2 * sigma_helix**2)
    ]))
    
    # Dimer mixture components
    dimer_components = [
        pm.Normal.dist(mu=mu_dimer[i], sigma=sigma_dimer[i]) for i in range(4)
    ]
    
    # Dimer mixture model
    dimer_obs_ = pm.Mixture(
        'dimer_obs',
        w=pi_dimer,
        comp_dists=dimer_components,
        observed=dimer_obs
    )
    
    # Sample from the posterior
    trace = pm.sample(5000, tune=2000, cores=4, random_seed=44, 
                      target_accept=0.9, return_inferencedata=True)
    
import arviz as az

# Summary of the posterior
summary = az.summary(trace, var_names=['pi', 'mu_helix', 'sigma_helix', 'mu_glob', 'sigma_glob', 'pi_dimer', 'delta'])

print(summary)

# Plot trace plots for key parameters
az.plot_trace(
    trace,
    var_names=['pi', 'mu_helix', 'sigma_helix', 'mu_glob', 'sigma_glob', 'delta'],
    compact=True
)
plt.show()

# Plot posterior distributions
az.plot_posterior(trace, var_names=['pi', 'mu_helix', 'sigma_helix', 'mu_glob', 'sigma_glob'], hdi_prob=0.95)
plt.show()

# Plot posterior distributions of interaction terms
az.plot_posterior(trace, var_names=['delta'], hdi_prob=0.95)
plt.show()

# Compare true values with posterior estimates
print(f"True pi: {pi_true}")
print(f"Estimated pi mean: {summary.loc['pi', 'mean']}")
print(f"True mu_helix: {mu_helix_true}")
print(f"Estimated mu_helix mean: {summary.loc['mu_helix', 'mean']}")
print(f"True mu_glob: {mu_glob_true}")
print(f"Estimated mu_glob mean: {summary.loc['mu_glob', 'mean']}")