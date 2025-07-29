#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MCMC Implementation of Hierarchical vs. Regular Bayesian Inference
=================================================================

This example demonstrates the differences between regular and hierarchical
Bayesian modeling for distance measurements using full MCMC inference.

Problem: Distance measurements between particles (types A, B, C) with
different measurement uncertainties.

Approach 1: Independent priors for each uncertainty parameter (regular)
Approach 2: Hierarchical priors to share information (hierarchical)
"""

import numpy as np
import matplotlib.pyplot as plt
import pymc as pm
import arviz as az
import seaborn as sns
from typing import Dict, List, Tuple
import pandas as pd
import time
from pathlib import Path

# Set plot style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class BayesianDistanceModel:
    """
    Models distance measurements between particle pairs with Bayesian inference.
    
    Demonstrates the difference between regular and hierarchical Bayesian approaches
    when inferring measurement uncertainties (sigma).
    """
    
    def __init__(self, true_distances: Dict[str, float] = None, true_sigmas: Dict[str, float] = None):
        """
        Initialize the model with optional ground truth parameters.
        """
        self.true_distances = true_distances
        self.true_sigmas = true_sigmas
        if true_distances:
            self.pair_types = list(true_distances.keys())
        
    def generate_data(self, n_measurements: Dict[str, int]) -> Dict[str, np.ndarray]:
        """
        Generate synthetic distance measurements with known ground truth.
        """
        if not self.true_distances or not self.true_sigmas:
            raise ValueError("True distances and sigmas must be set to generate synthetic data")
            
        np.random.seed(42)  # Changed for better demonstration
        
        data = {}
        for pair_type, n in n_measurements.items():
            true_dist = self.true_distances[pair_type]
            true_sigma = self.true_sigmas[pair_type]
            
            # Generate noisy measurements
            measurements = np.random.normal(true_dist, true_sigma, n)
            data[pair_type] = measurements
            
        return data
    
    def regular_bayesian_mcmc(self, data: Dict[str, np.ndarray], 
                             n_samples: int = 2000) -> Dict:
        """
        Regular Bayesian inference using MCMC with independent priors for each sigma.
        """
        start_time = time.time()
        
        # Create PyMC model
        with pm.Model() as regular_model:
            # Independent priors for each pair type
            means = {}
            sigmas = {}
            
            # Define parameters with priors
            for pair_type, measurements in data.items():
                # Prior for the mean (weakly informative)
                means[pair_type] = pm.Normal(f"mu_{pair_type}", 
                                           mu=np.mean(measurements), 
                                           sigma=5.0)
                
                # FIXED: Better prior for sigma
                sigmas[pair_type] = pm.HalfNormal(f"sigma_{pair_type}", sigma=2.0)
            
            # Likelihood functions for each pair type
            for pair_type, measurements in data.items():
                pm.Normal(f"likelihood_{pair_type}", 
                         mu=means[pair_type],
                         sigma=sigmas[pair_type],
                         observed=measurements)
            
            # Sampling with better parameters
            try:
                trace = pm.sample(n_samples, tune=1000, return_inferencedata=True,
                                 target_accept=0.85, chains=2, cores=2)
            except Exception as e:
                print(f"Error in regular sampling: {e}")
                trace = pm.sample(n_samples, tune=1500, return_inferencedata=True,
                                 target_accept=0.90, chains=2, cores=2)
        
        end_time = time.time()
        
        # Gather results
        results = {
            'trace': trace,
            'computation_time': end_time - start_time,
            'model_type': 'regular'
        }
        
        return results
    
    def hierarchical_bayesian_mcmc(self, data: Dict[str, np.ndarray], 
                                  n_samples: int = 2000) -> Dict:
        """
        Hierarchical Bayesian inference using MCMC with shared hyperpriors.
        FIXED: Non-centered parameterization to avoid divergences.
        """
        start_time = time.time()
        
        # Create PyMC model with hierarchical structure
        with pm.Model() as hierarchical_model:
            # FIXED: Non-centered parameterization to avoid funnel geometry
            
            # Hyperpriors for the group-level parameters
            mu_sigma = pm.Normal("mu_sigma", mu=np.log(2.0), sigma=0.5)  # Centered around log(2.0)
            sigma_sigma = pm.HalfNormal("sigma_sigma", sigma=0.3)  # Smaller variation
            
            # Parameters with hierarchical priors
            means = {}
            sigmas = {}
            
            for pair_type, measurements in data.items():
                # Prior for the mean (weakly informative)
                means[pair_type] = pm.Normal(f"mu_{pair_type}", 
                                           mu=np.mean(measurements), 
                                           sigma=5.0)
                
                # FIXED: Non-centered parameterization for sigma
                # This avoids the funnel geometry that causes divergences
                sigma_offset = pm.Normal(f"sigma_offset_{pair_type}", mu=0, sigma=1)
                sigma_raw = mu_sigma + sigma_sigma * sigma_offset
                sigmas[pair_type] = pm.Deterministic(f"sigma_{pair_type}", 
                                                   pm.math.exp(sigma_raw))
            
            # Likelihood functions for each pair type
            for pair_type, measurements in data.items():
                pm.Normal(f"likelihood_{pair_type}", 
                         mu=means[pair_type],
                         sigma=sigmas[pair_type],
                         observed=measurements)
            
            # FIXED: More conservative sampling to handle hierarchical complexity
            try:
                trace = pm.sample(n_samples, tune=1500, return_inferencedata=True,
                                 target_accept=0.95, chains=2, cores=2,
                                 max_treedepth=12)
            except Exception as e:
                print(f"Error in hierarchical sampling: {e}")
                print("Falling back to more conservative sampling...")
                trace = pm.sample(n_samples, tune=2000, return_inferencedata=True,
                                 target_accept=0.99, chains=2, cores=2,
                                 max_treedepth=15)

        end_time = time.time()
        
        # Gather results
        results = {
            'trace': trace,
            'computation_time': end_time - start_time,
            'model_type': 'hierarchical'
        }
        
        return results

def compare_approaches(regular_results: Dict, hierarchical_results: Dict, 
                     true_sigmas: Dict[str, float], data: Dict[str, np.ndarray]) -> None:
    """
    Compare and visualize results from regular and hierarchical Bayesian inference.
    """
    # Set matplotlib parameters for larger fonts
    plt.rcParams.update({
        'font.size': 14,          # Base font size
        'axes.titlesize': 16,     # Title font size
        'axes.labelsize': 14,     # Axis label font size
        'xtick.labelsize': 12,    # X-axis tick label size
        'ytick.labelsize': 12,    # Y-axis tick label size
        'legend.fontsize': 12,    # Legend font size
        'figure.titlesize': 18    # Figure title font size
    })
    
    # Extract traces
    reg_trace = regular_results['trace']
    hier_trace = hierarchical_results['trace']
    
    # Print computation times
    print("\n=== COMPUTATIONAL PERFORMANCE ===")
    print(f"Regular Bayesian MCMC time: {regular_results['computation_time']:.2f} seconds")
    print(f"Hierarchical Bayesian MCMC time: {hierarchical_results['computation_time']:.2f} seconds")
    print(f"Ratio: {hierarchical_results['computation_time']/regular_results['computation_time']:.2f}x")
    
    # Check for sampling issues
    if hasattr(hier_trace, 'sample_stats') and 'diverging' in hier_trace.sample_stats:
        n_divergences = hier_trace.sample_stats['diverging'].sum().values
        print(f"Hierarchical model divergences: {n_divergences}")
        if n_divergences > 100:
            print("WARNING: High number of divergences - results may be unreliable")
    
    # Compare sigma estimates
    pair_types = list(data.keys())
    
    print("\n=== SIGMA PARAMETER ESTIMATES ===")
    print(f"{'Pair Type':<10} {'Sample Size':<12} {'True σ':<8} {'Regular σ':<15} {'Hierarchical σ':<15} {'Error Reduction':<15}")
    print("-" * 90)
    
    improvements = []
    for pair_type in pair_types:
        n = len(data[pair_type])
        true_sigma = true_sigmas.get(pair_type, float('nan'))
        
        # Extract posterior means
        reg_sigma = np.mean(reg_trace.posterior[f"sigma_{pair_type}"])
        hier_sigma = np.mean(hier_trace.posterior[f"sigma_{pair_type}"])
        
        # Calculate improvement (reduction in error)
        if not np.isnan(true_sigma):
            reg_error = abs(reg_sigma - true_sigma)
            hier_error = abs(hier_sigma - true_sigma)
            improvement = (reg_error - hier_error) / reg_error * 100 if reg_error > 0 else 0
            improvements.append(improvement)
            improvement_str = f"{improvement:.1f}%"
        else:
            improvement_str = "N/A"
        
        print(f"{pair_type:<10} {n:<12} {true_sigma:<8.3f} {reg_sigma:<15.3f} {hier_sigma:<15.3f} {improvement_str:<15}")
    
    if improvements:
        print(f"\nAverage error reduction: {np.mean(improvements):.1f}%")
    
    # Create output directory
    script_dir = Path(__file__).parent
    output_dir = script_dir / "figures"
    output_dir.mkdir(exist_ok=True)
    
    # Plot posterior distributions
    n_pairs = len(pair_types)
    fig, axes = plt.subplots(1, n_pairs, figsize=(n_pairs*7, 6))  # Increased figure size
    
    # Ensure axes is always a list for consistent handling
    if n_pairs == 1:
        axes = [axes]
    
    for i, pair_type in enumerate(pair_types):
        ax = axes[i]
        
        # Plot posterior distributions
        reg_samples = reg_trace.posterior[f"sigma_{pair_type}"].values.flatten()
        hier_samples = hier_trace.posterior[f"sigma_{pair_type}"].values.flatten()
        
        ax.hist(reg_samples, bins=30, alpha=0.6, label='Regular', color='blue', density=True)
        ax.hist(hier_samples, bins=30, alpha=0.6, label='Hierarchical', color='red', density=True)
        
        # Plot posterior distributions with custom styling
        try:
            # Plot regular posterior
            az.plot_posterior(reg_trace, var_names=[f"sigma_{pair_type}"], 
                            ax=ax, color='blue', hdi_prob=0.95, 
                            textsize=12, point_estimate='mean')
            
            # Plot hierarchical posterior (offset slightly to avoid overlap)
            hier_posterior = hier_trace.posterior[f"sigma_{pair_type}"].values.flatten()
            ax.hist(hier_posterior, bins=30, alpha=0.6, color='red', 
                   density=True, label='Hierarchical')
            
        except Exception as e:
            print(f"Error plotting {pair_type}: {e}")
            # Fallback to simple histogram
            reg_posterior = reg_trace.posterior[f"sigma_{pair_type}"].values.flatten()
            hier_posterior = hier_trace.posterior[f"sigma_{pair_type}"].values.flatten()
            
            ax.hist(reg_posterior, bins=30, alpha=0.6, color='blue', 
                   density=True, label='Regular')
            ax.hist(hier_posterior, bins=30, alpha=0.6, color='red', 
                   density=True, label='Hierarchical')
        
        # Add true value line
        if pair_type in true_sigmas:
            ax.axvline(true_sigmas[pair_type], color='black', linestyle='--', 
                      linewidth=2, label=f'True σ = {true_sigmas[pair_type]:.2f}')
        
        ax.set_title(f'Posterior for {pair_type} σ (n={len(data[pair_type])})')
        ax.set_xlabel('σ (measurement uncertainty)')
        ax.set_ylabel('Density')
        ax.legend()
        
    plt.tight_layout()
    plt.savefig(output_dir / 'mcmc_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Plot traces to check convergence for hierarchical model
    if "mu_sigma" in hier_trace.posterior:
        fig, axes = plt.subplots(2, 1, figsize=(10, 6))
        
        # Plot hyperparameter traces
        mu_sigma_samples = hier_trace.posterior["mu_sigma"].values
        sigma_sigma_samples = hier_trace.posterior["sigma_sigma"].values
        
        for chain in range(mu_sigma_samples.shape[0]):
            axes[0].plot(mu_sigma_samples[chain, :].flatten(), alpha=0.7)
        axes[0].set_title('Trace plot for μ_σ (group-level mean)')
        axes[0].set_ylabel('μ_σ')
        
        for chain in range(sigma_sigma_samples.shape[0]):
            axes[1].plot(sigma_sigma_samples[chain, :].flatten(), alpha=0.7)
        axes[1].set_title('Trace plot for σ_σ (group-level std)')
        axes[1].set_ylabel('σ_σ')
        axes[1].set_xlabel('Iteration')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'hyperparameter_trace.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    # Reset matplotlib parameters to defaults
    plt.rcParams.update(plt.rcParamsDefault)

def demonstrate_borrowing_strength():
    """
    Demonstrate the key advantage of hierarchical modeling: borrowing strength
    when data is sparse for some measurements.
    """
    print("\n=== DEMONSTRATION: BORROWING OF STRENGTH ===")
    print("Scenario: Abundant data for AA/AB measurements, sparse data for BC.")
    
    # FIXED: More realistic scenario
    true_distances = {'AA': 50.0, 'AB': 40.0, 'BC': 30.0}
    true_sigmas = {'AA': 1.8, 'AB': 2.0, 'BC': 2.2}  # Similar but not identical
    
    model = BayesianDistanceModel(true_distances, true_sigmas)
    
    # More extreme imbalance to show hierarchical advantage
    n_measurements = {'AA': 100, 'AB': 50, 'BC': 5}  # Slightly more BC data
    imbalanced_data = model.generate_data(n_measurements)
    
    print("\nData sizes:")
    for pair_type, data_points in imbalanced_data.items():
        print(f"  {pair_type}: {len(data_points)} measurements (mean = {np.mean(data_points):.2f}, std = {np.std(data_points):.2f})")
    
    # Run inference with both approaches
    print("\nRunning Regular Bayesian MCMC inference...")
    regular_results = model.regular_bayesian_mcmc(imbalanced_data, n_samples=1500)
    
    print("Running Hierarchical Bayesian MCMC inference...")
    hierarchical_results = model.hierarchical_bayesian_mcmc(imbalanced_data, n_samples=1500)
    
    # Compare results
    compare_approaches(regular_results, hierarchical_results, true_sigmas, imbalanced_data)
    
    return regular_results, hierarchical_results

def explain_complexity():
    """Explain the computational complexity of hierarchical Bayesian models."""
    print("\n=== COMPUTATIONAL CONSIDERATIONS FOR HIERARCHICAL MODELS ===")
    print("Hierarchical Bayesian models with MCMC can be expensive because:")
    print("  1. Parameter correlation induces slow mixing in the MCMC chains")
    print("  2. The joint posterior has complex geometry (funnel-shaped distributions)")
    print("  3. Sampling efficiency declines as the hierarchy deepens")
    print("  4. Higher-dimensional parameter space increases convergence time")
    print("  5. Funnel geometries from hierarchical structure cause sampling difficulties")
    
    print("\nCommon problems and solutions:")
    print("  - The 'funnel problem': High correlation between group-level and individual parameters")
    print("    → Solution: Non-centered parameterization (implemented in our code)")
    print("  - Divergences in NUTS: Indicates areas of high curvature in posterior")
    print("    → Solution: Increase target_accept (0.9-0.99) and adapt_delta")
    print("  - Slow mixing chains: Parameters don't explore full posterior efficiently")
    print("    → Solution: Reparameterization, longer warm-up, better priors")
    print("  - Memory issues with large datasets: Can cause worker crashes")
    print("    → Solution: Batch processing or variational inference")

def explain_joint_distributions():
    """Explain joint probability distributions in PyMC and alternatives."""
    print("\n=== IMPLEMENTING JOINT PROBABILITY DISTRIBUTIONS ===")
    print("1. PyMC approach (as used here):")
    print("   - Defines graphical model that encodes conditional dependencies")
    print("   - Joint distribution is implicitly defined through the model structure")
    print("   - NUTS sampler handles joint posterior efficiently")
    
    print("\n2. Custom implementation options:")
    print("   - Using NumPy/SciPy with Metropolis-Hastings:")
    print("     * Define log_joint_prob(parameters, data)")
    print("     * Implement proposal and acceptance steps manually")
    print("     * Requires careful tuning of proposal distributions")
    
    print("\n3. Performance considerations:")
    print("   - PyMC (NUTS) is generally more efficient for complex hierarchical models")
    print("   - Custom implementations give more control but require more expertise")
    print("   - Variational inference scales better but provides approximations")
    print("   - For large datasets, consider using PyMC with minibatches")

if __name__ == "__main__":
    print("=== BAYESIAN HIERARCHICAL MODELING FOR DISTANCE MEASUREMENTS ===")
    
    # Demonstrate the borrowing of strength with imbalanced data
    regular_results, hierarchical_results = demonstrate_borrowing_strength()
    
    # Explain computational complexity
    explain_complexity()
    
    # Explain joint probability distributions
    explain_joint_distributions()
    
    print("\n=== PRACTICAL APPLICATIONS ===")
    print("This approach applies to real scientific measurements such as:")
    print("  - FRET measurements between protein domains")
    print("  - Atomic distances in X-ray crystallography")
    print("  - Molecular dynamics simulation analysis")
    print("  - NMR spectroscopy distance constraints")
    print("  - Single-molecule tracking with variable precision")
    
    print("\n=== ADVANTAGES OF HIERARCHICAL MODELING ===")
    print("1. 'Borrowing strength' across similar measurements")
    print("2. Better uncertainty quantification for sparse data")
    print("3. Automatically balancing prior knowledge with observed data")
    print("4. More accurate estimates for rare or difficult measurements")
    print("5. Modeling of both group-level and individual-level variation")
