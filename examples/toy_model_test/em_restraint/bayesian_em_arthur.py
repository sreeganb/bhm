#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bayesian EM Restraint Scoring Function
Converts cross-correlation coefficient to negative log-likelihood for MCMC sampling
"""

import numpy as np
import math
from typing import Optional, Tuple
import logging

class BayesianEMRestraint:
    """
    A Bayesian restraint that converts cross-correlation coefficient to negative log-likelihood
    for use in MCMC simulations.
    """
    
    def __init__(self, 
                 target_density_map,
                 resolution: float,
                 sigma_noise: float = 0.1,
                 correlation_prior_mean: float = 0.8,
                 correlation_prior_std: float = 0.2,
                 backend: str = 'cpu'):
        """
        Initialize the Bayesian EM restraint.
        
        Args:
            target_density_map: Pre-loaded MRC file object
            resolution: Resolution for Gaussian blurring (Angstroms)
            sigma_noise: Noise level in the correlation measurement
            correlation_prior_mean: Prior mean for expected correlation
            correlation_prior_std: Prior standard deviation for correlation
            backend: 'cpu' or 'gpu'
        """
        self.target_density_map = target_density_map
        self.resolution = resolution
        self.sigma_noise = sigma_noise
        self.correlation_prior_mean = correlation_prior_mean
        self.correlation_prior_std = correlation_prior_std
        self.backend = backend
        
        # Pre-compute the experimental data normalization
        self._precompute_experimental_normalization()
        
    def _precompute_experimental_normalization(self):
        """Pre-compute the experimental data mean and variance for efficiency."""
        exp_data = self.target_density_map.data.flatten()
        self.exp_mean = np.mean(exp_data)
        self.exp_centered = exp_data - self.exp_mean
        self.exp_norm_sq = np.sum(self.exp_centered**2)
        
    def calculate_correlation_coefficient(self, sphere_coords: np.ndarray, 
                                        sphere_radii: np.ndarray) -> float:
        """
        Calculate the cross-correlation coefficient between model and experimental data.
        
        Args:
            sphere_coords: Array of shape (N, 3) for sphere coordinates
            sphere_radii: Array of shape (N,) for sphere radii
            
        Returns:
            Cross-correlation coefficient
        """
        from working_example import calculate_ccc_score
        
        try:
            ccc = calculate_ccc_score(
                sphere_coords, 
                sphere_radii, 
                self.target_density_map, 
                self.resolution, 
                backend=self.backend
            )
            return ccc
        except Exception as e:
            logging.warning(f"Error calculating CCC: {e}, returning low correlation")
            return -1.0
    
    def correlation_to_likelihood(self, correlation: float) -> float:
        """
        Convert correlation coefficient to likelihood.
        
        Uses a truncated normal distribution to model the likelihood of observing
        the correlation given the true underlying correlation.
        
        Args:
            correlation: Observed correlation coefficient
            
        Returns:
            Likelihood value
        """
        # Clamp correlation to valid range [-1, 1]
        correlation = np.clip(correlation, -1.0, 1.0)
        
        # Model: observed_correlation ~ Normal(true_correlation, sigma_noise)
        # For simplicity, assume true_correlation = correlation_prior_mean
        likelihood = np.exp(-0.5 * ((correlation - self.correlation_prior_mean) / self.sigma_noise)**2)
        
        # Normalize by the noise standard deviation
        likelihood /= (self.sigma_noise * np.sqrt(2 * np.pi))
        
        return likelihood
    
    def correlation_prior_log_prob(self, correlation: float) -> float:
        """
        Calculate the log prior probability for the correlation coefficient.
        
        Args:
            correlation: Correlation coefficient
            
        Returns:
            Log prior probability
        """
        # Clamp correlation to valid range
        correlation = np.clip(correlation, -1.0, 1.0)
        
        # Use a truncated normal prior
        log_prior = -0.5 * ((correlation - self.correlation_prior_mean) / self.correlation_prior_std)**2
        log_prior -= np.log(self.correlation_prior_std * np.sqrt(2 * np.pi))
        
        return log_prior
    
    def calculate_negative_log_likelihood(self, sphere_coords: np.ndarray, 
                                        sphere_radii: np.ndarray,
                                        include_prior: bool = True) -> float:
        """
        Calculate the negative log-likelihood for the EM restraint.
        
        This is the main function to use in your MCMC sampler.
        
        Args:
            sphere_coords: Array of shape (N, 3) for sphere coordinates
            sphere_radii: Array of shape (N,) for sphere radii  
            include_prior: Whether to include the correlation prior
            
        Returns:
            Negative log-likelihood value (lower is better)
        """
        # Calculate correlation coefficient
        correlation = self.calculate_correlation_coefficient(sphere_coords, sphere_radii)
        
        # Calculate likelihood
        likelihood = self.correlation_to_likelihood(correlation)
        
        # Avoid log(0) by adding a small epsilon
        epsilon = 1e-10
        neg_log_likelihood = -np.log(likelihood + epsilon)
        
        # Add prior term if requested
        if include_prior:
            log_prior = self.correlation_prior_log_prob(correlation)
            neg_log_likelihood -= log_prior
        
        return neg_log_likelihood
    
    def calculate_log_probability(self, sphere_coords: np.ndarray,
                                sphere_radii: np.ndarray) -> Tuple[float, dict]:
        """
        Calculate the log probability and return detailed information.
        
        Args:
            sphere_coords: Array of shape (N, 3) for sphere coordinates
            sphere_radii: Array of shape (N,) for sphere radii
            
        Returns:
            Tuple of (log_probability, info_dict)
        """
        correlation = self.calculate_correlation_coefficient(sphere_coords, sphere_radii)
        likelihood = self.correlation_to_likelihood(correlation)
        log_prior = self.correlation_prior_log_prob(correlation)
        
        log_probability = np.log(likelihood + 1e-10) + log_prior
        
        info = {
            'correlation': correlation,
            'likelihood': likelihood,
            'log_prior': log_prior,
            'log_likelihood': np.log(likelihood + 1e-10),
            'log_probability': log_probability
        }
        
        return log_probability, info


# =====================================================================
# EXAMPLE USAGE AND TESTING
# =====================================================================
def test_bayesian_em_restraint():
    """Test the Bayesian EM restraint with the toy model."""
    
    # Import your working example functions
    from working_example import (
        create_dummy_map_from_model, parse_density
    )
    
    print("=== Testing Bayesian EM Restraint ===")
    
    # 1. Set up toy model data (same as your working example)
    array_A = np.array([[ 63.,0.,0.],[ 44.55,44.55,0.],[0.,63.,0.],[-44.55,44.55,0.],[-63.,0.,0.],[-44.55,-44.55,0.],[ -0.,-63.,0.],[ 44.55,-44.55,0.]])
    array_B = np.array([[ 63.,0.,-38.5],[ 44.55,44.55,-38.5],[0.,63.,-38.5],[-44.55,44.55,-38.5],[-63.,0.,-38.5],[-44.55,-44.55,-38.5],[-0.,-63.,-38.5],[ 44.55,-44.55,-38.5]]) 
    array_C = np.array([[ 47.,0.,-68.5],[ 79.,0.,-68.5],[ 55.86,55.86,-68.5],[ 33.23,33.23,-68.5],[0.,47.,-68.5],[0.,79.,-68.5],[-55.86,55.86,-68.5],[-33.23,33.23,-68.5],[-47.,0.,-68.5],[-79.,0.,-68.5],[-55.86,-55.86,-68.5],[-33.23,-33.23,-68.5],[0.,-47.,-68.5],[0.,-79.,-68.5],[ 55.86,-55.86,-68.5],[ 33.23,-33.23,-68.5]])

    radius_A, radius_B, radius_C = 24.0, 14.0, 16.0
    copy_numbers = {'A': 8, 'B': 8, 'C': 16}

    ideal_coords = np.vstack([array_A, array_B, array_C])
    ideal_radii = np.repeat([radius_A, radius_B, radius_C], 
                           repeats=[copy_numbers['A'], copy_numbers['B'], copy_numbers['C']])
    
    # 2. Create and load target map
    RESOLUTION = 20.0
    VOXEL_SIZE = 5.0
    BOX_SIZE = 250.0
    TARGET_MAP_FILE = "test_target_map.mrc"
    
    create_dummy_map_from_model(ideal_coords, ideal_radii, RESOLUTION, 
                               VOXEL_SIZE, BOX_SIZE, filename=TARGET_MAP_FILE)
    target_map = parse_density(TARGET_MAP_FILE)
    
    # 3. Initialize Bayesian restraint
    restraint = BayesianEMRestraint(
        target_density_map=target_map,
        resolution=RESOLUTION,
        sigma_noise=0.05,  # Low noise for synthetic data
        correlation_prior_mean=0.9,  # Expect high correlation
        correlation_prior_std=0.1,   # Fairly confident
        backend='cpu'
    )
    
    # 4. Test with ideal model (should have high probability, low NLL)
    print("\n--- Testing ideal model ---")
    nll_ideal = restraint.calculate_negative_log_likelihood(ideal_coords, ideal_radii)
    log_prob_ideal, info_ideal = restraint.calculate_log_probability(ideal_coords, ideal_radii)
    
    print(f"Ideal model:")
    print(f"  Correlation: {info_ideal['correlation']:.6f}")
    print(f"  Negative Log-Likelihood: {nll_ideal:.6f}")
    print(f"  Log Probability: {log_prob_ideal:.6f}")
    
    # 5. Test with perturbed model (should have lower probability, higher NLL)
    print("\n--- Testing perturbed model ---")
    np.random.seed(42)  # For reproducibility
    perturbed_coords = ideal_coords + np.random.normal(0, 5.0, ideal_coords.shape)
    
    nll_perturbed = restraint.calculate_negative_log_likelihood(perturbed_coords, ideal_radii)
    log_prob_perturbed, info_perturbed = restraint.calculate_log_probability(perturbed_coords, ideal_radii)
    
    print(f"Perturbed model:")
    print(f"  Correlation: {info_perturbed['correlation']:.6f}")
    print(f"  Negative Log-Likelihood: {nll_perturbed:.6f}")
    print(f"  Log Probability: {log_prob_perturbed:.6f}")
    
    # 6. Demonstrate the scoring function behavior
    print(f"\n--- Scoring Summary ---")
    print(f"Ideal model NLL: {nll_ideal:.6f} (lower is better)")
    print(f"Perturbed model NLL: {nll_perturbed:.6f}")
    print(f"Difference: {nll_perturbed - nll_ideal:.6f}")
    print(f"Ratio: {nll_perturbed / nll_ideal:.2f}x")
    
    target_map.close()
    
    return restraint


if __name__ == "__main__":
    # Run the test
    restraint = test_bayesian_em_restraint()
    
    print("\n=== Usage in MCMC ===")
    print("In your MCMC sampler, call:")
    print("nll = restraint.calculate_negative_log_likelihood(coords, radii)")
    print("Then add this NLL to your total energy function.")