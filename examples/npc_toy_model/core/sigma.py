import os
import re
import json
import random
import numpy as np
from typing import Dict, List, Tuple

class GMMSigmaProvider:
    """
    Provides sigma values from GMM fits and calculates negative log priors.
    Works with PairSampler, TetramerSampler, and OctetSampler outputs.
    """
    # Define sampler source mappings
    SAMPLER_SOURCES = {
        "PairSampler": "pairsampler",
        "TetramerSampler": "pairsampler",
        "OctetSampler": "tetramersampler"
    }
    
    # Default sigma values if GMM not available
    DEFAULT_SIGMA = {
        "AA": 2.0, 
        "AB": 2.0, 
        "BC": 2.0,
        "CC": 2.0
    }
    
    def __init__(self, sampler_name: str = "PairSampler"):
        """
        Initialize the sigma provider for a specific sampler type.
        
        Args:
            sampler_name: Name of the sampler ('PairSampler', 'TetramerSampler', 'OctetSampler')
        """
        self.sampler_name = sampler_name
        self.pair_types = ["AA", "AB", "BC"]
        
        # Determine source folder for GMM files
        source_sampler = self.SAMPLER_SOURCES.get(sampler_name, "pairsampler")
        self.output_dir = os.path.join(os.getcwd(), f"output_analysis/{source_sampler}_results")
        
        # Load GMM parameters
        self.gmm_params = self._load_gmm_parameters()
        
    def _load_gmm_parameters(self) -> Dict:
        """Load GMM parameters from JSON files for a randomly selected chain."""
        gmm_params = {sigma_type: None for sigma_type in self.pair_types}
        
        if not os.path.exists(self.output_dir):
            print(f"Output directory not found: {self.output_dir}")
            return gmm_params
        
        # Find all GMM files
        pattern = re.compile(r"gmm_fit_(AA|AB|BC)_chain_(\d+)\.json")
        chain_files = {}
        
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                sigma_type, chain = match.group(1), int(match.group(2))
                if chain not in chain_files:
                    chain_files[chain] = []
                chain_files[chain].append(sigma_type)
        
        if not chain_files:
            print(f"No GMM fit files found in {self.output_dir}")
            return gmm_params
        
        # Select random chain with the most sigma types
        chains_by_completeness = sorted(
            chain_files.keys(), 
            key=lambda c: len(chain_files[c]), 
            reverse=True
        )
        selected_chain = chains_by_completeness[0]
        print(f"Using GMM parameters from chain {selected_chain}")
        
        # Load GMM files for the selected chain
        for sigma_type in self.pair_types:
            file_path = os.path.join(self.output_dir, f"gmm_fit_{sigma_type}_chain_{selected_chain}.json")
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_params[sigma_type] = json.load(f)
                        print(f"Loaded GMM parameters for {sigma_type}")
                except Exception as e:
                    print(f"Error loading {file_path}: {e}")
        
        return gmm_params

    def sample_sigma_values(self) -> Dict[str, float]:
        """Sample sigma values from the loaded GMMs or use defaults."""
        sigma = {}
        
        for pair_type in self.pair_types:
            gmm_info = self.gmm_params.get(pair_type)
            
            # Use default if no valid GMM data
            if not gmm_info or not all(k in gmm_info for k in ['n_components', 'means', 'covariances', 'weights']):
                sigma[pair_type] = self.DEFAULT_SIGMA[pair_type]
                continue
                
            try:
                # Extract GMM parameters
                n_components = gmm_info['n_components']
                means = np.asarray(gmm_info['means']).reshape(n_components, 1)
                covariances = np.asarray(gmm_info['covariances'])
                weights = np.asarray(gmm_info['weights'])
                
                # Normalize weights
                weights = weights / np.sum(weights)
                
                # Sample from GMM
                component = np.random.choice(n_components, p=weights)
                mean_value = means[component].flatten()
                cov_value = float(covariances[component])
                
                # Sample from normal distribution
                sampled_value = np.random.normal(loc=mean_value, scale=np.sqrt(cov_value))
                sigma[pair_type] = float(sampled_value)
            except Exception as e:
                print(f"Error sampling GMM for {pair_type}: {e}. Using default.")
                sigma[pair_type] = self.DEFAULT_SIGMA[pair_type]
        
        return sigma
    
    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior across all sigma values."""
        if not sigma:
            return np.inf
            
        log_priors = []
        for pair_type, value in sigma.items():
            log_prob = self._calculate_gmm_log_prob(value, pair_type)
            if log_prob > -np.inf:
                log_priors.append(log_prob)
        
        if not log_priors:
            return 100000.0  # Large penalty for all invalid priors
            
        return -np.sum(log_priors)
    
    def _calculate_gmm_log_prob(self, sigma_value: float, pair_type: str) -> float:
        """Calculate log probability for a sigma value given a GMM."""
        if pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
            return -np.inf
            
        gmm_info = self.gmm_params[pair_type]
        
        try:
            n_components = gmm_info['n_components']
            means = np.asarray(gmm_info['means']).reshape(n_components, 1)
            covariances = np.asarray(gmm_info['covariances']).reshape(n_components, 1, 1) 
            weights = np.asarray(gmm_info['weights'])
            
            # Ensure numerical stability
            covariances = np.maximum(covariances, 1e-10)
            
            # Vectorized component log probability calculation
            x = np.array(sigma_value).reshape(1)
            diff = x - means.flatten()
            exponents = -0.5 * (diff**2 / covariances.flatten())
            norms = np.log(weights) - 0.5 * np.log(2 * np.pi * covariances.flatten())
            component_log_probs = norms + exponents
            
            # Numerically stable log-sum-exp
            max_log_prob = np.max(component_log_probs)
            return max_log_prob + np.log(np.sum(np.exp(component_log_probs - max_log_prob)))
            
        except Exception as e:
            print(f"Error calculating GMM log probability for {pair_type}: {e}")
            return -np.inf

# Simple function for sampling sigma values from a specified sampler
def get_sigma_values(sampler_name="PairSampler"):
    """
    Get sigma values for the specified sampler.
    
    Args:
        sampler_name: 'PairSampler', 'TetramerSampler', or 'OctetSampler'
        
    Returns:
        Tuple: (sigma_values, negative_log_prior)
    """
    provider = GMMSigmaProvider(sampler_name)
    sigma_values = provider.sample_sigma_values()
    neg_log_prior = provider.calculate_negative_log_prior(sigma_values)
    
    print(f"\nSampled sigma values for {sampler_name}:")
    for key, value in sigma_values.items():
        print(f"  {key}: {value:.4f}")
    print(f"Negative log prior: {neg_log_prior:.4f}")
    
    return sigma_values, neg_log_prior

# Example usage
if __name__ == "__main__":
    import sys
    sampler = "PairSampler" if len(sys.argv) < 2 else sys.argv[1]
    get_sigma_values(sampler)