import os
import re
import json
import random
import numpy as np
from typing import Dict, Optional, List, Tuple

class GMMSigmaProvider:
    """
    Provides sigma values from GMM fits and calculates negative log priors.
    Works with both PairSampler and TetramerSampler outputs.
    """
    def __init__(self, sampler_name: str = "PairSampler", output_dir: str = None, 
                 specific_chain: int = None, sigma_ranges: Dict[str, Tuple[float, float]] = None):
        """
        Initialize the sigma provider for a specific sampler type.
        
        Args:
            sampler_name: Name of the sampler ('PairSampler', 'TetramerSampler', 'OctetSampler')
            output_dir: Directory containing GMM files (default: output-analysis/{source_sampler}_results)
            specific_chain: Load from a specific chain number (default: random selection)
            sigma_ranges: Dictionary of (min, max) ranges for each sigma type (default: [0.5, 10.0] for all)
        """
        self.sampler_name = sampler_name
        self.pair_types = ["AA", "AB", "BC"]
        
        # Set output directory
        if output_dir is None:
            if sampler_name.lower() == "tetramersampler":
                source_sampler = "pairsampler"
            elif sampler_name.lower() == "octetsampler":
                source_sampler = "tetramersampler"
            else:
                source_sampler = "pairsampler"
            self.output_dir = os.path.join(os.getcwd(), f"output_analysis/{source_sampler}_results")
        else:
            self.output_dir = output_dir
            
        self.specific_chain = specific_chain
        self.gmm_params = self._load_gmm_parameters()
        
        # Load default values from SystemParameters (simplified version)
        self.default_sigma = {
            "AA": 2.0,
            "AB": 2.0,
            "BC": 2.0,
        }
        
        # Set sigma ranges
        if sigma_ranges is None:
            self.sigma_ranges = {pt: (0.5, 10.0) for pt in self.pair_types}
        else:
            self.sigma_ranges = sigma_ranges
        
    def _load_gmm_parameters(self) -> Dict:
        """Load GMM parameters from JSON files."""
        gmm_params = {sigma_type: None for sigma_type in self.pair_types}
        
        if not os.path.exists(self.output_dir):
            print(f"Output directory not found: {self.output_dir}")
            return gmm_params
        else:
            print("output_directory is: ", self.output_dir)
        
        # Find all GMM files for the sampler
        pattern = re.compile(f"gmm_fit_(AA|AB|BC|CC)_chain_(\\d+)\\.json")
        matches = []
        
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                matches.append((match.group(1), int(match.group(2))))
        
        if not matches:
            print(f"No GMM fit files found in {self.output_dir}")
            return gmm_params
        
        # Get all available chain numbers
        chain_numbers = set(chain_num for _, chain_num in matches)
        if not chain_numbers:
            return gmm_params
        
        # Select chain: either specified or random
        selected_chain = self.specific_chain if self.specific_chain is not None else random.choice(list(chain_numbers))
        print(f"Using GMM parameters from chain {selected_chain}")
        
        # Load GMM files for the selected chain
        file_paths = {
            sigma_type: os.path.join(self.output_dir, f"gmm_fit_{sigma_type}_chain_{selected_chain}.json")
            for sigma_type in gmm_params
        }
        
        for sigma_type, file_path in file_paths.items():
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_params[sigma_type] = json.load(f)
                        print(f"Loaded GMM parameters for {sigma_type}")
                except Exception as e:
                    print(f"Error loading {file_path}: {e}")
        
        return gmm_params

    def sample_sigma_values(self, max_attempts: int = 100) -> Dict[str, float]:
        """Sample sigma values from the loaded GMMs with bounds checking."""
        sigma = {}
        
        for pair_type in self.pair_types:
            gmm_info = self.gmm_params.get(pair_type)
            sampled_value = None
            
            if (
                gmm_info
                and all(k in gmm_info for k in ['n_components', 'means', 'covariances', 'weights'])
            ):
                # Try multiple times to get a valid sample
                for attempt in range(max_attempts):
                    try:
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
                        
                        # Limit the standard deviation to prevent extreme values
                        std_value = min(np.sqrt(abs(cov_value)), mean_value[0] * 0.3)  # Cap at 30% of mean
                        
                        # Sample from normal distribution
                        candidate = np.random.normal(loc=mean_value[0], scale=std_value)
                        candidate = float(candidate)
                        
                        # Check if within valid range
                        if pair_type in self.sigma_ranges:
                            min_val, max_val = self.sigma_ranges[pair_type]
                            if min_val <= candidate <= max_val:
                                sampled_value = candidate
                                break
                        
                    except Exception as e:
                        continue
                
                if sampled_value is None:
                    print(f"Failed to sample valid {pair_type} after {max_attempts} attempts. Using default.")
                    sampled_value = self.default_sigma[pair_type]
            else:
                sampled_value = self.default_sigma[pair_type]
            
            sigma[pair_type] = sampled_value
        
        # Final validation
        for pair_type, value in sigma.items():
            if not (0.1 <= value <= 15.0):  # Hard bounds
                print(f"Warning: {pair_type}={value:.3f} is outside safe bounds. Using default.")
                sigma[pair_type] = self.default_sigma[pair_type]
        
        return sigma
    
    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior across all sigma values with range enforcement."""
        if not sigma:
            return np.inf
        
        # Check if sigma values are within their specified ranges
        for pair_type, value in sigma.items():
            if pair_type in self.sigma_ranges:
                min_val, max_val = self.sigma_ranges[pair_type]
                if value < min_val or value > max_val:
                    return np.inf  # Reject out-of-range values
        
        log_priors = [self._calculate_gmm_log_prob(value, pair_type) 
                      for pair_type, value in sigma.items()]
        
        # Filter out -inf values to avoid returning inf unnecessarily
        valid_priors = [lp for lp in log_priors if lp > -np.inf]
        
        if not valid_priors:
            return 100000.0  # Large but finite penalty for all invalid priors
            
        return -np.sum(valid_priors)
    
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
            
        except (KeyError, ValueError, np.linalg.LinAlgError) as e:
            print(f"Error calculating GMM log probability for {pair_type}: {e}")
            return -np.inf
    
    def get_available_chains(self) -> List[int]:
        """Get list of available chains with GMM fits."""
        if not os.path.exists(self.output_dir):
            return []
            
        # Change the pattern to match any sigma type (AA, AB, BC, or CC)
        pattern = re.compile(r"gmm_fit_.*_chain_(\d+)\.json")
        chain_numbers = []
        
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                chain_numbers.append(int(match.group(1)))
                
        return sorted(set(chain_numbers))

# Example usage function
def sample_sigmas_from_gmm(sampler_name="PairSampler", output_dir=None, chain=None):
    """
    Sample sigma values from GMM fits for a specific sampler.
    
    Args:
        sampler_name: Name of the sampler ('PairSampler', 'TetramerSampler')
        output_dir: Directory containing GMM files
        chain: Specific chain to use (default: random selection)
        
    Returns:
        dict: Sampled sigma values
    """
    provider = GMMSigmaProvider(sampler_name, output_dir, chain)
    sigma_values = provider.sample_sigma_values()
    
    # Print the sampled values
    print(f"\nSampled sigma values for {sampler_name}:")
    for key, value in sigma_values.items():
        print(f"  {key}: {value:.4f}")
    
    # Calculate and print the negative log prior
    neg_log_prior = provider.calculate_negative_log_prior(sigma_values)
    print(f"Negative log prior: {neg_log_prior:.4f}")
    
    return sigma_values, neg_log_prior

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Sample sigma values from GMM fits')
    parser.add_argument('--sampler', type=str, default='PairSampler',
                        choices=['PairSampler', 'TetramerSampler', 'OctetSampler'],
                        help='Sampler name (PairSampler, TetramerSampler, or OctetSampler)')
    parser.add_argument('--dir', type=str, default=None,
                        help='Directory containing GMM files')
    parser.add_argument('--chain', type=int, default=None,
                        help='Specific chain to use (default: random)')
    parser.add_argument('--list-chains', action='store_true',
                        help='List available chains and exit')
    
    args = parser.parse_args()
    
    if args.list_chains:
        provider = GMMSigmaProvider(args.sampler, args.dir)
        chains = provider.get_available_chains()
        if chains:
            print(f"Available chains for {args.sampler}:")
            for chain in chains:
                print(f"  Chain {chain}")
        else:
            print(f"No chains found for {args.sampler}")
    else:
        sigma_values, neg_log_prior = sample_sigmas_from_gmm(args.sampler, args.dir, args.chain)