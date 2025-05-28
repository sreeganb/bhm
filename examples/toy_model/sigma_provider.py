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
        self._print_gmm_stats()
        
        # Load default values
        self.default_sigma = {"AA": 2.0, "AB": 2.0, "BC": 2.0}
        
        # Set sigma ranges
        self.sigma_ranges = sigma_ranges if sigma_ranges is not None else {pt: (0.5, 10.0) for pt in self.pair_types}
        
        # Configurable standard deviation cap factor (default 0.3)
        self.std_cap_factor = 0.3

    def _load_gmm_parameters(self) -> Dict:
        """Load GMM parameters from JSON files."""
        gmm_params = {sigma_type: None for sigma_type in self.pair_types}
        
        if not os.path.exists(self.output_dir):
            print(f"Output directory not found: {self.output_dir}")
            return gmm_params
        print("output_directory is: ", self.output_dir)
        
        pattern = re.compile(f"gmm_fit_(AA|AB|BC|CC)_chain_(\\d+)\\.json")
        matches = [(match.group(1), int(match.group(2))) 
                   for filename in os.listdir(self.output_dir) 
                   if (match := pattern.match(filename))]
        
        if not matches:
            print(f"No GMM fit files found in {self.output_dir}")
            return gmm_params
        
        chain_numbers = set(chain_num for _, chain_num in matches)
        selected_chain = self.specific_chain if self.specific_chain is not None else random.choice(list(chain_numbers))
        print(f"Using GMM parameters from chain {selected_chain}")
        
        file_paths = {sigma_type: os.path.join(self.output_dir, f"gmm_fit_{sigma_type}_chain_{selected_chain}.json")
                      for sigma_type in gmm_params}
        
        for sigma_type, file_path in file_paths.items():
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_params[sigma_type] = json.load(f)
                        print(f"Loaded GMM parameters for {sigma_type}")
                except Exception as e:
                    print(f"Error loading {file_path}: {e}")
        
        return gmm_params

    def _print_gmm_stats(self):
        """Print mean and standard deviation of loaded GMMs for diagnostics."""
        for pair_type, gmm_info in self.gmm_params.items():
            if gmm_info and all(k in gmm_info for k in ['means', 'covariances', 'weights']):
                try:
                    means = np.array(gmm_info['means'])
                    covariances = np.array(gmm_info['covariances'])
                    weights = np.array(gmm_info['weights']) / np.sum(gmm_info['weights'])
                    mean = np.sum(means * weights)
                    variance = np.sum(weights * (means - mean)**2) + np.sum(weights * covariances)
                    std = np.sqrt(variance)
                    print(f"{pair_type} GMM: mean={mean:.3f}, std={std:.3f}")
                except Exception as e:
                    print(f"Error computing stats for {pair_type}: {e}")

    def sample_sigma_values(self, max_attempts: int = 100) -> Dict[str, float]:
        """Sample sigma values from GMMs with bounds checking."""
        sigma = {}
        
        for pair_type in self.pair_types:
            gmm_info = self.gmm_params.get(pair_type)
            sampled_value = None
            
            if gmm_info and all(k in gmm_info for k in ['n_components', 'means', 'covariances', 'weights']):
                min_val, max_val = self.sigma_ranges[pair_type]
                for attempt in range(max_attempts):
                    try:
                        n_components = gmm_info['n_components']
                        means = np.asarray(gmm_info['means']).reshape(n_components, 1)
                        covariances = np.asarray(gmm_info['covariances'])
                        weights = np.asarray(gmm_info['weights']) / np.sum(gmm_info['weights'])
                        
                        component = np.random.choice(n_components, p=weights)
                        mean_value = means[component].item()
                        cov_value = float(covariances[component])
                        std_value = min(np.sqrt(abs(cov_value)), mean_value * self.std_cap_factor)
                        
                        candidate = np.random.normal(loc=mean_value, scale=std_value)
                        if min_val <= candidate <= max_val:
                            sampled_value = candidate
                            break
                        if attempt == max_attempts - 1:
                            print(f"{pair_type}: Failed after {max_attempts} attempts. Last candidate={candidate:.3f}")
                    except Exception as e:
                        continue
                
                sampled_value = sampled_value if sampled_value is not None else self.default_sigma[pair_type]
            else:
                sampled_value = self.default_sigma[pair_type]
            
            sigma[pair_type] = sampled_value
        
        # Enforce sigma_ranges as hard bounds (previously [0.1, 15.0])
        for pair_type, value in sigma.items():
            min_val, max_val = self.sigma_ranges[pair_type]
            if not (min_val <= value <= max_val):
                print(f"{pair_type}={value:.3f} outside {pair_type} range [{min_val}, {max_val}]. Using default.")
                sigma[pair_type] = self.default_sigma[pair_type]
        
        return sigma
    
    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior with strict range enforcement."""
        if not sigma:
            return np.inf
        
        for pair_type, value in sigma.items():
            if pair_type in self.sigma_ranges:
                min_val, max_val = self.sigma_ranges[pair_type]
                if value < min_val or value > max_val:
                    print(f"{pair_type}={value:.3f} outside range [{min_val}, {max_val}]")
                    return np.inf
        
        log_priors = [self._calculate_gmm_log_prob(value, pair_type) for pair_type, value in sigma.items()]
        valid_priors = [lp for lp in log_priors if lp > -np.inf]
        
        if not valid_priors:
            print(f"All sigmas have zero prior probability: {sigma}")
            return np.inf  # Strict rejection
        
        total_neg_log_prior = -np.sum(valid_priors)
        if total_neg_log_prior > 1000:  # Threshold for monitoring
            print(f"Large negative log prior: {total_neg_log_prior:.3f} for {sigma}")
        return total_neg_log_prior
    
    def _calculate_gmm_log_prob(self, sigma_value: float, pair_type: str) -> float:
        """Calculate log probability for a sigma value under its GMM."""
        if pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
            return -np.inf
        
        gmm_info = self.gmm_params[pair_type]
        try:
            n_components = gmm_info['n_components']
            means = np.asarray(gmm_info['means']).flatten()
            covariances = np.asarray(gmm_info['covariances']).flatten()
            weights = np.asarray(gmm_info['weights']) / np.sum(gmm_info['weights'])
            
            covariances = np.maximum(covariances, 1e-10)
            diff = sigma_value - means
            exponents = -0.5 * (diff**2 / covariances)
            norms = np.log(weights) - 0.5 * np.log(2 * np.pi * covariances)
            component_log_probs = norms + exponents
            
            max_log_prob = np.max(component_log_probs)
            log_sum_exp = max_log_prob + np.log(np.sum(np.exp(component_log_probs - max_log_prob)))
            return log_sum_exp if np.isfinite(log_sum_exp) else -np.inf
        except Exception as e:
            print(f"Error in GMM log prob for {pair_type}: {e}")
            return -np.inf
    
    def get_available_chains(self) -> List[int]:
        """Get list of available chains with GMM fits."""
        if not os.path.exists(self.output_dir):
            return []
        
        pattern = re.compile(r"gmm_fit_.*_chain_(\d+)\.json")
        chain_numbers = [int(match.group(1)) for filename in os.listdir(self.output_dir) 
                         if (match := pattern.match(filename))]
        return sorted(set(chain_numbers))

def sample_sigmas_from_gmm(sampler_name="PairSampler", output_dir=None, chain=None):
    """Sample sigma values from GMM fits."""
    provider = GMMSigmaProvider(sampler_name, output_dir, chain)
    sigma_values = provider.sample_sigma_values()
    
    print(f"\nSampled sigma values for {sampler_name}:")
    for key, value in sigma_values.items():
        print(f"  {key}: {value:.4f}")
    
    neg_log_prior = provider.calculate_negative_log_prior(sigma_values)
    print(f"Negative log prior: {neg_log_prior:.4f}")
    
    return sigma_values, neg_log_prior

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Sample sigma values from GMM fits')
    parser.add_argument('--sampler', type=str, default='PairSampler',
                        choices=['PairSampler', 'TetramerSampler', 'OctetSampler'])
    parser.add_argument('--dir', type=str, default=None)
    parser.add_argument('--chain', type=int, default=None)
    parser.add_argument('--list-chains', action='store_true')
    
    args = parser.parse_args()
    
    if args.list_chains:
        provider = GMMSigmaProvider(args.sampler, args.dir)
        chains = provider.get_available_chains()
        print(f"Available chains for {args.sampler}: {'None' if not chains else ', '.join(map(str, chains))}")
    else:
        sample_sigmas_from_gmm(args.sampler, args.dir, args.chain)