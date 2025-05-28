import os
import re
import json
import random
import numpy as np
from typing import Dict, Optional, List, Tuple

class SamplerSequenceManager:
    """Manages sampler sequence information and directory mapping."""
    
    def __init__(self, sampler_sequence: List[str]):
        self.sampler_sequence = sampler_sequence
        self.sampler_counts = self._count_samplers()
        self.current_indices = {name: 0 for name in set(sampler_sequence)}
    
    def _count_samplers(self) -> Dict[str, int]:
        """Count total occurrences of each sampler in the sequence."""
        counts = {}
        for sampler in self.sampler_sequence:
            counts[sampler] = counts.get(sampler, 0) + 1
        return counts
    
    def get_sampler_info(self, sequence_idx: int) -> Tuple[str, int, int]:
        """
        Get sampler information for a given position in the sequence.
        
        Returns:
            Tuple of (sampler_name, current_index, total_count)
        """
        sampler_name = self.sampler_sequence[sequence_idx]
        self.current_indices[sampler_name] += 1
        current_idx = self.current_indices[sampler_name]
        total_count = self.sampler_counts[sampler_name]
        
        return sampler_name, current_idx, total_count
    
    def get_directory_name(self, sampler_name: str, current_idx: int) -> str:
        """Generate directory name for the sampler."""
        return f"{sampler_name}sampler_results_{current_idx}"
    
    def get_previous_sampler_directory(self, sequence_idx: int) -> Optional[str]:
        """Get the directory name of the previous sampler in the sequence."""
        if sequence_idx == 0:
            return None
        
        # Get previous sampler info
        prev_sampler = self.sampler_sequence[sequence_idx - 1]
        prev_idx = sum(1 for i in range(sequence_idx) if self.sampler_sequence[i] == prev_sampler)
        
        return self.get_directory_name(prev_sampler, prev_idx)

class GMMSigmaProvider:
    """
    Provides sigma values from GMM fits and calculates negative log priors.
    Works with sampler sequences, using uniform/Jeffreys priors for first sampler
    and GMM priors from previous samplers for subsequent ones.
    """
    
    def __init__(self, sampler_sequence: List[str], sequence_idx: int, 
                 base_output_dir: str = None, specific_chain: int = None, 
                 sigma_ranges: Dict[str, Tuple[float, float]] = None,
                 prior_type: str = "uniform"):
        """
        Initialize the sigma provider for a specific position in sampler sequence.
        
        Args:
            sampler_sequence: List of sampler names in execution order
            sequence_idx: Current position in the sampler sequence
            base_output_dir: Base directory containing analysis results (default: output_analysis)
            specific_chain: Load from a specific chain number (default: random selection)
            sigma_ranges: Dictionary of (min, max) ranges for each sigma type (default: [0.5, 10.0] for all)
            prior_type: Type of prior for first sampler ("uniform" or "jeffreys")
        """
        self.sequence_manager = SamplerSequenceManager(sampler_sequence)
        self.sequence_idx = sequence_idx
        self.pair_types = ["AA", "AB", "BC"]
        self.prior_type = prior_type
        
        # Get current sampler info
        self.sampler_name, self.current_idx, self.total_count = self.sequence_manager.get_sampler_info(sequence_idx)
        
        # Set base output directory
        if base_output_dir is None:
            self.base_output_dir = os.path.join(os.getcwd(), "output_analysis")
        else:
            self.base_output_dir = base_output_dir
        
        # Determine if this is the first sampler and set up accordingly
        self.is_first_sampler = (sequence_idx == 0)
        
        if self.is_first_sampler:
            print(f"First sampler in sequence: {self.sampler_name}")
            self.gmm_params = None
            self.output_dir = None
        else:
            # Get previous sampler's directory
            prev_dir = self.sequence_manager.get_previous_sampler_directory(sequence_idx)
            self.output_dir = os.path.join(self.base_output_dir, prev_dir)
            print(f"Loading GMM from previous sampler: {prev_dir}")
            self.gmm_params = self._load_gmm_parameters()
            self._print_gmm_stats()
        
        self.specific_chain = specific_chain
        
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
        if self.gmm_params is None:
            return
            
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
        """Sample sigma values from GMMs or uniform distribution for first sampler."""
        sigma = {}
        
        if self.is_first_sampler:
            # Sample from uniform distribution within ranges
            for pair_type in self.pair_types:
                min_val, max_val = self.sigma_ranges[pair_type]
                sigma[pair_type] = np.random.uniform(min_val, max_val)
            print(f"Sampled sigma values uniformly for first sampler: {sigma}")
            return sigma
        
        # Sample from GMMs for subsequent samplers
        for pair_type in self.pair_types:
            gmm_info = self.gmm_params.get(pair_type) if self.gmm_params else None
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
        
        # Enforce sigma_ranges as hard bounds
        for pair_type, value in sigma.items():
            min_val, max_val = self.sigma_ranges[pair_type]
            if not (min_val <= value <= max_val):
                print(f"{pair_type}={value:.3f} outside {pair_type} range [{min_val}, {max_val}]. Using default.")
                sigma[pair_type] = self.default_sigma[pair_type]
        
        return sigma
    
    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior with uniform/Jeffreys prior for first sampler, GMM for others."""
        if not sigma:
            return np.inf
        
        # Check bounds for all samplers
        for pair_type, value in sigma.items():
            if pair_type in self.sigma_ranges:
                min_val, max_val = self.sigma_ranges[pair_type]
                if value < min_val or value > max_val:
                    print(f"{pair_type}={value:.3f} outside range [{min_val}, {max_val}]")
                    return np.inf
        
        if self.is_first_sampler:
            # Use uniform or Jeffreys prior for first sampler
            log_prior = 0.0
            for pair_type, value in sigma.items():
                min_val, max_val = self.sigma_ranges[pair_type]
                if self.prior_type == "jeffreys":
                    # Jeffreys prior: log p(σ) = -log(σ)
                    log_prior += -np.log(value)
                else:
                    # Uniform prior: log p(σ) = -log(max - min)
                    log_prior += -np.log(max_val - min_val)
            
            return -log_prior
        
        # Use GMM prior for subsequent samplers
        log_priors = [self._calculate_gmm_log_prob(value, pair_type) for pair_type, value in sigma.items()]
        valid_priors = [lp for lp in log_priors if lp > -np.inf]
        
        if not valid_priors:
            print(f"All sigmas have zero prior probability: {sigma}")
            return np.inf
        
        total_neg_log_prior = -np.sum(valid_priors)
        if total_neg_log_prior > 1000:
            print(f"Large negative log prior: {total_neg_log_prior:.3f} for {sigma}")
        return total_neg_log_prior
    
    def _calculate_gmm_log_prob(self, sigma_value: float, pair_type: str) -> float:
        """Calculate log probability for a sigma value under its GMM."""
        if self.gmm_params is None or pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
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
        if self.is_first_sampler or not os.path.exists(self.output_dir):
            return []
        
        pattern = re.compile(r"gmm_fit_.*_chain_(\d+)\.json")
        chain_numbers = [int(match.group(1)) for filename in os.listdir(self.output_dir) 
                         if (match := pattern.match(filename))]
        return sorted(set(chain_numbers))

def sample_sigmas_from_sequence(sampler_sequence: List[str], sequence_idx: int, 
                               base_output_dir: str = None, chain: int = None,
                               prior_type: str = "uniform"):
    """Sample sigma values for a specific position in sampler sequence."""
    provider = GMMSigmaProvider(sampler_sequence, sequence_idx, base_output_dir, chain, prior_type=prior_type)
    sigma_values = provider.sample_sigma_values()
    
    sampler_name = sampler_sequence[sequence_idx]
    print(f"\nSampled sigma values for {sampler_name} (position {sequence_idx}):")
    for key, value in sigma_values.items():
        print(f"  {key}: {value:.4f}")
    
    neg_log_prior = provider.calculate_negative_log_prior(sigma_values)
    print(f"Negative log prior: {neg_log_prior:.4f}")
    
    return sigma_values, neg_log_prior

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Sample sigma values from sampler sequence')
    parser.add_argument('--sequence', type=str, nargs='+', 
                        default=['PairSampler', 'TetramerSampler', 'OctetSampler'],
                        help='Sampler sequence')
    parser.add_argument('--index', type=int, default=0,
                        help='Position in sequence to sample for')
    parser.add_argument('--dir', type=str, default=None,
                        help='Base output directory')
    parser.add_argument('--chain', type=int, default=None)
    parser.add_argument('--prior', type=str, default='uniform', 
                        choices=['uniform', 'jeffreys'],
                        help='Prior type for first sampler')
    
    args = parser.parse_args()
    
    sample_sigmas_from_sequence(args.sequence, args.index, args.dir, args.chain, args.prior)