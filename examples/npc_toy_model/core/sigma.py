import os
import re
import json
import random
import numpy as np
from typing import Dict, List, Tuple, Union
from scipy import stats

from core.state import SystemState

class GMMSigmaProvider:
    """
    Provides sigma values from GMM fits and calculates negative log priors.
    Works with PairSampler, TetramerSampler, and OctetSampler outputs.
    """
    # Define sampler source mappings
    SAMPLER_SOURCES = {
        "PairSampler": "pairsampler",
        "TetramerSampler": "pairsampler",
        "OctetSampler": "tetramersampler",
        "FullSampler": "fullsampler"
    }
    
    # Default sigma values if GMM not available
    DEFAULT_SIGMA = {
        "AA": 2.0, 
        "AB": 2.0, 
        "BC": 2.0,
        "CC": 2.0
    }
    
    # Sigma range for simple priors
    SIGMA_MIN = 0.0
    SIGMA_MAX = 20.0
    
    # Available simple priors
    SIMPLE_PRIORS = {
        'uniform': lambda x, a=SIGMA_MIN, b=SIGMA_MAX: stats.uniform(loc=a, scale=b-a),
        'jeffreys': lambda x, a=SIGMA_MIN, b=SIGMA_MAX: stats.reciprocal(a=max(a, 1e-6), b=b),
        'half_cauchy': lambda x, scale=10.0: stats.halfcauchy(scale=scale),
        'gamma': lambda x, a=2.0, scale=5.0: stats.gamma(a=a, scale=scale),
        'half_normal': lambda x, scale=10.0: stats.halfnorm(scale=scale),
        'exponential': lambda x, scale=10.0: stats.expon(scale=scale)
    }
    
    def __init__(self, sampler_name: str = "PairSampler", prior_type: str = "uniform", 
                 sequence_position: int = 0):
        """
        Initialize the sigma provider for a specific sampler type.
        
        Args:
            sampler_name: Name of the sampler ('PairSampler', 'TetramerSampler', 'OctetSampler')
            prior_type: Type of simple prior to use ('uniform', 'jeffreys', 'half_cauchy', 'gamma', 'half_normal', 'exponential')
            sequence_position: Position in the sampler sequence (0-indexed). Use simple priors if 0, GMM if available otherwise
        """
        self.sampler_name = sampler_name
        self.prior_type = prior_type
        self.sequence_position = sequence_position
        self.pair_types = ["AA", "AB", "BC"]
        
        # Determine source folder for GMM files
        source_sampler = self.SAMPLER_SOURCES.get(sampler_name, "pairsampler")
        self.output_dir = os.path.join(os.getcwd(), f"output_analysis/{source_sampler}_results")
        
        # Load GMM parameters (only if not first in sequence)
        self.gmm_params = {}
        self.use_gmm = False
        
        if sequence_position > 0:
            self.gmm_params = self._load_gmm_parameters()
            # Check if we have valid GMM data
            self.use_gmm = any(
                gmm and all(k in gmm for k in ['n_components', 'means', 'covariances', 'weights'])
                for gmm in self.gmm_params.values()
            )
        
        if sequence_position == 0 or not self.use_gmm:
            print(f"Using simple '{prior_type}' prior for sigma values (sequence position: {sequence_position})")
        else:
            print(f"Using GMM prior from previous sampler results (sequence position: {sequence_position})")
        
    def _load_gmm_parameters(self) -> Dict:
        """Load GMM parameters from JSON files for a randomly selected chain."""
        gmm_params = {sigma_type: None for sigma_type in self.pair_types}
        
        if not os.path.exists(self.output_dir):
            print(f"Output directory not found: {self.output_dir}")
            print(f"Will use simple priors.")
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
            print(f"Will use simple priors.")
            return gmm_params
        
        # Select random chain with the most sigma types
        chains_by_completeness = sorted(
            chain_files.keys(), 
            key=lambda c: len(chain_files[c]), 
            reverse=True
        )
        selected_chain = chains_by_completeness[0]
        print(f"Loading GMM parameters from chain {selected_chain}")
        
        # Load GMM files for the selected chain
        for sigma_type in self.pair_types:
            file_path = os.path.join(self.output_dir, f"gmm_fit_{sigma_type}_chain_{selected_chain}.json")
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_params[sigma_type] = json.load(f)
                        print(f"  Loaded GMM for {sigma_type}")
                except Exception as e:
                    print(f"  Error loading {file_path}: {e}")
        
        return gmm_params

    def initialize_sigma(
        self,
        state: SystemState,
        sigma_source: str = "gmm",
        sampler_name: str = "PairSampler",
        spread_chains: bool = True,
        chain_id: int = 0,
        n_chains: int = 1
    ) -> None:
        """
        Initialize sigma values for a state, with option to spread chains
        
        Args:
            state: SystemState to initialize
            sigma_source: "gmm" or "uniform"
            sampler_name: Which sampler (affects which pairs are used)
            spread_chains: If True, spread initial values across parameter space
            chain_id: Chain identifier (0 to n_chains-1)
            n_chains: Total number of chains
        """
        # Get pair types for this sampler
        if sampler_name == "PairSampler":
            pair_types = ["AA", "AB", "BC"]
        elif sampler_name == "TetramerSampler":
            pair_types = ["AA", "AB", "BC"]  # Adjust as needed
        else:
            pair_types = list(state.sigma.keys())
        
        state.sigma = {}
        state.sigma_range = {}
        
        for pair_type in pair_types:
            # Define reasonable bounds
            # Adjust these based on your target distances
            if pair_type == "AA":
                low, high = 1.0, 5.0
            elif pair_type == "AB":
                low, high = 1.5, 8.0
            elif pair_type == "BC":
                low, high = 2.0, 12.0
            else:
                low, high = 1.0, 10.0
            
            state.sigma_range[pair_type] = (low, high)
            
            if sigma_source == "gmm":
                # Sample from GMM
                state.sigma[pair_type] = self.sample_from_gmm(pair_type)
            else:
                # Uniform initialization
                if spread_chains and n_chains > 1:
                    # Spread chains across parameter space in log-space
                    frac = chain_id / max(1, n_chains - 1)
                    log_val = np.log(low) + frac * (np.log(high) - np.log(low))
                    state.sigma[pair_type] = float(np.exp(log_val))
                else:
                    # Random initialization
                    state.sigma[pair_type] = float(
                        np.exp(np.random.uniform(np.log(low), np.log(high)))
                    )
        
        return state

    def _sample_from_simple_prior(self) -> Dict[str, float]:
        """Sample sigma values from simple prior distributions."""
        sigma = {}
        
        if self.prior_type not in self.SIMPLE_PRIORS:
            print(f"Warning: Unknown prior type '{self.prior_type}'. Using 'uniform'.")
            self.prior_type = 'uniform'
        
        prior_dist = self.SIMPLE_PRIORS[self.prior_type]
        
        for pair_type in self.pair_types:
            try:
                # Create distribution
                dist = prior_dist(None)
                
                # Sample and clip to valid range
                sampled_value = dist.rvs()
                sampled_value = np.clip(sampled_value, self.SIGMA_MIN, self.SIGMA_MAX)
                sigma[pair_type] = float(sampled_value)
                
            except Exception as e:
                print(f"Error sampling simple prior for {pair_type}: {e}. Using default.")
                sigma[pair_type] = self.DEFAULT_SIGMA[pair_type]
        
        return sigma

    def sample_sigma_values(self) -> Dict[str, float]:
        """Sample sigma values from the loaded GMMs or use simple priors."""
        if not self.use_gmm:
            return self._sample_from_simple_prior()
        
        sigma = {}
        
        for pair_type in self.pair_types:
            gmm_info = self.gmm_params.get(pair_type)
            
            # Use simple prior if no valid GMM data
            if not gmm_info or not all(k in gmm_info for k in ['n_components', 'means', 'covariances', 'weights']):
                # Sample from simple prior instead
                dist = self.SIMPLE_PRIORS[self.prior_type](None)
                sampled_value = dist.rvs()
                sigma[pair_type] = float(np.clip(sampled_value, self.SIGMA_MIN, self.SIGMA_MAX))
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
                
                # Sample from normal distribution and clip
                sampled_value = np.random.normal(loc=mean_value, scale=np.sqrt(cov_value))
                sigma[pair_type] = float(np.clip(sampled_value, self.SIGMA_MIN, self.SIGMA_MAX))
                
            except Exception as e:
                print(f"Error sampling GMM for {pair_type}: {e}. Using simple prior.")
                dist = self.SIMPLE_PRIORS[self.prior_type](None)
                sampled_value = dist.rvs()
                sigma[pair_type] = float(np.clip(sampled_value, self.SIGMA_MIN, self.SIGMA_MAX))
        
        return sigma
    
    def calculate_negative_log_prior(self, sigma: Union[Dict[str, float], 'SystemState']) -> float:
        """
        Calculate negative log prior across all sigma values.
        
        Args:
            sigma: Either a dictionary of sigma values or a SystemState object
            
        Returns:
            Negative log prior value
        """
        # Handle if a SystemState object is passed instead of dict
        if hasattr(sigma, 'sigma'):
            sigma = sigma.sigma
        
        if not sigma or not isinstance(sigma, dict):
            return np.inf
        
        # Use GMM prior if available, otherwise use simple prior
        if self.use_gmm:
            return self._calculate_gmm_log_prior(sigma)
        else:
            return self._calculate_simple_log_prior(sigma)
    
    def _calculate_simple_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior using simple prior distributions."""
        if self.prior_type not in self.SIMPLE_PRIORS:
            return 100000.0
        
        log_priors = []
        prior_dist = self.SIMPLE_PRIORS[self.prior_type]
        
        for pair_type, value in sigma.items():
            if not (self.SIGMA_MIN <= value <= self.SIGMA_MAX):
                log_priors.append(-np.inf)
                continue
            
            try:
                dist = prior_dist(None)
                log_prob = dist.logpdf(value)
                
                if np.isfinite(log_prob):
                    log_priors.append(log_prob)
                else:
                    log_priors.append(-np.inf)
                    
            except Exception as e:
                print(f"Error calculating simple prior for {pair_type}: {e}")
                log_priors.append(-np.inf)
        
        if not log_priors or all(lp == -np.inf for lp in log_priors):
            return 100000.0
        
        # Filter out -inf values and sum
        valid_priors = [lp for lp in log_priors if lp > -np.inf]
        if not valid_priors:
            return 100000.0
            
        return -np.sum(valid_priors)
    
    def _calculate_gmm_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior using GMM distributions."""
        log_priors = []
        
        for pair_type, value in sigma.items():
            log_prob = self._calculate_gmm_log_prob(value, pair_type)
            if log_prob > -np.inf:
                log_priors.append(log_prob)
        
        if not log_priors:
            return 100000.0
            
        return -np.sum(log_priors)
    
    def _calculate_gmm_log_prob(self, sigma_value: float, pair_type: str) -> float:
        """Calculate log probability for a sigma value given a GMM."""
        if pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
            # Fall back to simple prior
            try:
                dist = self.SIMPLE_PRIORS[self.prior_type](None)
                return dist.logpdf(sigma_value)
            except:
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
            # Fall back to simple prior
            try:
                dist = self.SIMPLE_PRIORS[self.prior_type](None)
                return dist.logpdf(sigma_value)
            except:
                return -np.inf

# Simple function for sampling sigma values from a specified sampler
def get_sigma_values(sampler_name="PairSampler", prior_type="uniform", sequence_position=0):
    """
    Get sigma values for the specified sampler.
    
    Args:
        sampler_name: 'PairSampler', 'TetramerSampler', or 'OctetSampler'
        prior_type: Type of simple prior ('uniform', 'jeffreys', 'half_cauchy', 'gamma', 'half_normal', 'exponential')
        sequence_position: Position in sequence (0 = use simple priors, >0 = try GMM)
        
    Returns:
        Tuple: (sigma_values, negative_log_prior)
    """
    provider = GMMSigmaProvider(sampler_name, prior_type, sequence_position)
    sigma_values = provider.sample_sigma_values()
    neg_log_prior = provider.calculate_negative_log_prior(sigma_values)
    
    print(f"\nSampled sigma values for {sampler_name} (position {sequence_position}):")
    for key, value in sigma_values.items():
        print(f"  {key}: {value:.4f}")
    print(f"Negative log prior: {neg_log_prior:.4f}")
    
    return sigma_values, neg_log_prior

# Example usage
if __name__ == "__main__":
    import sys
    
    sampler = "PairSampler" if len(sys.argv) < 2 else sys.argv[1]
    prior = "uniform" if len(sys.argv) < 3 else sys.argv[2]
    position = 0 if len(sys.argv) < 4 else int(sys.argv[3])
    
    print(f"\nTesting with sampler={sampler}, prior={prior}, position={position}")
    get_sigma_values(sampler, prior, position)
