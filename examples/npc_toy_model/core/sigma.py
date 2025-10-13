import os
import json
import numpy as np
from typing import Dict, Optional
from scipy import stats
from core.state import SystemState
from typing import List, Tuple, Any

class SigmaPrior:
    """Prior distribution for sigma parameters"""
    
    # Define available prior distributions as class attribute
    PRIOR_DISTRIBUTIONS = {
        'uniform': lambda low, high: stats.uniform(loc=low, scale=high-low),
        'jeffreys': lambda low, high: stats.loguniform(a=low, b=high),
        'log_uniform': lambda low, high: stats.loguniform(a=low, b=high),
        'gamma': lambda shape, scale: stats.gamma(a=shape, scale=scale),
        'cauchy': lambda loc, scale: stats.cauchy(loc=loc, scale=scale)
    }
    
    def __init__(self, 
                 pair_types: List[str],
                 sigma_ranges: Dict[str, Tuple[float, float]],
                 prior_type: str = "uniform",
                 use_gmm: bool = False,
                 gmm_params: Optional[Dict[str, Dict[str, Any]]] = None):
        """
        Initialize sigma prior.
        
        Args:
            pair_types: List of pair type names (e.g., ['AA', 'AB', 'BC'])
            sigma_ranges: Dict mapping pair types to (low, high) bounds
            prior_type: Type of prior ("uniform", "jeffreys", "log_uniform")
            use_gmm: Whether to use GMM distribution
            gmm_params: GMM parameters if use_gmm=True
        """
        self.pair_types = pair_types
        self.sigma_ranges = sigma_ranges
        self.prior_type = prior_type
        self.use_gmm = use_gmm
        self.gmm_params = gmm_params or {}
        
        # Validate prior type
        if not use_gmm and prior_type not in self.PRIOR_DISTRIBUTIONS:
            raise ValueError(f"Unknown prior type: {prior_type}. Choose from {list(self.PRIOR_DISTRIBUTIONS.keys())}")
        
        # Validate ranges
        for pt in pair_types:
            if pt not in sigma_ranges:
                raise ValueError(f"Missing sigma range for pair type: {pt}")
        
    def _load_gmm(self, gmm_file: str) -> Optional[Dict]:
        """Load GMM parameters from JSON file."""
        if not os.path.exists(gmm_file):
            print(f"Warning: GMM file not found: {gmm_file}")
            return None
            
        try:
            with open(gmm_file, 'r') as f:
                gmm_data = json.load(f)
            
            # Validate GMM structure
            for pair_type in self.pair_types:
                if pair_type not in gmm_data:
                    print(f"Warning: Missing GMM data for {pair_type}")
                    return None
                    
                gmm = gmm_data[pair_type]
                required_keys = ['weights', 'means', 'covariances', 'n_components']
                if not all(k in gmm for k in required_keys):
                    print(f"Warning: Invalid GMM structure for {pair_type}")
                    return None
            
            print(f"Loaded GMM prior from {gmm_file}")
            return gmm_data
            
        except Exception as e:
            print(f"Error loading GMM file: {e}")
            return None
    
    def initialize_sigma(self, rng: np.random.Generator = None) -> Dict[str, float]:
        """
        Initialize sigma values randomly within specified ranges.
        
        Args:
            rng: NumPy random generator (uses default if None)
            
        Returns:
            Dictionary of initialized sigma values
        """
        if rng is None:
            rng = np.random.default_rng()
        
        sigma = {}
        for pair_type in self.pair_types:
            low, high = self.sigma_ranges[pair_type]
            # Initialize in log-space for better coverage
            sigma[pair_type] = float(np.exp(rng.uniform(np.log(low), np.log(high))))
        
        return sigma
    
    def log_prior(self, sigma: Dict[str, float]) -> float:
        """
        Calculate log prior probability for given sigma values.
        
        Args:
            sigma: Dictionary of sigma values for each pair type
            
        Returns:
            Log prior probability (NOT negative log prior)
        """
        if self.use_gmm:
            return self._gmm_log_prior(sigma)
        else:
            return self._simple_log_prior(sigma)
    
    def _simple_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate log prior using simple distributions."""
        if self.prior_type not in self.PRIOR_DISTRIBUTIONS:
            raise ValueError(f"Unknown prior type: {self.prior_type}")
        
        log_prob = 0.0
        
        for pair_type, value in sigma.items():
            low, high = self.sigma_ranges[pair_type]
            
            # Check bounds
            if not (low <= value <= high):
                return -np.inf
            
            # Get distribution and evaluate
            dist = self.PRIOR_DISTRIBUTIONS[self.prior_type](low, high)
            log_prob += dist.logpdf(value)
        
        return log_prob

    @classmethod
    def from_gmm_file(cls, gmm_file: str, 
                    prior_type: str = "uniform",
                    sigma_ranges: Optional[Dict[str, Tuple[float, float]]] = None) -> 'SigmaPrior':
        """
        Load SigmaPrior from GMM parameters file.
        
        Args:
            gmm_file: Path to JSON file with GMM parameters
            prior_type: Fallback prior type if GMM fails
            sigma_ranges: Bounds for each pair type
            
        Returns:
            SigmaPrior configured with GMM
        """
        import json
        
        try:
            with open(gmm_file, 'r') as f:
                gmm_params = json.load(f)
            
            pair_types = list(gmm_params.keys())
            
            # Use provided ranges or extract from GMM means
            if sigma_ranges is None:
                sigma_ranges = {}
                for pt in pair_types:
                    means = np.array(gmm_params[pt]['means']).flatten()
                    stds = np.sqrt(np.array(gmm_params[pt]['covariances']).flatten())
                    low = max(0.1, np.min(means - 3*stds))
                    high = np.max(means + 3*stds)
                    sigma_ranges[pt] = (float(low), float(high))
            
            return cls(
                pair_types=pair_types,
                sigma_ranges=sigma_ranges,
                prior_type=prior_type,
                use_gmm=True,
                gmm_params=gmm_params
            )
            
        except Exception as e:
            print(f"Warning: Failed to load GMM from {gmm_file}: {e}")
            print(f"Falling back to simple {prior_type} prior")
            
            # Fallback to simple prior
            pair_types = list(sigma_ranges.keys()) if sigma_ranges else ['AA', 'AB', 'BC']
            default_ranges = {
                'AA': (1.0, 5.0),
                'AB': (1.5, 8.0),
                'BC': (2.0, 12.0)
            }
            if sigma_ranges is None:
                sigma_ranges = default_ranges
            
            return cls(
                pair_types=pair_types,
                sigma_ranges=sigma_ranges,
                prior_type=prior_type,
                use_gmm=False
            )    
    def _gmm_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate log prior using GMM from previous sampler."""
        log_prob = 0.0
        
        for pair_type, value in sigma.items():
            # Check bounds
            low, high = self.sigma_ranges[pair_type]
            if not (low <= value <= high):
                return -np.inf
            
            # Get GMM parameters
            gmm = self.gmm_params[pair_type]
            weights = np.array(gmm['weights'])
            means = np.array(gmm['means']).flatten()
            covs = np.array(gmm['covariances']).flatten()
            
            # Ensure numerical stability
            covs = np.maximum(covs, 1e-10)
            weights = weights / np.sum(weights)
            
            # Calculate log probability for each component
            stds = np.sqrt(covs)
            component_logprobs = (
                np.log(weights) 
                - 0.5 * np.log(2 * np.pi * covs)
                - 0.5 * ((value - means) / stds) ** 2
            )
            
            # Log-sum-exp for numerical stability
            max_logprob = np.max(component_logprobs)
            log_prob += max_logprob + np.log(np.sum(np.exp(component_logprobs - max_logprob)))
        
        return log_prob


def create_sigma_prior(state: SystemState, 
                      gmm_file: Optional[str] = None,
                      prior_type: str = "uniform") -> SigmaPrior:
    """
    Factory function to create a SigmaPrior from state or GMM file.
    
    Args:
        state: SystemState containing sigma dict with pair types
        gmm_file: Optional path to GMM parameters JSON file
        prior_type: Type of prior distribution ("uniform", "jeffreys", "log_uniform")
        
    Returns:
        Configured SigmaPrior object
    """
    # Extract pair types from state
    pair_types = list(state.sigma.keys())
    
    # Define default sigma ranges if not in state
    default_ranges = {
        'AA': (1.0, 5.0),
        'AB': (1.5, 8.0),
        'BC': (2.0, 12.0)
    }
    
    # Get ranges from state or use defaults
    if hasattr(state, 'sigma_range') and state.sigma_range:
        sigma_ranges = state.sigma_range
    else:
        sigma_ranges = {pt: default_ranges.get(pt, (1.0, 10.0)) for pt in pair_types}
    
    if gmm_file is None:
        # Use simple prior (uniform, jeffreys, etc.)
        return SigmaPrior(
            pair_types=pair_types,
            sigma_ranges=sigma_ranges,  # ← FIX: Pass sigma_ranges
            prior_type=prior_type,
            use_gmm=False
        )
    else:
        # Load GMM from file for subsequent stages
        return SigmaPrior.from_gmm_file(
            gmm_file, 
            prior_type,
            sigma_ranges=sigma_ranges  # ← FIX: Pass sigma_ranges
        )


# Example usage
if __name__ == "__main__":
    # First sampler in sequence (no GMM)
    pair_types = ['AA', 'AB', 'BC']
    sigma_ranges = {'AA': (1.0, 5.0), 'AB': (1.5, 8.0), 'BC': (2.0, 12.0)}

    prior = SigmaPrior(pair_types, sigma_ranges, gmm_file=None, prior_type='uniform')

    # Initialize sigma
    sigma = prior.initialize_sigma()
    print("Initialized sigma:", sigma)
    
    # Calculate log prior
    log_p = prior.log_prior(sigma)
    print(f"Log prior: {log_p:.4f}")
    print(f"Negative log prior: {-log_p:.4f}")
    
    # Subsequent sampler (with GMM)
    # prior2 = SigmaPrior(pair_types, sigma_ranges, 
    #                     gmm_file='output_analysis/pairsampler_results/gmm_all_chains.json',
    #                     prior_type='uniform')