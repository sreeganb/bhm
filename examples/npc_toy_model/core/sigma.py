import os
import numpy as np
import pickle
from typing import Dict, Tuple, Optional, Union, List

class SigmaProvider:
    """Base class for providing sigma values."""
    
    def __init__(self, prior_type: str = 'jeffreys'):
        """
        Initialize sigma provider with a prior type.
        
        Args:
            prior_type: Type of prior to use ('uniform', 'jeffreys', 'halfcauchy')
        """
        self.prior_type = prior_type
        
        # Default sigma values and ranges for different pair types
        self.default_sigma = {
            'AA': 1.5,
            'AB': 1.2,
            'BC': 1.1,
            'CC': 1.3,
        }
        
        self.default_ranges = {
            'AA': (0.1, 12.0),
            'AB': (0.1, 12.0),
            'BC': (0.1, 12.0),
            'CC': (0.1, 12.0),
        }
        
    def get_sigma_values(self) -> Dict[str, float]:
        """Get default sigma values."""
        return self.default_sigma.copy()
    
    def get_sigma_ranges(self) -> Dict[str, Tuple[float, float]]:
        """Get default sigma ranges."""
        return self.default_ranges.copy()
    
    def sample_sigma(self, pair_type: str = None) -> Union[float, Dict[str, float]]:
        """
        Sample a sigma value or dict of values from the prior.
        
        Args:
            pair_type: If specified, sample for just this pair type; otherwise all
            
        Returns:
            Single sigma value or dictionary of values
        """
        if pair_type is not None:
            # Sample for a single pair type
            sigma_min, sigma_max = self.default_ranges[pair_type]
            return np.random.uniform(sigma_min, sigma_max)
        else:
            # Sample for all pair types
            return {
                pair_type: np.random.uniform(sigma_range[0], sigma_range[1])
                for pair_type, sigma_range in self.default_ranges.items()
            }
    
    def calculate_neg_log_prior(
        self, 
        sigma: Dict[str, float], 
        sigma_range: Optional[Dict[str, Tuple[float, float]]] = None
    ) -> float:
        """
        Calculate negative log prior for given sigma values.
        
        Args:
            sigma: Dictionary of sigma values
            sigma_range: Optional dictionary of sigma ranges, or use defaults
            
        Returns:
            Total negative log prior
        """
        #if sigma_range is None:
        #    sigma_range = self.default_ranges
        sigma_range = self.default_ranges
            
        total_prior = 0.0
        for pair_type, value in sigma.items():
            if pair_type not in sigma_range:
                print(f"Warning: {pair_type} not in sigma range dictionary")
                continue
                
            range_min, range_max = sigma_range[pair_type]
            
            # Check if sigma is within range
            if value <= range_min or value >= range_max:
                return float('inf')
            
            # Calculate prior based on type
            if self.prior_type == 'uniform':
                # Uniform prior has constant density
                total_prior += 0.0
            elif self.prior_type == 'jeffreys':
                # Jeffreys prior is proportional to 1/sigma
                total_prior += np.log(value)
            elif self.prior_type == 'halfcauchy':
                # Half-Cauchy prior with scale=1.0
                scale = 1.0
                total_prior += np.log(1 + (value / scale) ** 2)
            else:
                raise ValueError(f"Unknown prior type: {self.prior_type}")
            
        return total_prior

class GMMSigmaProvider(SigmaProvider):
    """GMM-based sigma provider."""
    
    def __init__(
        self, 
        sampler_name: str, 
        prior_type: str = 'uniform',
        gmm_dir: str = 'output_analysis/gmm'
    ):
        """
        Initialize GMM-based sigma provider.
        
        Args:
            sampler_name: Name of sampler to load GMMs for
            prior_type: Prior type to use for calculations
            gmm_dir: Directory where GMM files are stored
        """
        super().__init__(prior_type)
        self.sampler_name = sampler_name
        self.gmm_dir = gmm_dir
        
        # Load GMM models
        self.gmm_models = self._load_gmm_models()
        
    def _load_gmm_models(self) -> Dict[str, any]:
        """
        Load GMM models from pickle files.
        
        Returns:
            Dictionary of GMM models by pair type
        """
        gmm_models = {}
        gmm_pattern = f"{self.sampler_name}_sigma_*.pkl"
        gmm_files = [f for f in os.listdir(self.gmm_dir) 
                     if f.startswith(f"{self.sampler_name}_sigma_") and f.endswith(".pkl")]
        
        if not gmm_files:
            print(f"Warning: No GMM files found for {self.sampler_name} in {self.gmm_dir}")
            return {}
        
        for gmm_file in gmm_files:
            # Extract pair type from filename (e.g., "PairSampler_sigma_AA.pkl" -> "AA")
            pair_type = gmm_file.replace(f"{self.sampler_name}_sigma_", "").replace(".pkl", "")
            
            try:
                with open(os.path.join(self.gmm_dir, gmm_file), 'rb') as f:
                    gmm_models[pair_type] = pickle.load(f)
                    print(f"Loaded GMM for {pair_type} from {gmm_file}")
            except Exception as e:
                print(f"Error loading GMM for {pair_type}: {e}")
                
        return gmm_models
    
    def get_sigma_values(self) -> Dict[str, float]:
        """
        Sample sigma values from GMM models when available, fall back to defaults.
        
        Returns:
            Dictionary of sigma values
        """
        sigma_values = {}
        
        # For each pair type, sample from GMM if available, otherwise use default
        for pair_type in self.default_sigma:
            if pair_type in self.gmm_models:
                # Sample from GMM and take the first sample
                samples = self.gmm_models[pair_type].sample(1)[0]
                sigma_values[pair_type] = float(samples[0])
            else:
                # Fall back to default
                sigma_values[pair_type] = self.default_sigma[pair_type]
                
        return sigma_values
    
    def sample_sigma(self, pair_type: str = None) -> Union[float, Dict[str, float]]:
        """
        Sample sigma values from GMM distributions.
        
        Args:
            pair_type: If specified, sample for just this pair type
            
        Returns:
            Sampled sigma value(s)
        """
        if pair_type is not None:
            if pair_type in self.gmm_models:
                samples = self.gmm_models[pair_type].sample(1)[0]
                return float(samples[0])
            else:
                return super().sample_sigma(pair_type)
        else:
            sigma_values = {}
            for pt in self.default_sigma:
                if pt in self.gmm_models:
                    samples = self.gmm_models[pt].sample(1)[0]
                    sigma_values[pt] = float(samples[0])
                else:
                    sigma_values[pt] = super().sample_sigma(pt)
            return sigma_values
    
    def calculate_neg_log_prior(
        self, 
        sigma: Dict[str, float], 
        sigma_range: Optional[Dict[str, Tuple[float, float]]] = None
    ) -> float:
        """
        Calculate negative log prior based on GMM density.
        
        Args:
            sigma: Dictionary of sigma values
            sigma_range: Optional sigma ranges (used for basic range checking)
            
        Returns:
            GMM-based negative log prior
        """
        # First do basic range checks
        for pair_type, value in sigma.items():
            if sigma_range and pair_type in sigma_range:
                range_min, range_max = sigma_range[pair_type]
                if value <= range_min or value >= range_max:
                    return float('inf')
        
        # Then calculate GMM-based prior if GMMs are available
        prior = 0.0
        for pair_type, value in sigma.items():
            if pair_type in self.gmm_models:
                # Reshape for sklearn GMM
                x = np.array([[value]])
                # Get log probability density
                log_density = self.gmm_models[pair_type].score_samples(x)[0]
                # Add negative log density to prior
                prior -= log_density
            else:
                # Fall back to parametric prior
                range_vals = sigma_range.get(pair_type, self.default_ranges.get(pair_type)) if sigma_range else self.default_ranges.get(pair_type)
                if range_vals:
                    # Add standard prior contribution
                    if self.prior_type == 'jeffreys':
                        prior += np.log(value)
                    elif self.prior_type == 'halfcauchy':
                        scale = 1.0
                        prior += np.log(1 + (value / scale) ** 2)
                        
        print(f"GMM-based negative log prior: {prior}")
        
        return prior


# Utility functions for direct use in other modules

def initialize_sigma(state, sampler_sequence, current_sampler_idx, sigma_source=None, prior_type='uniform'):
    """
    Initialize sigma values for a system state.
    
    Args:
        state: SystemState object
        sampler_sequence: List of sampler names
        current_sampler_idx: Index of current sampler
        sigma_source: 'default' or 'gmm'
        prior_type: Type of prior to use
        
    Returns:
        Dictionary of initialized sigma values
    """
    current_sampler = sampler_sequence[current_sampler_idx]
    
    # Auto-determine sigma source
    if sigma_source is None:
        if current_sampler_idx == 0:
            sigma_source = "default"
        else:
            sigma_source = "gmm"
    
    if sigma_source == "gmm" and current_sampler_idx > 0:
        # Use GMM from previous sampler
        prev_sampler = sampler_sequence[current_sampler_idx - 1]
        print(f"Initializing sigma from GMM based on previous sampler: {prev_sampler}")
        sigma_provider = GMMSigmaProvider(prev_sampler, prior_type=prior_type)
    else:
        print(f"Initializing sigma with default provider")
        sigma_provider = SigmaProvider(prior_type=prior_type)
    
    # Sample and update sigma values
    sigma_values = sigma_provider.get_sigma_values()
    if hasattr(state, 'update_sigma'):
        state.update_sigma(sigma_values)
    
    return sigma_values

def calculate_negative_log_prior(state, provider_type='default', sampler_name=None, prior_type='jeffreys'):
    """
    Calculate negative log prior for a system state.
    
    Args:
        state: SystemState with sigma values
        provider_type: 'default' or 'gmm'
        sampler_name: Required if provider_type is 'gmm'
        prior_type: Type of prior to use
        
    Returns:
        Negative log prior value
    """
    if not hasattr(state, 'sigma'):
        raise ValueError("State object must have sigma attribute")
    
    if provider_type == 'gmm' and sampler_name:
        provider = GMMSigmaProvider(sampler_name, prior_type=prior_type)
    else:
        provider = SigmaProvider(prior_type=prior_type)
        
    # Get sigma ranges from the state if available
    sigma_range = getattr(state, 'sigma_range', None)
    
    return provider.calculate_neg_log_prior(state.sigma, state.sigma_range)
