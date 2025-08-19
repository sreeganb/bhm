import os
import re
import json
import random
import logging
import numpy as np
from typing import Dict, Optional, List, Tuple
from parameters import SystemParameters  # Assumes an external file "parameters.py" defines SystemParameters

###############################################################################
# Debug Logger Setup
###############################################################################

def setup_logger(
    logger_name: str,
    debug_file: Optional[str] = None,
    debug_level: int = logging.INFO
) -> logging.Logger:
    """
    Create and configure a logger.

    Args:
        logger_name: The name of this logger.
        debug_file: Optional path to a log file.
        debug_level: Debug level (e.g. logging.DEBUG, logging.INFO).
    Returns:
        A configured logger.
    """
    # Create logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(debug_level)

    # Clear existing handlers so repeated setup_logger calls won't multiply logs
    logger.handlers.clear()

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(debug_level)
    formatter = logging.Formatter('[%(levelname)s] %(name)s: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Optional file handler
    if debug_file:
        fh = logging.FileHandler(debug_file, mode='w')
        fh.setLevel(debug_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger

###############################################################################
# SamplerSequenceManager
###############################################################################

class SamplerSequenceManager:
    """
    Manages sampler sequence information and directory mapping.
    Helps retrieve sampler names and indexes as you iterate through them.
    """
    def __init__(self, sampler_sequence: List[str]):
        self.logger = logging.getLogger(__name__)
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
        Get sampler information (name, current index, total count) for a given 
        position in the sequence.
        """
        sampler_name = self.sampler_sequence[sequence_idx]
        self.current_indices[sampler_name] += 1
        current_idx = self.current_indices[sampler_name]
        total_count = self.sampler_counts[sampler_name]

        self.logger.debug(
            f"Sampler info for sequence_idx={sequence_idx}: "
            f"{sampler_name}, current={current_idx}, total={total_count}"
        )
        return sampler_name, current_idx, total_count

    def get_directory_name(self, sampler_name: str, current_idx: int) -> str:
        """Generate directory name for the sampler."""
        return f"{sampler_name}sampler_results_{current_idx}"

    def get_previous_sampler_directory(self, sequence_idx: int) -> Optional[str]:
        """Get the directory name of the previous sampler in the sequence."""
        if sequence_idx == 0:
            return None

        prev_sampler = self.sampler_sequence[sequence_idx - 1]
        # The previous index is how many times we've encountered prev_sampler thus far
        prev_idx = sum(1 for i in range(sequence_idx) if self.sampler_sequence[i] == prev_sampler)
        return self.get_directory_name(prev_sampler, prev_idx)

###############################################################################
# SigmaProvider
###############################################################################

import os
import json
import random
import logging
import re
import numpy as np
from typing import Dict, List, Optional, Tuple
from scipy.stats import gamma


class GMMSigmaProvider:
    """
    Provides sigma values from GMM fits and calculates negative log priors.
    Works with sampler sequences: 
      - For the first sampler, uses uniform/Jeffreys/gamma priors.
      - For subsequent samplers, uses GMM parameters from a previous sampler's fit.
    """

    def __init__(
        self,
        sampler_sequence: List[str],
        sequence_idx: int,
        base_output_dir: Optional[str] = None,
        specific_chain: Optional[int] = None,
        sigma_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
        prior_type: str = "uniform",
        debug_logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the SigmaProvider for a specific position in a sampler sequence.
        """
        # Setup logging
        self.logger = debug_logger or logging.getLogger(__name__)

        # Re-seed the RNG using process ID to ensure different seeds per process
        process_seed = os.getpid() + int(np.random.random() * 100000)
        np.random.seed(process_seed)
        random.seed(process_seed)

        self.sequence_manager = SamplerSequenceManager(sampler_sequence)
        self.sequence_idx = sequence_idx
        self.prior_type = prior_type
        self.specific_chain = specific_chain
        self.params = SystemParameters()
        
        # Get all pair types from params
        self.pair_types = list(self.params.pair_distances.keys())
        
        # Get sampler info
        self.sampler_name, self.current_idx, self.total_count = self.sequence_manager.get_sampler_info(sequence_idx)

        # Set base output directory
        self.base_output_dir = base_output_dir or os.path.join(os.getcwd(), "output_analysis")

        # Is this the first sampler in the sequence?
        self.is_first_sampler = (sequence_idx == 0)

        # Configure sigma ranges BEFORE loading GMM
        self.sigma_ranges = self._setup_sigma_ranges(sigma_ranges)
        self.default_sigma = self._get_default_sigmas()

        # If not first sampler, load GMM from previous sampler
        if self.is_first_sampler:
            self.logger.info(f"First sampler in sequence: {self.sampler_name}")
            self.gmm_params = None
            self.output_dir = None
        else:
            prev_dir = self.sequence_manager.get_previous_sampler_directory(sequence_idx)
            self.output_dir = os.path.join(self.base_output_dir, prev_dir) if prev_dir else None
            if self.output_dir:
                self.logger.info(f"Loading GMM from previous sampler directory: {prev_dir}")
            self.gmm_params = self._load_gmm_parameters()
            self.print_gmm_stats()
        
        # Print the assigned ranges
        self.logger.debug(f"Sigma ranges for {self.sampler_name}:")
        for pair_type, (mn, mx) in self.sigma_ranges.items():
            self.logger.debug(f"  {pair_type}: [{mn:.3f}, {mx:.3f}]")

    ###########################################################################
    # Helper Setup
    ###########################################################################

    def _setup_sigma_ranges(
        self, user_ranges: Optional[Dict[str, Tuple[float, float]]]
    ) -> Dict[str, Tuple[float, float]]:
        """
        Determine final sigma ranges.
        """
        default_ranges = {}
        for pair_type in self.params.pair_distances.keys():
            sum_radii = self.params.radii[pair_type[0]] + self.params.radii[pair_type[1]]
            lower_bound = 0.02 * sum_radii
            upper_bound = 0.2 * sum_radii
            default_ranges[pair_type] = (lower_bound, upper_bound)

        if user_ranges:
            default_ranges.update(user_ranges)
        return default_ranges

    def _get_default_sigmas(self) -> Dict[str, float]:
        """
        Choose a default sigma for each pair type (geometric midpoint).
        """
        return {
            pt: np.sqrt(min_val * max_val)
            for pt, (min_val, max_val) in self.sigma_ranges.items()
        }

    ###########################################################################
    # GMM Loading and Diagnostics
    ###########################################################################

    def _load_gmm_parameters(self) -> Dict[str, Optional[dict]]:
        """Load GMM parameters for each pair type from JSON files."""
        gmm_params = {sigma_type: None for sigma_type in self.pair_types}
        
        if not self.output_dir or not os.path.exists(self.output_dir):
            self.logger.warning(f"Output directory not found: {self.output_dir}")
            return gmm_params

        # Find all GMM files and extract chain numbers
        pattern = re.compile(r"gmm_fit_([A-Z]+)_chain_(\d+)\.json")
        available_files = {}
        
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                pair_type, chain_num = match.group(1), int(match.group(2))
                if pair_type not in available_files:
                    available_files[pair_type] = []
                available_files[pair_type].append(chain_num)

        if not available_files:
            self.logger.warning(f"No GMM fit files found in {self.output_dir}")
            return gmm_params

        # Choose chain number
        all_chains = set()
        for chains in available_files.values():
            all_chains.update(chains)
        
        if self.specific_chain and self.specific_chain in all_chains:
            selected_chain = self.specific_chain
        elif all_chains:
            selected_chain = random.choice(list(all_chains))
        else:
            self.logger.warning("No valid chains found")
            return gmm_params
            
        self.logger.info(f"Using GMM parameters from chain {selected_chain}")

        # Load each GMM file
        for pair_type in self.pair_types:
            file_path = os.path.join(self.output_dir, f"gmm_fit_{pair_type}_chain_{selected_chain}.json")
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_data = json.load(f)
                        # Validate and normalize GMM data
                        if self._validate_gmm_data(gmm_data):
                            gmm_params[pair_type] = gmm_data
                            self.logger.debug(f"Loaded GMM parameters for {pair_type} from {file_path}")
                        else:
                            self.logger.warning(f"Invalid GMM data for {pair_type}")
                except Exception as e:
                    self.logger.error(f"Error loading {file_path}: {e}")
            else:
                self.logger.debug(f"GMM file not found: {file_path}")
                
        return gmm_params

    def _validate_gmm_data(self, gmm_data: dict) -> bool:
        """Validate GMM data structure."""
        required_keys = ['n_components', 'means', 'weights']
        if not all(key in gmm_data for key in required_keys):
            return False
        
        # Check for variance information
        variance_keys = ['variances', 'variances_', 'std', 'precisions']
        if not any(key in gmm_data for key in variance_keys):
            self.logger.warning("No variance information in GMM data")
        
        return True

    def print_gmm_stats(self):
        """Print mean and std deviation of loaded GMMs for diagnostics."""
        if not self.gmm_params:
            return

        for pair_type, gmm_info in self.gmm_params.items():
            if not gmm_info:
                continue
                
            try:
                means = np.array(gmm_info['means']).flatten()
                weights = np.array(gmm_info['weights']).flatten()
                weights = weights / np.sum(weights)
                
                # Get variances
                variances = self._extract_variances(gmm_info)
                
                # Calculate mixture statistics
                mixture_mean = np.sum(means * weights)
                mixture_var = np.sum(weights * (variances + (means - mixture_mean)**2))
                mixture_std = np.sqrt(mixture_var)

                self.logger.info(
                    f"{pair_type} GMM -> mean: {mixture_mean:.3f}, std: {mixture_std:.3f}, "
                    f"n_components: {gmm_info.get('n_components', len(means))}"
                )
            except Exception as e:
                self.logger.warning(f"Error computing stats for {pair_type}: {e}")

    def _extract_variances(self, gmm_info: dict) -> np.ndarray:
        """Extract variances from GMM info, handling different formats."""
        n_components = int(gmm_info.get('n_components', 1))
        
        if 'variances' in gmm_info:
            return np.array(gmm_info['variances']).flatten()
        elif 'variances_' in gmm_info:
            return np.array(gmm_info['variances_']).flatten()
        elif 'std' in gmm_info:
            std_devs = np.array(gmm_info['std']).flatten()
            return std_devs ** 2
        elif 'precisions' in gmm_info:
            precisions = np.array(gmm_info['precisions']).flatten()
            return 1.0 / np.maximum(precisions, 1e-12)
        else:
            # Default variance
            self.logger.debug(f"Using default variance for GMM")
            return np.ones(n_components) * 0.1

    ###########################################################################
    # Sigma Sampling
    ###########################################################################

    def sample_sigma_values(self, max_attempts: int = 100) -> Dict[str, float]:
        """
        Sample sigma values from GMMs or priors.
        """
        if self.is_first_sampler:
            return self._sample_sigma_first_sampler()
        else:
            return self._sample_sigma_gmm(max_attempts)

    def _sample_sigma_first_sampler(self) -> Dict[str, float]:
        """Sample sigma for first sampler based on prior type."""
        sigma = {}
        
        for pt in self.pair_types:
            min_val, max_val = self.sigma_ranges[pt]
            
            if self.prior_type == "jeffreys":
                # Sample from log-uniform distribution
                log_min, log_max = np.log(max(min_val, 1e-6)), np.log(max_val)
                sigma[pt] = np.exp(np.random.uniform(log_min, log_max))
            elif self.prior_type == "inverse_gamma":
                # Sample from inverse gamma
                alpha, beta = 2.0, 0.5
                variance = 1.0 / np.random.gamma(alpha, 1.0/beta)
                sigma[pt] = np.sqrt(variance)
            elif self.prior_type == "gamma":
                # Sample from gamma distribution
                shape, scale = 3.0, (min_val + max_val) / 6.0
                sigma[pt] = np.random.gamma(shape, scale)
            else:  # uniform
                sigma[pt] = np.random.uniform(min_val, max_val)
            
            # Ensure within bounds
            sigma[pt] = np.clip(sigma[pt], min_val, max_val)

        self.logger.info(f"Sampled sigma values for first sampler ({self.prior_type}): {sigma}")
        return sigma

    def _sample_sigma_gmm(self, max_attempts: int) -> Dict[str, float]:
        """Sample sigma values from GMM fits."""
        sigma = {}
        
        for pt in self.pair_types:
            if self.gmm_params and pt in self.gmm_params and self.gmm_params[pt]:
                sampled_value = self._sample_from_gmm(pt, self.gmm_params[pt], max_attempts)
                if sampled_value is not None:
                    sigma[pt] = sampled_value
                else:
                    # Fallback to default if sampling fails
                    sigma[pt] = self.default_sigma[pt]
                    self.logger.warning(f"Using default sigma for {pt}: {sigma[pt]:.3f}")
            else:
                # No GMM available, use default
                sigma[pt] = self.default_sigma[pt]
                self.logger.debug(f"No GMM for {pt}, using default: {sigma[pt]:.3f}")
        
        return sigma

    def _sample_from_gmm(self, pair_type: str, gmm_info: dict, max_attempts: int) -> Optional[float]:
        """Sample from GMM with optional bounds checking."""
        try:
            n_components = int(gmm_info['n_components'])
            means = np.array(gmm_info['means']).flatten()
            weights = np.array(gmm_info['weights']).flatten()
            weights = weights / np.sum(weights)
            variances = self._extract_variances(gmm_info)
            
            min_val, max_val = self.sigma_ranges[pair_type]
            
            # Try to sample within bounds
            for _ in range(max_attempts):
                component = np.random.choice(n_components, p=weights)
                mean = means[component]
                std = np.sqrt(max(variances[component], 1e-12))
                
                sample = np.random.normal(mean, std)
                
                # Accept if within reasonable range (3x bounds for flexibility)
                if min_val/3 <= sample <= max_val*3:
                    return np.clip(sample, min_val, max_val)
            
            # If all attempts fail, sample from truncated distribution
            # Use the component with highest weight
            best_component = np.argmax(weights)
            return np.clip(means[best_component], min_val, max_val)
            
        except Exception as e:
            self.logger.error(f"Error sampling from GMM for {pair_type}: {e}")
            return None

    ###########################################################################
    # Negative Log Priors
    ###########################################################################

    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior without hard boundaries."""
        # First sampler uses specified prior
        if self.is_first_sampler:
            return self._calc_first_sampler_log_prior(sigma)

        # Subsequent samplers use GMM-based priors
        total_log_prob = 0.0
        
        for pt, val in sigma.items():
            if val <= 0:  # Physical constraint
                return np.inf
            
            # Add soft penalty for values far outside expected range
            min_val, max_val = self.sigma_ranges[pt]
            penalty = 0.0
            if val < min_val/3 or val > max_val*3:
                # Quadratic penalty for extreme values
                if val < min_val/3:
                    penalty = ((min_val/3 - val) / min_val) ** 2
                else:
                    penalty = ((val - max_val*3) / max_val) ** 2
                penalty *= 10  # Scale factor
            
            lp = self._calculate_gmm_log_prob(val, pt)
            if not np.isfinite(lp):
                self.logger.debug(f"Log prob for {pt} was -inf -> using penalty")
                lp = -20 - penalty  # Large but finite penalty
            else:
                lp -= penalty
            
            total_log_prob += lp

        return -total_log_prob

    def _calc_first_sampler_log_prior(self, sigma: Dict[str, float]) -> float:
        """Compute the negative log prior for the first sampler."""
        log_prior = 0.0
        
        for pt, val in sigma.items():
            if val <= 0:
                return np.inf
            
            if self.prior_type == "jeffreys":
                log_prior -= np.log(val)
            elif self.prior_type == "inverse_gamma":
                alpha, beta = 2.0, 0.5
                variance = val ** 2
                # Inverse gamma on variance
                log_prior += (alpha * np.log(beta) - np.log(np.math.gamma(alpha))
                             - (alpha + 1) * np.log(variance) - beta / variance
                             + np.log(2 * val))  # Jacobian
            elif self.prior_type == "gamma":
                shape, scale = 3.0, 1.0
                log_prior += gamma.logpdf(val, a=shape, scale=scale)
            else:  # uniform
                min_val, max_val = self.sigma_ranges[pt]
                if min_val <= val <= max_val:
                    log_prior -= np.log(max_val - min_val)
                else:
                    return np.inf
        
        return -log_prior

    def _calculate_gmm_log_prob(self, val: float, pair_type: str) -> float:
        """Calculate log probability for a sigma under the GMM."""
        if not self.gmm_params or pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
            return -np.inf

        gmm = self.gmm_params[pair_type]
        
        try:
            n_components = int(gmm['n_components'])
            means = np.array(gmm['means']).flatten()
            weights = np.array(gmm['weights']).flatten()
            weights = weights / np.sum(weights)
            variances = self._extract_variances(gmm)
            
            # Ensure arrays have correct size
            if len(means) != n_components or len(weights) != n_components:
                self.logger.error(f"Size mismatch in GMM for {pair_type}")
                return -np.inf
            
            # Calculate log probability for each component
            variances = np.maximum(variances, 1e-12)
            diff_squared = (val - means) ** 2
            
            # Log of Gaussian PDF for each component
            log_probs = (np.log(weights) 
                        - 0.5 * np.log(2 * np.pi * variances)
                        - 0.5 * diff_squared / variances)
            
            # Use log-sum-exp trick for numerical stability
            max_log_prob = np.max(log_probs)
            if not np.isfinite(max_log_prob):
                return -np.inf
            
            log_sum = max_log_prob + np.log(np.sum(np.exp(log_probs - max_log_prob)))
            
            return log_sum
            
        except Exception as e:
            self.logger.error(f"Error in GMM log prob for {pair_type}: {e}")
            return -np.inf

    ###########################################################################
    # Utility
    ###########################################################################

    def get_available_chains(self) -> List[int]:
        """Get list of available chains with GMM fits."""
        if self.is_first_sampler or not self.output_dir or not os.path.exists(self.output_dir):
            return []

        pattern = re.compile(r"gmm_fit_[A-Z]+_chain_(\d+)\.json")
        chains = set()
        
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                chains.add(int(match.group(1)))
        
        return sorted(chains)

###############################################################################
# Convenience Function
###############################################################################

def sample_sigmas_from_sequence(
    sampler_sequence: List[str],
    sequence_idx: int,
    base_output_dir: Optional[str] = None,
    chain: Optional[int] = None,
    prior_type: str = "uniform",
    logger: Optional[logging.Logger] = None
) -> Tuple[Dict[str, float], float]:
    """
    Sample sigma values for a specific position in a sampler sequence.

    Args:
        sampler_sequence: The overall sampler sequence.
        sequence_idx: Position in the sequence to sample for.
        base_output_dir: Base directory for reading GMM fits.
        chain: Specific chain number (random if None).
        prior_type: "uniform" or "jeffreys" for the first sampler's prior.
        logger: Optional logger for debugging.

    Returns:
        (sigma_values, negative_log_prior)
    """
    # Create the GMMSigmaProvider
    provider = GMMSigmaProvider(
        sampler_sequence,
        sequence_idx,
        base_output_dir,
        chain,
        sigma_ranges=None,  # or provide a custom dict
        prior_type=prior_type,
        debug_logger=logger
    )

    # Sample sigmas
    sigma_values = provider.sample_sigma_values()
    # Compute negative log prior
    neg_log_prior = provider.calculate_negative_log_prior(sigma_values)

    # Logging / Print summary
    if logger:
        logger.info(f"Sampler: {provider.sampler_name} (idx={sequence_idx}) -> Sigma: {sigma_values}")
        logger.info(f"Negative log prior: {neg_log_prior:.4f}")
    else:
        print(f"Sampler: {provider.sampler_name} (idx={sequence_idx}) -> Sigma: {sigma_values}")
        print(f"Negative log prior: {neg_log_prior:.4f}")

    return sigma_values, neg_log_prior

###############################################################################
# Main Entry Point (CLI Usage)
###############################################################################

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Sample sigma values from a sampler sequence."
    )
    parser.add_argument(
        "--sequence",
        type=str,
        nargs="+",
        default=["PairSampler", "TetramerSampler", "OctetSampler"],
        help="List of sampler names in execution order."
    )
    parser.add_argument(
        "--index",
        type=int,
        default=0,
        help="Position in the sampler sequence to sample."
    )
    parser.add_argument(
        "--dir",
        type=str,
        default=None,
        help="Base output directory containing GMM analysis results."
    )
    parser.add_argument(
        "--chain",
        type=int,
        default=None,
        help="Specific chain number for GMM loading."
    )
    parser.add_argument(
        "--prior",
        type=str,
        default="gamma",
        choices=["uniform", "jeffreys", "gamma"],
        help="Prior type for the first sampler."
    )
    parser.add_argument(
        "--debug_file",
        type=str,
        default=None,
        help="Optional path to a debug log file."
    )
    parser.add_argument(
        "--debug_level",
        type=str,
        default="INFO",
        help="Log level (DEBUG, INFO, WARNING, ERROR)."
    )

    args = parser.parse_args()

    # Convert string to logging level
    numeric_level = getattr(logging, args.debug_level.upper(), logging.INFO)
    logger = setup_logger(__name__, debug_file=args.debug_file, debug_level=numeric_level)

    # Sample sigma
    sample_sigmas_from_sequence(
        sampler_sequence=args.sequence,
        sequence_idx=args.index,
        base_output_dir=args.dir,
        chain=args.chain,
        prior_type=args.prior,
        logger=logger
    )
