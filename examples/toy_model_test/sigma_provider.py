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

class GMMSigmaProvider:
    """
    Provides sigma values from GMM fits and calculates negative log priors.
    Works with sampler sequences: 
      - For the first sampler, uses uniform/Jeffreys priors.
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

        Args:
            sampler_sequence: List of sampler names in execution order.
            sequence_idx: Current position in the sampler sequence.
            base_output_dir: Base directory containing analysis results 
                             (default: current_dir/output_analysis).
            specific_chain: Use a specific chain number if available (random if None).
            sigma_ranges: Custom dictionary of (min, max) ranges for each sigma type.
            prior_type: Type of prior for the first sampler ("uniform" or "jeffreys").
            debug_logger: An optional logger for debug messages.
        """
        # Setup logging
        self.logger = debug_logger or logging.getLogger(__name__)

        # Re-seed the RNG using process ID to ensure different seeds per process
        process_seed = os.getpid() + int(np.random.random() * 100000)
        np.random.seed(process_seed)
        random.seed(process_seed)

        self.sequence_manager = SamplerSequenceManager(sampler_sequence)
        self.sequence_idx = sequence_idx
        self.pair_types = ["AA", "AB", "BC"]  # Example, adjust as needed
        self.prior_type = prior_type
        self.specific_chain = specific_chain
        self.params = SystemParameters()

        # Get sampler info
        self.sampler_name, self.current_idx, self.total_count = self.sequence_manager.get_sampler_info(sequence_idx)

        # Set base output directory
        if base_output_dir is None:
            self.base_output_dir = os.path.join(os.getcwd(), "output_analysis")
        else:
            self.base_output_dir = base_output_dir

        # Is this the first sampler in the sequence?
        self.is_first_sampler = (sequence_idx == 0)

        # If not first sampler, load GMM from previous sampler
        if self.is_first_sampler:
            self.logger.info(f"First sampler in sequence: {self.sampler_name}")
            self.gmm_params = None
            self.output_dir = None
        else:
            prev_dir = self.sequence_manager.get_previous_sampler_directory(sequence_idx)
            self.output_dir = os.path.join(self.base_output_dir, prev_dir) if prev_dir else None
            if self.output_dir is not None:
                self.logger.info(f"Loading GMM from previous sampler directory: {prev_dir}")
            self.gmm_params = self._load_gmm_parameters()
            self.print_gmm_stats()  # Moved from old _print_gmm_stats

        # Configure sigma ranges
        self.sigma_ranges = self._setup_sigma_ranges(sigma_ranges)
        self.default_sigma = self._get_default_sigmas()

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
        If user provided a dictionary, use it; otherwise compute defaults.
        """
        # Precompute default min/max for each pair type.
        # Adjust as needed for your system's expected distances.
        default_ranges = {}
        for pair_type in self.params.pair_distances.keys():
            sum_radii = self.params.radii[pair_type[0]] + self.params.radii[pair_type[1]]
            lower_bound = 0.02 * sum_radii
            upper_bound = 0.2 * sum_radii
            default_ranges[pair_type] = (lower_bound, upper_bound)

        if user_ranges is not None:
            # Merge user ranges if given
            for key in user_ranges:
                default_ranges[key] = user_ranges[key]
        return default_ranges

    def _get_default_sigmas(self) -> Dict[str, float]:
        """
        Choose a default sigma for each pair type (e.g. geometric midpoint 
        of the min & max).
        """
        sig = {}
        for pt in self.params.pair_distances.keys():
            min_val, max_val = self.sigma_ranges[pt]
            sigma_val = np.exp((np.log(min_val) + np.log(max_val)) / 2)
            sig[pt] = sigma_val
        return sig

    ###########################################################################
    # GMM Loading and Diagnostics
    ###########################################################################

    def _load_gmm_parameters(self) -> Dict[str, Optional[dict]]:
        """Load GMM parameters for each pair type (AA, AB, BC, etc.) from JSON files."""
        gmm_params = {sigma_type: None for sigma_type in self.pair_types}
        if not self.output_dir or not os.path.exists(self.output_dir):
            self.logger.warning(f"Output directory not found: {self.output_dir}")
            return gmm_params

        pattern = re.compile(r"gmm_fit_(AA|AB|BC|CC)_chain_(\d+)\.json")
        matches = []
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                matches.append((match.group(1), int(match.group(2))))

        if not matches:
            self.logger.warning(f"No GMM fit files found in {self.output_dir}")
            return gmm_params

        # Choose the chain number
        chain_numbers = list({chain_num for _, chain_num in matches})
        selected_chain = (self.specific_chain if self.specific_chain is not None
                          else random.choice(chain_numbers))
        self.logger.info(f"Using GMM parameters from chain {selected_chain}")

        # Load each type of GMM
        for stype in gmm_params.keys():
            file_path = os.path.join(self.output_dir, f"gmm_fit_{stype}_chain_{selected_chain}.json")
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_params[stype] = json.load(f)
                    self.logger.debug(f"Loaded GMM parameters for {stype} from {file_path}")
                except Exception as e:
                    self.logger.error(f"Error loading {file_path}: {e}")
        return gmm_params

    def print_gmm_stats(self):
        """Print mean and std deviation of loaded GMMs for diagnostics."""
        if not self.gmm_params:
            return

        for pair_type, gmm_info in self.gmm_params.items():
            if gmm_info and all(k in gmm_info for k in ['means', 'covariances', 'weights']):
                try:
                    means = np.array(gmm_info['means'])
                    covariances = np.array(gmm_info['covariances'])
                    weights = np.array(gmm_info['weights'])
                    normalized_weights = weights / np.sum(weights)

                    mean_val = np.sum(means * normalized_weights)
                    variance = np.sum(normalized_weights * (means - mean_val)**2) + np.sum(normalized_weights * covariances)
                    std_val = np.sqrt(variance)

                    self.logger.info(f"{pair_type} GMM -> mean: {mean_val:.3f}, std: {std_val:.3f}")
                except Exception as e:
                    self.logger.warning(f"Error computing stats for {pair_type}: {e}")

    ###########################################################################
    # Sigma Sampling
    ###########################################################################

    def sample_sigma_values(self, max_attempts: int = 100) -> Dict[str, float]:
        """
        Sample sigma values from GMMs (if not the first sampler).
        For the first sampler, sample from a uniform or Jeffreys prior.
        """
        if self.is_first_sampler:
            sig = self._sample_sigma_first_sampler()
        else:
            sig = self._sample_sigma_gmm(max_attempts)

        # Enforce bounds strictly
        for pair_type, val in sig.items():
            mn, mx = self.sigma_ranges[pair_type]
            if not (mn <= val <= mx):
                self.logger.debug(
                    f"{pair_type}={val:.3f} outside [{mn:.3f}, {mx:.3f}]. Using default."
                )
                sig[pair_type] = self.default_sigma[pair_type]

        return sig

    def _sample_sigma_first_sampler(self) -> Dict[str, float]:
        """
        Sample sigma values for the first sampler, either from a uniform or a
        Jeffreys prior.
        """
        sigma = {}
        for pt in self.pair_types:
            min_val, max_val = self.sigma_ranges[pt]
            # Sample in log space for better coverage
            if self.prior_type == "jeffreys":
                # Jeffreys prior for sigma: p(sigma) ~ 1/sigma
                # We can approximate by sampling uniformly in log space
                log_min, log_max = np.log(min_val), np.log(max_val)
                log_sigma = np.random.uniform(log_min, log_max)
                sigma[pt] = np.exp(log_sigma)
            elif self.prior_type == "inverse_gamma":
                # Sample from inverse gamma distribution
                # First sample from uniform in log space
                log_min, log_max = np.log(min_val), np.log(max_val)
                log_sigma = np.random.uniform(log_min, log_max)
                sigma_val = np.exp(log_sigma)
                
                # If needed, you could apply rejection sampling here
                # to properly sample from inverse gamma distribution
                
                # For now, just use the uniform sample
                sigma[pt] = sigma_val
            else:
                # "uniform" in log space
                log_min, log_max = np.log(min_val), np.log(max_val)
                log_sigma = np.random.uniform(log_min, log_max)
                sigma[pt] = np.exp(log_sigma)

        self.logger.info(f"Sampled sigma values for first sampler: {sigma}")
        return sigma

    def _sample_sigma_gmm(self, max_attempts: int) -> Dict[str, float]:
        """
        Sample sigma values from GMM fits for subsequent samplers.
        """
        sigma = {}
        for pt in self.pair_types:
            gmm_info = self.gmm_params.get(pt) if self.gmm_params else None

            if gmm_info and all(k in gmm_info for k in ['n_components', 'means', 'covariances', 'weights']):
                sigma[pt] = self._sample_from_gmm(pt, gmm_info, max_attempts)
            else:
                # Fallback to the default sigma if no GMM is loaded
                sigma[pt] = self.default_sigma[pt]
        return sigma

    def _sample_from_gmm(self, pair_type: str, gmm_info: dict, max_attempts: int) -> float:
        """Attempt to sample sigma from a GMM for a specific pair type."""
        min_val, max_val = self.sigma_ranges[pair_type]
        n_components = int(gmm_info['n_components'])
        means = np.asarray(gmm_info['means']).flatten()
        covariances = np.asarray(gmm_info['covariances']).flatten()
        weights = np.asarray(gmm_info['weights']) / np.sum(gmm_info['weights'])

        # Keep attempts limited
        for attempt in range(max_attempts):
            component = np.random.choice(n_components, p=weights)
            mean_value = means[component]
            cov_value = np.maximum(covariances[component], 1e-12)
            std_value = np.sqrt(cov_value)

            candidate = np.random.normal(loc=mean_value, scale=std_value)
            if min_val <= candidate <= max_val:
                return candidate

        # If all attempts fail, fallback
        fallback = self.default_sigma[pair_type]
        self.logger.warning(
            f"{pair_type}: Could not sample within [{min_val:.3f}, {max_val:.3f}] "
            f"after {max_attempts} attempts. Using fallback {fallback:.3f}"
        )
        return fallback

    ###########################################################################
    # Negative Log Priors
    ###########################################################################

    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """
        Calculate negative log prior for the sampled sigma values.
        Uses:
          - uniform/Jeffreys prior for first sampler
          - GMM prior for subsequent samplers
        """
        # Check bounds
        for pt, val in sigma.items():
            mn, mx = self.sigma_ranges[pt]
            if val < mn or val > mx:
                self.logger.debug(f"{pt}={val:.3f} outside [{mn:.3f}, {mx:.3f}] -> prior=inf")
                return np.inf

        # First sampler uses uniform or Jeffreys (with possible penalty)
        if self.is_first_sampler:
            return self._calc_first_sampler_log_prior(sigma)

        # Subsequent samplers use GMM-based priors
        total_log_prob = 0.0
        for pt, val in sigma.items():
            lp = self._calculate_gmm_log_prob(val, pt)
            if not np.isfinite(lp):  # If any one is -inf, entire prior is 0
                self.logger.debug(f"Log prob for {pt} was -inf -> total prior=inf")
                return np.inf
            total_log_prob += lp

        neg_log_prior = -total_log_prob
        return neg_log_prior

    def _calc_first_sampler_log_prior(self, sigma: Dict[str, float]) -> float:
        """Compute the negative log prior for the first sampler (uniform/Jeffreys)."""
        # Just an example structure — you can refine as needed.
        # For a Jeffreys prior: log p(sigma) ~ -log(sigma).
        log_prior = 0.0
        for pt, val in sigma.items():
            # simple example
            if self.prior_type == "jeffreys":
                log_prior += -np.log(val)
            elif self.prior_type == "inverse_gamma":
                # the distribution is defined by a shape and scale parameter
                # the mathematical form is 
                # p(sigma|alpha,beta) = beta**alpha / gamma(alpha) * sigma**(-alpha-1) * exp(-beta/sigma)
                beta = 0.001  # scale parameter
                alpha = 0.001 # shape parameter
                log_prior += (alpha * np.log(beta) - np.log(np.math.gamma(alpha))
                              - (alpha + 1) * np.log(val) - beta / val)
            else:
                log_prior += 0.0  # uniform in log space
        return -log_prior

    def _calculate_gmm_log_prob(self, val: float, pair_type: str) -> float:
        """Calculate log probability for a sigma under the loaded GMM of a given pair type."""
        if not self.gmm_params or pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
            return -np.inf

        gmm = self.gmm_params[pair_type]
        try:
            n_components = gmm['n_components']
            means = np.array(gmm['means']).flatten()
            covariances = np.array(gmm['covariances']).flatten()
            weights = np.array(gmm['weights']) / np.sum(gmm['weights'])

            # For each component, compute univariate Gaussian log prob
            covariances = np.maximum(covariances, 1e-12)  # Avoid zero variance
            diff = val - means
            exponents = -0.5 * (diff**2 / covariances)
            norms = np.log(weights) - 0.5 * np.log(2 * np.pi * covariances)
            component_log_probs = norms + exponents

            # Sum in log space
            max_lp = np.max(component_log_probs)
            log_sum_exp = max_lp + np.log(np.sum(np.exp(component_log_probs - max_lp)))
            return log_sum_exp
        except Exception as e:
            self.logger.error(f"Error in GMM log prob for {pair_type}: {e}")
            return -np.inf

    ###########################################################################
    # Utility
    ###########################################################################

    def get_available_chains(self) -> List[int]:
        """Get list of available chains with GMM fits in the output directory."""
        if self.is_first_sampler or not self.output_dir or not os.path.exists(self.output_dir):
            return []

        pattern = re.compile(r"gmm_fit_.*_chain_(\d+)\.json")
        chain_numbers = []
        for filename in os.listdir(self.output_dir):
            match = pattern.match(filename)
            if match:
                chain_numbers.append(int(match.group(1)))
        return sorted(set(chain_numbers))

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
        default="inverse_gamma",
        choices=["uniform", "jeffreys", "inverse_gamma"],
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