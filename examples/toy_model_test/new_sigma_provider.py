"""
Lightweight sigma provider.

* For the first sampler in a sequence, draw sigmas from simple, non-informative
  priors (uniform / Jeffreys / gamma / half-Cauchy).
* For every later sampler, read Gaussian-mixture fits produced by `fit_gmm.py`,
  sample sigmas from those mixtures, and evaluate the corresponding log prior.

The public surface matches the previous module (`GMMSigmaProvider`,
`sample_sigmas_from_sequence`) so other code keeps working.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import cauchy, gamma

from parameters import SystemParameters  # unchanged dependency

# --------------------------------------------------------------------------- #
# Helper utilities
# --------------------------------------------------------------------------- #


def geometric_midpoint(low: float, high: float) -> float:
    """Return the geometric midpoint between two positive numbers."""
    return float(np.sqrt(low * high))


def occurrence_count(sequence: Sequence[str], name: str, up_to_idx: int) -> int:
    """Number of times `name` appears in sequence[0:up_to_idx+1]."""
    return sum(1 for i in range(up_to_idx + 1) if sequence[i] == name)


def sampler_directory(sequence: Sequence[str], idx: int) -> str:
    """
    Name of the folder that contains results for sampler at position `idx`.

    Matches the convention used by `fit_gmm.py`:
        "{sampler_name}sampler_results_{occurrence}"
    """
    name = sequence[idx]
    occurrence = occurrence_count(sequence, name, idx)
    return f"{name}sampler_results_{occurrence}"


def previous_sampler_directory(
    sequence: Sequence[str], idx: int
) -> Optional[str]:
    """Directory for the sampler immediately preceding position `idx`."""
    if idx == 0:
        return None
    prev_idx = idx - 1
    return sampler_directory(sequence, prev_idx)


def default_sigma_ranges(system_params: SystemParameters) -> Dict[str, Tuple[float, float]]:
    """
    Conservative sigma bounds (2% – 20% of the sum of particle radii).
    These limits only apply to the first sampler and as loose clipping later.
    """
    ranges: Dict[str, Tuple[float, float]] = {}
    for pair in system_params.pair_distances:
        radius_sum = system_params.radii[pair[0]] + system_params.radii[pair[1]]
        ranges[pair] = (0.02 * radius_sum, 0.20 * radius_sum)
    return ranges


@dataclass
class GMMPayload:
    """Container for a single sigma-type GMM read from JSON."""
    means: np.ndarray
    variances: np.ndarray
    weights: np.ndarray

    @property
    def n_components(self) -> int:
        return len(self.means)

    def sample(self) -> float:
        """One draw from the stored mixture."""
        idx = np.random.choice(self.n_components, p=self.weights)
        mean = self.means[idx]
        std = np.sqrt(max(self.variances[idx], 1e-12))
        return float(np.random.normal(mean, std))

    def logpdf(self, value: float) -> float:
        """Log probability density of the mixture at `value`."""
        precisions = 1.0 / np.maximum(self.variances, 1e-12)
        log_components = (
            np.log(self.weights)
            - 0.5 * np.log(2.0 * np.pi / precisions)
            - 0.5 * precisions * (value - self.means) ** 2
        )
        # log-sum-exp
        m = log_components.max()
        return float(m + np.log(np.sum(np.exp(log_components - m))))


# --------------------------------------------------------------------------- #
# Main provider
# --------------------------------------------------------------------------- #


class GMMSigmaProvider:
    """
    Minimal sigma provider.

    * First sampler (sequence_idx == 0):
        Draw sigmas from the requested prior.
    * Later samplers:
        Load GMM JSON files from the previous sampler's analysis directory and
        sample from those mixtures.
    """

    def __init__(
        self,
        sampler_sequence: List[str],
        sequence_idx: int,
        base_output_dir: Optional[str] = None,
        specific_chain: Optional[str] = None,
        sigma_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
        prior_type: str = "uniform",
        debug_logger: Optional[logging.Logger] = None,
    ):
        self.sequence = sampler_sequence
        self.sequence_idx = sequence_idx
        self.prior_type = prior_type.lower()
        self.base_output_dir = base_output_dir or os.path.join(os.getcwd(), "output_analysis")
        self.specific_chain = specific_chain  # optional chain filter
        self.logger = debug_logger or logging.getLogger(__name__)

        self.params = SystemParameters()
        self.pair_types = list(self.params.pair_distances.keys())

        # Determine whether we are the first sampler.
        self.is_first_sampler = (sequence_idx == 0)

        # Bounds and defaults used for clipping.
        self.sigma_ranges = sigma_ranges or default_sigma_ranges(self.params)
        self.default_sigma = {
            pt: geometric_midpoint(*self.sigma_ranges[pt]) for pt in self.pair_types
        }

        # Pre-load mixture models if needed.
        self.gmm_models: Dict[str, GMMPayload] = {}
        if not self.is_first_sampler:
            self._load_gmms()

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    def sample_sigma_values(self) -> Dict[str, float]:
        """
        Draw one sigma per pair type.

        First sampler -> prior draw.
        Others        -> sample from stored GMM; fall back to default if missing.
        """
        if self.is_first_sampler:
            return self._sample_from_prior()

        samples: Dict[str, float] = {}
        for pair in self.pair_types:
            model = self.gmm_models.get(pair)
            if model is None:
                samples[pair] = self.default_sigma[pair]
                continue

            value = model.sample()
            low, high = self.sigma_ranges[pair]
            samples[pair] = float(np.clip(value, low, high))
        return samples

    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """
        Negative log prior at the provided sigma dictionary.
        """
        if self.is_first_sampler:
            return self._neg_log_prior_first_sampler(sigma)

        log_prob = 0.0
        for pair, value in sigma.items():
            if value <= 0:
                return np.inf
            model = self.gmm_models.get(pair)
            if model is None:
                # prefer loose quadratic penalty around default if no GMM
                default = self.default_sigma[pair]
                log_prob += -0.5 * ((value - default) / (0.5 * default)) ** 2
            else:
                log_prob += model.logpdf(value)
        return -log_prob

    def get_available_chains(self) -> List[str]:
        """
        Which chain identifiers have GMM JSON files for the previous sampler.
        """
        if self.is_first_sampler:
            return []
        prev_dir = previous_sampler_directory(self.sequence, self.sequence_idx)
        if prev_dir is None:
            return []

        directory = os.path.join(self.base_output_dir, prev_dir)
        if not os.path.isdir(directory):
            return []

        pattern = re.compile(r"gmm_fit_[A-Z]+_(.+)\.json")
        chains = {
            match.group(1)
            for fname in os.listdir(directory)
            for match in [pattern.match(fname)]
            if match
        }
        return sorted(chains)

    # ------------------------------------------------------------------ #
    # Internal helpers
    # ------------------------------------------------------------------ #

    def _sample_from_prior(self) -> Dict[str, float]:
        """Draw sigmas for the first sampler using the chosen non-informative prior."""
        sigma: Dict[str, float] = {}
        for pair in self.pair_types:
            low, high = self.sigma_ranges[pair]
            if self.prior_type == "uniform":
                value = np.random.uniform(low, high)
            elif self.prior_type == "jeffreys":
                value = np.exp(np.random.uniform(np.log(low), np.log(high)))
            elif self.prior_type == "gamma":
                # Shape/scale chosen so mean lies midway between bounds
                shape = 2.0
                scale = (low + high) / (2.0 * shape)
                value = gamma.rvs(shape, scale=scale)
            elif self.prior_type == "cauchy":
                # Half-Cauchy centred at geometric midpoint with broad scale
                location = geometric_midpoint(low, high)
                scale = (high - low) / 2.0
                value = abs(cauchy.rvs(loc=location, scale=scale))
            else:
                raise ValueError(f"Unknown prior type: {self.prior_type}")

            sigma[pair] = float(np.clip(value, low, high))
        return sigma

    def _neg_log_prior_first_sampler(self, sigma: Dict[str, float]) -> float:
        """Closed-form negative log prior for the first sampler."""
        log_prob = 0.0
        for pair, value in sigma.items():
            low, high = self.sigma_ranges[pair]
            if not (low <= value <= high) or value <= 0:
                return np.inf

            if self.prior_type == "uniform":
                log_prob += -np.log(high - low)
            elif self.prior_type == "jeffreys":
                log_prob += -np.log(value)
            elif self.prior_type == "gamma":
                shape = 2.0
                scale = (low + high) / (2.0 * shape)
                log_prob += gamma.logpdf(value, a=shape, scale=scale)
            elif self.prior_type == "cauchy":
                location = geometric_midpoint(low, high)
                scale = (high - low) / 2.0
                log_prob += np.log(2.0) + cauchy.logpdf(value, loc=location, scale=scale)
        return -log_prob

    def _load_gmms(self) -> None:
        """
        Read GMM JSON files from the previous sampler's directory.
        If `specific_chain` is provided, only use files that contain that token.
        """
        prev_dir_name = previous_sampler_directory(self.sequence, self.sequence_idx)
        if prev_dir_name is None:
            return

        directory = os.path.join(self.base_output_dir, prev_dir_name)
        if not os.path.isdir(directory):
            self.logger.warning("GMM directory not found: %s", directory)
            return

        for pair in self.pair_types:
            pattern = re.compile(rf"gmm_fit_{pair}_(.+)\.json")
            candidates = []
            for fname in os.listdir(directory):
                match = pattern.match(fname)
                if not match:
                    continue
                chain_id = match.group(1)
                if self.specific_chain and self.specific_chain not in chain_id:
                    continue
                candidates.append((chain_id, os.path.join(directory, fname)))

            if not candidates:
                continue

            # Prefer a deterministic choice if specific_chain was given, otherwise pick one at random.
            chain_id, filepath = (
                candidates[0] if self.specific_chain else random.choice(candidates)
            )
            with open(filepath, "r", encoding="utf-8") as fh:
                payload = json.load(fh)

            means = np.array(payload.get("means", []), dtype=float)
            variances = np.array(payload.get("variances", []), dtype=float)
            weights = np.array(payload.get("weights", []), dtype=float)

            if not len(means) or len(means) != len(variances) or len(means) != len(weights):
                self.logger.warning("Malformed GMM file: %s", filepath)
                continue

            weights = np.maximum(weights, 1e-12)
            weights /= weights.sum()

            self.gmm_models[pair] = GMMPayload(means=means, variances=variances, weights=weights)


# --------------------------------------------------------------------------- #
# Convenience wrapper (kept for compatibility with existing scripts)
# --------------------------------------------------------------------------- #


def sample_sigmas_from_sequence(
    sampler_sequence: List[str],
    sequence_idx: int,
    base_output_dir: Optional[str] = None,
    chain: Optional[str] = None,
    prior_type: str = "uniform",
    logger: Optional[logging.Logger] = None,
) -> Tuple[Dict[str, float], float]:
    """
    Wrapper used by command-line tools.

    Returns (sigma_dict, negative_log_prior).
    """
    provider = GMMSigmaProvider(
        sampler_sequence=sampler_sequence,
        sequence_idx=sequence_idx,
        base_output_dir=base_output_dir,
        specific_chain=chain,
        prior_type=prior_type,
        debug_logger=logger,
    )
    sigma = provider.sample_sigma_values()
    nlp = provider.calculate_negative_log_prior(sigma)
    return sigma, nlp


# --------------------------------------------------------------------------- #
# Command-line entry point (kept for parity with original file)
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Sample sigma values for a sampler step.")
    parser.add_argument("--sequence", type=str, nargs="+", required=True, help="Sampler names in order.")
    parser.add_argument("--index", type=int, required=True, help="Position within the sequence.")
    parser.add_argument("--dir", type=str, default=None, help="Base directory for GMM JSON files.")
    parser.add_argument("--chain", type=str, default=None, help="Specific chain identifier to load.")
    parser.add_argument(
        "--prior",
        type=str,
        default="uniform",
        choices=["uniform", "jeffreys", "gamma", "cauchy"],
        help="Prior for the first sampler.",
    )
    parser.add_argument("--log-level", type=str, default="INFO", help="Logging level (e.g. DEBUG).")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    sigma, nlp = sample_sigmas_from_sequence(
        sampler_sequence=args.sequence,
        sequence_idx=args.index,
        base_output_dir=args.dir,
        chain=args.chain,
        prior_type=args.prior,
    )

    print(f"Sigmas: {json.dumps(sigma, indent=2)}")
    print(f"Negative log prior: {nlp:.6f}")


if __name__ == "__main__":
    main()