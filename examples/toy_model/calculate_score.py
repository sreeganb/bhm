import numpy as np
from typing import Dict, Tuple
from scipy.spatial.distance import cdist
from parameters import SystemParameters

# Assuming SystemParameters is defined elsewhere and imported, for example:
# from system_parameters import SystemParameters

class ScoreCalculator:
    def __init__(self,
                 exclusion_weight: float = 1.0,
                 pair_weight: float = 1.0,
                 sigma: Dict[str, float] = None,
                 sigma_range: Dict[str, Tuple[float, float]] = None,
                 prior_type: str = 'jeffreys'):
        """
        Initialize the ScoreCalculator.

        Parameters:
          params: A SystemParameters object containing radii, pair distances, etc.
          exclusion_weight: Weight for the excluded volume contribution.
          pair_weight: Weight for the pairwise interaction contribution.
          sigma: Dictionary mapping pair types (e.g., 'AA', 'AB', etc.) to sigma values.
                 Defaults to a fixed value if not provided.
          sigma_range: Dictionary mapping pair types to (lb, ub) for sigma.
                       Defaults to (0.1, 2.0) if not provided.
          prior_type: Type of prior to use ('uniform', 'jeffreys', or 'halfcauchy').
        """
        self.params = SystemParameters()  # Initialize with default parameters.
        self.exclusion_weight = exclusion_weight
        self.pair_weight = pair_weight
        self.prior_type = prior_type

        # Default sigma values if none provided.
        self.sigma = sigma or {'AA': 0.9, 'AB': 0.6, 'BC': 0.4, 'CC': 0.6}
        # Default sigma ranges.
        self.sigma_range = sigma_range or {'AA': (0.1, 2.0),
                                           'AB': (0.1, 2.0),
                                           'BC': (0.1, 2.0),
                                           'CC': (0.1, 2.0)}
        # Define pair types to consider.
        self.pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]

    def excluded_volume_nll(self, positions: Dict[str, np.ndarray], step: int = 0) -> float:
        """
        Compute the negative log likelihood (NLL) for the excluded volume penalty,
        normalized for a Gaussian distribution. Uses self.params.radii to determine
        the minimum allowed distance between atom types.
        
        Parameters:
          positions: Dictionary mapping atom types (e.g., 'A', 'B', etc.) to their
                     coordinates (NumPy array of shape [n_atoms, 3]).
          step: Optional step parameter (not used here but kept for compatibility).
        
        Returns:
          The total excluded volume penalty score.
        """
        score = 0.0
        n_violations = 0  # Count the number of violations
        exvol_penalty = 0.1  # Excluded volume sigma penalty strength

        # Loop over all pair combinations
        for type1 in positions:
            for type2 in positions:
                # Only count each unique pair once
                if type1 <= type2:
                    pos1, pos2 = positions[type1], positions[type2]
                    # Minimum allowed distance based on the sum of the radii
                    min_dist = self.params.radii[type1] + self.params.radii[type2]
                    distances = cdist(pos1, pos2)
                    if type1 == type2:
                        mask = np.triu(np.ones_like(distances), k=1)
                        viol_mask = (distances < min_dist) & (mask > 0)
                    else:
                        viol_mask = distances < min_dist

                    if np.any(viol_mask):
                        # Calculate squared overlap where violations occur.
                        overlap = (min_dist - distances[viol_mask]) ** 2
                        n_violations += overlap.size
                        score += np.sum(overlap / (exvol_penalty ** 2))
        return score

    def calculate_pair_scores_matrix(self,
                                     pos1: np.ndarray,
                                     pos2: np.ndarray,
                                     target_dist: float,
                                     sigma_value: float) -> np.ndarray:
        """
        Calculate the pairwise score matrix using a Gaussian model.
        
        For each pair of atoms, the score is computed as:
            (distance - target_dist)^2 / (2 * sigma_value^2)
        
        Parameters:
          pos1, pos2: Arrays of atom positions (shape: [n_atoms, 3]).
          target_dist: The target distance for the pair type.
          sigma_value: The sigma value for the Gaussian likelihood.
        
        Returns:
          A matrix of scores for each atom pair.
        """
        diff = pos1[:, None, :] - pos2[None, :, :]
        dists = np.linalg.norm(diff, axis=2)
        score_matrix = ((dists - target_dist) ** 2) / (2 * sigma_value ** 2) + np.log(2 * np.pi * sigma_value ** 2)
        return score_matrix

    def log_prior_penalty(self, sigma_value: float, sigma_range: Tuple[float, float]) -> float:
        """
        Compute the log prior penalty for a given sigma value.
        
        Parameters:
          sigma_value: The sigma value to evaluate.
          sigma_range: A tuple (lb, ub) specifying the allowed range for sigma.
        
        Returns:
          The log prior penalty.
        """
        lb, ub = sigma_range
        if self.prior_type == 'uniform':
            if sigma_value <= lb or sigma_value >= ub:
                return np.inf
            return 0.0
        elif self.prior_type == 'jeffreys':
            if sigma_value <= lb or sigma_value >= ub:
                return np.inf
            return np.log(sigma_value)
        elif self.prior_type == 'halfcauchy':
            scale = 1.0
            if sigma_value <= 0:
                return np.inf
            return np.log(1 + (sigma_value / scale) ** 2)
        else:
            raise ValueError("Unknown prior type: " + self.prior_type)

    def calculate_score(self, positions: Dict[str, np.ndarray]) -> Tuple[float, float, float, float]:
        """
        Calculate the total negative log posterior score for the given configuration.
        This includes:
          1) The excluded volume penalty.
          2) The pairwise interaction penalty.
          3) The prior penalty for the sigma parameters.
        
        Parameters:
          positions: Dictionary mapping atom types (e.g., 'A', 'B', 'C') to their
                     coordinates (NumPy array of shape [n_atoms, 3]).
        
        Returns:
          A tuple (total_score, exclusion_score, pairwise_score, prior_penalty).
        """
        # 1) Excluded volume contribution.
        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(positions)
        
        # 2) Pairwise interaction contribution.
        pairwise_score = 0.0
        for type1, type2 in self.pair_types:
            pair_key = f"{type1}{type2}"
            if pair_key in self.params.pair_distances:
                target_dist = self.params.pair_distances[pair_key]
                sigma_value = self.sigma[pair_key]
                if type1 not in positions or type2 not in positions:
                    continue
                pos1 = positions[type1]
                pos2 = positions[type2]
                score_matrix = self.calculate_pair_scores_matrix(pos1, pos2, target_dist, sigma_value)
                
                # For same-type interactions, mask out self-interactions.
                if type1 == type2:
                    mask = np.triu(np.ones_like(score_matrix), k=1)
                    score_matrix = score_matrix * mask
                
                # For each matrix, choose the minimal score per row and per column.
                row_indices = np.argmin(score_matrix, axis=1)
                col_indices = np.argmin(score_matrix, axis=0)
                row_pairs = {(i, row_indices[i]) for i in range(len(row_indices))}
                col_pairs = {(col_indices[j], j) for j in range(len(col_indices))}
                unique_pairs = row_pairs.union(col_pairs)
                for i, j in unique_pairs:
                    pairwise_score += self.pair_weight * score_matrix[i, j]
        
        # 3) Prior penalty contribution.
        prior_penalty = 0.0
        for pair_key, sigma_val in self.sigma.items():
            prior_penalty += self.log_prior_penalty(sigma_val, self.sigma_range[pair_key])
        
        total_score = exclusion_score + pairwise_score + prior_penalty
        return total_score, exclusion_score, pairwise_score, prior_penalty


# ===== Example Usage =====
if __name__ == '__main__':
    # Create an instance of the ScoreCalculator with the parameters.
    scorer = ScoreCalculator()

    # Example positions dictionary.
    positions = {
        'A': np.array([[0.0, 0.0, 0.0],
                       [30.0, 0.0, 0.0]]),
        'B': np.array([[15.0, 15.0, 0.0]]),
        'C': np.array([[0.0, 30.0, 0.0],
                       [30.0, 30.0, 0.0]])
    }

    total, excl, pair, prior = scorer.calculate_score(positions)
    print("Total Score:", total)
    print("Exclusion Score:", excl)
    print("Pairwise Score:", pair)
    print("Prior Penalty:", prior)
