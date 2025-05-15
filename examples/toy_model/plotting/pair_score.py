import numpy as np
from typing import Dict, Tuple, Optional, Set, List, Union
from scipy.spatial.distance import cdist
import os

# Assuming parameters.py exists in the same directory or is in PYTHONPATH
from parameters import SystemParameters

class Priors:
    '''Defines different priors for the sigma parameter.'''
    def __init__(self, prior_type: str):
        self.prior_type = prior_type
    
    def neg_log_prior(self, sigma: Dict[str, float], sigma_range: Dict[str, Tuple[float, float]]) -> float:
        # for all the pair types loop over and find the total negative log prior
        total_prior = 0.0
        for pair_type in sigma.keys():
            if sigma[pair_type] <= sigma_range[pair_type][0] or sigma[pair_type] >= sigma_range[pair_type][1]:
                return np.inf
            
            if self.prior_type == 'uniform':
                total_prior += 0.0
            elif self.prior_type == 'jeffreys':
                total_prior += np.log(sigma[pair_type])
            elif self.prior_type == 'halfcauchy':
                scale = 1.0
                total_prior += np.log(1 + (sigma[pair_type] / scale) ** 2)
            else:
                raise ValueError(f"Unknown prior type: {self.prior_type}")
        return total_prior

class ScoringSystem:
    def __init__(self,
                 system_params: SystemParameters,
                 exclusion_weight: float = 1.0,
                 pair_weight: float = 1.0,
                 prior_type: str = 'jeffreys'):
        self.params = system_params
        self.exclusion_weight = exclusion_weight
        self.pair_weight = pair_weight
        self.priors = Priors(prior_type)
        self.pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')] # Define standard pair types

    def excluded_volume_nll(self, positions: Dict[str, np.ndarray]) -> float:
        """
        Compute the negative log likelihood of the excluded volume penalty.
        """
        score = 0.0
        exvol_sigma = 0.1  # Excluded volume sigma penalty strength
        
        particle_types = list(positions.keys())

        for i in range(len(particle_types)):
            for j in range(i, len(particle_types)):
                type1 = particle_types[i]
                type2 = particle_types[j]
                
                pos1 = positions[type1]
                pos2 = positions[type2]

                if len(pos1) == 0 or len(pos2) == 0:
                    continue

                min_dist = self.params.radii[type1] + self.params.radii[type2]
                
                if type1 == type2:
                    if len(pos1) < 2: # Not enough particles for a pair
                        continue
                    distances = cdist(pos1, pos1)
                    # Use upper triangle mask to avoid counting same particle twice and self-interaction
                    mask = np.triu(np.ones_like(distances, dtype=bool), k=1)
                    relevant_distances = distances[mask]
                else:
                    distances = cdist(pos1, pos2)
                    relevant_distances = distances.flatten()
                
                viol_mask = relevant_distances < min_dist
                
                if np.any(viol_mask):
                    overlaps = min_dist - relevant_distances[viol_mask]
                    score += np.sum((overlaps**2) / (exvol_sigma**2))
        return score
    
    def calculate_pair_scores_matrix(self, pos1: np.ndarray, pos2: np.ndarray, target_dist: float, sigma: float) -> np.ndarray:
        """Vectorized calculation of pair scores between two sets of positions."""
        distances = cdist(pos1, pos2)
        #print(f"Distances between {pos1.shape[0]} and {pos2.shape[0]} particles:\n{distances}")
        return ((distances - target_dist) ** 2) / (2 * sigma**2) + np.log(2 * np.pi * sigma**2)
    
    def calculate_score(
            self,
            pos: Dict[str, np.ndarray],  # input positions
            sig: Dict[str, float],       # input sigma values
            sig_range: Dict[str, Tuple[float, float]],  # input sigma ranges
            excluded_pairs=None,
            use_sigma_distribution=False,
            prior_penalty_from_distribution=0.0,
            debug=False,  # Whether to log detailed pair scoring info
            debug_file="pair_score_debug.csv"  # File to write debug info to
        ) -> Tuple[float, float, float, float]:
            """Calculate the log posterior for the pair-level interactions."""
            # Setup debug file if requested
            if debug:
                debug_fh = open(debug_file, 'w')
                debug_fh.write("Pair Type, Particle1 Type, Particle1 Index, Particle2 Type, Particle2 Index, Distance, Target Distance, Sigma, Score\n")
            else:
                debug_fh = None
                
            # 1) Excluded volume contribution (from BaseMCSampler)
            exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)

            # 2) Pairwise negative log-likelihood (pair-specific)
            pairwise_score = 0.0
            #pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
            pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C')]
            for type1, type2 in pair_types:
                pair_key = f"{type1}{type2}"
                if pair_key in self.params.pair_distances:
                    target_dist = self.params.pair_distances[pair_key]
                    sigma_value = sig[pair_key]  # Use the passed sigma

                    # Calculate pairwise score matrix
                    score_matrix = self.calculate_pair_scores_matrix(
                        pos[type1],
                        pos[type2],
                        target_dist,
                        sigma_value
                    )
                    
                    # For debugging, calculate the actual distances between particles
                    if debug_fh:
                        distance_matrix = cdist(pos[type1], pos[type2])
                        
                    # Mask the lower triangle and diagonal to avoid self-interactions and duplicates
                    if type1 == type2:
                        #just fill inf on the diagonal
                        #score_matrix[np.diag_indices_from(score_matrix)] = np.inf
                        score_matrix[np.tril_indices_from(score_matrix)] = np.inf
                        # print the full score matrix for debugging
                        # write it out in a nice format
                        np.set_printoptions(precision=1, suppress=True, linewidth=120)
                        #print(f"Score matrix for {pair_key}:\n{score_matrix}")

                    # Zero out excluded pairs if provided
                    if excluded_pairs:
                        for i in range(len(pos[type1])):
                            for j in range(len(pos[type2])):
                                if (type1, i, type2, j) in excluded_pairs:
                                    score_matrix[i, j] = np.inf

                    # Find the minimum elements for each row and column
                    row_indices = np.argmin(score_matrix, axis=1)
                    col_indices = np.argmin(score_matrix, axis=0)
                    #print(f"Row indices: {row_indices}, Column indices: {col_indices}")

                    # Create sets of unique pairs from rows and columns
                    row_pairs = {(i, row_indices[i]) for i in range(len(row_indices)) if np.isfinite(score_matrix[i, row_indices[i]])}
                    col_pairs = {(col_indices[j], j) for j in range(len(col_indices)) if np.isfinite(score_matrix[col_indices[j], j])}
                    unique_pairs = row_pairs.union(col_pairs)

                    # Log debug information for each unique pair that contributes to the score
                    if debug_fh:
                        for i, j in unique_pairs:
                            distance = distance_matrix[i, j]
                            pair_score = self.pair_weight * score_matrix[i, j]
                            debug_fh.write(f"{pair_key},{type1},{i},{type2},{j},{distance:.6f},{target_dist:.6f},{sigma_value:.6f},{pair_score:.6f}\n")

                    # Accumulate scores
                    for i, j in unique_pairs:
                        pairwise_score += self.pair_weight * score_matrix[i, j]

            # 3) Prior penalty: if not using the sigma distribution, compute using Priors
            if not use_sigma_distribution:
                prior_penalty = self.priors.neg_log_prior(sig, sig_range)
            else:
                prior_penalty = prior_penalty_from_distribution

            # Clean up debug file if opened
            if debug_fh:
                debug_fh.close()
                print(f"Particle pairing debug information written to {debug_file}")

            total_score = exclusion_score + pairwise_score + prior_penalty
            return total_score, exclusion_score, pairwise_score, prior_penalty
    

if __name__ == "__main__":
    # Example Usage:
    # 1. Create SystemParameters instance (ensure parameters.py is accessible)
    try:
        sp = SystemParameters()
        print("SystemParameters loaded successfully.")
    except ImportError:
        print("Error: Could not import SystemParameters. Make sure parameters.py is in the correct path.")
        print("Using dummy SystemParameters for demonstration.")
        class DummySystemParameters:
            def __init__(self):
                self.radii = {"A": 0.5, "B": 0.6, "C": 0.7}
                self.pair_distances = {"AA": 2.0, "AB": 2.5, "BC": 3.0, "CC": 2.2}
                self.box_size = 10.0
        sp = DummySystemParameters()

    # 2. Create ScoringSystem instance
    scorer = ScoringSystem(system_params=sp, exclusion_weight=1.0, pair_weight=1.0, prior_type='jeffreys')
    print(f"ScoringSystem initialized with prior type: {scorer.priors.prior_type}")
    example_sig = {"AA": 3.0, "AB": 2.0, "BC": 1.0, "CC": 1.5}
    example_sig_range = {"AA": (0.01, 10.0), "AB": (0.01, 10.0), "BC": (0.01, 10.0), "CC": (0.01, 10.0)}
    # 6. Calculate score
    print("\nCalculating score with debug_pairs=True...")
    
    id_coords = {'A': np.array([[ 63.  ,   0.  ,   0.  ],
       [ 44.55,  44.55,   0.  ],
       [  0.  ,  63.  ,   0.  ],
       [-44.55,  44.55,   0.  ],
       [-63.  ,   0.  ,   0.  ],
       [-44.55, -44.55,   0.  ],
       [ -0.  , -63.  ,   0.  ],
       [ 44.55, -44.55,   0.  ]]), 'B': np.array([[ 63.  ,   0.  , -38.5 ],
       [ 44.55,  44.55, -38.5 ],
       [  0.  ,  63.  , -38.5 ],
       [-44.55,  44.55, -38.5 ],
       [-63.  ,   0.  , -38.5 ],
       [-44.55, -44.55, -38.5 ],
       [ -0.  , -63.  , -38.5 ],
       [ 44.55, -44.55, -38.5 ]]), 
     'C': np.array([
       [ 47.00,   0.00, -68.50],
       [ 79.00,   0.00, -68.50],
       
       [ 55.86,  55.86, -68.50],
       [ 33.23,  33.23, -68.50],
       
       [  0.00,  47.00, -68.50],
       [  0.00,  79.00, -68.50],
       
       [-55.86,  55.86, -68.50],
       [-33.23,  33.23, -68.50],
       
       [-47.00,   0.00, -68.50],
       [-79.00,   0.00, -68.50],
       
       [-55.86, -55.86, -68.50],
       [-33.23, -33.23, -68.50],
       
       [  0.00, -47.00, -68.50],
       [  0.00, -79.00, -68.50],
       
       [ 55.86, -55.86, -68.50],
       [ 33.23, -33.23, -68.50],
])}


    total_score, exclusion_s, pairwise_s, prior_s = scorer.calculate_score(
        #pos=sp.ideal_coordinates,
        pos=id_coords,
        sig=example_sig,
        sig_range=example_sig_range
    )

    print("\n--- Final Score Components ---")
    print(f"Total Score: {total_score:.4f}")
    print(f"Exclusion Score: {exclusion_s:.4f}")
    print(f"Pairwise Score: {pairwise_s:.4f}")
    print(f"Prior Penalty: {prior_s:.4f}")
