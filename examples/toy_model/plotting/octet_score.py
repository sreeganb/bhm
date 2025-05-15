import numpy as np
from typing import Dict, Tuple, List
from scipy.spatial.distance import cdist
import pandas as pd
from parameters import SystemParameters
from pair_score import ScoringSystem
from tetramer_score import TetramerScorer

class OctetScorer:
    def __init__(self):
        self.params = SystemParameters()
        self.ps = ScoringSystem(self.params)
        self.ts = TetramerScorer()
        self.sigma_range = {"AA": (0.01, 10.0), "AB": (0.01, 10.0), "BC": (0.01, 10.0), "CC": (0.01, 10.0)}
        self.use_sigma_distribution = False
        self.sigma = {
            'AA': 1.0,  # defaults, user can override
            'AB': 1.0,
            'BC': 1.0,
            'CC': 1.0
        }
    
    def get_octets(self, positions: Dict[str, np.ndarray], temp: float = 0.99999) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
        """
        Group tetramers into octets with temperature-based selection, ignoring periodic boundaries.
        """
 
        tetramers = self.ts.get_tetramers(positions)

        if len(tetramers) < 2:
            return []

        # Compute each tetramer's geometric center
        centers = np.zeros((len(tetramers), 3), dtype=np.float64)
        for i, (a_idx, b_idx, c_idx1, c_idx2) in enumerate(tetramers):
            coords = np.vstack([
                positions['A'][a_idx],
                positions['B'][b_idx],
                positions['C'][c_idx1],
                positions['C'][c_idx2]
            ])
            centers[i] = np.mean(coords, axis=0)

        octets = []
        available = list(range(len(tetramers)))

        while len(available) >= 2:
            idx1 = np.random.choice(available)
            available.remove(idx1)

            indices = np.array(available)
            center1 = centers[idx1]
            deltas = centers[indices] - center1
            distances = np.linalg.norm(deltas, axis=1)

            # Temperature-based probabilities
            probs = np.exp(-distances / temp)
            probs_sum = probs.sum()
            if probs_sum > 1e-10:
                probs /= probs_sum
                idx2_rel = np.random.choice(len(available), p=probs)
                idx2 = available[idx2_rel]
            else:
                idx2 = np.random.choice(available)

            octets.append((tetramers[idx1], tetramers[idx2]))
            available.remove(idx2)

        return octets, tetramers 

    def calculate_octet_scores_batch(self, positions, octets, sig, debug_logging=True):
        """Calculate scores for all octets (pairs of adjacent tetramers) with optional debugging."""
        if not octets:
            return np.array([], dtype=np.float32)
        
        n_octets = len(octets)
        
        # for the two chosen tetramers use the ts.
        
        # Extract indices for A particles in each tetramer pair
        a1_indices = np.array([tet1[0] for tet1, tet2 in octets], dtype=np.int32)
        a2_indices = np.array([tet2[0] for tet1, tet2 in octets], dtype=np.int32)
        
        # Get A particle positions
        pos_a1 = positions['A'][a1_indices]
        pos_a2 = positions['A'][a2_indices]
        
        # Calculate distances between A1 and A2 in each octet
        aa_dists = np.sqrt(np.sum((pos_a1 - pos_a2)**2, axis=1))
        
        # Define target distance for A-A between adjacent tetramers
        aa_inter_tetramer_target = self.params.pair_distances['AA']
        
        # Calculate scores
        aa_scores = ((aa_dists - aa_inter_tetramer_target)**2) / (2 * sig['AA']**2) + np.log(2 * np.pi * sig['AA'])
        
        # Total octet scores (currently just A-A; add more specific pairs if needed)
        scores = aa_scores
        
        # sum the scores for all octets
        final_score = np.sum(scores)
        
        return final_score

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
        octet_weight: float = 1.0,
        debug: bool = False
    ) -> Tuple[float, float, float, float, float]:
        sigma = sig if sig is not None else self.sigma
        
        #excluded_pairs = tetramer_pairs.union(octet_pairs)
        excluded_pairs = None
        
        # Pair-level score
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions, sigma, self.sigma_range, excluded_pairs, False, 0.0
        )
        pair_score *= pair_weight
        ex_score *= exclusion_weight
        
        # Tetramer score
        # Compute the tetramers and octamers in the system
        octets, tetramers = self.get_octets(positions)
        
        tetramer_score = 0.0
        tet_score = self.ts.calculate_tetramer_scores_batch(positions, tetramers, sigma, debug_logging=debug)
        tetramer_score = tetramer_weight * tet_score.sum()
        
        # Octet score
        octet_scores = self.calculate_octet_scores_batch(positions, octets, sigma, debug)
        octet_score = octet_weight * octet_scores.sum()
        
        # Total score
        total_score = score + tetramer_score + octet_score
        
        if debug:
            print("\n===== SCORE SUMMARY =====")
            print(f"Exclusion Score: {ex_score:.2f}")
            print(f"Pair Score (non-tetramer/octet): {pair_score:.2f}")
            print(f"Tetramer Score: {tetramer_score:.2f}")
            print(f"Octet Score: {octet_score:.2f}")
            print(f"Total Score: {total_score:.2f}")
        
        return total_score, ex_score, pair_score, tetramer_score, octet_score