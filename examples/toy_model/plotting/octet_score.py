import numpy as np
from typing import Dict, Tuple, List
from scipy.spatial.distance import cdist
import pandas as pd
from parameters import SystemParameters
from pair_score import ScoringSystem
from tetramer_score import TetramerScorer
import networkx as nx

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

    def get_octets(self, positions: Dict[str, np.ndarray]) -> Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], List[Tuple[int, ...]]]:
        tetramers = self.ts.get_tetramers(positions)
        print("tetramers: ", tetramers)

        if len(tetramers) < 2:
            return [], tetramers

        # 1) Compute geometric centers for each tetramer
        centers = np.zeros((len(tetramers), 3), dtype=np.float64)
        for i, (a_idx, b_idx, c_idx1, c_idx2) in enumerate(tetramers):
            coords = np.vstack([
                positions['A'][a_idx],
                positions['B'][b_idx],
                positions['C'][c_idx1],
                positions['C'][c_idx2]
            ])
            centers[i] = np.mean(coords, axis=0)

        # 2) If only one tetramer or fewer, no pairs possible
        if len(tetramers) < 2:
            return [], tetramers

        # 3) Build a graph of tetramers (nodes) with edge weights = distances
        G = nx.Graph()
        for i_t in range(len(tetramers)):
            G.add_node(i_t)
        for i_t in range(len(tetramers)):
            for j_t in range(i_t + 1, len(tetramers)):
                dist_ij = np.linalg.norm(centers[i_t] - centers[j_t])
                G.add_edge(i_t, j_t, weight=dist_ij)

        # 4) Negate the weights to convert min-weight to max-weight problem
        for u, v, d in G.edges(data=True):
            d['weight'] = -d['weight']

        # 5) Compute the maximum-weight perfect matching (which minimizes original distances)
        matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)

        # 6) Convert the matching (set of edges) into a list of octets
        octets = []
        for i_t, j_t in matching:
            # Sort the node IDs for consistency
            i_t, j_t = sorted([i_t, j_t])
            octets.append((tetramers[i_t], tetramers[j_t]))

        print(f"Octets: {octets}")
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
        print(f"AA distances: {aa_dists}")
        
        # Define target distance for A-A between adjacent tetramers
        aa_inter_tetramer_target = self.params.pair_distances['AA']
        
        # Calculate scores
        aa_scores = ((aa_dists - aa_inter_tetramer_target)**2) / (2 * sig['AA']**2) + np.log(2 * np.pi * sig['AA'])
        
        # Total octet scores (currently just A-A; add more specific pairs if needed)
        scores = aa_scores
        
        return scores

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
        octet_weight: float = 1.0,
        debug: bool = True
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
        tet_score = self.ts.calculate_tetramer_scores_batch(positions, tetramers, 
                                                            sigma, debug_logging=debug)
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
    
if __name__ == "__main__":
    os = OctetScorer()
    sp = SystemParameters()
    ideal_coords = sp.ideal_coordinates
    sigma_fixed = {
        'AA': 1.0,
        'AB': 1.0,
        'BC': 1.0,
        'CC': 1.0
    }
    os.neg_log_posterior(
        ideal_coords,
        prior_penalty_from_distribution=0.0,
        sig=sigma_fixed,
        exclusion_weight=1.0,
        pair_weight=1.0,
        tetramer_weight=1.0,
        octet_weight=1.0,
        debug=True
    )