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
        
    def calculate_pairwise_distances(self, positions: dict, output_file: str = 'pairwise_distances.csv'):
        """
        Calculate all possible pairwise distances for AA, AB, BC, and CC particle types
        and save the results to a CSV file.

        Parameters:
        - positions: dict with keys 'A', 'B', 'C' and values as NumPy arrays of shape (n_particles, 3)
        - output_file: str, path to the output CSV file (default: 'pairwise_distances.csv')

        The CSV file will have columns: 'Type1', 'Index1', 'Type2', 'Index2', 'Distance'
        """
        # Validate input
        required_keys = ['A', 'B', 'C']
        for key in required_keys:
            if key not in positions or not isinstance(positions[key], np.ndarray) or positions[key].ndim != 2 or positions[key].shape[1] != 3:
                raise ValueError(f"positions must contain '{key}' with a NumPy array of shape (n_particles, 3)")

        # Initialize list to hold all pair data
        pair_data = []

        # Helper function to add pairs to the list
        def add_pairs(type1, type2, indices1, indices2, distances):
            for idx1, idx2, dist in zip(indices1, indices2, distances):
                pair_data.append({
                    'Type1': type1,
                    'Index1': idx1,
                    'Type2': type2,
                    'Index2': idx2,
                    'Distance': dist
                })

        # Calculate AA distances (upper triangle only)
        if 'A' in positions:
            pos_A = positions['A']
            dist_AA = cdist(pos_A, pos_A)
            # Get upper triangle indices
            i, j = np.triu_indices(pos_A.shape[0], k=1)
            add_pairs('A', 'A', i, j, dist_AA[i, j])
            min_dist_AA = np.min(dist_AA[i, j])

        # Calculate AB distances
        if 'A' in positions and 'B' in positions:
            pos_A = positions['A']
            pos_B = positions['B']
            dist_AB = cdist(pos_A, pos_B)
            # Get all pairs
            i, j = np.indices(dist_AB.shape)
            add_pairs('A', 'B', i.flatten(), j.flatten(), dist_AB.flatten())
            min_dist_AB = np.min(dist_AB)

        # Calculate BC distances
        if 'B' in positions and 'C' in positions:
            pos_B = positions['B']
            pos_C = positions['C']
            dist_BC = cdist(pos_B, pos_C)
            # Get all pairs
            i, j = np.indices(dist_BC.shape)
            add_pairs('B', 'C', i.flatten(), j.flatten(), dist_BC.flatten())
            min_dist_BC = np.min(dist_BC)

        # Calculate CC distances (upper triangle only)
        if 'C' in positions:
            pos_C = positions['C']
            dist_CC = cdist(pos_C, pos_C)
            # Get upper triangle indices
            i, j = np.triu_indices(pos_C.shape[0], k=1)
            add_pairs('C', 'C', i, j, dist_CC[i, j])
            min_dist_CC = np.min(dist_CC[i, j])

        # Convert to DataFrame
        df = pd.DataFrame(pair_data)

        # Sort the DataFrame for better readability (optional)
        df = df.sort_values(by=['Type1', 'Type2', 'Index1', 'Index2']).reset_index(drop=True)

        # Save to CSV
        df.to_csv(output_file, index=False)
        print(f"Pairwise distances saved to {output_file}")
        
        # Write the minimum distances to the same file
        with open(output_file, 'a') as f:
            f.write("\nMinimum Distances:\n")
            if 'A' in positions:
                f.write(f"AA: {min_dist_AA:.4f}\n")
            if 'A' in positions and 'B' in positions:
                f.write(f"AB: {min_dist_AB:.4f}\n")
            if 'B' in positions and 'C' in positions:
                f.write(f"BC: {min_dist_BC:.4f}\n")
            if 'C' in positions:
                f.write(f"CC: {min_dist_CC:.4f}\n")
    
    def get_octets(self, positions: Dict[str, np.ndarray], tetramers=None, temp: float = 0.99999) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
        """
        Group tetramers into octets with temperature-based selection, ignoring periodic boundaries.
        """
        if tetramers is None:
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

        return octets

    def calculate_octet_scores_batch(self, positions, octets, sig, debug_logging=True):
        """Calculate scores for all octets (pairs of adjacent tetramers) with optional debugging."""
        if not octets:
            return np.array([], dtype=np.float32)
        
        n_octets = len(octets)
        
        # Extract indices for A particles in each tetramer pair
        a1_indices = np.array([tet1[0] for tet1, tet2 in octets], dtype=np.int32)
        a2_indices = np.array([tet2[0] for tet1, tet2 in octets], dtype=np.int32)
        
        # Get A particle positions
        pos_a1 = positions['A'][a1_indices]
        pos_a2 = positions['A'][a2_indices]
        
        # Calculate distances between A1 and A2 in each octet
        aa_dists = np.sqrt(np.sum((pos_a1 - pos_a2)**2, axis=1))
        
        # Define target distance for A-A between adjacent tetramers
        # For an octagon with radius 80.0: 2 * 80 * sin(pi/8) ≈ 61.24
        aa_inter_tetramer_target = 2 * 80.0 * np.sin(np.pi / 8)
        
        # Calculate scores
        aa_scores = ((aa_dists - aa_inter_tetramer_target)**2) / (2 * sig['AA']**2) + np.log(2 * np.pi * sig['AA'])
        
        # Total octet scores (currently just A-A; add more specific pairs if needed)
        scores = aa_scores
        
        # Debug logging
        if debug_logging:
            with open("logfile_octet.txt", "a") as log_file:
                log_file.write("\n" + "="*80 + "\n")
                log_file.write(f"OCTET SCORING LOG - {n_octets} octets\n")
                log_file.write("="*80 + "\n\n")
                log_file.write(f"Sigma AA: {sig['AA']:.4f}\n")
                log_file.write(f"Target A-A distance (inter-tetramer): {aa_inter_tetramer_target:.4f}\n\n")
                log_file.write(f"{'Octet':^10} | {'A1 idx':^6} | {'A2 idx':^6} | {'A-A dist':^8} | {'Score':^8}\n")
                log_file.write("-"*80 + "\n")
                for i in range(n_octets):
                    log_file.write(f"{i:^10} | {a1_indices[i]:^6} | {a2_indices[i]:^6} | {aa_dists[i]:8.3f} | {scores[i]:8.3f}\n")
                log_file.write("\nSummary:\n")
                log_file.write(f"  Mean A-A distance: {np.mean(aa_dists):.4f}\n")
                log_file.write(f"  Mean score: {np.mean(scores):.4f}\n")
                log_file.write(f"  Total octet score: {np.sum(scores):.4f}\n")
                log_file.write("="*80 + "\n")
        
        return scores

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        tetramers: List[Tuple[int, ...]],
        octets: List[Tuple[Tuple[int, ...], Tuple[int, ...]]],  # List of (tet1, tet2) pairs
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
        octet_weight: float = 1.0,
        debug: bool = False
    ) -> Tuple[float, float, float, float, float]:
        sigma = sig if sig is not None else self.sigma
        
        # Exclude tetramer pairs from pair scoring
        tetramer_pairs = set()
        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
            tetramer_pairs.add(('A', a_idx, 'B', b_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c1_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c2_idx))
            tetramer_pairs.add(('C', c1_idx, 'C', c2_idx))
        
        # Exclude octet A-A pairs from pair scoring
        octet_pairs = set()
        for tet1, tet2 in octets:
            octet_pairs.add(('A', tet1[0], 'A', tet2[0]))
        
        excluded_pairs = tetramer_pairs.union(octet_pairs)
        
        # Pair-level score
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions, sigma, self.sigma_range, excluded_pairs, False, 0.0
        )
        
        # Tetramer score
        #tet_scores = self.calculate_tetramer_scores_batch(positions, tetramers, sigma, debug)
        #tetramer_score = tetramer_weight * tet_scores.sum()
        tetramer_score = 0.0
        
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