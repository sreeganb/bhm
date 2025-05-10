# Calculate tetramer score for a given geometry and sigma values 

import numpy as np
from typing import Dict, Tuple, Optional, Set, List, Union
from scipy.spatial.distance import cdist
from parameters import SystemParameters
from pair_score import ScoringSystem

class TetramerScorer: 
    def __init__(self):
        self.params = SystemParameters()
        self.ps = ScoringSystem(self.params)
        self.sigma_range = {"AA": (0.01, 10.0), "AB": (0.01, 10.0), "BC": (0.01, 10.0), "CC": (0.01, 10.0)}
        self.use_sigma_distribution = False
        
    def get_tetramers(self, positions: Dict[str, np.ndarray], temp: float = 0.999) -> List[Tuple[int, ...]]:
        """Generate tetramers with particle exclusivity and distance-weighted selection."""
        try:
            # Quick validation
            if not all(k in positions and len(positions[k]) > 0 for k in ['A', 'B', 'C']) or len(positions['C']) < 2:
                return []
            
            a_pos, b_pos, c_pos = positions['A'], positions['B'], positions['C']
            used_a, used_b, used_c = set(), set(), set()
            tetramers = []
            
            # A-B distance calculation
            dist_AB = cdist(a_pos, b_pos) / max(0.1, temp)
            
            # Process A particles in order of increasing minimum distance to any B
            for a_idx in np.argsort(np.min(dist_AB, axis=1)):
                if a_idx in used_a:
                    continue
                    
                # Find available B particles and their selection probabilities
                b_mask = np.ones(len(b_pos), dtype=bool)
                for idx in used_b:
                    b_mask[idx] = False
                    
                if not np.any(b_mask):
                    continue
                    
                # Calculate B selection probabilities
                b_dists = dist_AB[a_idx].copy()
                b_dists[~b_mask] = np.inf
                min_b_dist = np.min(b_dists)
                if np.isinf(min_b_dist):
                    continue
                    
                b_probs = np.exp(-(b_dists - min_b_dist))
                b_probs = b_probs / np.sum(b_probs)
                
                # Select B particle
                available_b = np.where(b_mask)[0]
                b_idx = np.random.choice(available_b, p=b_probs[b_mask])
                
                # Find available C particles
                c_mask = np.ones(len(c_pos), dtype=bool)
                for idx in used_c:
                    c_mask[idx] = False
                    
                if np.sum(c_mask) < 2:  # Need at least 2 C particles
                    continue
                    
                # Calculate C selection probabilities
                c_dists = cdist(b_pos[b_idx].reshape(1, -1), c_pos)[0] / max(0.1, temp)
                c_dists[~c_mask] = np.inf
                min_c_dist = np.min(c_dists)
                if np.isinf(min_c_dist):
                    continue
                    
                c_probs = np.exp(-(c_dists - min_c_dist))
                c_probs = c_probs / np.sum(c_probs)
                
                # Select two C particles
                available_c = np.where(c_mask)[0]
                c_indices = np.random.choice(
                    available_c, 
                    size=2, 
                    replace=False, 
                    p=c_probs[c_mask] / np.sum(c_probs[c_mask])
                )
                
                # Add tetramer and mark particles as used
                tetramers.append((a_idx, b_idx, c_indices[0], c_indices[1]))
                used_a.add(a_idx)
                used_b.add(b_idx)
                used_c.update(c_indices)
                
                # Stop if we have enough tetramers
                if len(tetramers) >= min(len(a_pos), len(b_pos), len(c_pos) // 2):
                    break
                    
            return tetramers
            
        except Exception as e:
            print(f"Error in tetramer generation: {e}")
            return []
        
    def calculate_tetramer_scores_batch(self, positions, tetramers, sig, debug_logging=True):
        """Calculate scores for all tetramers with optional debugging."""
        if not tetramers:
            return np.array([], dtype=np.float32)
            
        n_tetramers = len(tetramers)
        
        # Extract all tetramer indices efficiently
        a_indices = np.array([t[0] for t in tetramers], dtype=np.int32)
        b_indices = np.array([t[1] for t in tetramers], dtype=np.int32)
        c1_indices = np.array([t[2] for t in tetramers], dtype=np.int32)
        c2_indices = np.array([t[3] for t in tetramers], dtype=np.int32)
        
        # Get all positions in single vectorized operations
        pos_a = positions['A'][a_indices]
        pos_b = positions['B'][b_indices]
        pos_c1 = positions['C'][c1_indices]
        pos_c2 = positions['C'][c2_indices]
        
        # Calculate all distances at once
        ab_dists = np.sqrt(np.sum((pos_a - pos_b)**2, axis=1))
        bc1_dists = np.sqrt(np.sum((pos_b - pos_c1)**2, axis=1))
        bc2_dists = np.sqrt(np.sum((pos_b - pos_c2)**2, axis=1))
        cc_dists = np.sqrt(np.sum((pos_c1 - pos_c2)**2, axis=1))
        
        # Cache target distances for performance
        ab_target = self.params.pair_distances['AB']
        bc_target = self.params.pair_distances['BC']
        cc_target = self.params.pair_distances['CC']
        
        # Calculate scores individually to enable detailed logging
        ab_scores = ((ab_dists - ab_target)**2)/(2*sig['AB']**2) + np.log(2 * np.pi *sig['AB'])
        bc1_scores = ((bc1_dists - bc_target)**2)/(2*sig['BC']**2) + np.log(2 * np.pi *sig['BC'])
        bc2_scores = ((bc2_dists - bc_target)**2)/(2*sig['BC']**2) + np.log(2 * np.pi *sig['BC'])
        cc_scores = ((cc_dists - cc_target)**2)/(2*sig['CC']**2) + np.log(2 * np.pi *sig['CC'])
        
        scores = ab_scores + bc1_scores + bc2_scores + cc_scores
        
        # Log all distances, sigmas, and scores if debug_logging is enabled
        if debug_logging:
            with open("logfile_tetramer.txt", "a") as log_file:
                # Write header
                log_file.write("\n" + "="*80 + "\n")
                log_file.write(f"TETRAMER SCORING LOG - {len(tetramers)} tetramers\n")
                log_file.write("="*80 + "\n\n")
                
                # Write sigma values
                log_file.write("Sigma values used:\n")
                log_file.write(f"  AB: {sig['AB']:.4f}\n")
                log_file.write(f"  BC: {sig['BC']:.4f}\n")
                log_file.write(f"  CC: {sig['CC']:.4f}\n\n")
                
                # Write target distances
                log_file.write("Target distances:\n")
                log_file.write(f"  AB: {ab_target:.4f}\n")
                log_file.write(f"  BC: {bc_target:.4f}\n")
                log_file.write(f"  CC: {cc_target:.4f}\n\n")
                
                # Write detailed information for each tetramer
                log_file.write("Detailed tetramer information:\n")
                log_file.write(f"{'Tetramer':^10} | {'A idx':^6} | {'B idx':^6} | {'C1 idx':^6} | {'C2 idx':^6} | " +
                            f"{'AB dist':^8} | {'BC1 dist':^8} | {'BC2 dist':^8} | {'CC dist':^8} | " +
                            f"{'AB score':^8} | {'BC1 score':^8} | {'BC2 score':^8} | {'CC score':^8} | {'Total':^8}\n")
                log_file.write("-"*150 + "\n")
                
                for i in range(len(tetramers)):
                    log_file.write(f"{i:^10} | {a_indices[i]:^6} | {b_indices[i]:^6} | {c1_indices[i]:^6} | {c2_indices[i]:^6} | " +
                                f"{ab_dists[i]:8.3f} | {bc1_dists[i]:8.3f} | {bc2_dists[i]:8.3f} | {cc_dists[i]:8.3f} | " +
                                f"{ab_scores[i]:8.3f} | {bc1_scores[i]:8.3f} | {bc2_scores[i]:8.3f} | {cc_scores[i]:8.3f} | {scores[i]:8.3f}\n")
                
                # Write summary statistics
                log_file.write("\nSummary Statistics:\n")
                log_file.write(f"  Mean AB distance: {np.mean(ab_dists):.4f} (target: {ab_target:.4f})\n")
                log_file.write(f"  Mean BC1 distance: {np.mean(bc1_dists):.4f} (target: {bc_target:.4f})\n")
                log_file.write(f"  Mean BC2 distance: {np.mean(bc2_dists):.4f} (target: {bc_target:.4f})\n")
                log_file.write(f"  Mean CC distance: {np.mean(cc_dists):.4f} (target: {cc_target:.4f})\n\n")
                
                log_file.write(f"  Mean AB score: {np.mean(ab_scores):.4f}\n")
                log_file.write(f"  Mean BC1 score: {np.mean(bc1_scores):.4f}\n")
                log_file.write(f"  Mean BC2 score: {np.mean(bc2_scores):.4f}\n")
                log_file.write(f"  Mean CC score: {np.mean(cc_scores):.4f}\n")
                log_file.write(f"  Mean total score: {np.mean(scores):.4f}\n")
                log_file.write(f"  Total tetramer score: {np.sum(scores):.4f}\n")
                log_file.write("="*80 + "\n")
        
        return scores
#-----------------------------------------------------------------------
    
    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        tetramers: List[Tuple[int, ...]],
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
        debug: bool = False
    ) -> Tuple[float, float, float, float]:
        """Calculate the negative log posterior with optional debugging."""
        # Use sigma if provided, otherwise use class sigma
        sigma = sig if sig is not None else self.sigma
        
        # 1) Build set of tetramer pairs
        tetramer_pairs = set()
        
        # Process all tetramers
        for t_idx, (a_idx, b_idx, c1_idx, c2_idx) in enumerate(tetramers):
            # Add pairs to excluded set (these are scored separately as tetramer pairs)
            ab_pair = ('A', a_idx, 'B', b_idx)
            bc1_pair = ('B', b_idx, 'C', c1_idx)
            bc2_pair = ('B', b_idx, 'C', c2_idx)
            cc_pair = ('C', c1_idx, 'C', c2_idx)
            
            tetramer_pairs.add(ab_pair)
            tetramer_pairs.add(bc1_pair)
            tetramer_pairs.add(bc2_pair)
            tetramer_pairs.add(cc_pair)

        # Log statistics about tetramer pairs only when debug is True
        unique_pairs = len(tetramer_pairs)
        expected_pairs = len(tetramers) * 4
        
        # 2) Calculate score excluding tetramer pairs
        result = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            tetramer_pairs, self.use_sigma_distribution,
            prior_penalty_from_distribution
        )
        
        if isinstance(result, tuple) and len(result) >= 3:
            score, ex_score, pair_score = result[:3]
        else:
            if debug:
                print("Error: Unexpected return value from calculate_score")
            return 0.0, 0.0, 0.0, 0.0

        # 3) Calculate tetramer score using vectorized batch method 
        # Pass debug flag to control debugging in the batch method
        scores_array = self.calculate_tetramer_scores_batch(positions, tetramers, sigma, debug_logging=debug)
        total_tet_score = scores_array.sum()
        
        # Apply weighting and return all score components
        weighted_tet_score = tetramer_weight * total_tet_score
        score += weighted_tet_score

        # Only log score summary when debug is True
        if debug:
            print("\n===== SCORE SUMMARY =====")
            print(f"Exclusion Score: {ex_score:.2f}")
            print(f"Pair Score (non-tetramer): {pair_score:.2f}")
            print(f"Tetramer Score: {weighted_tet_score:.2f}")
            print(f"Total Score: {score:.2f}")

        return score, ex_score, pair_score, weighted_tet_score
    

# Example usage
if __name__ == "__main__":
    # Create an instance of the scorer
    scorer = TetramerScorer()
    
    params = SystemParameters()

    ideal_coords = params.ideal_coordinates

    # Example ideal coordinates    
    tets = scorer.get_tetramers(ideal_coords)
    sample_sigma = {
        'AA': 6.0,
        'AB': 1.0,
        'BC': 1.0,
        'CC': 6.0
    }
    
    
    score, _, _, _ = scorer.neg_log_posterior(
        positions=ideal_coords,
        tetramers=tets,
        sig=sample_sigma,
        exclusion_weight=1.0,
        pair_weight=1.0,
        tetramer_weight=1.0,
        debug=True
    )
    
    