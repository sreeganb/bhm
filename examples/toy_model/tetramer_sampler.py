#----------------------------------------------------------------------
# TetramerSampler class which will use the base sampler functions 
# and run the MCMC simulations for the tetramers sampling and 
# track sigma values
#----------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import json
import os
import random
import re
import pickle
import h5py
from typing import List, Dict, Tuple, Optional
from base_sampler import BaseMCSampler, Priors  # Import BaseMCSampler and related classes
from parameters import SystemParameters
from visualization import visualize_3d_configuration
from scipy.spatial.distance import cdist
from typing import Dict, List, Tuple
from scipy.stats import multivariate_normal
from pair_sampler import PairSampler
import cProfile
import pstats
import io
from sigma_provider import GMMSigmaProvider
#import numba as nb
#-----------------------------------------------------------------------
class TetramerSampler(BaseMCSampler):
    """
    Sampler for Tetramer-level interactions, inheriting from BaseMCSampler.
    Includes run_mc method.
    """
    def __init__(self, use_sigma_distribution=False, positions_ts=None, 
                sig_passed=None, sig_range_passed=None):
        super().__init__()
        self.params = SystemParameters()
        self.use_sigma_distribution = use_sigma_distribution
        
        # Initialize features for tetramer-specific sampling
        self.tetramer_trans_step = 0.02
        self.tetramer_rot_step = 0.02
        self.target_acceptance = 0.5
        self.tet_trans_acc_rate = self.target_acceptance
        self.sigma_prior_dist = {}
        self.sig_provider = GMMSigmaProvider('TetramerSampler')
        
        # Validate required pair types
        required_pairs = {"AA", "AB", "BC", "CC"}
        if not required_pairs.issubset(self.params.pair_distances.keys()):
            raise ValueError("Missing required pair types for tetramer sampling.")
        
        # Handle positions
        if positions_ts is None:
            self.positions_ts = self.get_positions()
            print("using the positions from the trajectory")
        else:
            self.positions_ts = positions_ts
            print("using the passed in positions")
        
        #visualize_3d_configuration(self.positions_ts, self.params.radii)
        
        # Initialize sigma values based on input parameters
        self._initialize_sigma_values(use_sigma_distribution, sig_passed, sig_range_passed)
        
        # Initialize PairSampler
        self.ps = PairSampler(
            use_def_sig_pos=False, 
            sig_passed=self.sigma, 
            sig_range_passed=self.sigma_range,
            pos_passed=self.positions_ts
        )

    def _initialize_sigma_values(self, use_sigma_distribution, sig_passed, sig_range_passed):
        """Initialize sigma values correctly, ensuring GMM data is actually used if available."""
        if use_sigma_distribution:
            print("using sigma distribution now")
            #self.gmm_params = self._load_gmm_parameters()
            #self.sigma = self.initialize_sigma_from_gmm()
            self.sigma = self.sig_provider.sample_sigma_values()
            self.sigma_range = None
        elif sig_passed is not None:
            self.sigma = sig_passed
            self.sigma_range = sig_range_passed
            print("using passed in sigma values")
        else:
            self.sigma, self.sigma_range = self.initialize_sigma()
            self.base_priors = Priors("jeffreys")
            print("Using default sigma initialization")

    def get_positions(self) -> Dict[str, np.ndarray]:
        """Simply load positions from the last frame of a trajectory file."""
        import pathlib
        
        traj_dir = pathlib.Path(os.getcwd()) / "output_analysis/pairsampler_results"
        
        try:
            # Get first trajectory file
            trajectory_files = list(traj_dir.glob("trajectory_chain_*.h5"))
            
            if not trajectory_files:
                raise FileNotFoundError(f"No trajectory files found in {traj_dir}")
                
            # Select random file
            filepath = random.choice(trajectory_files)
            chain_num = int(filepath.stem.split('_')[-1])
            print(f"Selected chain: {chain_num}")
                
            with h5py.File(filepath, 'r') as f:
                if 'trajectory' not in f:
                    raise KeyError("Invalid trajectory file format: missing 'trajectory' group")
                    
                traj_grp = f['trajectory']
                keys = sorted(traj_grp.keys())
                
                if not keys:
                    raise ValueError("Empty trajectory file")
                    
                # Get last frame
                last_key = keys[-1]
                print(f"Using last frame: {last_key}")
                
                # Read positions
                positions = {}
                pos_grp = traj_grp[last_key]['positions']
                
                for type_name in pos_grp:
                    positions[type_name] = pos_grp[type_name][:].copy()
                
            return positions
            
        except Exception as e:
            print(f"Error loading trajectory: {e}")
            print("Falling back to initialized positions")
            return self.initialize_positions()

    def get_tetramers(self, positions: Dict[str, np.ndarray], temp: float = 0.90) -> List[Tuple[int, ...]]:
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
#-----------------------------------------------------------------------  
    def run_mc(self, n_steps=50000, save_freq=1000, output_dir="output_analysis/tetramersampler_results/"):
        """
        Monte Carlo sampling with position, sigma, and tetramer moves.
        Simplified for better performance and readability.
        """
        # Setup output directory
        os.makedirs(output_dir, exist_ok=True)
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        with h5py.File(trajectory_file, 'w') as f:
            pass  # Create empty file
        
        # Initialize tracking variables
        best_positions = None
        best_score = float('inf')
        sigma_history = {key: np.zeros(n_steps // save_freq + 1) for key in self.sigma}
        accepts = {'position': 0, 'sigma': 0, 'tetramer': 0}
        attempts = {'position': 0, 'sigma': 0, 'tetramer': 0}
        
        # Initialize state
        current_tetramers = self.get_tetramers(self.positions_ts)
        prior_penalty = (self.sig_provider.calculate_negative_log_prior(self.sigma) if self.use_sigma_distribution 
                        else self.base_priors.neg_log_prior(self.sigma, self.sigma_range))
#        current_score, curr_ex, curr_pair, curr_tet = self.neg_log_posterior(
#            self.positions_ts, current_tetramers, prior_penalty, self.sigma)
        # And also for the initial score calculation at the beginning:
        current_score, curr_ex, curr_pair, curr_tet = self.neg_log_posterior(
            self.positions_ts, current_tetramers, prior_penalty, self.sigma, 
            debug=False  # Debug the initial state
        )
        
        # Store initial sigma values
        for key in self.sigma:
            sigma_history[key][0] = self.sigma[key]
        
        # Main MCMC loop
        move_types = ['position', 'sigma', 'tetramer']
        move_probs = [0.4, 0.1, 0.5]  # position, sigma, tetramer
        
        # Simple cooling schedule
        temp_start, temp_end = 6.0, 2.0
        temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
        
        print(f"Starting MCMC sampling for {n_steps} steps...")
        
        accepted_moves = 0
        total_moves = 0
        max_iterations = n_steps * 30  # Safety cap
        
        while accepted_moves < n_steps and total_moves < max_iterations:
            total_moves += 1
            
            # Temperature schedule
            temp = temp_start * (temp_decay ** accepted_moves)
            
            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1
            
            # Propose move
            if move_type == 'position':
                proposed_positions = self.propose_position_move(self.positions_ts, accepts['position'] / max(1, accepted_moves))
                proposed_sigma = self.sigma
                proposed_tetramers = self.get_tetramers(proposed_positions)
            elif move_type == 'sigma':
                proposed_positions = self.positions_ts
                proposed_sigma, _ = self.propose_sigma_move(self.sigma, accepts['sigma'] / max(1, accepted_moves))
                proposed_tetramers = current_tetramers  # Reuse
            else:  # tetramer move
                proposed_positions = self.propose_tetramer_move(self.positions_ts, accepts['tetramer'] / max(1, accepted_moves))
                proposed_sigma = self.sigma
                proposed_tetramers = self.get_tetramers(proposed_positions)
            
            # Calculate prior
            new_prior = (self.sig_provider.calculate_negative_log_prior(proposed_sigma) if self.use_sigma_distribution 
                        else self.base_priors.neg_log_prior(proposed_sigma, self.sigma_range))
            
            # Calculate new score
#            proposed_score, prop_ex, prop_pair, prop_tet = self.neg_log_posterior(
#                proposed_positions, proposed_tetramers, new_prior, proposed_sigma)
            proposed_score, prop_ex, prop_pair, prop_tet = self.neg_log_posterior(
                proposed_positions, proposed_tetramers, new_prior, proposed_sigma, 
                debug=False  # Enable debug when score is too high
            )
            
            # Metropolis criterion
            delta = proposed_score - current_score
            accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
            
            if accept:
                # Update state
                self.positions_ts = proposed_positions
                self.sigma = proposed_sigma
                current_tetramers = proposed_tetramers
                current_score = proposed_score
                curr_ex, curr_pair, curr_tet = prop_ex, prop_pair, prop_tet
                
                # Track best state
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_ts.items()}
                
                accepts[move_type] += 1
                accepted_moves += 1
                
                # Save state periodically
                if accepted_moves % save_freq == 0:
                    # Store sigma history
                    save_idx = accepted_moves // save_freq
                    if save_idx < len(sigma_history[list(sigma_history.keys())[0]]):
                        for key in self.sigma:
                            sigma_history[key][save_idx] = self.sigma[key]
                    
                    # Save to disk using existing method
                    self.save_state_to_disk(
                        accepted_moves, self.positions_ts, self.sigma, current_score,
                        prior_score=new_prior, pair_score=curr_pair,
                        exvol_score=curr_ex, tet_score=curr_tet,
                        traj_file=trajectory_file
                    )
                    
                    # Print progress
                    acceptance_rate = accepted_moves / total_moves
                    print(f"Step {accepted_moves}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={acceptance_rate:.2f}")
        
        # Save sigma history
        import pandas as pd
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        # Print final statistics
        print("\nSampling complete:")
        for move_type in move_types:
            rate = accepts[move_type] / max(1, attempts[move_type])
            print(f"- {move_type}: {rate:.2f} acceptance ({accepts[move_type]}/{attempts[move_type]})")
        
        return best_positions, trajectory_file
#-----------------------------------------------------------------------
    def propose_tetramer_move(self, positions: Dict[str, np.ndarray], acceptance_rate: float) -> Dict[str, np.ndarray]:
        """
        Tetramer move: translate and rotate a randomly selected tetramer using adaptive step sizes.
        Ensures symmetric proposals (Gaussian translations + isotropic rotations).
        """
        new_pos = {k: v.copy() for k, v in positions.items()}

        # Get tetramers
        tetramers = self.get_tetramers(positions)
        if not tetramers:
            return new_pos

        # Select a random tetramer
        tetramer = tetramers[np.random.randint(len(tetramers))]
        a_idx, b_idx, c_idx1, c_idx2 = tetramer
        particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]

        # Get coordinates of the selected tetramer
        coords = np.array([new_pos[part][idx] for part, idx in particles])
        centroid = np.mean(coords, axis=0)

        # Smaller base step sizes for large radii/distances
        base_trans_step = 0.02  # Example: reduce from 0.1
        base_rot_step   = 0.02  # Example: reduce from 0.1

        # Adaptive step factor with narrower clipping
        factor = np.clip(1.0 + 1.5 * (acceptance_rate - self.target_acceptance), 0.6, 1.6)

        # Symmetric translation (Gaussian around 0)
        trans_step = base_trans_step * factor
        displacement = np.random.normal(0.0, trans_step, 3)

        # Symmetric rotation: pick random axis, angle ~ N(0, rot_step)
        rot_step = base_rot_step * factor
        rotation_axis = np.random.randn(3)
        rotation_axis /= np.linalg.norm(rotation_axis) + 1e-10
        rotation_angle = np.random.normal(0.0, rot_step)

        # Build rotation matrix via quaternion
        half_angle = rotation_angle / 2.0
        qw = np.cos(half_angle)
        qx = rotation_axis[0] * np.sin(half_angle)
        qy = rotation_axis[1] * np.sin(half_angle)
        qz = rotation_axis[2] * np.sin(half_angle)
        rot_matrix = np.array([
            [1 - 2*(qy**2 + qz**2),   2*(qx*qy - qw*qz),     2*(qx*qz + qw*qy)],
            [2*(qx*qy + qw*qz),       1 - 2*(qx**2 + qz**2), 2*(qy*qz - qw*qx)],
            [2*(qx*qz - qw*qy),       2*(qy*qz + qw*qx),     1 - 2*(qx**2 + qy**2)]
        ])

        # Apply rotation + translation about centroid
        for i, (part, idx) in enumerate(particles):
            vec = coords[i] - centroid
            rotated = rot_matrix @ vec
            final_pos = centroid + rotated + displacement
            new_pos[part][idx] = np.mod(final_pos, self.params.box_size)

        return new_pos
#-----------------------------------------------------------------------
    def calculate_tetramer_scores_batch(self, positions, tetramers, sig, debug_logging=False):
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
#        
        return scores

#    def neg_log_posterior(
#        self,
#        positions: Dict[str, np.ndarray],
#        tetramers: List[Tuple[int, ...]],
#        prior_penalty_from_distribution: float = 0.0,
#        sig: Dict[str, float] = None,
#        exclusion_weight: float = 1.0,
#        pair_weight: float = 1.0,
#        tetramer_weight: float = 1.0,
#    ) -> Tuple[float, float, float, float]:
#        """
#        Efficiently calculate the total score for a tetramer system with vectorized operations.
#        
#        Returns:
#            Tuple of (total_score, exclusion_score, pair_score, tetramer_score)
#        """
#        # Use sigma if provided, otherwise use class sigma
#        sigma = sig if sig is not None else self.sigma
#        
#        # Quick return if no tetramers
#        if not tetramers:
#            # Calculate score with no excluded pairs
#            score, ex_score, pair_score, _ = self.ps.calculate_score(
#                positions, sigma, self.sigma_range, 
#                set(), self.use_sigma_distribution, 
#                prior_penalty_from_distribution
#            )
#            return score, ex_score, pair_score, 0.0
#        
#        # 1) Efficiently build set of tetramer pairs
#        tetramer_pairs = set()
#        tetramer_pairs_add = tetramer_pairs.add  # Local reference for faster calls
#        
#        # Build set of pairs with optimized batch processing
#        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
#            tetramer_pairs_add(('A', a_idx, 'B', b_idx))
#            tetramer_pairs_add(('B', b_idx, 'C', c1_idx))
#            tetramer_pairs_add(('B', b_idx, 'C', c2_idx))
#            tetramer_pairs_add(('C', c1_idx, 'C', c2_idx))
#
#        # 2) Calculate score excluding tetramer pairs
#        score, ex_score, pair_score, _ = self.ps.calculate_score(
#            positions, sigma, self.sigma_range,
#            tetramer_pairs, self.use_sigma_distribution,
#            prior_penalty_from_distribution
#        )
#
#        # 3) Calculate tetramer score using the vectorized batch method
#        scores_array = self.calculate_tetramer_scores_batch(positions, tetramers, sigma)
#        total_tet_score = scores_array.sum()
#        
#        # Apply weighting and return all score components
#        weighted_tet_score = tetramer_weight * total_tet_score
#        score += weighted_tet_score
#
#        return score, ex_score, pair_score, weighted_tet_score
    
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
        
        # Only create/write to log file when debug is True
        if debug:
            tetramer_pairs_log = "output_analysis/tetramer_pairs.txt"
            with open(tetramer_pairs_log, "a") as log_file:
                log_file.write(f"\n===== MCMC STEP - {len(tetramers)} TETRAMERS =====\n")
        
        # Quick return if no tetramers
        if not tetramers:
            # Calculate score with no excluded pairs
            result = self.ps.calculate_score(
                positions, sigma, self.sigma_range, 
                set(), self.use_sigma_distribution, 
                prior_penalty_from_distribution
            )
            
            if isinstance(result, tuple) and len(result) >= 3:
                score, ex_score, pair_score = result[:3]
            else:
                if debug:
                    print("Error: Unexpected return value from calculate_score")
                return 0.0, 0.0, 0.0, 0.0
                
            if debug:
                print("No tetramers - using all pairs for scoring")
                with open(tetramer_pairs_log, "a") as log_file:
                    log_file.write("No tetramers found - no intra-tetramer pairs to score\n")
                    
            return score, ex_score, pair_score, 0.0
        
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
            
            # Only log when debug is True
            if debug:
                # Log to file
                with open(tetramer_pairs_log, "a") as log_file:
                    log_file.write(f"Tetramer {t_idx}:\n")
                    log_file.write(f"  A({a_idx})-B({b_idx})\n")
                    log_file.write(f"  B({b_idx})-C({c1_idx})\n")
                    log_file.write(f"  B({b_idx})-C({c2_idx})\n")
                    log_file.write(f"  C({c1_idx})-C({c2_idx})\n")
                
                # Log to console
                print(f"Tetramer {t_idx}:")
                print(f"  A({a_idx})-B({b_idx})")
                print(f"  B({b_idx})-C({c1_idx})")
                print(f"  B({b_idx})-C({c2_idx})")
                print(f"  C({c1_idx})-C({c2_idx})")

        # Log statistics about tetramer pairs only when debug is True
        unique_pairs = len(tetramer_pairs)
        expected_pairs = len(tetramers) * 4
        
        if debug:
            with open(tetramer_pairs_log, "a") as log_file:
                log_file.write(f"Total unique tetramer pairs: {unique_pairs} (expected {expected_pairs})\n")
                if unique_pairs != expected_pairs:
                    log_file.write("WARNING: Some tetramers share components - potential source of problems!\n")
            
            print(f"Total unique tetramer pairs: {unique_pairs} (expected {expected_pairs})")
            if unique_pairs != expected_pairs:
                print("WARNING: Some tetramers share components - potential source of problems!")
        
        # 2) Calculate score excluding tetramer pairs
        result = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            tetramer_pairs, self.use_sigma_distribution,
            prior_penalty_from_distribution, debug_pairs = False
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
    #-----------------------------------------------------------------------
