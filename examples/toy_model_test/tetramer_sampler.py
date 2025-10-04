#----------------------------------------------------------------------
# TetramerSampler class which will use the base sampler functions 
# and run the MCMC simulations for the tetramers sampling and 
# track sigma values
#----------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import h5py
import random
from typing import List, Dict, Tuple
from base_sampler import BaseMCSampler
from parameters import SystemParameters
from scipy.spatial.distance import cdist
from pair_sampler import PairSampler
from new_sigma_provider import GMMSigmaProvider

#-----------------------------------------------------------------------
class TetramerSampler(BaseMCSampler):
    """
    Sampler for Tetramer-level interactions, inheriting from BaseMCSampler.
    """
    def __init__(self, sampler_sequence: List[str], sequence_idx: int,
                 base_output_dir: str = "output_analysis", positions_ts=None,
                 specific_chain: int = None, sigma_ranges: Dict[str, Tuple[float, float]] = None,
                 prior_type: str = "gamma"):
        """
        Initialize TetramerSampler with sampler sequence information.
        
        Args:
            sampler_sequence: List of sampler names in execution order
            sequence_idx: Current position in the sampler sequence
            base_output_dir: Base directory for analysis results
            positions_ts: Optional positions to use instead of loading
            specific_chain: Specific chain to load from (if not specified, random)
            sigma_ranges: Ranges for sigma values
            prior_type: Prior type for first sampler ("uniform" or "jeffreys")
        """
        super().__init__()
        self.params = SystemParameters()
        self.sampler_sequence = sampler_sequence
        self.sequence_idx = sequence_idx
        self.specific_chain = specific_chain
        
        # Initialize features for tetramer-specific sampling
        self.target_acceptance = 0.5
        
        # Validate required pair types
        required_pairs = {"AA", "AB", "BC"}
        if not required_pairs.issubset(self.params.pair_distances.keys()):
            raise ValueError("Missing required pair types for tetramer sampling.")
        
        # Handle positions
        if positions_ts is None:
            self.positions_ts = self.get_positions()
            print("Using positions from previous trajectory")
        else:
            self.positions_ts = positions_ts
            print("Using passed in positions")
        
        # Initialize sigma provider
        self.sigma_provider = GMMSigmaProvider(
            sampler_sequence=sampler_sequence,
            sequence_idx=sequence_idx,
            base_output_dir=base_output_dir,
            specific_chain=specific_chain,
            sigma_ranges=sigma_ranges,
            prior_type=prior_type
        )
        
        # Initialize sigma values and ranges from provider
        self.sigma = self.sigma_provider.sample_sigma_values()
        self.sigma_range = self.sigma_provider.sigma_ranges
        
        # Initialize PairSampler for scoring
        self.ps = PairSampler(
            sampler_sequence=sampler_sequence,
            sequence_idx=sequence_idx,
            base_output_dir=base_output_dir,
            specific_chain=specific_chain,
            sigma_ranges=sigma_ranges,
            prior_type=prior_type,
            pos_passed=self.positions_ts
        )
        
    def get_positions(self) -> Dict[str, np.ndarray]:
        """
        Load positions from the previous sampler in the sequence or initialize if first.
        This method uses the class attributes and calls the standalone function.
        """
        return self.get_positions_from_previous_sampler(
            sampler_sequence=self.sampler_sequence,
            current_sequence_idx=self.sequence_idx,
            specific_chain=getattr(self, 'specific_chain', None),
            base_output_dir="output_analysis",
            params=self.params
        )

    # Alternative: Hungarian algorithm approach for optimal matching
    def get_tetramers(self, positions: Dict[str, np.ndarray]) -> List[Tuple[int, ...]]:
        """Use Hungarian algorithm for optimal A-B matching, then greedy C selection."""
        try:
            from scipy.optimize import linear_sum_assignment
            
            # Validate input
            if not all(k in positions and len(positions[k]) > 0 for k in ['A', 'B', 'C']) or len(positions['C']) < 2:
                return []
                
            a_pos, b_pos, c_pos = positions['A'], positions['B'], positions['C']
            
            # Get target distances
            ab_target = self.params.pair_distances['AB']
            bc_target = self.params.pair_distances['BC']
            
            # Calculate AB cost matrix (deviation from target distance)
            dist_AB = cdist(a_pos, b_pos)
            cost_matrix = np.abs(dist_AB - ab_target)
            
            # Solve optimal assignment problem
            a_indices, b_indices = linear_sum_assignment(cost_matrix)
            
            # Pre-calculate BC distances
            dist_BC = cdist(b_pos, c_pos)
            
            # Now assign C particles greedily based on B assignments
            c_used = set()
            tetramers = []
            
            # Sort A-B pairs by their cost (best matches first)
            pair_costs = cost_matrix[a_indices, b_indices]
            sorted_pairs = np.argsort(pair_costs)
            
            for pair_idx in sorted_pairs:
                a_idx = a_indices[pair_idx]
                b_idx = b_indices[pair_idx]
                
                # Find available C particles
                available_c = [i for i in range(len(c_pos)) if i not in c_used]
                if len(available_c) < 2:
                    break
                
                # Get best C pair for this B
                bc_dists = dist_BC[b_idx, available_c]
                c_scores = np.abs(bc_dists - bc_target)
                
                best_c_local = np.argsort(c_scores)[:2]
                best_c_indices = [available_c[i] for i in best_c_local]
                
                # Form tetramer
                tetramers.append((a_idx, b_idx, best_c_indices[0], best_c_indices[1]))
                c_used.update(best_c_indices)
                
                if len(tetramers) >= min(len(a_pos), len(b_pos), len(c_pos) // 2):
                    break
            
            return tetramers
            
        except Exception as e:
            print(f"Error in Hungarian tetramer generation: {e}")
            return []
    #=======================================================================
    def run_mc(self, n_steps=50000, save_freq=1000, output_dir="output_analysis/tetramersampler_results/", debug=False):
        """Monte Carlo sampling with position, sigma, and tetramer moves."""
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
        accepted_moves = 0
        
        # Initialize state
        current_tetramers = self.get_tetramers(self.positions_ts)
        prior_penalty = self.sigma_provider.calculate_negative_log_prior(self.sigma)
        current_score, curr_ex, curr_pair, curr_tet = self.neg_log_posterior(
            self.positions_ts, current_tetramers, prior_penalty, self.sigma, debug=False)
        
        # Store initial sigma values
        for key in self.sigma:
            sigma_history[key][0] = self.sigma[key]
        
        # Main MCMC loop
        move_types = ['position', 'sigma', 'tetramer']
        move_probs = [0.4, 0.1, 0.5]
        
        # Cooling schedule
        temp_start, temp_end = 5.0, 0.10
        temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
        
        print(f"Starting MCMC sampling for {n_steps} total steps...")
        
        for step in range(1, n_steps + 1):
            temp = temp_start * (temp_decay ** step)
            
            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1
            
            # Propose move
            pair_type = None
            if move_type == 'position':
                proposed_positions = self.propose_position_move(self.positions_ts, accepts['position'] / max(1, step))
                proposed_sigma = self.sigma
                proposed_tetramers = self.get_tetramers(proposed_positions)
            elif move_type == 'sigma':
                proposed_positions = self.positions_ts
                proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)
                proposed_tetramers = current_tetramers
            else:  # tetramer move
                proposed_positions = self.propose_tetramer_move(self.positions_ts, accepts['tetramer'] / max(1, step))
                proposed_sigma = self.sigma
                proposed_tetramers = self.get_tetramers(proposed_positions)
            
            # Calculate prior and score
            new_prior = self.sigma_provider.calculate_negative_log_prior(proposed_sigma)
            
            proposed_score, prop_ex, prop_pair, prop_tet = self.neg_log_posterior(
                proposed_positions, proposed_tetramers, new_prior, proposed_sigma, debug=False)
            
            # Metropolis criterion with Jacobian correction for sigma moves
            delta = proposed_score - current_score
            
            # Accept/reject
            accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
            
            if accept:
                self.positions_ts = proposed_positions
                self.sigma = proposed_sigma
                current_tetramers = proposed_tetramers
                current_score = proposed_score
                curr_ex, curr_pair, curr_tet = prop_ex, prop_pair, prop_tet
                
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_ts.items()}
                
                accepts[move_type] += 1
                accepted_moves += 1
            
            # Save state
            if step % save_freq == 0:
                save_idx = step // save_freq
                if save_idx < len(sigma_history[list(sigma_history.keys())[0]]):
                    for key in self.sigma:
                        sigma_history[key][save_idx] = self.sigma[key]
                
                self.save_state_to_disk(
                    step, self.positions_ts, self.sigma, current_score,
                    prior_score=new_prior, pair_score=curr_pair,
                    exvol_score=curr_ex, tet_score=curr_tet,
                    traj_file=trajectory_file
                )
                
                acceptance_rate = accepted_moves / step
                print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={acceptance_rate:.2f}")
        
        # Save results
        import pandas as pd
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        # Print final statistics
        print("\nSampling complete:")
        for move_type in move_types:
            rate = accepts[move_type] / max(1, attempts[move_type])
            print(f"- {move_type}: {rate:.2f} acceptance ({accepts[move_type]}/{attempts[move_type]})")
        print(f"- Best score: {best_score:.2f}")
        
        return best_positions, trajectory_file

    def propose_tetramer_move(self, positions: Dict[str, np.ndarray], acceptance_rate: float) -> Dict[str, np.ndarray]:
        """
        Optimized tetramer move proposal with decoupled translation/rotation.
        60% probability for translation, 30% for rotation, 10% for mixed moves.
        Clips coordinates to stay within box boundaries minus particle radius.
        """
        new_pos = {k: v.copy() for k, v in positions.items()}
        tetramers = self.get_tetramers(positions)
        
        if not tetramers:
            return new_pos
        
        # Get box boundaries with particle radius buffer
        box_size = getattr(self.params, 'box_size', 800.0)
        particle_radius = getattr(self.params, 'particle_radius', 5.0)  # Default radius
        half_box = box_size / 2.0
        max_coord = half_box - 2 * particle_radius
        min_coord = -max_coord
        
        # Select random tetramer
        tetramer = tetramers[np.random.randint(len(tetramers))]
        a_idx, b_idx, c_idx1, c_idx2 = tetramer
        
        # Pre-extract tetramer particle information efficiently
        particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]
        
        # Get coordinates in a vectorized way
        coords = np.array([new_pos[part][idx] for part, idx in particles])
        centroid = np.mean(coords, axis=0)
        
        # Calculate tetramer size for adaptive scaling
        distances_from_center = np.linalg.norm(coords - centroid, axis=1)
        tetramer_radius = np.max(distances_from_center)
        size_factor = max(0.5, min(2.0, tetramer_radius))
        
        # Fixed step sizes
        trans_step = 0.25
        rot_step = 0.2
        
        # Choose move type: 60% translation, 30% rotation, 10% mixed
        rand_val = np.random.random()
        
        if rand_val < 0.6:
            # --- TRANSLATION MOVE ---
            # Generate displacement with size-adaptive scaling
            displacement = np.random.normal(0.0, trans_step * size_factor, 3)
            
            # Apply translation to all particles in the tetramer
            for i, (part, idx) in enumerate(particles):
                final_pos = coords[i] + displacement
                # Clip coordinates to stay within bounds
                final_pos = np.clip(final_pos, min_coord, max_coord)
                new_pos[part][idx] = final_pos
        
        elif rand_val < 0.9:  # 0.6 to 0.9 = 30% probability
            # --- ROTATION MOVE ---
            # Generate random rotation using Marsaglia method for unit vector
            while True:
                x1, x2 = np.random.uniform(-1, 1, 2)
                if x1*x1 + x2*x2 < 1:
                    break
            
            sqrt_term = np.sqrt(1 - x1*x1 - x2*x2)
            rotation_axis = np.array([2*x1*sqrt_term, 2*x2*sqrt_term, 1 - 2*(x1*x1 + x2*x2)])
            
            # Generate rotation angle
            rotation_angle = np.random.normal(0.0, rot_step)
            
            # Build rotation matrix using quaternion (numerically stable)
            half_angle = rotation_angle / 2.0
            qw = np.cos(half_angle)
            qx = rotation_axis[0] * np.sin(half_angle)
            qy = rotation_axis[1] * np.sin(half_angle)
            qz = rotation_axis[2] * np.sin(half_angle)
            
            rot_matrix = np.array([
                [1 - 2*(qy**2 + qz**2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                [2*(qx*qy + qw*qz), 1 - 2*(qx**2 + qz**2), 2*(qy*qz - qw*qx)],
                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx**2 + qy**2)]
            ])
            
            # Apply rotation around centroid to all particles
            for i, (part, idx) in enumerate(particles):
                vec = coords[i] - centroid
                rotated = rot_matrix @ vec
                final_pos = centroid + rotated
                # Clip coordinates to stay within bounds
                final_pos = np.clip(final_pos, min_coord, max_coord)
                new_pos[part][idx] = final_pos
        
        else:  # 0.9 to 1.0 = 10% probability
            # --- MIXED MOVE (Translation + Rotation) ---
            # Use smaller step sizes for combined moves to maintain reasonable acceptance
            trans_scale = 0.7
            rot_scale = 0.7
            
            # Generate displacement with reduced step size
            displacement = np.random.normal(0.0, trans_step * size_factor * trans_scale, 3)
            
            # Generate rotation with reduced step size
            while True:
                x1, x2 = np.random.uniform(-1, 1, 2)
                if x1*x1 + x2*x2 < 1:
                    break
            
            sqrt_term = np.sqrt(1 - x1*x1 - x2*x2)
            rotation_axis = np.array([2*x1*sqrt_term, 2*x2*sqrt_term, 1 - 2*(x1*x1 + x2*x2)])
            
            rotation_angle = np.random.normal(0.0, rot_step * rot_scale)
            
            # Build rotation matrix
            half_angle = rotation_angle / 2.0
            qw = np.cos(half_angle)
            qx = rotation_axis[0] * np.sin(half_angle)
            qy = rotation_axis[1] * np.sin(half_angle)
            qz = rotation_axis[2] * np.sin(half_angle)
            
            rot_matrix = np.array([
                [1 - 2*(qy**2 + qz**2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                [2*(qx*qy + qw*qz), 1 - 2*(qx**2 + qz**2), 2*(qy*qz - qw*qx)],
                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx**2 + qy**2)]
            ])
            
            # Apply combined transformation: rotate around centroid, then translate
            for i, (part, idx) in enumerate(particles):
                vec = coords[i] - centroid
                rotated = rot_matrix @ vec
                final_pos = centroid + rotated + displacement
                # Clip coordinates to stay within bounds
                final_pos = np.clip(final_pos, min_coord, max_coord)
                new_pos[part][idx] = final_pos
        
        return new_pos
#------------------------------------------------------------------------------
    def calculate_tetramer_scores_batch(self, positions, tetramers, sig, debug_logging=False):
        """Calculate scores for all tetramers."""
        if not tetramers:
            return np.array([], dtype=np.float32)
        
        # Extract indices
        a_indices = np.array([t[0] for t in tetramers], dtype=np.int32)
        b_indices = np.array([t[1] for t in tetramers], dtype=np.int32)
        c1_indices = np.array([t[2] for t in tetramers], dtype=np.int32)
        c2_indices = np.array([t[3] for t in tetramers], dtype=np.int32)
        
        # Get positions
        pos_a = positions['A'][a_indices]
        pos_b = positions['B'][b_indices]
        pos_c1 = positions['C'][c1_indices]
        pos_c2 = positions['C'][c2_indices]
        
        # Calculate distances
        ab_dists = np.sqrt(np.sum((pos_a - pos_b)**2, axis=1))
        bc1_dists = np.sqrt(np.sum((pos_b - pos_c1)**2, axis=1))
        bc2_dists = np.sqrt(np.sum((pos_b - pos_c2)**2, axis=1))
        
        # Target distances
        ab_target = self.params.pair_distances['AB']
        bc_target = self.params.pair_distances['BC']
        
        # Calculate scores
        ab_scores = ((ab_dists - ab_target)**2)/(2*sig['AB']**2) + 0.5 * np.log(2 * np.pi * sig['AB']**2)
        bc1_scores = ((bc1_dists - bc_target)**2)/(2*sig['BC']**2) + 0.5 * np.log(2 * np.pi * sig['BC']**2)
        bc2_scores = ((bc2_dists - bc_target)**2)/(2*sig['BC']**2) + 0.5 * np.log(2 * np.pi * sig['BC']**2)
        
        scores = ab_scores + bc1_scores + bc2_scores
        return scores

    def neg_log_posterior(self, positions: Dict[str, np.ndarray], tetramers: List[Tuple[int, ...]],
                         prior_penalty_from_distribution: float = 0.0, sig: Dict[str, float] = None,
                         exclusion_weight: float = 1.0, pair_weight: float = 1.0,
                         tetramer_weight: float = 1.0, debug: bool = False) -> Tuple[float, float, float, float]:
        """Calculate the negative log posterior."""
        sigma = sig if sig is not None else self.sigma
        
        if not tetramers:
            result = self.ps.calculate_score(
                positions, sigma, self.sigma_range,
                set(), prior_penalty_from_distribution
            )
            
            if isinstance(result, tuple) and len(result) >= 3:
                score, ex_score, pair_score = result[:3]
                return score, ex_score, pair_score, 0.0
            else:
                return 0.0, 0.0, 0.0, 0.0
        
        # Calculate score excluding tetramer pairs (empty set for now)
        tetramer_pairs = set()
        result = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            tetramer_pairs, prior_penalty_from_distribution
        )
        
        if isinstance(result, tuple) and len(result) >= 3:
            score, ex_score, pair_score = result[:3]
        else:
            return 0.0, 0.0, 0.0, 0.0

        # Calculate tetramer score
        scores_array = self.calculate_tetramer_scores_batch(positions, tetramers, sigma, debug_logging=debug)
        total_tet_score = scores_array.sum()
        weighted_tet_score = tetramer_weight * total_tet_score
        score += weighted_tet_score

        if debug:
            print(f"Exclusion: {ex_score:.2f}, Pair: {pair_score:.2f}, Tetramer: {weighted_tet_score:.2f}, Total: {score:.2f}")

        return score, ex_score, pair_score, weighted_tet_score
