#----------------------------------------------------------------------
# OctetSampler class which will use the base sampler functions 
# and run the MCMC simulations for the octamers sampling and 
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
from typing import List, Dict, Tuple, Optional
from base_sampler import BaseMCSampler, Priors  # Import BaseMCSampler and related classes
from parameters import SystemParameters
from tetramer_sampler import TetramerSampler
from visualization import visualize_3d_configuration
from scipy.spatial.distance import cdist
from typing import Dict, List, Tuple
from scipy.stats import multivariate_normal
from pair_sampler import PairSampler
import h5py  # For reading HDF5 files
#-----------------------------------------------------------------------
class OctetSampler(BaseMCSampler):
    """
    Sampler for Octamer-level interactions, inheriting from BaseMCSampler.
    Includes run_mc method.
    """
    def __init__(self, use_sigma_distribution=False, sig_passed=None, sig_range_passed=None):
        super().__init__()  # Call BaseMCSampler constructor
        self.use_sigma_distribution = use_sigma_distribution
        self.params = SystemParameters()  # Initialize system parameters
        self._initialize_step_sizes()

        # Basic sampler parameters
        self.octet_trans_step = 0.01
        self.octet_rot_step = 0.01
        self.octet_trans_acc_rate = 0.5
        self.target_acceptance = 0.5

        if use_sigma_distribution:
            # use positions from trajectory files
            self.positions_os = self.get_positions()
        else:
            # Initialize positions using uniform prior
            self.positions_os = self.initialize_positions()
        
        # Initialize sigma values using sigma_provider if needed
        self._initialize_sigma_values(use_sigma_distribution, sig_passed, sig_range_passed)
        
        # Initialize samplers for lower levels
        self.ts = TetramerSampler(
            use_sigma_distribution=use_sigma_distribution, 
            positions_ts=self.positions_os,
            sig_passed=self.sigma, 
            sig_range_passed=self.sigma_range
        )
        
        self.ps = PairSampler(
            use_def_sig_pos=False, 
            sig_passed=self.sigma, 
            sig_range_passed=self.sigma_range,
            pos_passed=self.positions_os
        )

    def _initialize_sigma_values(self, use_sigma_distribution, sig_passed, sig_range_passed):
        """Initialize sigma values correctly using sigma_provider if available."""
        if use_sigma_distribution:
            # Use sigma_provider to get values and handle GMM calculations
            from sigma_provider import GMMSigmaProvider
            self.sig_provider = GMMSigmaProvider(sampler_name="OctetSampler")
            self.sigma = self.sig_provider.sample_sigma_values()
            self.sigma_range = None
            print("OctetSampler initialized with GMM-based sigma values.")
        elif sig_passed is not None:
            # Use passed-in values
            self.sigma = sig_passed
            self.sigma_range = sig_range_passed
            print("OctetSampler initialized with passed sigma values.")
        else:
            # Use default initialization
            self.sigma, self.sigma_range = self.initialize_sigma()
            self.base_priors = Priors("jeffreys")
            print("OctetSampler initialized with default sigma values.")
    
    def get_positions(self) -> Dict[str, np.ndarray]:
        """Simply load positions from the last frame of a trajectory file."""
        import pathlib
        
        traj_dir = pathlib.Path(os.getcwd()) / "output_analysis/tetramersampler_results"
        
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
        
    
    def propose_octet_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Propose a move for a randomly selected octet by applying a small translation and rotation.
        This version does not apply periodic boundary conditions.
        """
        # Create a deep copy of positions to modify
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # Get tetramers and then form octets
        tetramers = self.ts.get_tetramers(new_pos) # Use new_pos in case get_tetramers modifies it
        octets = self.get_octets(new_pos, tetramers)
        
        if not octets:
            return new_pos # No octets to move, return original positions
        
        # Select a random octet
        octet_idx = np.random.randint(len(octets))
        tetramer1, tetramer2 = octets[octet_idx]
        
        # Collect all particle types and indices belonging to the selected octet
        octet_particles_info = []
        for tet in [tetramer1, tetramer2]:
            # tet is (a_idx, b_idx, c_idx1, c_idx2)
            octet_particles_info.extend([
                ('A', tet[0]), ('B', tet[1]), ('C', tet[2]), ('C', tet[3])
            ])
        
        # Group indices by particle type for efficient transformation
        grouped_indices: Dict[str, List[int]] = {}
        for particle_type, particle_idx in octet_particles_info:
            grouped_indices.setdefault(particle_type, []).append(particle_idx)
        
        # Collect current coordinates of all particles in the octet to calculate the centroid
        octet_coords_list = []
        for particle_type, indices_list in grouped_indices.items():
            octet_coords_list.append(new_pos[particle_type][indices_list])
        
        if not octet_coords_list: # Should not happen if octets were found
            return new_pos 
            
        all_octet_coords = np.vstack(octet_coords_list)
        centroid = np.mean(all_octet_coords, axis=0)
        
        # Define step sizes for translation and rotation
        # Using the adaptive step sizes from the class
        trans_step = self.octet_trans_step 
        rot_step = self.octet_rot_step
        
        # Generate a random translation vector
        displacement = np.random.normal(0, trans_step, 3)
        
        # Generate a random rotation axis and angle
        axis = self._random_unit_vector()
        angle = np.random.normal(0, rot_step)
        rotation_matrix = self._rotation_matrix(axis, angle)
        
        # Apply the transformation to each particle in the octet
        for particle_type, indices_list in grouped_indices.items():
            # Get current positions of particles of this type in the octet
            current_particle_positions = new_pos[particle_type][indices_list]
            
            # Translate to origin (centroid), rotate, then translate back and apply displacement
            transformed_positions = (current_particle_positions - centroid) @ rotation_matrix.T + centroid + displacement
            
            # Update positions in the new_pos dictionary
            new_pos[particle_type][indices_list] = transformed_positions
            
        return new_pos
    
    def _initialize_step_sizes(self):
        """Set much smaller initial step sizes for better acceptance rates"""
        self.octet_trans_step = 0.04  # Start with very small steps
        self.octet_rot_step = 0.02
        self.step_adaptation_factor = 0.95
        self.min_step = 0.001  # Even smaller minimum
        self.max_step = 0.2    # More conservative maximum
        self.octet_trans_acc_rate = 0.1  # Start with a more realistic value

    def adapt_step_sizes(self, acceptance_rate):
        """
        Adapt step sizes based on acceptance rate to target ~30% acceptance
        """
        target_rate = 0.3
        
        if acceptance_rate < target_rate:
            # Decrease step sizes if acceptance is too low
            self.octet_trans_step *= self.step_adaptation_factor
            self.octet_rot_step *= self.step_adaptation_factor
        else:
            # Increase step sizes if acceptance is too high
            self.octet_trans_step /= self.step_adaptation_factor
            self.octet_rot_step /= self.step_adaptation_factor
        
        # Enforce bounds
        self.octet_trans_step = max(self.min_step, min(self.max_step, self.octet_trans_step))
        self.octet_rot_step = max(self.min_step, min(self.max_step, self.octet_rot_step))

    def get_octets(self, positions: Dict[str, np.ndarray], tetramers=None, temp=0.9) -> List[Tuple]:
        """
        Group tetramers into octets (pairs of tetramers) with temperature-based selection.
        Vectorized implementation for better performance.
        """
        if tetramers is None:
            tetramers = self.ts.get_tetramers(positions)
        
        if len(tetramers) < 2:
            return []
        
        # Calculate tetramer centers using vectorized operations
        centers = np.zeros((len(tetramers), 3))
        for i, (a_idx, b_idx, c_idx1, c_idx2) in enumerate(tetramers):
            coords = np.vstack([
                positions['A'][a_idx],
                positions['B'][b_idx],
                positions['C'][c_idx1],
                positions['C'][c_idx2]
            ])
            centers[i] = np.mean(coords, axis=0)
        
        # Form octets by pairing tetramers
        octets = []
        available = list(range(len(tetramers)))
        box_size = self.params.box_size
        
        while len(available) >= 2:
            # Pick first tetramer randomly
            idx1 = np.random.choice(available)
            available.remove(idx1)
            
            # Calculate all distances at once using vectorized operations
            indices = np.array(available)
            center1 = centers[idx1]
            deltas = centers[indices] - center1
            
            # Periodic boundary correction (vectorized)
            mask = np.abs(deltas) > box_size/2
            deltas[mask] -= np.sign(deltas[mask]) * box_size
            
            # Calculate distances (vectorized)
            distances = np.linalg.norm(deltas, axis=1)
            
            # Calculate selection probabilities
            probs = np.exp(-distances / temp)
            probs_sum = probs.sum()
            
            # Select second tetramer
            if probs_sum > 1e-10:
                probs = probs / probs_sum
                idx2_rel = np.random.choice(len(available), p=probs)
                idx2 = available[idx2_rel]
            else:
                idx2 = np.random.choice(available)
            
            # Add the pair to octets and remove from available
            octets.append((tetramers[idx1], tetramers[idx2]))
            available.remove(idx2)
        
        return octets
    
    def _random_unit_vector(self):
        """Generate a random unit vector."""
        vec = np.random.randn(3)
        vec /= np.linalg.norm(vec) + 1e-10  # Add small epsilon to avoid division by zero
        return vec

    def _rotation_matrix(self, axis, theta):
        """
        Create a 3D rotation matrix using Rodrigues' formula.
        More numerically stable implementation.
        """
        axis = np.asarray(axis)
        axis = axis / np.linalg.norm(axis)
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        
        return np.array([
            [a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
            [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
            [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]
        ])
        
    def run_mc(self, n_steps=50000, save_freq=1000, output_dir="output_analysis/octetsampler_results/"):
        """
        Monte Carlo sampling with position, sigma, tetramer, and octet moves.
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
        accepts = {'position': 0, 'sigma': 0, 'tetramer': 0, 'octet': 0}
        attempts = {'position': 0, 'sigma': 0, 'tetramer': 0, 'octet': 0}
        
        # Initialize state
        tetramers = self.ts.get_tetramers(self.positions_os)
        octets = self.get_octets(self.positions_os, tetramers)
        
        prior_penalty = (self.sig_provider.calculate_negative_log_prior(self.sigma) if self.use_sigma_distribution 
                        else self.base_priors.neg_log_prior(self.sigma, self.sigma_range))
        
        current_score, curr_ex, curr_pair, curr_oct = self.neg_log_posterior(
            self.positions_os, tetramers, octets, prior_penalty, self.sigma)
        
        # Store initial sigma values
        for key in self.sigma:
            sigma_history[key][0] = self.sigma[key]
        
        # Main MCMC loop parameters
        move_types = ['position', 'sigma', 'tetramer', 'octet']
        move_probs = [0.2, 0.1, 0.3, 0.4]  # position, sigma, tetramer, octet
        
        # Simple cooling schedule
        temp_start, temp_end = 40.0, 1.0
        temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
        
        print(f"Starting MCMC sampling for {n_steps} steps...")
        
        accepted_moves = 0
        total_moves = 0
        max_iterations = n_steps * 20  # Safety cap
        
        while accepted_moves < n_steps and total_moves < max_iterations:
            total_moves += 1
            
            # Temperature schedule
            temp = temp_start * (temp_decay ** accepted_moves)
            
            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1
            
            # Propose move
            if move_type == 'position':
                # also pass the acceptance rate for position moves
                proposed_positions = self.propose_position_move(self.positions_os, 
                                                                accepts['position'] / max(1, attempts['position']))
                proposed_sigma = self.sigma
            elif move_type == 'sigma':
                proposed_positions = self.positions_os
                proposed_sigma, _ = self.propose_sigma_move(self.sigma, 
                                                            accepts['sigma'] / max(1, attempts['sigma']))
            elif move_type == 'tetramer':
                proposed_positions = self.ts.propose_tetramer_move(self.positions_os, 
                                                                   accepts['tetramer'] / max(1, attempts['tetramer']))
                proposed_sigma = self.sigma
            else:  # octet move
                proposed_positions = self.propose_octet_move(self.positions_os)
                proposed_sigma = self.sigma
            
            # Get tetramers and octets for proposed configuration
            proposed_tetramers = self.ts.get_tetramers(proposed_positions)
            proposed_octets = self.get_octets(proposed_positions, proposed_tetramers)
            
            # Calculate prior
            new_prior = (self.sig_provider.calculate_negative_log_prior(proposed_sigma) if self.use_sigma_distribution 
                        else self.base_priors.neg_log_prior(proposed_sigma, self.sigma_range))
            
            # Calculate new score
            proposed_score, prop_ex, prop_pair, prop_oct = self.neg_log_posterior(
                proposed_positions, proposed_tetramers, proposed_octets, new_prior, proposed_sigma)
            
            # Check for non-finite scores
            if not np.isfinite(proposed_score):
                print(f"WARNING: Non-finite score {proposed_score} detected!")
                print(f"Components: ex={prop_ex}, pair={prop_pair}, oct={prop_oct}, prior={new_prior}")
                # Skip this iteration or handle accordingly
                continue
            
            # Metropolis criterion
            delta = proposed_score - current_score
            accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
            
            debug = True  # Uncomment for detailed debug output
            
            if accept:
                # Update state
                self.positions_os = proposed_positions
                self.sigma = proposed_sigma
                current_score = proposed_score
                curr_ex, curr_pair, curr_oct = prop_ex, prop_pair, prop_oct
                
                if debug:
                    # open a file and write the score components to it 
                    # overwrite if it exists
                    # (this is for debugging purposes)
                    with open(os.path.join(output_dir, "debug_scores.txt"), 'w') as f:
                        f.write(f"Step {accepted_moves}: Score={current_score:.2f}, "
                                f"Exclusion={curr_ex:.2f}, Pair={curr_pair:.2f}, Octet={curr_oct:.2f}, Prior={new_prior:.2f}\n")
                
                # Track best state
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_os.items()}
                
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
                        accepted_moves, self.positions_os, self.sigma, current_score,
                        prior_score=new_prior, pair_score=curr_pair,
                        exvol_score=curr_ex, oct_score=curr_oct,
                        traj_file=trajectory_file
                    )
                    
                    # Print progress
                    acceptance_rate = accepted_moves / total_moves
                    print(f"Step {accepted_moves}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={acceptance_rate:.2f}")

                    # Adapt step sizes periodically
                    if accepted_moves % 500 == 0:
                        for move_type in move_types:
                            rate = accepts[move_type] / max(1, attempts[move_type])
                            if move_type == 'octet':
                                self.octet_trans_acc_rate = rate  # Update the acc_rate for propose_octet_move
                                self.adapt_step_sizes(rate)
                                print(f"Adapting octet step sizes: trans={self.octet_trans_step:.4f}, rot={self.octet_rot_step:.4f}")

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

    def calculate_octet_scores_batch(self, positions, octets, sig):
        """
        Efficiently calculate octet scores with minimal memory allocations and vectorized operations.
        
        Args:
            positions: Dictionary of particle positions
            octets: List of pairs of tetramers [(tet1, tet2), ...]
            sig: Dictionary of sigma values
        
        Returns:
            np.ndarray: Array of octet scores (one per octet)
        """
        if not octets:
            return np.array([])  # Empty array with default float64 precision
        
        n_octets = len(octets)
        
        # Extract all pairs in one pass to minimize object creation
        a1, b1, c1a, c1b = [], [], [], []
        a2, b2, c2a, c2b = [], [], [], []
        
        for tet1, tet2 in octets:
            a1.append(tet1[0])
            b1.append(tet1[1])
            c1a.append(tet1[2])
            c1b.append(tet1[3])
            
            a2.append(tet2[0])
            b2.append(tet2[1])
            c2a.append(tet2[2])
            c2b.append(tet2[3])
        
        # Convert to arrays (int is more efficient than int32/int64 here)
        a1, b1, c1a, c1b = map(np.array, [a1, b1, c1a, c1b])
        a2, b2, c2a, c2b = map(np.array, [a2, b2, c2a, c2b])
        
        # Get positions (single lookup per type)
        posA1 = positions['A'][a1]
        posB1 = positions['B'][b1]
        posC1a = positions['C'][c1a]
        posC1b = positions['C'][c1b]
        
        posA2 = positions['A'][a2]
        posB2 = positions['B'][b2]
        posC2a = positions['C'][c2a]
        posC2b = positions['C'][c2b]
        
        # Preallocate the results array
        octet_scores = np.zeros(n_octets)
        
        # Compute all distances in one vectorized pass per pair type
        # Correct Gaussian NLL with proper sigma² term
        def score_pairs(pos1, pos2, target_dist, sigma_val):
            delta = pos1 - pos2
            dist = np.sqrt(np.sum(delta**2, axis=1))
            return ((dist - target_dist)**2)/(2.0*sigma_val**2) + 0.5*np.log(2.0*np.pi*sigma_val**2)
        
        # 1. A-A interactions
        octet_scores += score_pairs(posA1, posA2, self.params.pair_distances['AA'], sig['AA'])
        
        # 2. A-B interactions (both directions)
        octet_scores += score_pairs(posA1, posB2, self.params.pair_distances['AB'], sig['AB'])
        octet_scores += score_pairs(posB1, posA2, self.params.pair_distances['AB'], sig['AB'])
        
        # 3. B-C interactions (four combinations, vectorized)
        bc_target = self.params.pair_distances['BC']
        bc_sigma = sig['BC']
        
        # B1-C2a, B1-C2b
        octet_scores += score_pairs(posB1, posC2a, bc_target, bc_sigma)
        octet_scores += score_pairs(posB1, posC2b, bc_target, bc_sigma)
        
        # C1a-B2, C1b-B2
        octet_scores += score_pairs(posC1a, posB2, bc_target, bc_sigma)
        octet_scores += score_pairs(posC1b, posB2, bc_target, bc_sigma)
        
        # 4. C-C interactions (all four combinations)
        cc_target = self.params.pair_distances['CC']
        cc_sigma = sig['CC']
        
        octet_scores += score_pairs(posC1a, posC2a, cc_target, cc_sigma)
        octet_scores += score_pairs(posC1a, posC2b, cc_target, cc_sigma)
        octet_scores += score_pairs(posC1b, posC2a, cc_target, cc_sigma)
        octet_scores += score_pairs(posC1b, posC2b, cc_target, cc_sigma)
        
        return octet_scores

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        tetramers: List[Tuple[int, ...]],
        octets: List[Tuple[Tuple[int, ...], Tuple[int, ...]]],
        prior_penalty: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        octet_weight: float = 1.0,
    ) -> Tuple[float, float, float, float]:
        """
        Calculate the total negative log-posterior score for an octet system.
        Optimized implementation with reduced redundancy and memory usage.
        """
        # Use the provided sigma if given, else use the sampler's
        sigma = sig if sig is not None else self.sigma

        # Fast path for no octets
        if not octets:
            score, ex_score, pair_score, _ = self.ps.calculate_score(
                positions, sigma, self.sigma_range,
                excluded_pairs=set(),
                use_sigma_distribution=self.use_sigma_distribution,
                prior_penalty_from_distribution=prior_penalty
            )
            return score, ex_score, pair_score, 0.0

        # Use a dict-based structure first (faster for insertions)
        # Format: {(type1, idx1, type2, idx2): True}
        excluded_dict = {}
        
        # Helper to add a pair to exclusions only once
        def add_pair(t1, i1, t2, i2):
            excluded_dict[(t1, i1, t2, i2)] = True
        
        # 1) Process tetramers (intra-tetramer pairs)
        for a, b, c1, c2 in tetramers:
            # A-B connections (both directions)
            add_pair('A', a, 'B', b)
            add_pair('B', b, 'A', a)
            
            # B-C connections 
            for c in [c1, c2]:
                add_pair('B', b, 'C', c)
                add_pair('C', c, 'B', b)
            
            # C1-C2 connection
            add_pair('C', c1, 'C', c2)
            add_pair('C', c2, 'C', c1)
        
        # 2) Process octets (inter-tetramer pairs)
        for (a1, b1, c1a, c1b), (a2, b2, c2a, c2b) in octets:
            # A-A connection
            add_pair('A', a1, 'A', a2)
            add_pair('A', a2, 'A', a1)
            
            # A-B cross-connections
            add_pair('A', a1, 'B', b2)
            add_pair('B', b2, 'A', a1)
            add_pair('A', a2, 'B', b1)
            add_pair('B', b1, 'A', a2)
            
            # B-C cross-connections
            for c2 in [c2a, c2b]:
                add_pair('B', b1, 'C', c2)
                add_pair('C', c2, 'B', b1)
            
            for c1 in [c1a, c1b]:
                add_pair('B', b2, 'C', c1)
                add_pair('C', c1, 'B', b2)
            
            # C-C cross-connections (4 pairs)
            for c1 in [c1a, c1b]:
                for c2 in [c2a, c2b]:
                    add_pair('C', c1, 'C', c2)
                    add_pair('C', c2, 'C', c1)
        
        # Convert to set only once at the end
        excluded_pairs = set(excluded_dict.keys())
        
        # 3) Calculate scores with minimum memory overhead
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            excluded_pairs=excluded_pairs,
            use_sigma_distribution=self.use_sigma_distribution,
            prior_penalty_from_distribution=prior_penalty
        )
        
        # 4) Calculate octet score efficiently
        octet_scores = self.calculate_octet_scores_batch(positions, octets, sigma)
        weighted_octet_score = octet_weight * np.sum(octet_scores)
        total_score = score + weighted_octet_score
        
        return total_score, ex_score, pair_score, weighted_octet_score