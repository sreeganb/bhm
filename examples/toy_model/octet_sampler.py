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
    Sampler for Tetramer-level interactions, inheriting from BaseMCSampler.
    Includes run_mc method.
    """
    def __init__(self, use_sigma_distribution=False, sig_passed=None, sig_range_passed=None):
        super().__init__()  # Call BaseMCSampler constructor
        self.use_sigma_distribution = use_sigma_distribution
        self.params = SystemParameters()  # Initialize system parameters

        # Basic sampler parameters
        self.octet_trans_step = 0.2
        self.octet_rot_step = 0.1
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
            print("PairSampler initialized with GMM-based sigma values.")
        elif sig_passed is not None:
            # Use passed-in values
            self.sigma = sig_passed
            self.sigma_range = sig_range_passed
            print("PairSampler initialized with passed sigma values.")
        else:
            # Use default initialization
            self.sigma, self.sigma_range = self.initialize_sigma()
            self.base_priors = Priors("jeffreys")
            print("PairSampler initialized with default sigma values.")
    
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
        Simple octet move without applying periodic boundaries. We retry the random
        displacement and rotation if any updated positions land outside the box. If,
        after several attempts, we still can't keep the octet in the box,
        we discard the move.
        """
        # Create a deep copy of positions
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # Get tetramers and form octets
        tetramers = self.ts.get_tetramers(positions)
        octets = self.get_octets(positions, tetramers)
        
        # If no octets found, return unchanged positions
        if not octets:
            return new_pos
        
        # Select a random octet
        octet_idx = np.random.randint(len(octets))
        tetramer1, tetramer2 = octets[octet_idx]
        
        # Collect all particles in the octet
        octet_particles = []
        for a_idx, b_idx, c_idx1, c_idx2 in [tetramer1, tetramer2]:
            octet_particles.extend([
                ('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)
            ])
        
        # Calculate centroid of the octet
        octet_coords = [new_pos[p][i] for p, i in octet_particles]
        centroid = np.mean(octet_coords, axis=0)
        
        # We'll allow a few tries to keep the octet within [0, box_size]
        max_tries = 10
        box_size = self.params.box_size
        
        for attempt in range(max_tries):
            # Copy positions so we can revert if needed
            temp_pos = {k: v.copy() for k, v in new_pos.items()}
            
            # 1. Apply translation
            displacement = np.random.normal(0, self.octet_trans_step, 3)
            
            # 2. Generate random rotation
            axis = self._random_unit_vector()
            angle = np.random.normal(0, self.octet_rot_step)
            rot_matrix = self._rotation_matrix(axis, angle)
            
            # Apply translation and rotation to all particles
            out_of_bounds = False
            for part, idx in octet_particles:
                vec = temp_pos[part][idx] - centroid
                rotated_pos = centroid + rot_matrix @ vec
                final_pos = rotated_pos + displacement
                
                # Check boundaries
                if np.any(final_pos < 0) or np.any(final_pos > box_size):
                    out_of_bounds = True
                    break
                
                # Update position
                temp_pos[part][idx] = final_pos
            
            # If everything is within bounds, accept and break
            if not out_of_bounds:
                new_pos = temp_pos
                break
        else:
            # If we exhaust max_tries with no success, return the original positions
            # (the move is effectively ignored)
            return positions
        
        return new_pos

    def get_octets(self, positions: Dict[str, np.ndarray], tetramers=None, temp=0.90) -> List[Tuple]:
        """
        Group tetramers into octets (pairs of tetramers) with temperature-based selection.
        
        Args:
            positions: Dictionary of particle positions
            tetramers: Pre-computed tetramers (optional)
            temp: Temperature parameter for proximity-based pairing
        
        Returns:
            List of octet tuples, each containing two tetramers
        """
        if tetramers is None:
            tetramers = self.ts.get_tetramers(positions)
        
        if len(tetramers) < 2:
            return []
        
        # Calculate centers of each tetramer
        centers = []
        for tetramer in tetramers:
            a_idx, b_idx, c_idx1, c_idx2 = tetramer
            coords = [
                positions['A'][a_idx],
                positions['B'][b_idx],
                positions['C'][c_idx1],
                positions['C'][c_idx2]
            ]
            centers.append(np.mean(coords, axis=0))
        
        # Form octets by pairing tetramers
        octets = []
        available = list(range(len(tetramers)))
        
        while len(available) >= 2:
            # Pick first tetramer randomly
            idx1 = np.random.choice(available)
            available.remove(idx1)
            
            # Calculate distances to all other tetramers with periodic boundaries
            distances = []
            for idx2 in available:
                delta = centers[idx1] - centers[idx2]
                # Periodic boundary correction
                delta = np.where(np.abs(delta) > self.params.box_size/2,
                            delta - np.sign(delta) * self.params.box_size,
                            delta)
                dist = np.linalg.norm(delta)
                distances.append(dist)
            
            # Convert distances to selection probabilities using temperature
            probs = np.exp(-np.array(distances) / temp)
            if probs.sum() > 0:
                probs = probs / probs.sum()  # Normalize
                idx2_rel = np.random.choice(len(available), p=probs)
                idx2 = available[idx2_rel]
            else:
                # If all probabilities are zero, choose randomly
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
        temp_start, temp_end = 5.0, 1.0
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
                proposed_positions = self.propose_position_move(self.positions_os, accepts['position'] / max(1, attempts['position']))
                proposed_sigma = self.sigma
            elif move_type == 'sigma':
                proposed_positions = self.positions_os
                proposed_sigma, _ = self.propose_sigma_move(self.sigma, accepts['sigma'] / max(1, attempts['sigma']))
            elif move_type == 'tetramer':
                proposed_positions = self.ts.propose_tetramer_move(self.positions_os, accepts['tetramer'] / max(1, attempts['tetramer']))
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
        Calculate octet scores in a vectorized operation.
        Each octet consists of two tetramers (Tet1 and Tet2),
        and we want to account for inter-tetramer pairs:
        A–A, A–B, B–B, B–C, and C–C
        so that the 'normal' pair sampler can exclude them.
        
        Args:
            positions: Dictionary of particle positions
            octets: List of pairs of tetramers [(tet1, tet2), ...]
            sig: Dictionary of sigma values
        
        Returns:
            numpy.ndarray: Array of octet scores (one per octet)
        """

        if not octets:
            return np.array([], dtype=np.float32)

        # Separate the two tetramers in each octet
        tetramers1 = []
        tetramers2 = []
        for tet1, tet2 in octets:
            tetramers1.append(tet1)
            tetramers2.append(tet2)

        # We can include the intra-tetramer terms by reusing the tetramer sampler
        # if you wish, but typically we handle intra-tetramer in a separate step.
        # For now, this function focuses on inter-tetramer interactions.

        # Unpack indices: (A, B, C1, C2)
        a1_1 = np.array([t[0] for t in tetramers1], dtype=np.int32)
        b1_1 = np.array([t[1] for t in tetramers1], dtype=np.int32)
        c1a_1 = np.array([t[2] for t in tetramers1], dtype=np.int32)
        c1b_1 = np.array([t[3] for t in tetramers1], dtype=np.int32)

        a1_2 = np.array([t[0] for t in tetramers2], dtype=np.int32)
        b1_2 = np.array([t[1] for t in tetramers2], dtype=np.int32)
        c1a_2 = np.array([t[2] for t in tetramers2], dtype=np.int32)
        c1b_2 = np.array([t[3] for t in tetramers2], dtype=np.int32)

        # We'll map them to their particle positions
        posA1_1 = positions['A'][a1_1]
        posA1_2 = positions['A'][a1_2]
        posB1_1 = positions['B'][b1_1]
        posB1_2 = positions['B'][b1_2]
        posC1a_1 = positions['C'][c1a_1]
        posC1a_2 = positions['C'][c1a_2]
        posC1b_1 = positions['C'][c1b_1]
        posC1b_2 = positions['C'][c1b_2]

        # We'll accumulate inter-tetramer scores for each octet in this array
        octet_scores = np.zeros(len(octets), dtype=np.float32)

        # Helper function to compute Gaussian negative log-likelihood for distances
        def gaussian_nll(dists, target, sigma_val):
            return ((dists - target)**2)/(2*sigma_val**2) + np.log(2*np.pi*sigma_val)

        # 1) A–A across tetramers
        delta_aa = posA1_1 - posA1_2
        dists_aa = np.sqrt(np.sum(delta_aa**2, axis=1))
        aa_target = self.params.pair_distances['AA']
        aa_score_arr = gaussian_nll(dists_aa, aa_target, sig['AA'])

        # 2) A–B across tetramers (two ways: A in Tet1 with B in Tet2, and B in Tet1 with A in Tet2)
        # A in Tet1 vs B in Tet2
        delta_ab_12 = posA1_1 - posB1_2
        dists_ab_12 = np.sqrt(np.sum(delta_ab_12**2, axis=1))
        ab_target = self.params.pair_distances['AB']
        ab_score_arr_12 = gaussian_nll(dists_ab_12, ab_target, sig['AB'])

        # B in Tet1 vs A in Tet2
        delta_ab_21 = posB1_1 - posA1_2
        dists_ab_21 = np.sqrt(np.sum(delta_ab_21**2, axis=1))
        # same AB target, same sigma
        ab_score_arr_21 = gaussian_nll(dists_ab_21, ab_target, sig['AB'])

        # 3) B–C across tetramers
        # We have 2 Cs in each tetramer. So we have B1_1 to (C1a_2, C1b_2) etc.
        # We'll just do this in a loop for clarity:
        bc_target = self.params.pair_distances['BC']
        bc_score_sum = np.zeros(len(octets), dtype=np.float32)
        for c_pos2 in [posC1a_2, posC1b_2]:
            # B in Tet1 vs C in Tet2
            delta_bc_12 = posB1_1 - c_pos2
            dists_bc_12 = np.sqrt(np.sum(delta_bc_12**2, axis=1))
            bc_score_sum += gaussian_nll(dists_bc_12, bc_target, sig['BC'])

        for c_pos1 in [posC1a_1, posC1b_1]:
            # C in Tet1 vs B in Tet2
            delta_bc_21 = c_pos1 - posB1_2
            dists_bc_21 = np.sqrt(np.sum(delta_bc_21**2, axis=1))
            bc_score_sum += gaussian_nll(dists_bc_21, bc_target, sig['BC'])

        # 4) C–C across tetramers (4 combos for each octet)
        cc_target = self.params.pair_distances['CC']
        cc_score_sum = np.zeros(len(octets), dtype=np.float32)
        c_pairs = [
            (posC1a_1, posC1a_2),
            (posC1a_1, posC1b_2),
            (posC1b_1, posC1a_2),
            (posC1b_1, posC1b_2),
        ]
        for p1, p2 in c_pairs:
            delta_cc = p1 - p2
            dists_cc = np.sqrt(np.sum(delta_cc**2, axis=1))
            cc_score_sum += gaussian_nll(dists_cc, cc_target, sig['CC'])

        # Combine them all
        # If you'd like to weight certain interactions, apply your weighting factors here
        octet_scores += aa_score_arr
        octet_scores += ab_score_arr_12
        octet_scores += ab_score_arr_21
        octet_scores += bc_score_sum
        octet_scores += cc_score_sum

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
        
        Steps:
        1) Exclude from the normal pair-sampler any pairs that belong to tetramers or octets,
            so we don't double-count them.
        2) Use the pair-sampler to compute the "background" pair cost for everything else.
        3) Calculate the tetramer scores (intra-tetramer) if needed separately
            (often done in calculate_tetramer_scores_batch).
        4) Calculate the new octet scores (all inter-tetramer interactions).
        5) Sum them, add the prior penalty, and return.

        Returns:
        (total_score, exclusion_score, pair_score, octet_score)
        """
        # Use the provided sigma if given, else use the sampler's
        sigma = sig if sig is not None else self.sigma

        # If no octets, just do the normal pair-sampler + tetramer logic and return
        if not octets:
            # The user wants no special octet logic if there are none
            score, ex_score, pair_score, _ = self.ps.calculate_score(
                positions, sigma, self.sigma_range,
                excluded_pairs=set(),  # no exclusions
                use_sigma_distribution=self.use_sigma_distribution,
                prior_penalty_from_distribution=prior_penalty
            )
            return score, ex_score, pair_score, 0.0

        # Build a set of all pairs belonging to:
        #   1) Each tetramer (intra-tetramer pairs: (A–B), (B–C), (C–C))
        #   2) Each octet (inter-tetramer pairs: A–A, A–B, B–C, C–C across tetra1 and tetra2)
        # so that these are excluded from the "normal" pair-sampler.
        excluded_pairs = set()

        # --- 1) Exclude intratetramer pairs ---
        for (a_idx, b_idx, c1_idx, c2_idx) in tetramers:
            # A–B
            excluded_pairs.add(('A', a_idx, 'B', b_idx))
            excluded_pairs.add(('B', b_idx, 'A', a_idx))
            # B–C
            excluded_pairs.add(('B', b_idx, 'C', c1_idx))
            excluded_pairs.add(('B', b_idx, 'C', c2_idx))
            excluded_pairs.add(('C', c1_idx, 'B', b_idx))
            excluded_pairs.add(('C', c2_idx, 'B', b_idx))
            # C–C
            excluded_pairs.add(('C', c1_idx, 'C', c2_idx))
            excluded_pairs.add(('C', c2_idx, 'C', c1_idx))

        # --- 2) Exclude inter-tetramer pairs for each octet ---
        for tet1, tet2 in octets:
            a1, b1, c1a1, c1b1 = tet1
            a2, b2, c2a1, c2b1 = tet2

            # A–A
            excluded_pairs.add(('A', a1, 'A', a2))
            excluded_pairs.add(('A', a2, 'A', a1))
            # A–B
            excluded_pairs.add(('A', a1, 'B', b2))
            excluded_pairs.add(('B', b2, 'A', a1))
            excluded_pairs.add(('B', b1, 'A', a2))
            excluded_pairs.add(('A', a2, 'B', b1))
            # B–C
            for c_idx2 in [c2a1, c2b1]:
                excluded_pairs.add(('B', b1, 'C', c_idx2))
                excluded_pairs.add(('C', c_idx2, 'B', b1))
            for c_idx1 in [c1a1, c1b1]:
                excluded_pairs.add(('C', c_idx1, 'B', b2))
                excluded_pairs.add(('B', b2, 'C', c_idx1))
            # C–C
            for c_idx1 in [c1a1, c1b1]:
                for c_idx2 in [c2a1, c2b1]:
                    excluded_pairs.add(('C', c_idx1, 'C', c_idx2))
                    excluded_pairs.add(('C', c_idx2, 'C', c_idx1))

        # --- 3) Let the pair sampler handle all non-excluded pairs ---
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            excluded_pairs=excluded_pairs,
            use_sigma_distribution=self.use_sigma_distribution,
            prior_penalty_from_distribution=prior_penalty
        )

        # --- 4) Calculate the octet score for all octets (A–A, A–B, B–C, C–C inter-tetramer) ---
        octet_scores = self.calculate_octet_scores_batch(positions, octets, sigma)
        total_octet_score = octet_scores.sum()

        # Combine everything
        weighted_octet_score = octet_weight * total_octet_score
        total_score = score + weighted_octet_score

        return total_score, ex_score, pair_score, weighted_octet_score