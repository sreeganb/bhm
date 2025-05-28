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
import networkx as nx  # For graph operations
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
        
        # Need a method to take in sampler sequence and find out at which sampler level we are and then 
        # initialize and read sigma values accordingly, for example if the sequence is 
        # ["pair", "tetramer", "octet", "pair", "octet", tetramer] and we are at element 4 in this sequence then the 
        # code should automatically take the sigma values from element 3 which is "pair" and read in the GMM values from this 

        # Basic sampler parameters
        self.octet_trans_step = 0.05
        self.octet_rot_step = 0.05
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
        
    def propose_octet_move(self, positions: Dict[str, np.ndarray], octets) -> Dict[str, np.ndarray]:
        """
        Propose a move for a randomly selected octet by applying a small translation and rotation.
        This version does not apply periodic boundary conditions.
        """
        # Create a deep copy of positions to modify
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # Get tetramers and then form octets
        #tetramers = self.ts.get_tetramers(new_pos) # Use new_pos in case get_tetramers modifies it
        #octets = self.get_octets(new_pos, tetramers)
        
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
    
    def get_octets(self, positions: Dict[str, np.ndarray]) -> Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], List[Tuple[int, ...]]]:
        tetramers = self.ts.get_tetramers(positions)
        #print("tetramers: ", tetramers)

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

        #print(f"Octets: {octets}")
        return octets, tetramers
    
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
        n_steps is the total number of moves (accepted + rejected).
        Includes Jacobian correction for sigma moves and detailed debugging logs.
        """
        # Setup output directory
        os.makedirs(output_dir, exist_ok=True)
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        # Debug log files
        general_debug_log_file = os.path.join(output_dir, "mcmc_debug_log.txt")
        sigma_debug_log_file = os.path.join(output_dir, "sigma_debug.txt")

        # Clear/Create debug log files at the start of the run
        with open(general_debug_log_file, 'w') as f:
            f.write("Step,MoveType,Accepted,Temp,CurrentScore,ProposedScore,Delta,Exclusion,Pair/Tet,Octet,Prior\n")
        with open(sigma_debug_log_file, 'w') as f:
            f.write("Step,CurrentSigma,ProposedSigma,CurrentPrior,ProposedPrior,CurrentScore,ProposedScore,Delta,AcceptProb,Accepted\n")
        # Create empty trajectory file
        with h5py.File(trajectory_file, 'w') as f:
            pass

        # Initialize tracking variables
        best_positions = {k: v.copy() for k, v in self.positions_os.items()}  # Init best pos
        best_score = float('inf')
        sigma_history = {key: np.zeros(n_steps // save_freq + 1) for key in self.sigma}
        accepts = {'position': 0, 'sigma': 0, 'tetramer': 0, 'octet': 0}
        attempts = {'position': 0, 'sigma': 0, 'tetramer': 0, 'octet': 0}

        # --- Initialize state ---
        current_positions = {k: v.copy() for k, v in self.positions_os.items()}  # Use a copy to avoid modifying original
        current_sigma = self.sigma.copy()
        current_octets, current_tetramers = self.get_octets(current_positions)

        # Calculate initial prior penalty
        current_prior = (self.sig_provider.calculate_negative_log_prior(current_sigma) if self.use_sigma_distribution
                        else self.base_priors.neg_log_prior(current_sigma, self.sigma_range))

        # Calculate initial score
        current_score, curr_ex, curr_pair, curr_tet, curr_oct = self.neg_log_posterior(
            current_positions, current_tetramers, current_octets, current_prior, current_sigma)

        if not np.isfinite(current_score):
            print(f"FATAL: Initial score is non-finite ({current_score}). Exiting.")
            print(f"Components: ex={curr_ex}, tet={curr_tet}, oct={curr_oct}, prior={current_prior}")
            return None, None  # Cannot proceed

        best_score = current_score  # Initialize best score

        # Store initial sigma values
        for key in current_sigma:
            sigma_history[key][0] = current_sigma[key]

        # Save initial state
        self.save_state_to_disk(
            0, current_positions, current_sigma, current_score,
            prior_score=current_prior,
            tet_score=curr_tet,
            exvol_score=curr_ex,
            oct_score=curr_oct,
            traj_file=trajectory_file
        )

        # --- Main MCMC loop parameters ---
        move_types = ['position', 'sigma', 'tetramer', 'octet']
        move_probs = [0.2, 0.1, 0.3, 0.4]  # Adjust probabilities if needed later
        temp_start, temp_end = 40.0, 1.0  # Cooling schedule
        temp_decay = (temp_end / temp_start) ** (1.0 / (n_steps - 1)) if n_steps > 1 else 1.0

        print(f"Starting MCMC sampling for {n_steps} total steps...")

        # --- Main MCMC loop ---
        for step in range(n_steps):
            # Temperature schedule based on total steps
            temp = max(temp_start * (temp_decay ** step), temp_end)

            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1

            # --- Propose move ---
            proposed_positions = {k: v.copy() for k, v in current_positions.items()}
            proposed_sigma = current_sigma.copy()
            pair_type = None  # Will track which sigma parameter was modified

            if move_type == 'position':
                proposed_positions = self.propose_position_move(current_positions,
                                                                accepts['position'] / max(1, attempts['position']))
            elif move_type == 'sigma':
                proposed_sigma, pair_type = self.propose_sigma_move(current_sigma,
                                                                    accepts['sigma'] / max(1, attempts['sigma']))
            elif move_type == 'tetramer':
                proposed_positions = self.ts.propose_tetramer_move(current_positions,
                                                                accepts['tetramer'] / max(1, attempts['tetramer']))
            else:  # octet
                proposed_positions = self.propose_octet_move(current_positions, current_octets)

            # --- Calculate components for the proposed state ---
            proposed_octets, proposed_tetramers = self.get_octets(proposed_positions)
            proposed_prior = (self.sig_provider.calculate_negative_log_prior(proposed_sigma) if self.use_sigma_distribution
                            else self.base_priors.neg_log_prior(proposed_sigma, self.sigma_range))
            proposed_score, prop_ex, prop_pair, prop_tet, prop_oct = self.neg_log_posterior(
                proposed_positions, proposed_tetramers, proposed_octets, proposed_prior, proposed_sigma)

            # --- Check for non-finite scores ---
            if not np.isfinite(proposed_score) or not np.isfinite(current_score):
                print(f"WARNING @ step {step}: Non-finite score detected! Proposed={proposed_score}, Current={current_score}. Skipping move.")
                with open(general_debug_log_file, 'a') as f:
                    f.write(f"{step},{move_type},REJECTED(InvalidScore),{temp:.3f},{current_score:.4f},{proposed_score:.4f},NaN,"
                            f"{prop_ex:.4f},{prop_tet:.4f},{prop_oct:.4f},{proposed_prior:.4f}\n")
                continue

            # --- Metropolis criterion with Jacobian correction ---
            delta = proposed_score - current_score
            accept_prob = 0.0

            if move_type == 'sigma' and pair_type is not None:
                # Jacobian correction for sigma moves proposed in log-space
                jacobian_term = np.log(proposed_sigma[pair_type] / current_sigma[pair_type])
                delta += jacobian_term
                if step % save_freq == 0:
                    print(f"  Sigma move: {pair_type} {current_sigma[pair_type]:.4f}->{proposed_sigma[pair_type]:.4f}, "
                        f"Jacobian term: {jacobian_term:.4f}")

            if delta < 0:
                accept_prob = 1.0
                accept = True
            else:
                accept_prob = np.exp(-delta / temp) if temp > 1e-9 else 0.0
                accept = np.random.random() < accept_prob

            # --- Detailed Logging ---
            if move_type == 'sigma':
                with open(sigma_debug_log_file, 'a') as f:
                    f.write(f"{step},"
                            f"\"{current_sigma}\",\"{proposed_sigma}\","
                            f"{current_prior:.4f},{proposed_prior:.4f},"
                            f"{current_score:.4f},{proposed_score:.4f},"
                            f"{delta:.4f},{accept_prob:.4g},"
                            f"{'ACCEPTED' if accept else 'REJECTED'}\n")

            log_ex, log_tet, log_oct, log_prior = (prop_ex, prop_tet, prop_oct, proposed_prior) if accept else (curr_ex, curr_tet, curr_oct, current_prior)
            with open(general_debug_log_file, 'a') as f:
                f.write(f"{step},{move_type},{'ACCEPTED' if accept else 'REJECTED'},"
                        f"{temp:.3f},{current_score:.4f},{proposed_score:.4f},{delta:.4f},"
                        f"{log_ex:.4f},{log_tet:.4f},{log_oct:.4f},{log_prior:.4f}\n")

            # --- Update state if accepted ---
            if accept:
                current_positions = proposed_positions
                current_sigma = proposed_sigma
                current_tetramers = proposed_tetramers
                current_octets = proposed_octets
                current_score = proposed_score
                curr_ex, curr_tet, curr_oct = prop_ex, prop_tet, prop_oct
                current_prior = proposed_prior
                accepts[move_type] += 1

                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in current_positions.items()}

            # --- Save state every save_freq steps ---
            if step % save_freq == 0 and step > 0:  # Skip step 0 since initial state is saved above
                save_idx = step // save_freq
                for key in current_sigma:
                    sigma_history[key][save_idx] = current_sigma[key]
                self.save_state_to_disk(
                    step, current_positions, current_sigma, current_score,
                    prior_score=current_prior,
                    tet_score=curr_tet,
                    exvol_score=curr_ex,
                    oct_score=curr_oct,
                    traj_file=trajectory_file
                )
                acceptance_rate = sum(accepts.values()) / (step + 1)
                print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, AcceptRate={acceptance_rate:.2f}")

            # --- Adapt step sizes every 500 steps ---
            if step % 500 == 0 and step > 0:
                pos_rate = accepts['position'] / max(1, attempts['position'])
                tet_rate = accepts['tetramer'] / max(1, attempts['tetramer'])
                oct_rate = accepts['octet'] / max(1, attempts['octet'])
                sigma_rate = accepts['sigma'] / max(1, attempts['sigma'])
                print(f"Step {step}: Acceptance rates - Pos: {pos_rate:.2f}, Sigma: {sigma_rate:.2f}, Tet: {tet_rate:.2f}, Oct: {oct_rate:.2f}")
                self.octet_trans_acc_rate = oct_rate  # Update rate for propose_octet_move

        # --- End of MCMC loop ---
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)

        # Print final statistics
        print("\nSampling complete:")
        for mv_type in move_types:
            rate = accepts[mv_type] / max(1, attempts[mv_type])
            print(f"- {mv_type}: {rate:.2f} acceptance ({accepts[mv_type]}/{attempts[mv_type]})")
        total_accepts = sum(accepts.values())
        print(f"- Total steps: {n_steps}")
        print(f"- Accepted moves: {total_accepts}")
        print(f"- Overall acceptance rate: {total_accepts / n_steps:.2f}")
        print(f"- Best score: {best_score:.4f}")
        print(f"- Final sigma values: {current_sigma}")

        return best_positions, trajectory_file
    
    def calculate_octet_scores_batch(self, positions, octets, sig):
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
        #print(f"AA distances: {aa_dists}")
        
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
        tetramers: List[Tuple[int, ...]],
        octets: List[Tuple[Tuple[int, ...], Tuple[int, ...]]],
        prior_penalty: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
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

        excluded_pairs = None
        
        # 3) Calculate pair scores
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            excluded_pairs=excluded_pairs,
            use_sigma_distribution=self.use_sigma_distribution,
            prior_penalty_from_distribution=prior_penalty
        )
        
        # 2) Calculate tetramer scores 
        tetramer_score = 0.0
        tet_score = self.ts.calculate_tetramer_scores_batch(positions, tetramers, 
                                                            sigma)
        tetramer_score = tetramer_weight * tet_score.sum()
        
        # 3) Octet score
        octet_score = 0.0
        octet_scores = self.calculate_octet_scores_batch(positions, octets, sigma)
        octet_score = octet_weight * octet_scores.sum()
        
        # 4) Total score
        total_score = score + tetramer_score + octet_score
            
        #return total_score, ex_score, pair_score, weighted_octet_score
        return total_score, ex_score, pair_score, tetramer_score, octet_score
