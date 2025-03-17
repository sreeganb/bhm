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
import numba as nb
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
        self.tetramer_trans_step = 0.1
        self.tetramer_rot_step = 0.1
        self.target_acceptance = 0.5
        self.tet_trans_acc_rate = self.target_acceptance
        self.sigma_prior_dist = {}
        
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
        
        visualize_3d_configuration(self.positions_ts, self.params.radii)
        
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
            self.gmm_params = self._load_gmm_parameters()
            self.sigma = self.initialize_sigma_from_gmm()
            self.sigma_range = None
        elif sig_passed is not None:
            self.sigma = sig_passed
            self.sigma_range = sig_range_passed
            print("using passed in sigma values")
        else:
            self.sigma, self.sigma_range = self.initialize_sigma()
            self.base_priors = Priors("jeffreys")
            print("Using default sigma initialization")


    def _load_gmm_parameters(self):
        """Load GMM parameters from JSON files without discarding valid data."""
        gmm_params = {sigma_type: None for sigma_type in ["AA", "AB", "BC", "CC"]}
        output_dir = os.path.join(os.getcwd(), "output_analysis/pairsampler_results")
        sampler_name = "PairSampler"
        
        if not os.path.exists(output_dir):
            print(f"Output directory not found: {output_dir}")
            return gmm_params
        
        pattern = re.compile(f"gmm_fit_(AA|AB|BC|CC)_{re.escape(sampler_name)}_chain_(\\d+)\\.json")
        matches = []
        for filename in os.listdir(output_dir):
            match = pattern.match(filename)
            if match:
                matches.append((match.group(1), int(match.group(2))))
        
        if not matches:
            print(f"No GMM fit files found in {output_dir}")
            return gmm_params
        
        chain_numbers = set(chain_num for _, chain_num in matches)
        if not chain_numbers:
            return gmm_params
        
        selected_chain = random.choice(list(chain_numbers))
        self.sel_chain_number = selected_chain
        
        file_paths = {
            sigma_type: os.path.join(output_dir, f"gmm_fit_{sigma_type}_{sampler_name}_chain_{selected_chain}.json")
            for sigma_type in gmm_params
        }
        
        for sigma_type, file_path in file_paths.items():
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'r') as f:
                        gmm_params[sigma_type] = json.load(f)
                except Exception as e:
                    print(f"Error loading {file_path}: {e}")
        
        print("the chosen gmm params are ", gmm_params)
        return gmm_params


    def initialize_sigma_from_gmm(self) -> Dict[str, float]:
        """Sample sigma values from the loaded GMMs if valid; otherwise use defaults."""
        sigma = {}
        defaults = self.params.pair_distances
        
        for pair_type in defaults:
            gmm_info = self.gmm_params.get(pair_type)
            
            # Only use GMM data if it exists and has the needed keys
            if (
                gmm_info
                and all(k in gmm_info for k in ['n_components', 'means', 'covariances', 'weights'])
            ):
                try:
                    n_components = gmm_info['n_components']
                    means = np.asarray(gmm_info['means']).reshape(n_components, 1)
                    covariances = np.asarray(gmm_info['covariances'])
                    weights = np.asarray(gmm_info['weights'])
                    
                    # Normalize weights if needed
                    wsum = np.sum(weights)
                    if abs(wsum - 1.0) > 1e-6:
                        weights = weights / wsum
                    
                    # Randomly pick a GMM component
                    component = np.random.choice(n_components, p=weights)
                    
                    # Reshape mean to (1,) and covariance to (1,1)
                    # so multivariate_normal sees a 2D, square covariance matrix
                    mean_value = means[component].flatten()
                    cov_value = float(covariances[component])
                    cov_matrix = np.array([[cov_value]])
                    
                    # Sample from the 1D GMM
                    sampled_value = np.random.multivariate_normal(mean=mean_value, cov=cov_matrix, size=1)
                    sigma[pair_type] = float(sampled_value)
                except (KeyError, ValueError, np.linalg.LinAlgError) as e:
                    print(f"Error sampling GMM for {pair_type}: {e}. Using default.")
                    sigma[pair_type] = defaults[pair_type]
            else:
                # No valid GMM data for this pair type
                sigma[pair_type] = defaults[pair_type]
        
        print("Final sigma values from GMM or defaults:", sigma)
        return sigma

    def get_positions(self) -> Dict[str, np.ndarray]:
        """Efficiently load positions from trajectory files with robust error handling."""
        # Use pathlib for more modern path handling
        import pathlib
        
        traj_dir = pathlib.Path(os.getcwd()) / "output_analysis/pairsampler_results"
        
        try:
            # Fast file matching using glob pattern
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
                    
                # Get last 10% of frames, ensuring at least one
                n_equil = max(1, int(len(keys) * 0.1))
                key = random.choice(keys[-n_equil:])
                print(f"Selected frame: {key}")
                
                # Efficiently extract position data directly using HDF5 dataset access
                positions = {}
                pos_grp = traj_grp[key]['positions']
                
                # Direct array copy is faster than iterative building
                for type_name in pos_grp:
                    positions[type_name] = pos_grp[type_name][:].copy()
                
            return positions
            
        except (FileNotFoundError, KeyError, ValueError, OSError) as e:
            print(f"Error loading trajectory: {e}")
            print("Falling back to initialized positions")
            return self.initialize_positions()

    def _calculate_gmm_log_prob(self, sigma_value: float, pair_type: str) -> float:
        """Vectorized GMM log probability calculation with numerical stability."""
        if not self.gmm_params or pair_type not in self.gmm_params or not self.gmm_params[pair_type]:
            return -np.inf
            
        gmm_info = self.gmm_params[pair_type]
        
        try:
            n_components = gmm_info['n_components']
            means = np.asarray(gmm_info['means']).reshape(n_components, 1)
            covariances = np.asarray(gmm_info['covariances']).reshape(n_components, 1, 1) 
            weights = np.asarray(gmm_info['weights'])
            
            # Ensure numerical stability
            covariances = np.maximum(covariances, 1e-10)
            
            # Vectorized component log probability calculation
            x = np.array(sigma_value).reshape(1)
            diff = x - means.flatten()
            exponents = -0.5 * (diff**2 / covariances.flatten())
            norms = np.log(weights) - 0.5 * np.log(2 * np.pi * covariances.flatten())
            component_log_probs = norms + exponents
            
            # Numerically stable log-sum-exp
            max_log_prob = np.max(component_log_probs)
            return max_log_prob + np.log(np.sum(np.exp(component_log_probs - max_log_prob)))
            
        except (KeyError, ValueError, np.linalg.LinAlgError) as e:
            return -np.inf

    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculate negative log prior across all sigma values efficiently."""
        if not sigma:
            return np.inf
            
        log_priors = [self._calculate_gmm_log_prob(value, pair_type) 
                    for pair_type, value in sigma.items()]
        
        # Filter out -inf values to avoid returning inf unnecessarily
        valid_priors = [lp for lp in log_priors if lp > -np.inf]
        
        if not valid_priors:
            return 1000.0  # Large but finite penalty for all invalid priors
            
        return -np.sum(valid_priors)

    def _compute_probabilities(self, dist_matrix: np.ndarray) -> np.ndarray:
        """
        Given a distance matrix, subtract the minimum distance along each row,
        exponentiate, and then normalize to get probability vectors.
        """
        min_vals = dist_matrix.min(axis=1, keepdims=True)
        exp_terms = np.exp(-(dist_matrix - min_vals))
        row_sums = exp_terms.sum(axis=1, keepdims=True) + 1e-10
        return exp_terms / row_sums

    def _random_indices_from_probs(self, prob_matrix: np.ndarray) -> np.ndarray:
        """
        For each row of prob_matrix, randomly sample a single index 
        according to the probabilities in that row.
        """
        n_rows, n_cols = prob_matrix.shape
        indices = np.empty(n_rows, dtype=int)
        for i in range(n_rows):
            indices[i] = np.random.choice(n_cols, p=prob_matrix[i])
        return indices

    def get_tetramers(self, positions: Dict[str, np.ndarray], temp: float = 0.9) -> List[Tuple[int, ...]]:
        """Simplified tetramer selection with clear helper routines."""
        try:
            # Basic checks
            for comp in ['A', 'B', 'C']:
                if comp not in positions or len(positions[comp]) == 0:
                    return []
            
            a_positions = positions['A']
            b_positions = positions['B']
            c_positions = positions['C']
            
            n_a, n_b, n_c = len(a_positions), len(b_positions), len(c_positions)
            if n_a == 0 or n_b == 0 or n_c < 2:
                return []
            
            # Distance from A to B
            dist_AB = cdist(a_positions, b_positions) / max(0.1, temp)
            probs_B = self._compute_probabilities(dist_AB)
            probs_B = np.nan_to_num(probs_B, nan=1.0/n_b)
            b_indices = self._random_indices_from_probs(probs_B)
            
            # Distance from chosen B to all C
            dist_BC = cdist(b_positions[b_indices], c_positions) / max(0.1, temp)
            probs_C = self._compute_probabilities(dist_BC)
            
            # Build the tetramers
            tetramers = []
            for a_idx in range(n_a):
                if np.sum(probs_C[a_idx]) < 1e-10 or n_c < 2:
                    c_indices = np.random.choice(n_c, size=2, replace=(n_c < 2))
                else:
                    p_norm = probs_C[a_idx] / (probs_C[a_idx].sum() + 1e-10)
                    c_indices = np.random.choice(n_c, size=2, replace=False, p=p_norm)
                tetramers.append((a_idx, b_indices[a_idx], c_indices[0], c_indices[1]))

            return tetramers

        except Exception as e:
            print(f"Error generating tetramers: {e}")
            return [] # Return empty list on error for graceful failure
#-----------------------------------------------------------------------  
    def run_mc(self, n_steps: int = 50000, save_freq: int = 100, 
            output_dir: str = "output_analysis/tetramersampler_results/") -> Tuple:
        """Monte Carlo sampling with tetramer moves and efficient caching.
        Runs until n_steps accepted moves (Markov chain steps) are completed.
        Overlap checking has been removed."""
        # Pre-allocate memory for results & tracking
        best_positions = None
        best_score = float('inf')
        trajectory = []
        sigma_history = {key: np.zeros(n_steps // save_freq + 1, dtype=np.float32) for key in self.sigma}
        
        # Setup output directory and file handles
        os.makedirs(output_dir, exist_ok=True)
        csv_log_file = os.path.join(output_dir, "all_info_mcmc_tetramer.csv")
        #debug_file = os.path.join(output_dir, "debug_mcmc.txt")
        csv_buffer = []  # Buffer for delayed CSV writes
        
        # Debug file setup 
#        with open(debug_file, "w") as f:
#            f.write("# TetramerSampler MCMC Debug Log\n")
#            f.write("# Overlap checking has been removed\n\n")
        
        # Initialize counters with pre-allocation
        moves_counts = {'position': 0, 'sigma': 0, 'tetramer': 0}
        accepts_counts = {'position': 0, 'sigma': 0, 'tetramer': 0}
        
        # Write CSV header once
        with open(csv_log_file, "w") as f:
            f.write("Step,Prior,Exvol_score,Pair_score,Tet_score,Score,Accepted\n")
        
        # Temperature schedule calculation (vectorized)
        # Changed to use accepted_moves rather than total iterations
        cooling_factor = -np.log(1.0/5.0) / n_steps
        get_temp = lambda accepted_moves: 5.0 * np.exp(-cooling_factor * accepted_moves)
        
        # Calculate initial state with caching
        current_tetramers = self.get_tetramers(self.positions_ts)
        prior_penalty = (self.calculate_negative_log_prior(self.sigma) if self.use_sigma_distribution 
                        else self.base_priors.neg_log_prior(self.sigma, self.sigma_range))
        current_score, curr_ex, curr_pair, curr_tet = self.neg_log_posterior(
            self.positions_ts, current_tetramers, prior_penalty, self.sigma)
        
        # Debug initial scores
#        with open(debug_file, "a") as f:
#            f.write(f"INITIAL STATE\n")
#            f.write(f"Initial Score: {current_score:.1f}, ExVol: {curr_ex:.3f}, ")
#            f.write(f"Pair: {curr_pair:.1f}, Tet: {curr_tet:.1f}, Prior: {prior_penalty:.1f}\n\n")
        
        # Store initial sigma values
        for i, (key, value) in enumerate(self.sigma.items()):
            sigma_history[key][0] = value
        
        # Main MCMC loop with optimizations
        move_probs = np.array([0.4, 0.1, 0.5])  # position, sigma, tetramer
        move_types = ['position', 'sigma', 'tetramer']
        
        # Progress reporting setup
        print_freq = max(1, min(n_steps // 20, 1000))
        
        # Modified loop to run until we have n_steps accepted moves
        accepted_moves = 0
        total_moves = 0
        max_iterations = n_steps * 50  # Safety cap
        
        while accepted_moves < n_steps and total_moves < max_iterations:
            total_moves += 1
            
            # Temperature based on accepted moves, not iterations
            temp = get_temp(accepted_moves)
            
            # Fast random move selection
            move_type = move_types[np.random.choice(3, p=move_probs)]
            moves_counts[move_type] += 1
            
            # Move proposal with caching
            cached_tetramers = None
            if move_type == 'position':
                proposed_positions = self.propose_position_move(self.positions_ts)
                proposed_sigma = self.sigma  # No change, just reference
            elif move_type == 'sigma':
                proposed_positions = self.positions_ts  # No change, just reference
                proposed_sigma, changed_pair_type = self.propose_sigma_move(self.sigma)
                cached_tetramers = current_tetramers  # Reuse tetramers for sigma-only moves
            else:  # tetramer move
                proposed_positions = self.propose_tetramer_move(self.positions_ts)
                proposed_sigma = self.sigma  # No change, just reference
                self.tet_trans_acc_rate = accepts_counts['tetramer'] / max(1, moves_counts['tetramer'])
            
            # Calculate prior efficiently
            new_prior_penalty = (self.calculate_negative_log_prior(proposed_sigma) if self.use_sigma_distribution 
                            else self.base_priors.neg_log_prior(proposed_sigma, self.sigma_range))
            
            # Optimize tetramer calculation with caching
            proposed_tetramers = cached_tetramers if cached_tetramers is not None else self.get_tetramers(proposed_positions)
            
            # Score calculation
            proposed_score, prop_ex, prop_pair, prop_tet = self.neg_log_posterior(
                proposed_positions, proposed_tetramers, new_prior_penalty, proposed_sigma)
            
            # Metropolis acceptance criterion (optimized)
            delta_e = proposed_score - current_score
            accept_move = delta_e < 0 or np.random.random() < np.exp(-delta_e / max(temp, 1e-10))
            
            # Debug logging (simplified, overlap info removed)
#            if total_moves % 1000 == 0:
#                with open(debug_file, "a") as f:
#                    f.write(f"\nSUMMARY at move {total_moves}:\n")
#                    f.write(f"- Accepted moves: {accepted_moves}/{n_steps}\n")
#                    f.write(f"- Overall acceptance rate: {accepted_moves/max(1,total_moves):.3f}\n")
#                    f.write("-" * 50 + "\n\n")
            
            # Update state if accepted
            if accept_move:
                self.positions_ts = proposed_positions
                self.sigma = proposed_sigma
                current_score, curr_ex, curr_pair, curr_tet = proposed_score, prop_ex, prop_pair, prop_tet
                current_tetramers = proposed_tetramers if cached_tetramers is None else current_tetramers
                accepts_counts[move_type] += 1
                accepted_moves += 1
                
                # Track best configuration
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_ts.items()}
            
                # Efficient batch logging based on accepted moves
                if accepted_moves % save_freq == 0:
                    save_idx = accepted_moves // save_freq
                    if save_idx < len(sigma_history[list(sigma_history.keys())[0]]):
                        for key in self.sigma:
                            sigma_history[key][save_idx] = self.sigma[key]
                    
                    trajectory.append(
                        self.save_state(
                            accepted_moves, self.positions_ts, self.sigma, current_score,
                            prior_score=new_prior_penalty, pair_score=curr_pair,
                            exvol_score=curr_ex, tet_score=curr_tet
                        )
                    )
                    
                    csv_buffer.append(
                        f"{accepted_moves},{new_prior_penalty:.1f},{curr_ex:.1f},{curr_pair:.1f},"
                        f"{curr_tet:.1f},{current_score:.1f},1\n"
                    )
                    
                    if len(csv_buffer) >= 10:
                        with open(csv_log_file, "a") as f:
                            f.writelines(csv_buffer)
                        csv_buffer = []
            
                if accepted_moves % print_freq == 0:
                    overall_accept = accepted_moves / total_moves
                    print(f"Accepted move {accepted_moves}/{n_steps} ({accepted_moves/n_steps*100:.1f}%), "
                        f"Score: {current_score:.1f}, Temp: {temp:.2f}, "
                        f"Accept: {overall_accept:.2f}")
        
        if csv_buffer:
            with open(csv_log_file, "a") as f:
                f.writelines(csv_buffer)
        
        sigma_history_df = pd.DataFrame({k: v[:n_steps//save_freq+1] for k, v in sigma_history.items()})
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        final_file = self.save_trajectory(trajectory, trajectory_file)
        
#        with open(debug_file, "a") as f:
#            f.write("\nFINAL SUMMARY:\n")
#            f.write(f"- Total iterations: {total_moves}\n")
#            f.write(f"- Accepted moves: {accepted_moves}\n")
#            f.write(f"- Overall acceptance rate: {accepted_moves/total_moves:.3f}\n")
#            f.write(f"- Final score: {current_score:.2f}\n")
#            f.write(f"- Best score: {best_score:.2f}\n")
        
        print(f"\nSampling complete: {accepted_moves} accepted moves out of {total_moves} iterations")
        for move_type in moves_counts:
            rate = accepts_counts[move_type] / max(1, moves_counts[move_type])
            print(f"{move_type.capitalize()} moves: {rate*100:.1f}% acceptance ({accepts_counts[move_type]}/{moves_counts[move_type]})")
        
        return best_positions, trajectory, final_file
#-----------------------------------------------------------------------
    def propose_tetramer_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Optimized proposal generation for tetramer moves with multiple strategies."""
        # Create a deep copy of positions - unavoidable for MCMC proposal
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # Get tetramers with temperature parameter that encourages exploration
        tetramers = self.get_tetramers(positions, temp=0.9)
        
        if not tetramers:
            return new_pos
            
        # Define move types and normalized probabilities (ensuring they sum to exactly 1.0)
        move_types = ['single_tetramer', 'coordinated_pair', 'radial', 'global_rotation', 'aggressive']
        #probs = np.array([0.99, 0.0015, 0.0025, 0.0015, 0.0045])
        probs = np.array([0.89, 0.0, 0.03, 0.07, 0.01])
        probs /= np.sum(probs)  # Normalize to ensure sum is exactly 1.0
        
        # Choose a move type with probability
        move_type = np.random.choice(move_types, p=probs)
        
        # Calculate adaptive step size once
        trans_base = self.tetramer_trans_step
        rot_base = self.tetramer_rot_step
        trans_adjust = np.clip(1.0 + 5.0 * (self.tet_trans_acc_rate - self.target_acceptance), 0.3, 2.0)
        
        # Use a strategy pattern with a dispatch dictionary
        return {
            'single_tetramer': lambda: self._move_single_tetramer(new_pos, tetramers, trans_adjust, trans_base, rot_base),
            'coordinated_pair': lambda: self._move_tetramer_pair(new_pos, tetramers, trans_adjust, trans_base, rot_base),
            'radial': lambda: self._move_radial_symmetry(new_pos, tetramers, trans_adjust, trans_base),
            'global_rotation': lambda: self._move_global_rotation(new_pos, tetramers, rot_base * trans_adjust),
            'aggressive': lambda: self._move_single_tetramer(new_pos, tetramers, trans_adjust * 2.0, trans_base, rot_base)
        }[move_type]()

    def _apply_boundary_conditions(self, positions):
        """Apply periodic boundary conditions efficiently."""
        # Vectorized boundary handling
        box_size = self.params.box_size
        return np.mod(positions, box_size)

    def _calculate_system_center(self, positions):
        """Calculate system center efficiently."""
        # Fast concatenation of all positions
        all_coords = np.vstack([pos for pos in positions.values()])
        return np.mean(all_coords, axis=0)

    def _move_single_tetramer(self, new_pos, tetramers, trans_adjust, trans_base, rot_base):
        """Optimized single tetramer move with vectorized operations."""
        # Select a random tetramer
        a_idx, b_idx, c_idx1, c_idx2 = tetramers[np.random.randint(len(tetramers))]
        particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]
        
        # Collect tetramer coordinates as a NumPy array for efficient operations
        coords_array = np.array([new_pos[part][idx] for part, idx in particles])
        
        # Generate a single 3D displacement vector for all particles
        trans_step = trans_base * trans_adjust * np.random.uniform(0.5, 1.5)
        displacement = np.random.normal(0, trans_step, 3)  # Fixed: Generate a 3D vector
            
        # Apply the same displacement to all particles
        for i, (part, idx) in enumerate(particles):
            new_pos[part][idx] = self._apply_boundary_conditions(coords_array[i] + displacement)
        
        # Apply rotation with probability
        if np.random.random() < 0.9:
            # Use pre-calculated centroid
            centroid = np.mean(coords_array, axis=0)
            
            # Generate rotation parameters
            rot_step = rot_base * trans_adjust * np.random.uniform(0.7, 1.3)
            rotation_axis = self._random_unit_vector()
            rotation_angle = np.random.normal(0, rot_step)
            rot_matrix = self._rotation_matrix(rotation_axis, rotation_angle)
            
            # Apply rotation to all particles
            for i, (part, idx) in enumerate(particles):
                vec = coords_array[i] - centroid
                new_pos[part][idx] = self._apply_boundary_conditions(centroid + rot_matrix @ vec)
        
        return new_pos

    def _move_tetramer_pair(self, new_pos, tetramers, trans_adjust, trans_base, rot_base):
        """Optimized pair move with reduced redundancy."""
        if len(tetramers) < 2:
            return self._move_single_tetramer(new_pos, tetramers, trans_adjust, trans_base, rot_base)
        
        # Select two tetramers
        indices = np.random.choice(len(tetramers), size=2, replace=False)
        
        # Create single displacement for both tetramers
        trans_step = trans_base * trans_adjust
        displacement = np.random.normal(0, trans_step, 3)
        
        # Apply to both tetramers efficiently
        for idx in indices:
            a_idx, b_idx, c_idx1, c_idx2 = tetramers[idx]
            for part, idx in [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]:
                new_pos[part][idx] = self._apply_boundary_conditions(new_pos[part][idx] + displacement)
        
        return new_pos

    def _move_radial_symmetry(self, new_pos, tetramers, trans_adjust, trans_base):
        """Efficient radial movement preserving symmetry."""
        if len(tetramers) < 3:
            return self._move_single_tetramer(new_pos, tetramers, trans_adjust, trans_base, self.tetramer_rot_step)
        
        # Calculate system center once
        system_center = self._calculate_system_center(new_pos)
        
        # Single radial adjustment for all tetramers
        radial_adjust = np.random.normal(0, trans_base * trans_adjust)
        
        # Process all tetramer particles
        particle_indices = []
        for tetramer in tetramers:
            a_idx, b_idx, c_idx1, c_idx2 = tetramer
            particle_indices.extend([('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)])
        
        # Apply radial adjustment to all particles
        for part, idx in particle_indices:
            vec = new_pos[part][idx] - system_center
            distance = np.linalg.norm(vec)
            if distance > 1e-10:  # Avoid division by zero with better threshold
                new_distance = max(0.1, distance + radial_adjust)
                new_pos[part][idx] = self._apply_boundary_conditions(
                    system_center + vec * (new_distance / distance)
                )
        
        return new_pos

    def _move_global_rotation(self, new_pos, tetramers, rot_step):
        """Highly optimized global rotation around system center."""
        # Calculate system center once
        system_center = self._calculate_system_center(new_pos)
        
        # Generate rotation parameters
        rotation_axis = self._random_unit_vector()
        rotation_angle = np.random.normal(0, rot_step)
        rotation_matrix = self._rotation_matrix(rotation_axis, rotation_angle)
        
        # Apply rotation to all particles at once for each component type
        for type_name, positions in new_pos.items():
            # Vectorized operation: translate to origin, rotate, translate back
            centered = positions - system_center
            # Apply rotation to all particles of this type at once
            rotated = np.dot(centered, rotation_matrix.T)
            new_pos[type_name] = self._apply_boundary_conditions(rotated + system_center)
        
        return new_pos

    def _random_unit_vector(self):
        """Generate a random unit vector with numerical stability."""
        vec = np.random.randn(3)
        norm = np.linalg.norm(vec)
        if norm < 1e-10:
            vec = np.array([0, 0, 1])  # Default vector if random is near zero
        else:
            vec /= norm
        return vec

    def _rotation_matrix(self, axis, theta):
        """
        Create a rotation matrix using quaternions - more efficient and stable
        than Rodrigues' formula for repeated calculations.
        """
        # Normalize axis
        axis = np.asarray(axis)
        axis = axis / np.linalg.norm(axis)
        
        # Quaternion representation (more stable)
        q = np.array([np.cos(theta/2), *(-axis*np.sin(theta/2))])
        
        # Fast quaternion-based rotation matrix construction
        return np.array([
            [1-2*(q[2]**2+q[3]**2), 2*(q[1]*q[2]-q[0]*q[3]), 2*(q[1]*q[3]+q[0]*q[2])],
            [2*(q[1]*q[2]+q[0]*q[3]), 1-2*(q[1]**2+q[3]**2), 2*(q[2]*q[3]-q[0]*q[1])],
            [2*(q[1]*q[3]-q[0]*q[2]), 2*(q[2]*q[3]+q[0]*q[1]), 1-2*(q[1]**2+q[2]**2)]
        ])
        
    def calculate_tetramer_score(self, positions: Dict[str, np.ndarray], 
                                tetramer: Tuple[int, ...], sig: Dict[str, float] = None) -> float:
        """Calculate score for a single tetramer with vectorized operations."""
        # Extract indices once
        a_idx, b_idx, c1_idx, c2_idx = tetramer
        
        # Pre-extract positions to avoid dict lookups in calculations
        pos_a = positions['A'][a_idx]
        pos_b = positions['B'][b_idx]
        pos_c1 = positions['C'][c1_idx]
        pos_c2 = positions['C'][c2_idx]
        
        # Calculate all scores in one batch
        score = self.pair_score_nll(pos_a, pos_b, self.params.pair_distances['AB'], sig['AB'])
        score += self.pair_score_nll(pos_b, pos_c1, self.params.pair_distances['BC'], sig['BC'])
        score += self.pair_score_nll(pos_b, pos_c2, self.params.pair_distances['BC'], sig['BC']) 
        score += self.pair_score_nll(pos_c1, pos_c2, self.params.pair_distances['CC'], sig['CC'])
        
        return score

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        tetramers: List[Tuple[int, ...]],
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
    ) -> Tuple[float, float, float, float]:
        """
        Efficiently calculate the total score for a tetramer system with vectorized operations.
        
        Returns:
            Tuple of (total_score, exclusion_score, pair_score, tetramer_score)
        """
        # Use sigma if provided, otherwise use class sigma
        sigma = sig if sig is not None else self.sigma
        
        # Quick return if no tetramers
        if not tetramers:
            # Calculate score with no excluded pairs
            score, ex_score, pair_score, _ = self.ps.calculate_score(
                positions, sigma, self.sigma_range, 
                set(), self.use_sigma_distribution, 
                prior_penalty_from_distribution
            )
            return score, ex_score, pair_score, 0.0
        
        # 1) Efficiently build set of tetramer pairs
        # Pre-allocate with expected size (4 pairs per tetramer)
        n_tetramers = len(tetramers)
        tetramer_pairs = set()
        tetramer_pairs_add = tetramer_pairs.add  # Local reference for faster calls
        
        # Build set of pairs with optimized batch processing
        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
            tetramer_pairs_add(('A', a_idx, 'B', b_idx))
            tetramer_pairs_add(('B', b_idx, 'C', c1_idx))
            tetramer_pairs_add(('B', b_idx, 'C', c2_idx))
            tetramer_pairs_add(('C', c1_idx, 'C', c2_idx))

        # 2) Calculate score excluding tetramer pairs
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            tetramer_pairs, self.use_sigma_distribution,
            prior_penalty_from_distribution
        )

        # 3) Calculate tetramer score - using a comprehension for better performance
        if not tetramers:
            return score, ex_score, pair_score, 0.0
            
        # Calculate all tetramer scores in one go
        total_tet_score = sum(self.calculate_tetramer_score(positions, tetramer, sigma) 
                            for tetramer in tetramers)
        
        # Apply weighting and return all score components
        weighted_tet_score = tetramer_weight * total_tet_score
        score += weighted_tet_score

        return score, ex_score, pair_score, weighted_tet_score
#-----------------------------------------------------------------------
