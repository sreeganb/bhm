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
#-----------------------------------------------------------------------
class OctetSampler(BaseMCSampler):
    """
    Sampler for Tetramer-level interactions, inheriting from BaseMCSampler.
    Includes run_mc method.
    """
    def __init__(self, use_sigma_distribution=False):
        super().__init__()  # Call BaseMCSampler constructor
        self.use_sigma_distribution = use_sigma_distribution
        self.params = SystemParameters()  # Initialize system parameters
        self.positions = self.initialize_positions() # Initialize positions
        self.ps = PairSampler() # Initialize PairSampler
        self.ts = TetramerSampler() # Initialize TetramerSampler

        # Additional initialization for tetramer-specific features
        self.octet_trans_step = 0.25
        self.octet_rot_step = 0.1
        self.sigma_prior_dist = {}
        #-------------------------------------
        # tracking acceptance rates
        self.octet_trans_acc_rate = 0.5
        self.target_acceptance = 0.5
        #-------------------------------------

        required_pairs = {"AA", "AB", "BC", "CC"}
        if not required_pairs.issubset(self.params.pair_distances.keys()):
            raise ValueError("Missing required pair types for tetramer sampling.")

        if use_sigma_distribution:
            print("using sigma distribution now: ")
            pwd = os.getcwd()
            sampler_output_dir = os.path.join(pwd, "output_analysis/tetramersampler_results") # hard coded here, needs to be changed later
            sampler_name = "TetramerSampler" # hard coded here, needs to be changed later
            if os.path.exists(sampler_output_dir):  # Check if the directory exists
                pattern = re.compile(r"gmm_fit_(AA|AB|BC|CC)_" + re.escape(sampler_name) + r"_chain_(\d+)\.json")
                chain_numbers = set()

                for filename in os.listdir(sampler_output_dir):
                    match = pattern.match(filename)
                    if match:
                        chain_numbers.add(int(match.group(2)))
                        print("chain numbers: ", chain_numbers)

                if chain_numbers:  # Check if any chains were found
                    selected_chain = random.choice(list(chain_numbers))
                    self.gmm_params = {}  # Store parameters in an instance variable
                    for sigma_type in ["AA", "AB", "BC", "CC"]:
                        json_path = os.path.join(sampler_output_dir, f"gmm_fit_{sigma_type}_{sampler_name}_chain_{selected_chain}.json")
                        if os.path.exists(json_path):  # Check file existence
                            try:
                                with open(json_path, 'r') as f:
                                    self.gmm_params[sigma_type] = json.load(f)
                            except (FileNotFoundError, json.JSONDecodeError) as e:
                                print(f"Error loading {json_path}: {e}")
                                # Decide how to handle the error, e.g., set a default value, raise, etc.
                                self.gmm_params[sigma_type] = None # or some other default value
                        else:  # File not found for this sigma_type
                            print(f"GMM fit file not found for {sigma_type}, chain {selected_chain}: {json_path}")
                            self.gmm_params[sigma_type] = None  # Optional: store None
                else:
                    print(f"No GMM fit JSON files found for {sampler_name} in {sampler_output_dir}")
                    self.gmm_params = None  # Or some other default behavior
            else:  # Directory doesn't exist
                print(f"Output directory not found: {sampler_output_dir}")
                self.gmm_params = None  # No GMM parameters available

            # Initialize self.sigma by sampling from the GMM
            self.sigma = self.initialize_sigma_from_gmm()
            self.sigma_range = None  # Not used in this case
        else:
            # Normal PairSampler initialization logic
            self.sigma, self.sigma_range = self.initialize_sigma()  # Initialize sigma values from BaseMCSampler
            self.base_priors = Priors("jeffreys")  # Initialize priors
            print("Using default sigma initialization (no GMM).")

    def initialize_sigma_from_gmm(self):
        """
        Initializes sigma values by sampling from the pre-fit GMMs.

        Returns:
            dict: A dictionary where keys are pair types (e.g., "AA", "AB") and
                values are the sampled sigma values (floats).
        """
        sigma: Dict[str, float] = {}
        for pair_type in self.params.pair_distances.keys():
            if self.gmm_params[pair_type] is not None:
                gmm_info = self.gmm_params[pair_type]
                if gmm_info is not None:
                    # Correctly sample from the GMM
                    n_components = gmm_info['n_components']
                    means = np.array(gmm_info['means']).reshape(n_components, 1)
                    covariances = np.array(gmm_info['covariances']).reshape(n_components, 1, 1)
                    weights = np.array(gmm_info['weights'])

                    # 1. Choose a component based on weights
                    component_choice = np.random.choice(n_components, p=weights)

                    # 2. Sample from the chosen component's Gaussian
                    sampled_value = np.random.multivariate_normal(
                        mean=means[component_choice].flatten(),
                        cov=covariances[component_choice]
                    )
                    sigma[pair_type] = float(sampled_value[0])  # Extract the single value

                else:
                    print(f"No GMM parameters available for {pair_type}. Using default sigma.")
                    sigma[pair_type] = self.params.pair_distances[pair_type]
            else:
                print(f"No GMM data for {pair_type}. Using default sigma.")
                sigma[pair_type] = self.params.pair_distances[pair_type]

        return sigma
    
    def _calculate_gmm_log_prob(self, sigma_value: float, pair_type: str) -> float:
        """
        Calculates the log probability density of a GMM for a *single* sigma value,
        for a specific pair type.
        """
        # If no GMM parameters or none for this pair type, return -inf or handle differently
        if self.gmm_params is None or self.gmm_params[pair_type] is None:
            return -np.inf

        gmm_info = self.gmm_params[pair_type]
        n_components = gmm_info['n_components']
        means = np.array(gmm_info['means']).reshape(n_components, 1)
        covariances = np.array(gmm_info['covariances']).reshape(n_components, 1, 1)
        weights = np.array(gmm_info['weights'])

        log_prob = -np.inf
        for i in range(n_components):
            component_log_prob = multivariate_normal.logpdf(
                sigma_value,
                mean=means[i].flatten(),
                cov=covariances[i].squeeze()
            )
            log_prob = np.logaddexp(log_prob, np.log(weights[i]) + component_log_prob)
        return log_prob

    def calculate_negative_log_prior(self, sigma: Dict[str, float]) -> float:
        """Calculates the negative log prior for a set of sigma values."""
        total_negative_log_prior = 0.0
        for pair_type, sigma_value in sigma.items():
            log_prior = self._calculate_gmm_log_prob(sigma_value, pair_type)
            total_negative_log_prior += -log_prior
        return total_negative_log_prior
    
    def get_initial_positions(self) -> Dict[str, np.ndarray]:
        """Read in the tetramer sampler trajectory file and get the positions from it"""
        cwd = os.getcwd()
        # Directory where the pairsampler trajectory files are stored
        traj_loc = os.path.join(cwd, "output_analysis/tetramersampler_results")
        # Collect chain numbers from filenames of the form "trajectory_PairSampler_chain_#.h5"
        chain_numbers = []
        for filename in os.listdir(traj_loc):
            if filename.startswith("trajectory_chain_") and filename.endswith(".h5"):
                chain_numbers.append(int(filename.split("_")[-1].split(".")[0]))
        chain_numbers.sort()
        # Choose a random chain
        cnum = random.choice(chain_numbers)
        print("selected chain: ", cnum)
        # Construct full path to file
        filename = os.path.join(traj_loc, f"trajectory_chain_{cnum}.h5")
        
        # Open the HDF5 file and pick a random frame from the last 10% of frames.
        with h5py.File(filename, 'r') as f:
            traj_grp = f['trajectory']
            keys = list(traj_grp.keys())
            keys.sort()  # Assumes keys are zero-padded (e.g., "state_00000")
            n = len(keys)
            n_equil = max(1, int(0.1 * n))  # Last 10% of frames; ensure at least one frame
            equil_keys = keys[-n_equil:]
            key = random.choice(equil_keys)
            print("selected frame: ", key)
            state_grp = traj_grp[key]
            
            # Build a dictionary of positions for parts A, B, C from the chosen frame
            positions = {}
            pos_grp = state_grp['positions']
            for type_name in pos_grp:
                positions[type_name] = pos_grp[type_name][:]
        
    
    def propose_octet_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Enhanced proposal for octet moves optimized for 8-fold symmetry systems.
        An octet consists of two tetramers (ABCC + ABCC).
        """
        # Create a deep copy of positions
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # First, get all tetramers using the existing method
        tetramers = self.ts.get_tetramers(positions, temp=1.5)  # Higher temp for more diversity
        
        # Then get octets (pairs of tetramers)
        octets = self.get_octets(positions, tetramers)
        
        if not octets:
            return new_pos
            
        # Choose a move type with probability
        move_type = np.random.choice([
            'single_octet',       # Move one octet (two tetramers together)
            'octet_pair',         # Move a pair of octets (maintaining relative orientation)
            'radial',             # Radial movement preserving symmetry
            'global_rotation',    # Rotate all octets around system center
            'aggressive'          # Larger steps for escaping local minima
        ], p=[0.96, 0.01, 0.01, 0.01, 0.01])
        
        # Dynamic step sizes based on acceptance rate
        trans_base = self.octet_trans_step
        rot_base = self.octet_rot_step
        
        # Adjust based on acceptance rate
        trans_adjust = np.clip(1.0 + 3.0 * (self.octet_trans_acc_rate - self.target_acceptance), 0.3, 3.0)
        
        if move_type == 'single_octet':
            return self._move_single_octet(new_pos, octets, trans_adjust, trans_base, rot_base)
            
        elif move_type == 'octet_pair':
            return self._move_octet_pair(new_pos, octets, trans_adjust, trans_base, rot_base)
            
        elif move_type == 'radial':
            return self._move_radial_symmetry(new_pos, octets, trans_adjust, trans_base)
            
        elif move_type == 'global_rotation':
            return self._move_global_rotation(new_pos, trans_adjust * rot_base)
            
        else:  # aggressive
            return self._move_single_octet(
                new_pos, octets, 
                trans_adjust * 2.5,  # More aggressive step size
                trans_base, rot_base
            )

    def get_octets(self, positions: Dict[str, np.ndarray], tetramers=None, temp=1.0) -> List[Tuple]:
        """
        Group tetramers into octets (pairs of ABCC tetramers).
        
        Args:
            positions: Dictionary of particle positions
            tetramers: Pre-computed tetramers (optional)
            temp: Temperature parameter for probabilistic pairing
        
        Returns:
            List of octet tuples, each containing two tetramers
        """
        if tetramers is None:
            tetramers = self.ts.get_tetramers(positions, temp=temp)
        
        if len(tetramers) < 8:  # Need at least 8 tetramers for 4 octets
            print(f"Warning: Only found {len(tetramers)} tetramers, need at least 8 for 4 octets")
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
        
        # Calculate distances between all tetramer centers (with periodic boundary conditions)
        dist_matrix = np.zeros((len(tetramers), len(tetramers)))
        for i in range(len(tetramers)):
            for j in range(i+1, len(tetramers)):
                # Calculate minimum distance considering periodic boundaries
                delta = centers[i] - centers[j]
                delta = np.where(np.abs(delta) > self.params.box_size/2,
                                delta - np.sign(delta) * self.params.box_size,
                                delta)
                dist = np.linalg.norm(delta)
                
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist
        
        # Apply temperature parameter to transform distances into probabilities
        prob_matrix = np.exp(-dist_matrix / temp)
        np.fill_diagonal(prob_matrix, 0)  # Cannot pair a tetramer with itself
        
        # Form octets by pairing tetramers
        octets = []
        available = set(range(len(tetramers)))
        
        # Try to make exactly 4 octets
        target_octets = 4
        
        while len(octets) < target_octets and len(available) >= 2:
            if len(available) == 0:
                break
                
            # Convert available set to list for indexing
            avail_list = list(available)
            
            if len(avail_list) == 2:
                # Only two tetramers left, pair them
                octets.append((tetramers[avail_list[0]], tetramers[avail_list[1]]))
                available.clear()
            else:
                # Choose first tetramer randomly from available
                idx1 = np.random.choice(avail_list)
                
                # Choose second tetramer based on probability matrix
                probs = prob_matrix[idx1, [i for i in avail_list if i != idx1]]
                if probs.sum() == 0:
                    # If all probabilities are zero, choose randomly
                    candidates = [i for i in avail_list if i != idx1]
                    idx2 = np.random.choice(candidates)
                else:
                    # Normalize probabilities
                    probs = probs / probs.sum()
                    idx2_rel = np.random.choice(len(probs), p=probs)
                    candidates = [i for i in avail_list if i != idx1]
                    idx2 = candidates[idx2_rel]
                
                # Add the pair to octets
                octets.append((tetramers[idx1], tetramers[idx2]))
                available.remove(idx1)
                available.remove(idx2)
        
        return octets

    def _move_single_octet(self, new_pos, octets, trans_adjust, trans_base, rot_base):
        """Move a single octet (a pair of tetramers) with translation and rotation."""
        # Select a random octet
        octet_idx = np.random.randint(len(octets))
        tetramer1, tetramer2 = octets[octet_idx]
        
        # Collect all particles in the octet
        octet_particles = []
        for a_idx, b_idx, c_idx1, c_idx2 in [tetramer1, tetramer2]:
            octet_particles.extend([
                ('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)
            ])
        
        # Collect octet coordinates for centroid calculation
        octet_coords = [new_pos[p][i] for p, i in octet_particles]
        centroid = np.mean(octet_coords, axis=0)
        
        # Apply translation with high probability
        if np.random.random() < 0.95:
            # Translation with jitter for exploration
            trans_step = trans_base * trans_adjust * np.random.uniform(0.7, 1.3)
            displacement = np.random.normal(0, trans_step, 3)
            
            # Apply translation to all particles in the octet
            for part, idx in octet_particles:
                new_pos[part][idx] += displacement
                new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        # Apply rotation with high probability
        if np.random.random() < 0.95:
            # Enhanced rotation with variable amplitude
            rot_step = rot_base * trans_adjust * np.random.uniform(0.7, 1.3)
            rotation_axis = self._random_unit_vector()
            rotation_angle = np.random.normal(0, rot_step)
            rot_matrix = self._rotation_matrix(rotation_axis, rotation_angle)
            
            # Apply rotation to all particles in the octet
            for part, idx in octet_particles:
                vec = new_pos[part][idx] - centroid
                new_pos[part][idx] = centroid + rot_matrix @ vec
                new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        return new_pos

    def _move_octet_pair(self, new_pos, octets, trans_adjust, trans_base, rot_base):
        """Move a pair of octets while maintaining their relative orientation."""
        if len(octets) < 2:
            return self._move_single_octet(new_pos, octets, trans_adjust, trans_base, rot_base)
        
        # Select two octets
        indices = np.random.choice(len(octets), size=2, replace=False)
        
        # Collect all particles in both octets
        all_particles = []
        for idx in indices:
            tetramer1, tetramer2 = octets[idx]
            for a_idx, b_idx, c_idx1, c_idx2 in [tetramer1, tetramer2]:
                all_particles.extend([
                    ('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)
                ])
        
        # Apply translation to both octets together
        trans_step = trans_base * trans_adjust
        displacement = np.random.normal(0, trans_step, 3)
        
        for part, idx in all_particles:
            new_pos[part][idx] += displacement
            new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        return new_pos

    def _move_radial_symmetry(self, new_pos, octets, trans_adjust, trans_base):
        """Move octets preserving radial symmetry."""
        if len(octets) < 2:
            return self._move_single_octet(new_pos, octets, trans_adjust, trans_base, self.octet_rot_step)
        
        # Calculate system center
        all_coords = []
        for type_name, positions_array in new_pos.items():
            all_coords.extend(positions_array)
        system_center = np.mean(all_coords, axis=0)
        
        # Randomly adjust the radial distance for all octets
        radial_adjust = np.random.normal(0, trans_base * trans_adjust)
        
        # Apply radial adjustment to all octets
        for tetramer1, tetramer2 in octets:
            octet_particles = []
            for a_idx, b_idx, c_idx1, c_idx2 in [tetramer1, tetramer2]:
                octet_particles.extend([
                    ('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)
                ])
            
            # Apply radial scaling to all particles in this octet
            for part, idx in octet_particles:
                vec = new_pos[part][idx] - system_center
                distance = np.linalg.norm(vec)
                if distance > 0:  # Avoid division by zero
                    new_distance = max(0.1, distance + radial_adjust)
                    scaling = new_distance / distance
                    new_pos[part][idx] = system_center + vec * scaling
                    new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        return new_pos

    def _move_global_rotation(self, new_pos, rot_step):
        """Apply a global rotation to all particles around the system center."""
        # Calculate system center
        all_coords = []
        for type_name, positions_array in new_pos.items():
            for pos in positions_array:
                all_coords.append(pos)
        system_center = np.mean(all_coords, axis=0)
        
        # Generate rotation
        rotation_axis = self._random_unit_vector()
        rotation_angle = np.random.normal(0, rot_step)
        rotation_matrix = self._rotation_matrix(rotation_axis, rotation_angle)
        
        # Apply to all particles
        for type_name in new_pos:
            for i in range(len(new_pos[type_name])):
                vec = new_pos[type_name][i] - system_center
                new_pos[type_name][i] = system_center + rotation_matrix @ vec
                new_pos[type_name][i] = np.mod(new_pos[type_name][i], self.params.box_size)
        
        return new_pos
    
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
        
    def run_mc(self, n_steps: int = 50000, save_freq: int = 100, output_dir: str = "output_analysis/octamersampler_results/") -> Tuple:
        """Monte Carlo sampling with octet moves."""
        best_positions = None
        trajectory = []
        sigma_history_octet = {key: [] for key in self.sigma}  # Pre-allocate sigma history

        # Calculate the prior for the initial sigma values
        if self.use_sigma_distribution:
            prior_penalty = self.calculate_negative_log_prior(self.sigma)
        else:
            prior_penalty = self.base_priors.neg_log_prior(self.sigma, self.sigma_range)

        # Get initial octets (pairs of tetramers)
        tetramers = self.ts.get_tetramers(self.positions)
        octets = self.get_octets(self.positions, tetramers)

        # Initial score calculation
        current_score, curr_ex, curr_pair, curr_oct = self.neg_log_posterior(
            self.positions,
            tetramers,
            octets,
            prior_penalty, 
            self.sigma
        )
        curr_prior = prior_penalty
        best_score = current_score  # Initialize with the initial score

        # Pre-calculate temperature schedule
        initial_temp = 5.0
        final_temp = 1.0
        cooling_factor = -np.log(final_temp / initial_temp) / n_steps
        temperatures = initial_temp * np.exp(-cooling_factor * np.arange(n_steps))

        # Prepare output directories and files
        os.makedirs(output_dir, exist_ok=True)
        csv_log_file = os.path.join(output_dir, "all_info_mcmc_octet.csv")

        # Write CSV header (include individual score components)
        with open(csv_log_file, "w") as f:
            f.write("Step,T,Prior,Exvol_score,Pair_score,Octet_score,Score,Accepted\n")
        
        oct_moves = 0
        oct_accepted = 0  # Counter for octet moves
        accepted_moves = 0
        
        for step in range(n_steps):
            temp = temperatures[step]
            move_type = np.random.choice(['position', 'sigma', 'tetramer', 'octet'],
                                        p=[0.2, 0.1, 0.3, 0.4])

            if move_type == 'position':
                proposed_positions = self.propose_position_move(self.positions)
                proposed_sigma = self.sigma
            elif move_type == 'sigma':
                proposed_positions = self.positions
                proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)
            elif move_type == 'tetramer':
                proposed_positions = self.ts.propose_tetramer_move(self.positions)
                proposed_sigma = self.sigma
            else:  # octet move
                oct_moves += 1
                self.oct_trans_acc_rate = oct_accepted / max(1, oct_moves)
                proposed_positions = self.propose_octet_move(self.positions)
                proposed_sigma = self.sigma
            
            # Get tetramers and octets for proposed configuration
            proposed_tetramers = self.ts.get_tetramers(proposed_positions)
            proposed_octets = self.get_octets(proposed_positions, proposed_tetramers)
            
            new_prior_penalty = 0.0
            # Calculate the prior for the proposed sigma values
            if self.use_sigma_distribution:
                new_prior_penalty = self.calculate_negative_log_prior(proposed_sigma)
            else:
                new_prior_penalty = self.base_priors.neg_log_prior(proposed_sigma, self.sigma_range)
                
            # Calculate proposed score    
            proposed_score, prop_ex, prop_pair, prop_oct = self.neg_log_posterior(
                proposed_positions,
                proposed_tetramers,
                proposed_octets,
                new_prior_penalty, 
                proposed_sigma
            )
            
            # Metropolis criterion
            delta_e = proposed_score - current_score
            acceptance = 0
            if delta_e < 0 or np.random.random() < np.exp(-delta_e / temp):
                self.positions = proposed_positions
                self.sigma = proposed_sigma
                current_score = proposed_score
                curr_ex = prop_ex
                curr_pair = prop_pair
                curr_oct = prop_oct
                curr_prior = new_prior_penalty
                accepted_moves += 1
                acceptance = 1

                # Update octet acceptance rate
                if move_type == 'octet':
                    oct_accepted += 1

                # Keep track of best positions
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions.items()}
            
            # Logging and saving trajectory at specified frequency
            if step % save_freq == 0:
                accept_rate = accepted_moves / (step + 1)
                # Update sigma history
                for key in sigma_history_octet:
                    sigma_history_octet[key].append(self.sigma[key])
                
                # Append current state with detailed score breakdown
                trajectory.append(
                    self.save_state(
                        step, self.positions, self.sigma, current_score,
                        prior_score=curr_prior,
                        pair_score=curr_pair,
                        exvol_score=curr_ex,
                        oct_score=curr_oct
                    )
                )
                
                # Write to CSV log
                with open(csv_log_file, "a") as f:
                    f.write(
                        f"{step},{temp:.1f},{curr_prior:.1f},{curr_ex:.1f},{curr_pair:.1f},"
                        f"{curr_oct:.1f},{current_score:.1f},{acceptance:.1f}\n"
                    )

                print(f"Step {step}, Score: {current_score:.2f}, "
                    f"Temp: {temp:.2f}, Accept: {accept_rate:.2f}")

        # Save sigma history
        sigma_history_oct_df = pd.DataFrame(sigma_history_octet)
        sigma_history_oct_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        # Save the trajectory in HDF5 format
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        final_file = self.save_trajectory(trajectory, trajectory_file)

        return best_positions, trajectory, final_file

    def calculate_octet_score(self, positions: Dict[str, np.ndarray], 
                            octet: Tuple[Tuple[int, ...], Tuple[int, ...]], 
                            sig: Dict[str, float] = None) -> float:
        """
        Calculate score for a single octet (a pair of tetramers).
        An octet consists of 2 A's, 2 B's, and 4 C's.
        
        Args:
            positions: Dictionary mapping atom types to position arrays
            octet: A tuple containing two tetramer tuples
            sig: Dictionary of sigma values
        
        Returns:
            float: The negative log likelihood score for the octet
        """
        tetramer1, tetramer2 = octet
        
        # First calculate individual tetramer scores
        score = self.ts.calculate_tetramer_score(positions, tetramer1, sig)
        score += self.ts.calculate_tetramer_score(positions, tetramer2, sig)
        
        # Now calculate inter-tetramer interactions within the octet
        a_idx1, b_idx1, c1_idx1, c2_idx1 = tetramer1
        a_idx2, b_idx2, c1_idx2, c2_idx2 = tetramer2
        
        # A-A interaction between tetramers
        score += self.pair_score_nll(
            positions['A'][a_idx1],
            positions['A'][a_idx2],
            self.params.pair_distances['AA'],
            sig['AA']
        )
        
        # B-B interaction between tetramers
#        score += self.pair_score_nll(
#            positions['B'][b_idx1],
#            positions['B'][b_idx2],
#            self.params.pair_distances['BB'],
#            sig['BB']
#        )
        
        # Cross-tetramer C-C interactions (can be customized based on your model)
        for c_idx1 in [c1_idx1, c2_idx1]:
            for c_idx2 in [c1_idx2, c2_idx2]:
                score += 0.5 * self.pair_score_nll(  # Weight factor of 0.5 to prevent overweighting
                    positions['C'][c_idx1],
                    positions['C'][c_idx2],
                    self.params.pair_distances['CC'],
                    sig['CC']
                )
        
        return score

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        tetramers: List[Tuple[int, ...]],
        octets: List[Tuple[Tuple[int, ...], Tuple[int, ...]]],
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        octet_weight: float = 1.0,
    ) -> Tuple[float, float, float, float]:
        """
        Calculate the total score for an octet system:
        1. Identify pairs within all tetramers
        2. Identify additional pairs that form octets
        3. Exclude those pairs from the parent PairSampler calculation
        4. Add additional terms for octet-specific scoring
        
        Returns:
            Tuple of (total_score, exclusion_score, pair_score, octet_score)
        """
        # 1) Identify pairs within tetramers
        tetramer_pairs = set()
        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
            tetramer_pairs.add(('A', a_idx, 'B', b_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c1_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c2_idx))
            tetramer_pairs.add(('C', c1_idx, 'C', c2_idx))
        
        # 2) Identify additional octet-specific pairs
        octet_pairs = set()
        for tetramer1, tetramer2 in octets:
            a_idx1, b_idx1, c1_idx1, c2_idx1 = tetramer1
            a_idx2, b_idx2, c1_idx2, c2_idx2 = tetramer2
            
            # Add A-A interaction between tetramers
            octet_pairs.add(('A', a_idx1, 'A', a_idx2))
            
#            # Add B-B interaction between tetramers
#            octet_pairs.add(('B', b_idx1, 'B', b_idx2))
            
            # Add cross-tetramer C-C interactions
            for c_idx1 in [c1_idx1, c2_idx1]:
                for c_idx2 in [c1_idx2, c2_idx2]:
                    octet_pairs.add(('C', c_idx1, 'C', c_idx2))
        
        # 3) Calculate score using the parent PairSampler, excluding all internal pairs
        excluded_pairs = tetramer_pairs.union(octet_pairs)
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions,
            sig,
            self.sigma_range,
            excluded_pairs,
            self.use_sigma_distribution,
            prior_penalty_from_distribution
        )

        # 4) Add octet-specific score
        total_oct_score = 0.0
        
        # Score for each octet
        for octet in octets:
            oct_score = self.calculate_octet_score(positions, octet, sig)
            total_oct_score += oct_score
        
        # Add the octet score to the total
        score += octet_weight * total_oct_score

        return score, ex_score, pair_score, octet_weight * total_oct_score