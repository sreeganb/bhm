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
from typing import List, Dict, Tuple, Optional
from base_sampler import BaseMCSampler, Priors  # Import BaseMCSampler and related classes
from parameters import SystemParameters
from visualization import visualize_3d_configuration
from scipy.spatial.distance import cdist
from typing import Dict, List, Tuple
from scipy.stats import multivariate_normal
from pair_sampler import PairSampler
#-----------------------------------------------------------------------
class TetramerSampler(BaseMCSampler):
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

        # Additional initialization for tetramer-specific features
        self.tetramer_trans_step = 0.4
        self.tetramer_rot_step = 0.1
        self.sigma_prior_dist = {}
        #-------------------------------------
        # tracking acceptance rates
        self.tet_trans_acc_rate = 0.5
        #self.tet_rot_acc_rate = 0.5
        self.target_acceptance = 0.5
        #-------------------------------------

        required_pairs = {"AA", "AB", "BC", "CC"}
        if not required_pairs.issubset(self.params.pair_distances.keys()):
            raise ValueError("Missing required pair types for tetramer sampling.")

        if use_sigma_distribution:
            print("using sigma distribution now: ")
            pwd = os.getcwd()
            sampler_output_dir = os.path.join(pwd, "output_analysis/pairsampler_results") # hard coded here, needs to be changed later
            sampler_name = "PairSampler" # hard coded here, needs to be changed later
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

    def get_tetramers(self, positions: Dict[str, np.ndarray], temp=1.0) -> List[Tuple[int, ...]]:
        """Efficient and robust tetramer selection. For this systems returns 8 tetramers."""
        
        a_positions = positions['A']
        b_positions = positions['B']
        c_positions = positions['C']
        
        n_a = len(a_positions)
        n_b = len(b_positions)
        n_c = len(c_positions)
        
        tetramers = []
        
        dist_AB = cdist(a_positions, b_positions) / temp
        probs_B = np.exp(-dist_AB)
        probs_B = np.nan_to_num(probs_B, nan=1 / n_b)
        probs_B /= probs_B.sum(axis=1, keepdims=True)
        
        b_indices = np.array([np.random.choice(n_b, p=probs_B[i]) for i in range(n_a)])
        
        dist_BC = cdist(b_positions[b_indices], c_positions) / temp
        probs_C = np.exp(-dist_BC)
        probs_C = np.nan_to_num(probs_C, nan=1 / n_c)
        
        for a_idx, b_idx in enumerate(b_indices):
            if probs_C[a_idx].sum() == 0 or n_c < 2:
                c_indices = np.random.choice(n_c, size=2, replace=(n_c < 2))
            else:
                probs_C_normalized = probs_C[a_idx] / probs_C[a_idx].sum()
                c_indices = np.random.choice(n_c, size=2, replace=False, p=probs_C_normalized)
            
            tetramers.append((a_idx, b_idx, c_indices[0], c_indices[1]))
        
        return tetramers
#-----------------------------------------------------------------------  
    def run_mc(self, n_steps: int = 50000, save_freq: int = 100) -> Tuple:
        """Monte Carlo sampling with tetramer moves."""
        best_positions = None
        trajectory = []
        sigma_history_tetramer = {key: [] for key in self.sigma}  # Pre-allocate sigma history

        # Calculate the prior for the initial sigma values
        if self.use_sigma_distribution:
            prior_penalty = self.calculate_negative_log_prior(self.sigma)
        else:
            prior_penalty = self.base_priors.neg_log_prior(self.sigma, self.sigma_range)

        # Initial score calculation
        current_score, curr_ex, curr_pair, curr_tet  = self.neg_log_posterior(self.positions,
                                            self.get_tetramers(self.positions),
                                            prior_penalty, self.sigma)
        curr_prior = prior_penalty
        best_score = current_score  # Initialize with the initial score

        # Pre-calculate temperature schedule
        initial_temp = 5.0
        final_temp = 1.0
        cooling_factor = -np.log(final_temp / initial_temp) / n_steps
        temperatures = initial_temp * np.exp(-cooling_factor * np.arange(n_steps))

        # Prepare output directories and files
        output_dir = os.path.join("output_analysis", "tetramersampler_results")
        os.makedirs(output_dir, exist_ok=True)
        csv_log_file = os.path.join(output_dir, "all_info_mcmc_tetramer.csv")

        # Write CSV header (include individual score components)
        with open(csv_log_file, "w") as f:
            f.write("Step,T,Prior,Exvol_score,Pair_score,Score,Accepted\n")
        
        tet_moves = 0
        tet_accepted = 0  # Single counter for tetramer moves (combined trans+rot)
        accepted_moves = 0
        
        for step in range(n_steps):
            temp = temperatures[step]
            move_type = np.random.choice(['position', 'sigma', 'tetramer'],
                                        p=[0.4, 0.1, 0.5])

            if move_type == 'position':
                proposed_positions = self.propose_position_move(self.positions)
                proposed_sigma = self.sigma
            elif move_type == 'sigma':
                proposed_positions = self.positions
                proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)
            else:  # tetramer move
                tet_moves += 1
                self.tet_trans_acc_rate = tet_accepted / max(1, tet_moves)
                proposed_positions = self.propose_tetramer_move(self.positions)
                proposed_sigma = self.sigma
            
            new_prior_penalty = 0.0
            # Calculate the prior for the initial sigma values
            if self.use_sigma_distribution:
                new_prior_penalty = self.calculate_negative_log_prior(proposed_sigma)
            else:
                new_prior_penalty = self.base_priors.neg_log_prior(proposed_sigma, self.sigma_range)
            # Calculate proposed score    
            proposed_score, prop_ex, prop_pair, prop_tet = self.neg_log_posterior(proposed_positions,
                                            self.get_tetramers(proposed_positions),
                                            new_prior_penalty, proposed_sigma)
            
            delta_e = proposed_score - current_score
            acceptance = 0
            if delta_e < 0 or np.random.random() < np.exp(-delta_e / temp):
                self.positions = proposed_positions
                self.sigma = proposed_sigma
                current_score = proposed_score
                curr_ex = prop_ex
                curr_pair = prop_pair
                curr_tet = prop_tet
                curr_prior = new_prior_penalty
                accepted_moves += 1
                acceptance = 1

                if move_type == 'tetramer':
                    tet_accepted += 1  # Increment tetramer acceptance counter

                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions.items()}
            
            # Batch logging and saving trajectory
            if step % save_freq == 0:
                accept_rate = accepted_moves / (step + 1)
                # update sigma history
                for key in sigma_history_tetramer:
                    sigma_history_tetramer[key].append(self.sigma[key])
                
                # Append current state with detailed score breakdown.
                trajectory.append(
                    self.save_state(
                        step, self.positions, self.sigma, current_score,
                        prior_score=curr_prior,
                        pair_score=curr_pair,
                        exvol_score=curr_ex,
                        tet_score=curr_tet
                    )
                )
                
                with open(csv_log_file, "a") as f:
                    f.write(
                        f"{step},{curr_prior:.1f},{curr_ex:.1f},{curr_pair:.1f},"
                        f"{curr_tet:.1f},{current_score:.1f},{acceptance:.1f}\n"
                    )

                print(f"Step {step}, Score: {current_score:.1f}, "
                      f"Temp: {temp:.1f}, Accept: {accept_rate:.1f}")

        sigma_history_tet_df = pd.DataFrame(sigma_history_tetramer)
        sigma_history_tet_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        # Save the trajectory in HDF5 format in the output folder.
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        final_file = self.save_trajectory(trajectory, trajectory_file)

        return best_positions, trajectory, final_file
#-----------------------------------------------------------------------
    
    def propose_tetramer_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Enhanced proposal for tetramer moves optimized for 8-fold symmetry systems.
        Implements several types of moves with varying aggressiveness to improve mixing.
        """
        # Create a deep copy of positions
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # Get tetramers with temperature parameter that encourages exploration
        # Higher temperature = more uniform selection, less dependent on current distances
        exploration_temp = 1.0  # Increased from 1.0 for better exploration
        tetramers = self.get_tetramers(positions, temp=exploration_temp)
        
        if not tetramers:
            return new_pos
            
        # Choose a move type with probability
        move_type = np.random.choice([
            'single_tetramer',     # Move one tetramer
            'coordinated_pair',    # Move a pair of tetramers (maintaining relative orientation)
            'radial',              # Radial movement preserving symmetry
            'global_rotation',     # Rotate all tetramers around system center
            'aggressive'           # Larger steps for escaping local minima
        ], p=[0.99, 0.0025, 0.0025, 0.0025, 0.0025])
        
        # Dynamic step sizes based on acceptance rate
        # More sophisticated adaptation - separate translation and rotation
        trans_base = self.tetramer_trans_step
        rot_base = self.tetramer_rot_step
        
        # Adjust based on acceptance rate, with wider range (0.3-3.0)
        trans_adjust = np.clip(1.0 + 5.0 * (self.tet_trans_acc_rate - self.target_acceptance), 0.3, 2.0)
        
        if move_type == 'single_tetramer':
            # Standard single tetramer move (enhanced)
            return self._move_single_tetramer(new_pos, tetramers, trans_adjust, trans_base, rot_base)
            
        elif move_type == 'coordinated_pair':
            # Move a pair of tetramers while maintaining their relative orientation
            return self._move_tetramer_pair(new_pos, tetramers, trans_adjust, trans_base, rot_base)
            
        elif move_type == 'radial':
            # Move preserving radial symmetry
            return self._move_radial_symmetry(new_pos, tetramers, trans_adjust, trans_base)
            
        elif move_type == 'global_rotation':
            # Global rotation around system center
            return self._move_global_rotation(new_pos, tetramers, rot_base * trans_adjust)
            
        else:  # aggressive
            # Aggressive move with larger steps - helps escape local minima
            return self._move_single_tetramer(
                new_pos, tetramers, 
                trans_adjust * 2.0,  # Double the step size
                trans_base, rot_base
            )

    def _move_single_tetramer(self, new_pos, tetramers, trans_adjust, trans_base, rot_base):
        """Move a single tetramer with translation and rotation."""
        # Select a random tetramer
        tetramer_idx = np.random.randint(len(tetramers))
        a_idx, b_idx, c_idx1, c_idx2 = tetramers[tetramer_idx]
        c_indices = [c_idx1, c_idx2]
        
        # Independent translation and rotation with probability
        if np.random.random() < 0.99:  # 80% chance of translation
            # More efficient translation step
            trans_step = trans_base * trans_adjust * np.random.uniform(0.7, 1.3)  # Add jitter
            displacement = np.random.normal(0, trans_step, 3)
            
            # Apply translation
            for part, idx in [('A', a_idx), ('B', b_idx)] + [('C', c_idx) for c_idx in c_indices]:
                new_pos[part][idx] += displacement
                # Improved boundary handling with periodic boundary conditions
                new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        if np.random.random() < 0.99:  # 80% chance of rotation
            # Collect tetramer coordinates
            tetramer_coords = [
                new_pos['A'][a_idx], 
                new_pos['B'][b_idx], 
                new_pos['C'][c_idx1], 
                new_pos['C'][c_idx2]
            ]
            centroid = np.mean(tetramer_coords, axis=0)
            
            # Enhanced rotation with variable amplitude
            rot_step = rot_base * trans_adjust * np.random.uniform(0.7, 1.3)  # Add jitter
            rotation_axis = self._random_unit_vector()
            rotation_angle = np.random.normal(0, rot_step)
            rot_matrix = self._rotation_matrix(rotation_axis, rotation_angle)
            
            # Apply rotation
            for part, idx in [('A', a_idx), ('B', b_idx)] + [('C', c_idx) for c_idx in c_indices]:
                vec = new_pos[part][idx] - centroid
                new_pos[part][idx] = centroid + rot_matrix @ vec
                # Apply boundary conditions
                new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        return new_pos

    def _move_tetramer_pair(self, new_pos, tetramers, trans_adjust, trans_base, rot_base):
        """Move a pair of tetramers while maintaining their relative orientation."""
        if len(tetramers) < 2:
            return self._move_single_tetramer(new_pos, tetramers, trans_adjust, trans_base, rot_base)
        
        # Select two tetramers
        indices = np.random.choice(len(tetramers), size=2, replace=False)
        
        # Calculate centers of the tetramers
        centers = []
        tetramer_particles = []
        
        for idx in indices:
            a_idx, b_idx, c_idx1, c_idx2 = tetramers[idx]
            particles = [
                ('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)
            ]
            tetramer_particles.append(particles)
            
            coords = np.array([new_pos[p][i] for p, i in particles])
            centers.append(np.mean(coords, axis=0))
        
        # Move both tetramers together
        trans_step = trans_base * trans_adjust
        displacement = np.random.normal(0, trans_step, 3)
        
        for particles in tetramer_particles:
            for part, idx in particles:
                new_pos[part][idx] += displacement
                new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        return new_pos

    def _move_radial_symmetry(self, new_pos, tetramers, trans_adjust, trans_base):
        """Move tetramers preserving radial symmetry."""
        if len(tetramers) < 3:  # Need multiple tetramers for symmetry moves
            return self._move_single_tetramer(new_pos, tetramers, trans_adjust, trans_base, self.tetramer_rot_step)
        
        # Calculate system center
        all_coords = []
        for type_name, positions_array in new_pos.items():
            all_coords.extend(positions_array)
        system_center = np.mean(all_coords, axis=0)
        
        # Randomly adjust the radial distance for all tetramers
        radial_adjust = np.random.normal(0, trans_base * trans_adjust)
        
        for tetramer in tetramers:
            a_idx, b_idx, c_idx1, c_idx2 = tetramer
            particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]
            
            for part, idx in particles:
                vec = new_pos[part][idx] - system_center
                distance = np.linalg.norm(vec)
                if distance > 0:  # Avoid division by zero
                    new_distance = max(0.1, distance + radial_adjust)
                    scaling = new_distance / distance
                    new_pos[part][idx] = system_center + vec * scaling
                    new_pos[part][idx] = np.mod(new_pos[part][idx], self.params.box_size)
        
        return new_pos

    def _move_global_rotation(self, new_pos, tetramers, rot_step):
        """Apply a global rotation to all tetramers around the system center."""
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
        
    def calculate_tetramer_score(self, positions: Dict[str, np.ndarray], 
                                tetramer: Tuple[int, ...], sig: Dict[str,float] = None) -> float:
        """Calculate score for a single tetramer."""
        a_idx, b_idx, c1_idx, c2_idx = tetramer
        score = 0.0
        score += self.pair_score_nll(
            positions['A'][a_idx],
            positions['B'][b_idx],
            self.params.pair_distances['AB'],
            sig['AB']
        )
        for c_idx in [c1_idx, c2_idx]:
            score += self.pair_score_nll(
                positions['B'][b_idx],
                positions['C'][c_idx],
                self.params.pair_distances['BC'],
                sig['BC']
            )
        score += self.pair_score_nll(
            positions['C'][c1_idx],
            positions['C'][c2_idx],
            self.params.pair_distances['CC'],
            sig['CC']
        )
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
    ) -> float:
        """
        Calculate the total score for a tetramer system:
        1. Identify the A-B-C-C tetramers and collect the pairs inside each tetramer.
        2. Exclude those pairs from the parent PairSampler calculation.
        3. Add an additional term for tetramer-specific scoring.
        """
        # 1) Identify pairs within tetramers
        #tetramers = self.get_tetramers(positions)  # Method that finds tuples like (a_idx, b_idx, c1_idx, c2_idx)
        tetramer_pairs = set()
        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
            tetramer_pairs.add(('A', a_idx, 'B', b_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c1_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c2_idx))
            tetramer_pairs.add(('C', c1_idx, 'C', c2_idx))

        # 2) Calculate score using the parent PairSampler, excluding tetramer-internal pairs
        #    Note the named argument `sigma=self.sigma` to ensure the parent sees a dictionary for sigma.
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions,
            sig,  # Pass the dictionary so parent can do sigma[pair_key]
            self.sigma_range,
            tetramer_pairs,
            self.use_sigma_distribution,
            prior_penalty_from_distribution
        )

        # 3) Add tetramer-specific score
        total_tet_score = 0.0
        for tetramer in tetramers:
            tet_score = self.calculate_tetramer_score(positions, tetramer, sig)
            total_tet_score += tet_score
        score += tetramer_weight * total_tet_score

        return score, ex_score, pair_score, tetramer_weight * total_tet_score
#-----------------------------------------------------------------------