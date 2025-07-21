import numpy as np
from typing import Dict, Tuple
import copy
from scipy.spatial.distance import cdist
import random
from parameters import SystemParameters
import h5py
import os
#---------------------------------------------------------------------------
# Modify Priors class to have a stronger uniform option
class Priors:
    """Calculate prior penalties for parameters."""
    
    def __init__(self, prior_type="uniform"):
        self.prior_type = prior_type
        
    def neg_log_prior(self, sigma: Dict[str, float], sigma_range: Dict[str, Tuple[float, float]]) -> float:
        """Calculate negative log prior for all sigma values."""
        penalty = 0.0
        
        for key, value in sigma.items():
            if key in sigma_range:
                min_val, max_val = sigma_range[key]
                
                if self.prior_type == "uniform":
                    # Hard uniform prior
                    if value < min_val or value > max_val:
                        penalty += 1e6  # Extreme penalty for out-of-bounds
                    # Add small penalty as sigma approaches bounds
                    elif value < min_val * 1.1:
                        penalty += 100 * (min_val * 1.1 - value) / min_val
                    elif value > max_val * 0.9:
                        penalty += 100 * (value - max_val * 0.9) / max_val
                elif self.prior_type == "jeffreys":
                    # Jeffreys prior (proportional to 1/sigma)
                    if value < min_val or value > max_val:
                        penalty += 1e6
                    else:
                        penalty += np.log(value)
                elif self.prior_type == "inverse_gamma":
                    # Inverse gamma prior (shape and scale parameters)
                    beta = 0.001  # scale parameter
                    alpha = 0.001  # shape parameter
                    if value < min_val or value > max_val:
                        penalty += 1e6
                    else:
                        penalty += -1.0 * (alpha * np.log(beta) - np.log(np.math.gamma(alpha))
                                          - (alpha + 1) * np.log(value) - beta / value)
                else:
                    raise ValueError(f"Unknown prior type: {self.prior_type}")
                        
        return penalty
#---------------------------------------------------------------------------
class BaseMCSampler:
    def __init__(self):
        self.params = SystemParameters()  # Initialize system parameters
        self.target_acceptance = 0.5  # Target acceptance rate for moves

    def initialize_sigma(self):
        """Initialize sigma values and priors based on pair distances and radii."""
        # Dictionaries for sigma parameters, ranges, and priors.
        sigma: Dict[str, float] = {}
        sigma_range: Dict[str, Tuple[float, float]] = {}

        for pair_type in self.params.pair_distances.keys():
            # Cache the sum of radii
            sum_radii = self.params.radii[pair_type[0]] + self.params.radii[pair_type[1]]
            lower_bound = 0.03 * sum_radii
            upper_bound = 0.4 * sum_radii
            # Propose sigma in log-space for a uniform proposal in that space.
            sigma_val = np.exp(np.random.uniform(np.log(lower_bound), np.log(upper_bound)))
            sigma[pair_type] = sigma_val
            sigma_range[pair_type] = (lower_bound, upper_bound)
        
        # return sigma, range
        return sigma, sigma_range 

    def initialize_positions(self, is_ideal: int = 0) -> Dict[str, np.ndarray]:
        """Initialize positions either as ideal or randomly within the box."""
        positions = {}
        if is_ideal:
            for type_name, _ in self.params.component_counts.items():
                # Use deep copy for ideal coordinates to avoid accidental modification.
                positions[type_name] = copy.deepcopy(self.params.ideal_coordinates[type_name])
        else:
            for type_name, count in self.params.component_counts.items():
                positions[type_name] = np.random.uniform(0, self.params.box_size, (count, 3))

        return positions

    def excluded_volume_nll(self, positions: Dict[str, np.ndarray]) -> float:
        """
        Compute the negative log likelihood of the excluded volume penalty.
        """
        score = 0.0
        exvol_sigma = 0.1  # Excluded volume sigma penalty strength
        
        # Process each pair of particle types
        for type1, pos1 in positions.items():
            for type2, pos2 in positions.items():
                # Only compute once for each unordered pair to avoid double-counting
                if type1 <= type2:
                    min_dist = self.params.radii[type1] + self.params.radii[type2]
                    distances = cdist(pos1, pos2)
                    
                    # Handle same-type particles (avoid self-interactions)
                    if type1 == type2:
                        # Use upper triangle mask to avoid counting same particle twice
                        mask = np.triu(np.ones_like(distances), k=1)
                        viol_mask = (distances < min_dist) & (mask > 0)
                    else:
                        viol_mask = distances < min_dist
                    
                    # Vectorized calculation of penalties
                    if np.any(viol_mask):
                        overlaps = min_dist - distances[viol_mask]
                        score += np.sum((overlaps**2) / (exvol_sigma**2))
                    
        return score

    def pair_score_nll(self, pos1: np.ndarray, pos2: np.ndarray, target_dist: float, sigma: float) -> float:
        """Compute the pairwise negative log likelihood score."""
        distance = np.linalg.norm(pos1 - pos2)
        return ((distance - target_dist) ** 2) / (2 * sigma**2) + np.log(2 * np.pi * sigma**2)

    def calculate_pair_scores_matrix(self, pos1: np.ndarray, pos2: np.ndarray, target_dist: float, sigma: float) -> np.ndarray:
        """Vectorized calculation of pair scores between two sets of positions."""
        distances = cdist(pos1, pos2)
        return ((distances - target_dist) ** 2) / (2 * sigma**2) + np.log(2 * np.pi * sigma**2)
    #---------------------------------------------------------------------------    
#    def propose_sigma_move(self, sigma: Dict[str, float], accept_rate: float = 0.5) -> Tuple[Dict[str, float], str]:
#        """Propose a move for one sigma parameter using improved log-space steps."""
#        import random
#        
#        # Randomly select which sigma parameter to update
#        pair_type = random.choice(list(sigma.keys()))
#        
#        # Current value in log-space
#        log_current = np.log(sigma[pair_type])
#        
#        log_step_size = 0.2  # ~40x larger than your current value
#        
#        # Symmetric proposal: add zero-mean Gaussian noise in log-space
#        log_proposed = log_current + np.random.normal(0, log_step_size)
#        
#        # Remove clipping to allow full exploration (your prior should handle extreme values)
#        # If needed, use less restrictive bounds: np.clip(log_proposed, -15, 15)
#        
#        # Create new sigma dictionary with proposed value
#        new_sigma = dict(sigma)
#        new_sigma[pair_type] = np.exp(log_proposed)
#        
#        return new_sigma, pair_type
    def propose_sigma_move(self, sigma: Dict[str, float], accept_rate: float = 0.5) -> Tuple[Dict[str, float], str]:
        """Propose a move for one sigma parameter using small direct-space steps."""
        import random
        
        # Randomly select which sigma parameter to update
        pair_type = random.choice(list(sigma.keys()))
        
        # Current value
        current_value = sigma[pair_type]
        
        # Step size as percentage of current value (1% here)
        relative_step_size = 0.005 # changed from 0.01 to 0.005
        step_size = relative_step_size * current_value
        
        # Symmetric proposal: add zero-mean Gaussian noise directly to value
        proposed_value = current_value + np.random.normal(0, step_size)
        
        # Ensure positivity
        proposed_value = max(1e-6, proposed_value)
        
        # Clipping to avoid extreme values
        proposed_value = np.clip(proposed_value, 0.0, 20.0)  # Adjust bounds as needed
        
        # Create new sigma dictionary with proposed value
        new_sigma = dict(sigma)
        new_sigma[pair_type] = proposed_value
        
        return new_sigma, pair_type
    #---------------------------------------------------------------------------
    def propose_position_move(self, positions: Dict[str, np.ndarray], accept_rate: float = 0.5) -> Dict[str, np.ndarray]:
        """
        MCMC position proposal with inverse radius scaling: larger particles move less.
        Uses simple clipping to stay within box boundaries.
        """
        import random
        
        # Deep copy to avoid modifying original
        new_positions = {key: np.copy(array) for key, array in positions.items()}
        
        # Randomly select particle type and index
        type_names = list(self.params.component_counts.keys())
        type_name = random.choice(type_names)
        idx = np.random.randint(self.params.component_counts[type_name])

        # Get min radius across all particle types for scaling
        min_radius = min(self.params.radii.values())

        # Inverse scaling: smaller radius = larger moves
        # Base step is normalized to the smallest particle
        base_step = 2.5  # Adjust this value as needed for optimal acceptance rate
        step_size = base_step * (min_radius / self.params.radii[type_name])

        # Current position
        current_pos = positions[type_name][idx]

        # Symmetric Gaussian proposal
        proposal = current_pos + np.random.normal(0, step_size, 3)

        # Simple clipping to box boundaries [0, box_size]
        proposal = np.clip(proposal, 0.0, self.params.box_size)

        # Update the selected particle
        new_positions[type_name][idx] = proposal

        return new_positions
    #---------------------------------------------------------------------------
    def save_state_to_disk(self, step, positions, sigmas, score, 
                        prior_score=0, pair_score=0, exvol_score=0, tet_score=0, oct_score=0,
                        types=None, bead_numbers=None, traj_file=None):
        """
        Save state directly to an HDF5 file in a memory-efficient manner.
        This version writes the same information as the previous in-memory trajectory:
        - Attributes: step, total_score, prior_score, pair_score, exvol_score, tet_score, oct_score.
        - A subgroup 'sigma' with sigma values stored as attributes.
        - A subgroup 'positions' with each component saved as a dataset (gzip-compressed).
        - Datasets 'types_keys', 'types_vals', 'bead_keys', and 'bead_vals'.
        
        Parameters:
        step: The current step number.
        positions: Dictionary mapping component names to position arrays.
        sigmas: Dictionary of sigma values.
        score: The total score (will be stored as 'total_score').
        prior_score, pair_score, exvol_score, tet_score, oct_score: Additional scores.
        types: Dictionary mapping bead indices to type names.
        bead_numbers: Dictionary mapping bead indices to bead numbers.
        traj_file: Path to the HDF5 file to write to.
        """
        # If types or bead_numbers are not provided, default to empty dictionaries.
        if types is None:
            types = {}
        if bead_numbers is None:
            bead_numbers = {}
        if traj_file is None:
            return

        try:
            import numpy as np
            import h5py

            with h5py.File(traj_file, 'a') as f:
                # Create or get the trajectory group.
                if 'trajectory' not in f:
                    traj_grp = f.create_group('trajectory')
                else:
                    traj_grp = f['trajectory']
                
                # Create a new state group named "state_XXXXX" where XXXXX is the step number zero-padded.
                state_name = f"state_{step:05d}"
                state_grp = traj_grp.create_group(state_name)
                
                # Save state attributes (same keys as before).
                state_grp.attrs["step"] = step
                state_grp.attrs["total_score"] = float(score)
                state_grp.attrs["prior_score"] = float(prior_score)
                state_grp.attrs["pair_score"] = float(pair_score)
                state_grp.attrs["exvol_score"] = float(exvol_score)
                state_grp.attrs["tet_score"] = float(tet_score)
                state_grp.attrs["oct_score"] = float(oct_score)
                
                # Save sigma as a subgroup.
                sigma_grp = state_grp.create_group("sigma")
                for key, value in sigmas.items():
                    sigma_grp.attrs[key] = float(value)
                
                # Save positions as datasets (with gzip compression).
                pos_grp = state_grp.create_group("positions")
                for comp, coords in positions.items():
                    pos_grp.create_dataset(comp, data=coords.astype(np.float32), compression="gzip")
                
                # Save types and bead_numbers as datasets.
                types_keys = list(types.keys())
                types_vals = [types[k] for k in types_keys]
                state_grp.create_dataset("types_keys", data=np.array(types_keys, dtype="S"))
                state_grp.create_dataset("types_vals", data=np.array(types_vals, dtype="S"))
                
                bead_keys = list(bead_numbers.keys())
                bead_vals = [bead_numbers[k] for k in bead_keys]
                state_grp.create_dataset("bead_keys", data=np.array(bead_keys))
                state_grp.create_dataset("bead_vals", data=np.array(bead_vals))
        
        except Exception as e:
            print(f"Warning: Failed to save state to HDF5: {e}")
