import copy
import numpy as np
from typing import Dict, Tuple
from scipy.spatial.distance import cdist
from parameters import SystemParameters

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
        """Initialize positions either as ideal or randomly within the box.
        
        Both methods (ideal and random) center particles around origin (0,0,0).
        - Ideal: Uses pre-defined coordinates from SystemParameters
        - Random: Distributes particles in range [-box_size/2, box_size/2]
        """
        positions = {}
        if is_ideal:
            for type_name, _ in self.params.component_counts.items():
                # Use deep copy for ideal coordinates to avoid accidental modification.
                positions[type_name] = copy.deepcopy(self.params.ideal_coordinates[type_name])
        else:
            # Generate positions centered on (0,0,0) instead of corner at origin
            half_box = self.params.box_size / 2
            for type_name, count in self.params.component_counts.items():
                positions[type_name] = np.random.uniform(-half_box, half_box, (count, 3))
        
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

    def propose_sigma_move(self, sigma: Dict[str, float], accept_rate: float = None) -> Tuple[Dict[str, float], str]:
        """
        Additive symmetric proposal in linear sigma (no Jacobian needed).
        accept_rate kept for backward compatibility (ignored).
        """
        import random
        
        # Randomly select which sigma parameter to update
        pair_type = random.choice(list(sigma.keys()))
        
        # Current value
        current_value = sigma[pair_type]
        
        # Step size as percentage of current value
        relative_step_size = 0.005
        step_size = relative_step_size * max(current_value, 1e-12)
        
        # Symmetric proposal: add zero-mean Gaussian noise directly to value
        proposed_value = current_value + np.random.normal(0, step_size)
        
        # Ensure positivity and reasonable bounds
        proposed_value = max(1e-6, min(20.0, proposed_value))
        
        # Create new sigma dictionary with proposed value
        new_sigma = dict(sigma)
        new_sigma[pair_type] = proposed_value
        
        return new_sigma, pair_type

    def propose_position_move(self, positions: Dict[str, np.ndarray], accept_rate: float = 0.5) -> Dict[str, np.ndarray]:
        """
        MCMC position proposal with inverse radius scaling: larger particles move less.
        Uses clipping to stay within centered box boundaries.
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
        base_step = 2.5  # Adjust this value as needed for optimal acceptance rate
        step_size = base_step * (min_radius / self.params.radii[type_name])

        # Current position
        current_pos = positions[type_name][idx]

        # Symmetric Gaussian proposal
        proposal = current_pos + np.random.normal(0, step_size, 3)

        # Clip to centered box boundaries [-box_size/2, box_size/2]
        half_box = self.params.box_size / 2
        proposal = np.clip(proposal, -half_box, half_box)

        # Update the selected particle
        new_positions[type_name][idx] = proposal

        return new_positions

    def save_state_to_disk(self, step, positions, sigmas, score, 
                          prior_score=0, pair_score=0, exvol_score=0, tet_score=0, oct_score=0,
                          types=None, bead_numbers=None, traj_file=None):
        """
        Save state directly to an HDF5 file in a memory-efficient manner.
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