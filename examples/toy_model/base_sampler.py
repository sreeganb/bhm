import numpy as np
from typing import Dict, Tuple
import copy
from scipy.spatial.distance import cdist
import random
from parameters import SystemParameters
import h5py
import os
#---------------------------------------------------------------------------
class Priors:
    '''Defines different priors for the sigma parameter.'''
    def __init__(self, prior_type: str):
        self.prior_type = prior_type
    
    def neg_log_prior(self, sigma: Dict[str, float], sigma_range: Dict[str, Tuple[float, float]]) -> float:
        # for all the pair types loop over and find the total negative log prior
        total_prior = 0.0
        for pair_type in sigma.keys():
            if sigma[pair_type] <= sigma_range[pair_type][0] or sigma[pair_type] >= sigma_range[pair_type][1]:
                return np.inf
            
            if self.prior_type == 'uniform':
                total_prior += 0.0
            elif self.prior_type == 'jeffreys':
                total_prior += np.log(sigma[pair_type])
            elif self.prior_type == 'halfcauchy':
                scale = 1.0
                total_prior += np.log(1 + (sigma[pair_type] / scale) ** 2)
            else:
                raise ValueError(f"Unknown prior type: {self.prior_type}")
        return total_prior
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
            lower_bound = 0.035 * sum_radii
            upper_bound = 0.35 * sum_radii
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
    def propose_sigma_move(self, sigma: Dict[str, float], accept_rate: float = 0.5) -> Tuple[Dict[str, float], str]:
        """
        Propose a move for one sigma parameter using symmetrical sampling in log-space.
        If the proposed sigma is outside [min_sigma, max_sigma], the prior penalty likely rejects it.
        """
        import random
        pair_type = random.choice(list(sigma.keys()))

        # Current value in log-space
        log_current = np.log(sigma[pair_type])

        # Smaller base step size and narrower factor
        base_step_size = 0.04 # switching from 0.002 to 0.004
        step_factor = (1.0 + 3.0 * np.clip(accept_rate - self.target_acceptance, -0.1, 0.3))
        step_size = base_step_size * step_factor

        log_proposed = log_current + np.random.normal(0, step_size)

        new_sigma = dict(sigma)
        new_sigma[pair_type] = np.exp(log_proposed)

        return new_sigma, pair_type
    #---------------------------------------------------------------------------
    def propose_position_move(self, positions: Dict[str, np.ndarray], accept_rate: float = 0.5) -> Dict[str, np.ndarray]:
        """
        Propose a move for a randomly selected component's position using a symmetric, 
        zero-mean Gaussian proposal with reflection at boundaries.
        """
        def reflect_in_box(coord: float, box_size: float) -> float:
            """
            Reflect a single coordinate if it goes out of [0, box_size].
            This ensures that proposals remain within the box boundaries symmetrically.
            """
            # A simple 'bounce' reflection: if we go below 0, reflect up;
            # if we go above box_size, reflect down.
            if coord < 0:
                coord = -coord
                if coord > box_size:
                    coord = 2 * box_size - coord
            elif coord > box_size:
                coord = 2 * box_size - coord
                if coord < 0:
                    coord = -coord
            return coord

        import random
        new_positions = positions.copy()
        type_names = list(self.params.component_counts.keys())
        type_name = random.choice(type_names)

        base_step = self.params.radii[type_name] * 0.25
        # Adaptive factor, with a zero-mean Gaussian ensuring symmetry
        adjustment = np.clip(1.0 + 5.0 * (accept_rate - self.target_acceptance), 0.1, 2.75)
        step_size = base_step * adjustment

        idx = np.random.randint(self.params.component_counts[type_name])
        # Copy only the array to be modified
        new_pos_array = np.copy(new_positions[type_name])

        # Generate a zero-mean proposal
        proposal = new_pos_array[idx] + np.random.normal(0, step_size, 3)
        # Reflect in all dimensions that exceed boundaries
        for i in range(3):
            proposal[i] = reflect_in_box(proposal[i], self.params.box_size)

        new_pos_array[idx] = proposal
        new_positions[type_name] = new_pos_array

        return new_positions
    #---------------------------------------------------------------------------
    def save_state(
        self, step: int, positions: Dict[str, np.ndarray], sigma: Dict[str, float], 
        total_score: float, prior_score: float, pair_score: float, exvol_score: float, 
        tet_score: float = 0.0, oct_score: float = 0.0) -> Dict:
        """Save the current state of the MCMC simulation with detailed score breakdown."""
        state = {
            "step": step,
            "positions": positions,
            "sigma": sigma,
            "total_score": total_score,
            "prior_score": prior_score,
            "pair_score": pair_score,
            "exvol_score": exvol_score,
            "tet_score": tet_score,
            "oct_score": oct_score,
            "types": {},
            "bead_numbers": {},
        }
        
        current_idx = 0
        for type_name, count in self.params.component_counts.items():
            for _ in range(count):
                state["types"][current_idx] = type_name
                state["bead_numbers"][current_idx] = current_idx + 1
                current_idx += 1

        return state
    
    def save_state_to_group(self, group: h5py.Group, state: Dict):
        """Save the state dictionary into an HDF5 group with all components."""
        group.attrs["step"] = state["step"]
        group.attrs["total_score"] = state["total_score"]
        group.attrs["prior_score"] = state["prior_score"]
        group.attrs["pair_score"] = state["pair_score"]
        group.attrs["exvol_score"] = state["exvol_score"]
        group.attrs["tet_score"] = state["tet_score"]
        group.attrs["oct_score"] = state["oct_score"]

        # Save sigma as a subgroup
        sigma_grp = group.create_group("sigma")
        for key, value in state["sigma"].items():
            sigma_grp.attrs[key] = value

        # Save positions as datasets
        pos_grp = group.create_group("positions")
        for type_name, array in state["positions"].items():
            pos_grp.create_dataset(type_name, data=array, compression="gzip")

        # Save types and bead_numbers as datasets
        types_keys = list(state["types"].keys())
        types_vals = [state["types"][k] for k in types_keys]
        group.create_dataset("types_keys", data=np.array(types_keys, dtype="S"))
        group.create_dataset("types_vals", data=np.array(types_vals, dtype="S"))

        bead_keys = list(state["bead_numbers"].keys())
        bead_vals = [state["bead_numbers"][k] for k in bead_keys]
        group.create_dataset("bead_keys", data=np.array(bead_keys))
        group.create_dataset("bead_vals", data=np.array(bead_vals))
        
    def save_trajectory(self, trajectory, filename: str = "trajectory.h5") -> str:
        """Save the trajectory of the MCMC simulation to an HDF5 file with all scores."""
        with h5py.File(filename, "w") as f:
            traj_grp = f.create_group("trajectory")
            for i, state in enumerate(trajectory):
                state_name = f"state_{i:05d}"
                state_grp = traj_grp.create_group(state_name)
                self.save_state_to_group(state_grp, state)
                
        return filename

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
