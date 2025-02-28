import numpy as np
from typing import Dict, Tuple
import copy
from scipy.spatial.distance import cdist
import random
from parameters import SystemParameters
import h5py
#---------------------------------------------------------------------------
class Priors:
    '''Defines different priors for the sigma parameter.'''
    def __init__(self, lb: float, ub: float, sigma: float, prior_type: str):
        self.lb = lb
        self.ub = ub
        self.sigma = sigma
        self.prior_type = prior_type
    
    def neg_log_prior(self) -> float:
        if self.sigma <= self.lb or self.sigma >= self.ub:
            return np.inf
        
        if self.prior_type == 'uniform':
            return 0.0
        elif self.prior_type == 'jeffreys':
            return np.log(self.sigma)
        elif self.prior_type == 'halfcauchy':
            scale = 1.0
            return np.log(1 + (self.sigma / scale) ** 2)
        else:
            raise ValueError(f"Unknown prior type: {self.prior_type}")
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
        priors: Dict[str, Priors] = {}

        for pair_type in self.params.pair_distances.keys():
            # Cache the sum of radii
            sum_radii = self.params.radii[pair_type[0]] + self.params.radii[pair_type[1]]
            lower_bound = 0.01 * sum_radii
            upper_bound = 0.25 * sum_radii
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

    def excluded_volume_nll(self, positions: Dict[str, np.ndarray], step: int = 0) -> float:
        """
        Compute the negative log likelihood of the excluded volume penalty.
        Uses vectorized operations over distance matrices.
        """
        score = 0.0
        exvol_sigma = 0.1  # Excluded volume sigma penalty strength

        for type1, pos1 in positions.items():
            for type2, pos2 in positions.items():
                # Only compute once for each unordered pair
                if type1 <= type2:
                    min_dist = self.params.radii[type1] + self.params.radii[type2]
                    distances = cdist(pos1, pos2)
                    if type1 == type2:
                        mask = np.triu(np.ones_like(distances), k=1)
                        viol_mask = (distances < min_dist) & (mask > 0)
                    else:
                        viol_mask = distances < min_dist

                    if np.any(viol_mask):
                        overlap = (min_dist - distances[viol_mask]) ** 2
                        score += np.sum(overlap / (exvol_sigma ** 2))
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
    
    def propose_sigma_move(self, sigma: Dict[str, float], accept_rate: float = 0.5 ) -> Tuple[Dict[str, float], str]:
        """
        Propose a move for one sigma parameter using a symmetric log-space proposal,
        with adaptive step size and reflection on the bounds.
        Takes the current acceptance rate as input to adjust the proposal step size.
        """
        new_sigma = sigma.copy()
        pair_types = list(sigma.keys())
        pair_type = random.choice(pair_types)
        
        # Cache physical bounds for the selected pair type.
        sum_radii = self.params.radii[pair_type[0]] + self.params.radii[pair_type[1]]
        min_sigma, max_sigma = 0.01 * sum_radii, 0.25 * sum_radii
        log_min, log_max = np.log(min_sigma), np.log(max_sigma)
        
        current_val = sigma[pair_type]
        base_step = min(0.1 * current_val, 0.02)
        # Adaptive adjustment based on acceptance rate.
        adjustment = np.clip(1.0 + 5.0 * (accept_rate - self.target_acceptance), 0.5, 2.0)
        step_size = base_step * adjustment
        
        log_current = np.log(current_val)
        log_proposed = log_current + np.random.normal(0, step_size)
        # Reflecting the proposal off the bounds.
        log_range = log_max - log_min
        log_proposed = log_min + np.abs(log_proposed - log_min) % (2 * log_range)
        log_proposed = log_proposed if log_proposed <= log_max else 2 * log_max - log_proposed
        
        new_sigma[pair_type] = np.exp(log_proposed)
        
        return new_sigma, pair_type  # Return the new sigma and the type for logging
    #---------------------------------------------------------------------------

    def propose_position_move(self, positions: Dict[str, np.ndarray], accept_rate: float = 0.5) -> Dict[str, np.ndarray]:
        """
        Propose a move for a randomly selected component's position with adaptive step size.
        Only the chosen type's positions are copied and updated.
        """
        new_positions = positions.copy()
        type_names = list(self.params.component_counts.keys())
        type_name = random.choice(type_names)
        
        base_step = self.params.radii[type_name] * 0.1
        adjustment = np.clip(1.0 + 5.0 * (accept_rate - self.target_acceptance), 0.5, 2.0)
        step_size = base_step * adjustment
        
        idx = np.random.randint(self.params.component_counts[type_name])
        # Copy only the array to be modified to avoid full deep copy.
        new_pos_array = np.copy(new_positions[type_name])
        new_pos_array[idx] = np.clip(new_pos_array[idx] + np.random.normal(0, step_size, 3),
                                     0, self.params.box_size)
        new_positions[type_name] = new_pos_array
        
        return new_positions
    
    def save_state(
        self, step: int, positions: Dict[str, np.ndarray], sigma: Dict[str, float], 
        total_score: float, prior_score: float, pair_score: float, exvol_score: float
    ) -> Dict:
        """Save the current state of the MCMC simulation with detailed score breakdown."""
        state = {
            "step": step,
            "positions": positions,
            "sigma": sigma,
            "total_score": total_score,
            "prior_score": prior_score,
            "pair_score": pair_score,
            "exvol_score": exvol_score,
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