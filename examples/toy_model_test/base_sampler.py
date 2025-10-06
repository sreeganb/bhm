import copy
import numpy as np
from typing import Dict, Tuple
from scipy.spatial.distance import cdist
from parameters import SystemParameters
from typing import List, Optional
import os
import pathlib
import random
import h5py

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
     
    def get_positions_from_previous_sampler(self,
        sampler_sequence: List[str], 
        current_sequence_idx: int, 
        specific_chain: Optional[int] = None,
        base_output_dir: str = "output_analysis",
        params: Optional[SystemParameters] = None
    ) -> Dict[str, np.ndarray]:
        """
        Load positions from the previous sampler in the sequence or initialize if first.
        
        Args:
            sampler_sequence: Complete sequence of samplers (e.g., ['pair', 'tetramer', 'pair', 'octet'])
            current_sequence_idx: Index of current sampler in the sequence (0-based)
            specific_chain: Optional chain number to load from (if None, selects random)
            base_output_dir: Base directory for output files
            params: SystemParameters instance (if None, creates default)
        
        Returns:
            Dictionary of positions by particle type
        """
        
        # Initialize params if not provided
        if params is None:
            params = SystemParameters()
        
        if current_sequence_idx == 0:
            print("First sampler in sequence - using initialized positions")
            return self.initialize_positions()
        
        # Get the previous sampler info
        previous_sampler = sampler_sequence[current_sequence_idx - 1]
        
        # Count occurrences of the previous sampler up to current position
        occurrence_count = 0
        for i in range(current_sequence_idx):
            if sampler_sequence[i] == previous_sampler:
                occurrence_count += 1
        
        # Construct directory name
        traj_dir = pathlib.Path(os.getcwd()) / f"{base_output_dir}/{previous_sampler}sampler_results_{occurrence_count}"
        
        try:
            # Get trajectory files
            trajectory_files = list(traj_dir.glob("trajectory_chain_*.h5"))
            
            if not trajectory_files:
                print(f"No trajectory files found in {traj_dir}")
                print("Falling back to initialized positions")
                return self.initialize_positions()
            
            # Select specific chain or random
            if specific_chain is not None:
                target_file = traj_dir / f"trajectory_chain_{specific_chain}.h5"
                if target_file.exists():
                    filepath = target_file
                    chain_num = specific_chain
                else:
                    print(f"Specified chain {specific_chain} not found, selecting random")
                    filepath = random.choice(trajectory_files)
                    chain_num = int(filepath.stem.split('_')[-1])
            else:
                filepath = random.choice(trajectory_files)
                chain_num = int(filepath.stem.split('_')[-1])
            
            print(f"Loading positions from {previous_sampler}sampler_results_{occurrence_count}, chain: {chain_num}")
            
            with h5py.File(filepath, 'r') as f:
                if 'trajectory' not in f:
                    raise KeyError("Invalid trajectory file format: missing 'trajectory' group")
                
                traj_grp = f['trajectory']
                keys = list(traj_grp.keys())
                
                if not keys:
                    raise ValueError("Empty trajectory file")
                
                # Sort keys properly (state_00000, state_00001, etc.)
                def extract_step_number(state_name):
                    try:
                        return int(state_name.split('_')[-1])
                    except (ValueError, IndexError):
                        return 0
                
                keys.sort(key=extract_step_number)
                
                # Get last frame
                last_key = keys[-1]
                print(f"Using last frame: {last_key}")
                
                # Debug: Print step number and total score of last frame
                last_state_grp = traj_grp[last_key]
                step_num = last_state_grp.attrs.get("step", 0)
                total_score = last_state_grp.attrs.get("total_score", 0.0)
                print(f"Last frame details: step={step_num}, score={total_score:.4f}")
                
                # Read positions
                positions = {}
                pos_grp = traj_grp[last_key]['positions']
                
                for type_name in pos_grp:
                    # read dataset to numpy array
                    positions[type_name] = pos_grp[type_name][:].copy()
                    print(f"Loaded {len(positions[type_name])} {type_name} particles")
                    
                    # Debug: Print first few positions to verify they're reasonable
                    if len(positions[type_name]) > 0:
                        print(f"  First {type_name} position: {positions[type_name][0]}")
                        if len(positions[type_name]) > 1:
                            print(f"  Second {type_name} position: {positions[type_name][1]}")
            
            # Additional validation / recentering: use centered box [-box_size/2, box_size/2]
            box_size = float(getattr(params, 'box_size', 100.0))
            half_box = box_size / 2.0
            tol = box_size * 1e-6  # small tolerance
            
            for type_name, pos_array in positions.items():
                if len(pos_array) == 0:
                    continue
                min_coords = np.min(pos_array, axis=0)
                max_coords = np.max(pos_array, axis=0)
                print(f"{type_name} position range: min={min_coords}, max={max_coords}")
                
                # Case A: positions already centered in [-half, half] (expected)
                if np.all(min_coords >= -half_box - tol) and np.all(max_coords <= half_box + tol):
                    # already correct coordinate system
                    print(f"{type_name} positions appear centered in [-{half_box}, {half_box}].")
                
                # Case B: positions in [0, box_size] -> recenter to [-half, half]
                elif np.all(min_coords >= -tol) and np.all(max_coords <= box_size + tol):
                    print(f"{type_name} positions appear in [0, {box_size}]; recentering to [-{half_box}, {half_box}].")
                    positions[type_name] = positions[type_name] - half_box
                    min_coords2 = np.min(positions[type_name], axis=0)
                    max_coords2 = np.max(positions[type_name], axis=0)
                    print(f"  After recenter: min={min_coords2}, max={max_coords2}")
                
                else:
                    # Unexpected coordinate range: warn but keep data
                    print(f"WARNING: {type_name} positions outside expected centered bounds [-{half_box}, {half_box}] "
                        f"and not clearly in [0,{box_size}]. Leaving as-is.")
            
            return positions

        except Exception as e:
            print(f"Error loading trajectory from {traj_dir}: {e}")
            print("Falling back to initialized positions")
            return self.initialize_positions()

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

#    def propose_sigma_move(
#        self,
#        sigma: Dict[str, float],
#        accept_rate: Optional[float] = None
#    ) -> Tuple[Dict[str, float], str]:
#        """
#        Log-space random walk for sigma:
#            log σ' = log σ + 𝒩(0, τ)
#        This keeps proposals strictly positive and scale-aware.
#
#        Returns a fresh sigma dict and the key that was perturbed.
#        """
#        pair_type = random.choice(list(sigma.keys()))
#        current_val = float(sigma[pair_type])
#        low, high = self.sigma_range.get(pair_type, (1e-6, 20.0))
#
#        # Base log-step (≈ multiplicative factor e^{±0.1} ≈ 1.1×)
#        tau = 0.10
#
#        # Simple Robbins–Monro adaptation toward target_acceptance
#        if accept_rate is not None:
#            diff = accept_rate - self.target_acceptance
#            tau *= np.clip(1.0 + 2.0 * diff, 0.25, 2.5)
#
#        # Propose in log space and exponentiate
#        current_log = np.log(current_val)
#        proposed_log = current_log + np.random.normal(0.0, tau)
#        proposed_val = float(np.exp(proposed_log))
#
#        # Soft clip into prior range (reflect if we hit the bounds)
#        if proposed_val < low:
#            proposed_val = low * low / max(proposed_val, 1e-12)
#        elif proposed_val > high:
#            proposed_val = high * high / proposed_val
#
#        new_sigma = dict(sigma)
#        new_sigma[pair_type] = proposed_val
#        return new_sigma, pair_type

    def propose_sigma_move(
        self,
        sigma: Dict[str, float],
        accept_rate: Optional[float] = None
    ) -> Tuple[Dict[str, float], str]:
        """
        Non-adaptive Metropolis proposal that preserves detailed balance.

        - Selects a single pair_type uniformly at random.
        - Uses an additive Gaussian step in linear sigma with constant scale
          per parameter (independent of the current value/state).
        - Applies exact reflective boundary conditions on [low, high],
          which keeps proposals unbiased and symmetric within bounds.

        Returns:
            (new_sigma, pair_type): a shallow copy with one updated entry.
        """
        import numpy as np
        import random

        # Choose which parameter to update (uniform)
        pair_type = random.choice(list(sigma.keys()))
        current_val = float(sigma[pair_type])

        # Bounds (strictly positive)
        low, high = self.sigma_range.get(pair_type, (1e-6, 20.0))
        # Guard invalid current values by snapping into bounds
        if not np.isfinite(current_val) or current_val <= 0.0:
            current_val = np.clip((low + high) * 0.5 if np.isfinite(current_val) else (low + high) * 0.5, low, high)

        # Constant, state-independent proposal width (ensures symmetry q(x->y)=q(y->x))
        width = max(high - low, 1e-9)
        step_sd = 0.15 * width  # tune as needed; constant for this pair_type

        # Symmetric additive Gaussian proposal in sigma-space
        proposed = current_val + np.random.normal(0.0, step_sd)

        # Reflective boundary conditions preserve symmetry on [low, high]
        def reflect(x: float, a: float, b: float) -> float:
            w = b - a
            if w <= 0.0:
                return float(np.clip(x, a, b))
            # Repeated reflections until inside [a, b]
            while x < a or x > b:
                if x < a:
                    x = a + (a - x)
                if x > b:
                    x = b - (x - b)
            return float(x)

        proposed_val = reflect(proposed, low, high)

        new_sigma = dict(sigma)
        new_sigma[pair_type] = proposed_val
        return new_sigma, pair_type
    
#    def propose_sigma_move(self, sigma: Dict[str, float], accept_rate: float = None) -> Tuple[Dict[str, float], str]:
#        """
#        Additive symmetric proposal in linear sigma (no Jacobian needed).
#        accept_rate kept for backward compatibility (ignored).
#        """
#        import random
#        
#        # Randomly select which sigma parameter to update
#        pair_type = random.choice(list(sigma.keys()))
#        
#        # Current value
#        current_value = sigma[pair_type]
#        
#        # Step size as percentage of current value
#        relative_step_size = 0.005  # changed from 0.005 to 0.005
#        step_size = relative_step_size * max(current_value, 1e-12)
#        
#        # Symmetric proposal: add zero-mean Gaussian noise directly to value
#        proposed_value = current_value + np.random.normal(0, step_size)
#        
#        # Ensure positivity and reasonable bounds
#        proposed_value = max(1e-6, min(20.0, proposed_value))
#        
#        # Create new sigma dictionary with proposed value
#        new_sigma = dict(sigma)
#        new_sigma[pair_type] = proposed_value
#        
#        return new_sigma, pair_type

#    def propose_position_move(self, positions: Dict[str, np.ndarray], accept_rate: float = 0.5) -> Dict[str, np.ndarray]:
#        """
#        MCMC position proposal with inverse radius scaling: larger particles move less.
#        Uses clipping to stay within centered box boundaries.
#        """
#        import random
#        
#        # Deep copy to avoid modifying original
#        new_positions = {key: np.copy(array) for key, array in positions.items()}
#        
#        # Randomly select particle type and index
#        type_names = list(self.params.component_counts.keys())
#        type_name = random.choice(type_names)
#        idx = np.random.randint(self.params.component_counts[type_name])
#
#        # Get min radius across all particle types for scaling
#        min_radius = min(self.params.radii.values())
#
#        # Inverse scaling: smaller radius = larger moves
#        base_step = 4.0  # Adjust this value as needed for optimal acceptance rate
#        step_size = base_step * (min_radius / self.params.radii[type_name])
#
#        # Current position
#        current_pos = positions[type_name][idx]
#
#        # Symmetric Gaussian proposal
#        proposal = current_pos + np.random.normal(0, step_size, 3)
#
#        # Clip to centered box boundaries [-box_size/2, box_size/2]
#        half_box = self.params.box_size / 2
#        proposal = np.clip(proposal, -half_box, half_box)
#
#        # Update the selected particle
#        new_positions[type_name][idx] = proposal
#
#        return new_positions

    def propose_position_move(
        self,
        positions: Dict[str, np.ndarray],
        accept_rate: float = 0.5
    ) -> Dict[str, np.ndarray]:
        """
        Single-particle Gaussian move with reflective walls.
        """
        new_positions = {k: v.copy() for k, v in positions.items()}
        box_half = self.params.box_size / 2.0

        type_names = list(self.params.component_counts.keys())
        counts = [self.params.component_counts[t] for t in type_names]
        cum_counts = np.cumsum(counts)
        flat_idx = np.random.randint(cum_counts[-1])

        for type_name, count, cum in zip(type_names, counts, cum_counts):
            if flat_idx < cum:
                local_idx = flat_idx - (cum - count)
                break

        radius = self.params.radii[type_name]
        max_radius = max(self.params.radii.values())
        step_sigma = 2.0 * (max_radius / radius)

        proposal = positions[type_name][local_idx] + np.random.normal(0.0, step_sigma, 3)

        for dim in range(3):
            while proposal[dim] > box_half or proposal[dim] < -box_half:
                if proposal[dim] > box_half:
                    proposal[dim] = 2 * box_half - proposal[dim]
                else:
                    proposal[dim] = -2 * box_half - proposal[dim]

        new_positions[type_name][local_idx] = proposal
        return new_positions

    def save_state_to_disk(self, step, positions, sigmas, score, 
                          prior_score=0, pair_score=0, exvol_score=0, tet_score=0, oct_score=0,
                          types=None, bead_numbers=None, traj_file=None):
        """
        Save state directly to an HDF5 file.
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
