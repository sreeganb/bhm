#----------------------------------------------------------------------
# TetramerSampler class which will use the base sampler functions 
# and run the MCMC simulations for the tetramers sampling and 
# track sigma values
#----------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import h5py
import random
from typing import List, Dict, Tuple
from base_sampler import BaseMCSampler
from parameters import SystemParameters
from scipy.spatial.distance import cdist
from pair_sampler import PairSampler
from sigma_provider import GMMSigmaProvider

#-----------------------------------------------------------------------
class TetramerSampler(BaseMCSampler):
    """
    Sampler for Tetramer-level interactions, inheriting from BaseMCSampler.
    """
    def __init__(self, sampler_sequence: List[str], sequence_idx: int,
                 base_output_dir: str = "output_analysis", positions_ts=None,
                 specific_chain: int = None, sigma_ranges: Dict[str, Tuple[float, float]] = None,
                 prior_type: str = "gamma"):
        """
        Initialize TetramerSampler with sampler sequence information.
        
        Args:
            sampler_sequence: List of sampler names in execution order
            sequence_idx: Current position in the sampler sequence
            base_output_dir: Base directory for analysis results
            positions_ts: Optional positions to use instead of loading
            specific_chain: Specific chain to load from (if not specified, random)
            sigma_ranges: Ranges for sigma values
            prior_type: Prior type for first sampler ("uniform" or "jeffreys")
        """
        super().__init__()
        self.params = SystemParameters()
        self.sampler_sequence = sampler_sequence
        self.sequence_idx = sequence_idx
        self.specific_chain = specific_chain
        
        # Initialize features for tetramer-specific sampling
        self.target_acceptance = 0.5
        
        # Validate required pair types
        required_pairs = {"AA", "AB", "BC"}
        if not required_pairs.issubset(self.params.pair_distances.keys()):
            raise ValueError("Missing required pair types for tetramer sampling.")
        
        # Handle positions
        if positions_ts is None:
            self.positions_ts = self.get_positions()
            print("Using positions from previous trajectory")
        else:
            self.positions_ts = positions_ts
            print("Using passed in positions")
        
        # Initialize sigma provider
        self.sigma_provider = GMMSigmaProvider(
            sampler_sequence=sampler_sequence,
            sequence_idx=sequence_idx,
            base_output_dir=base_output_dir,
            specific_chain=specific_chain,
            sigma_ranges=sigma_ranges,
            prior_type=prior_type
        )
        
        # Initialize sigma values and ranges from provider
        self.sigma = self.sigma_provider.sample_sigma_values()
        self.sigma_range = self.sigma_provider.sigma_ranges
        
        # Initialize PairSampler for scoring
        self.ps = PairSampler(
            sampler_sequence=sampler_sequence,
            sequence_idx=sequence_idx,
            base_output_dir=base_output_dir,
            specific_chain=specific_chain,
            sigma_ranges=sigma_ranges,
            prior_type=prior_type,
            pos_passed=self.positions_ts
        )

    def get_positions(self) -> Dict[str, np.ndarray]:
        """
        Load positions from the previous sampler in the sequence.
        """
        import pathlib
        
        if self.sequence_idx == 0:
            print("First sampler in sequence - using initialized positions")
            return self.initialize_positions()
        
        # Get the previous sampler info
        previous_sampler = self.sampler_sequence[self.sequence_idx - 1]
        
        # Count occurrences of the previous sampler up to current position
        occurrence_count = 0
        for i in range(self.sequence_idx):
            if self.sampler_sequence[i] == previous_sampler:
                occurrence_count += 1
        
        # Construct directory name
        traj_dir = pathlib.Path(os.getcwd()) / f"output_analysis/{previous_sampler}sampler_results_{occurrence_count}"
        
        try:
            # Get trajectory files
            trajectory_files = list(traj_dir.glob("trajectory_chain_*.h5"))
            
            if not trajectory_files:
                print(f"No trajectory files found in {traj_dir}")
                print("Falling back to initialized positions")
                return self.initialize_positions()
            
            # Select specific chain or random
            if self.specific_chain is not None:
                target_file = traj_dir / f"trajectory_chain_{self.specific_chain}.h5"
                if target_file.exists():
                    filepath = target_file
                    chain_num = self.specific_chain
                else:
                    print(f"Specified chain {self.specific_chain} not found, selecting random")
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
                    positions[type_name] = pos_grp[type_name][:].copy()
                    print(f"Loaded {len(positions[type_name])} {type_name} particles")
                    
                    # Debug: Print first few positions to verify they're reasonable
                    if len(positions[type_name]) > 0:
                        print(f"  First {type_name} position: {positions[type_name][0]}")
                        if len(positions[type_name]) > 1:
                            print(f"  Second {type_name} position: {positions[type_name][1]}")
            
            # Additional validation: Check if positions are within expected bounds
            box_size = getattr(self.params, 'box_size', 100.0)  # Default fallback
            for type_name, pos_array in positions.items():
                if len(pos_array) > 0:
                    min_coords = np.min(pos_array, axis=0)
                    max_coords = np.max(pos_array, axis=0)
                    print(f"{type_name} position range: min={min_coords}, max={max_coords}")
                    
                    # Check if any coordinates are outside expected bounds
                    if np.any(min_coords < 0) or np.any(max_coords > box_size):
                        print(f"WARNING: {type_name} positions outside expected bounds [0, {box_size}]")
            
            return positions

        except Exception as e:
            print(f"Error loading trajectory from {traj_dir}: {e}")
            print("Falling back to initialized positions")
            return self.initialize_positions()

#    def get_tetramers(self, positions: Dict[str, np.ndarray], temp: float = 0.9) -> List[Tuple[int, ...]]:
#        """Generate tetramers with particle exclusivity and distance-weighted selection."""
#        try:
#            # Quick validation
#            if not all(k in positions and len(positions[k]) > 0 for k in ['A', 'B', 'C']) or len(positions['C']) < 2:
#                return []
#            
#            a_pos, b_pos, c_pos = positions['A'], positions['B'], positions['C']
#            used_a, used_b, used_c = set(), set(), set()
#            tetramers = []
#            
#            # A-B distance calculation
#            dist_AB = cdist(a_pos, b_pos) / max(0.1, temp)
#            
#            # Process A particles in order of increasing minimum distance to any B
#            for a_idx in np.argsort(np.min(dist_AB, axis=1)):
#                if a_idx in used_a:
#                    continue
#                
#                # Find available B particles
#                b_mask = np.ones(len(b_pos), dtype=bool)
#                for idx in used_b:
#                    b_mask[idx] = False
#                
#                if not np.any(b_mask):
#                    continue
#                
#                # Calculate B selection probabilities
#                b_dists = dist_AB[a_idx].copy()
#                b_dists[~b_mask] = np.inf
#                min_b_dist = np.min(b_dists)
#                if np.isinf(min_b_dist):
#                    continue
#                
#                b_probs = np.exp(-(b_dists - min_b_dist))
#                b_probs = b_probs / np.sum(b_probs)
#                
#                # Select B particle
#                available_b = np.where(b_mask)[0]
#                b_idx = np.random.choice(available_b, p=b_probs[b_mask])
#                
#                # Find available C particles
#                c_mask = np.ones(len(c_pos), dtype=bool)
#                for idx in used_c:
#                    c_mask[idx] = False
#                
#                if np.sum(c_mask) < 2:  # Need at least 2 C particles
#                    continue
#                
#                # Calculate C selection probabilities
#                c_dists = cdist(b_pos[b_idx].reshape(1, -1), c_pos)[0] / max(0.1, temp)
#                c_dists[~c_mask] = np.inf
#                min_c_dist = np.min(c_dists)
#                if np.isinf(min_c_dist):
#                    continue
#                
#                c_probs = np.exp(-(c_dists - min_c_dist))
#                c_probs = c_probs / np.sum(c_probs)
#                
#                # Select two C particles
#                available_c = np.where(c_mask)[0]
#                c_indices = np.random.choice(
#                    available_c, size=2, replace=False,
#                    p=c_probs[c_mask] / np.sum(c_probs[c_mask])
#                )
#                
#                # Add tetramer and mark particles as used
#                tetramers.append((a_idx, b_idx, c_indices[0], c_indices[1]))
#                used_a.add(a_idx)
#                used_b.add(b_idx)
#                used_c.update(c_indices)
#                
#                # Stop if we have enough tetramers
#                if len(tetramers) >= min(len(a_pos), len(b_pos), len(c_pos) // 2):
#                    break
#            
#            return tetramers
#            
#        except Exception as e:
#            print(f"Error in tetramer generation: {e}")
#            return []
    #=======================================================================
    # New tetramer generation method using biophysical constraints
    #=======================================================================
    def get_tetramers(self, positions: Dict[str, np.ndarray]) -> List[Tuple[int, ...]]:
        """Generate tetramers based on optimal distance matching."""
        try:
            # Validate input
            if not all(k in positions and len(positions[k]) > 0 for k in ['A', 'B', 'C']) or len(positions['C']) < 2:
                return []
            
            a_pos, b_pos, c_pos = positions['A'], positions['B'], positions['C']
            
            # Get target distances from parameters
            ab_target = self.params.pair_distances['AB']
            bc_target = self.params.pair_distances['BC']
            
            # Calculate distance matrices
            dist_AB = cdist(a_pos, b_pos)
            
            # Score each A-B pair based on closeness to target distance
            ab_scores = np.abs(dist_AB - ab_target)
            
            # Create arrays to track used particles
            a_used = np.zeros(len(a_pos), dtype=bool)
            b_used = np.zeros(len(b_pos), dtype=bool)
            c_used = np.zeros(len(c_pos), dtype=bool)
            
            tetramers = []
            
            # Process A-B pairs in order of increasing score (closest to target)
            flat_indices = np.argsort(ab_scores.flatten())
            
            for flat_idx in flat_indices:
                a_idx = flat_idx // len(b_pos)
                b_idx = flat_idx % len(b_pos)
                
                # Skip if either particle is used
                if a_used[a_idx] or b_used[b_idx]:
                    continue
                
                # Find available C particles
                available_c = np.where(~c_used)[0]
                if len(available_c) < 2:
                    break
                    
                # Calculate distances from B to all available C particles
                b_c_dists = cdist(b_pos[b_idx].reshape(1, -1), c_pos[available_c])[0]
                
                # Score C particles by distance to target
                c_scores = np.abs(b_c_dists - bc_target)
                
                # Get the two best C particles
                best_c_indices = available_c[np.argsort(c_scores)[:2]]
                
                # Form tetramer
                tetramers.append((a_idx, b_idx, best_c_indices[0], best_c_indices[1]))
                
                # Mark particles as used
                a_used[a_idx] = True
                b_used[b_idx] = True
                c_used[best_c_indices] = True
                
                # Stop if we have enough tetramers
                if len(tetramers) >= min(len(a_pos), len(b_pos), len(c_pos) // 2):
                    break
            
            return tetramers
            
        except Exception as e:
            print(f"Error in tetramer generation: {e}")
            return []
    #=======================================================================
    def run_mc(self, n_steps=50000, save_freq=1000, output_dir="output_analysis/tetramersampler_results/", debug=False):
        """Monte Carlo sampling with position, sigma, and tetramer moves."""
        # Setup output directory
        os.makedirs(output_dir, exist_ok=True)
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        with h5py.File(trajectory_file, 'w') as f:
            pass  # Create empty file
        
        # Initialize tracking variables
        best_positions = None
        best_score = float('inf')
        sigma_history = {key: np.zeros(n_steps // save_freq + 1) for key in self.sigma}
        accepts = {'position': 0, 'sigma': 0, 'tetramer': 0}
        attempts = {'position': 0, 'sigma': 0, 'tetramer': 0}
        accepted_moves = 0
        
        # Initialize state
        current_tetramers = self.get_tetramers(self.positions_ts)
        prior_penalty = self.sigma_provider.calculate_negative_log_prior(self.sigma)
        current_score, curr_ex, curr_pair, curr_tet = self.neg_log_posterior(
            self.positions_ts, current_tetramers, prior_penalty, self.sigma, debug=False)
        
        # Store initial sigma values
        for key in self.sigma:
            sigma_history[key][0] = self.sigma[key]
        
        # Main MCMC loop
        move_types = ['position', 'sigma', 'tetramer']
        move_probs = [0.4, 0.1, 0.5]
        
        # Cooling schedule
        temp_start, temp_end = 5.0, 0.10
        temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
        
        print(f"Starting MCMC sampling for {n_steps} total steps...")
        
        for step in range(1, n_steps + 1):
            temp = temp_start * (temp_decay ** step)
            
            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1
            
            # Propose move
            pair_type = None
            if move_type == 'position':
                proposed_positions = self.propose_position_move(self.positions_ts, accepts['position'] / max(1, step))
                proposed_sigma = self.sigma
                proposed_tetramers = self.get_tetramers(proposed_positions)
            elif move_type == 'sigma':
                proposed_positions = self.positions_ts
                proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)
                proposed_tetramers = current_tetramers
            else:  # tetramer move
                proposed_positions = self.propose_tetramer_move(self.positions_ts, accepts['tetramer'] / max(1, step))
                proposed_sigma = self.sigma
                proposed_tetramers = self.get_tetramers(proposed_positions)
            
            # Calculate prior and score
            new_prior = self.sigma_provider.calculate_negative_log_prior(proposed_sigma)
            
            proposed_score, prop_ex, prop_pair, prop_tet = self.neg_log_posterior(
                proposed_positions, proposed_tetramers, new_prior, proposed_sigma, debug=False)
            
            # Metropolis criterion with Jacobian correction for sigma moves
            delta = proposed_score - current_score
            
            if move_type == 'sigma' and pair_type is not None:
                jacobian_term = np.log(proposed_sigma[pair_type] / self.sigma[pair_type])
                delta += jacobian_term
            
            # Accept/reject
            accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
            
            if accept:
                self.positions_ts = proposed_positions
                self.sigma = proposed_sigma
                current_tetramers = proposed_tetramers
                current_score = proposed_score
                curr_ex, curr_pair, curr_tet = prop_ex, prop_pair, prop_tet
                
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_ts.items()}
                
                accepts[move_type] += 1
                accepted_moves += 1
            
            # Save state
            if step % save_freq == 0:
                save_idx = step // save_freq
                if save_idx < len(sigma_history[list(sigma_history.keys())[0]]):
                    for key in self.sigma:
                        sigma_history[key][save_idx] = self.sigma[key]
                
                self.save_state_to_disk(
                    step, self.positions_ts, self.sigma, current_score,
                    prior_score=new_prior, pair_score=curr_pair,
                    exvol_score=curr_ex, tet_score=curr_tet,
                    traj_file=trajectory_file
                )
                
                acceptance_rate = accepted_moves / step
                print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={acceptance_rate:.2f}")
        
        # Save results
        import pandas as pd
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        # Print final statistics
        print("\nSampling complete:")
        for move_type in move_types:
            rate = accepts[move_type] / max(1, attempts[move_type])
            print(f"- {move_type}: {rate:.2f} acceptance ({accepts[move_type]}/{attempts[move_type]})")
        print(f"- Best score: {best_score:.2f}")
        
        return best_positions, trajectory_file

    def propose_tetramer_move(self, positions: Dict[str, np.ndarray], acceptance_rate: float) -> Dict[str, np.ndarray]:
        """Tetramer move: translate and rotate a randomly selected tetramer."""
        new_pos = {k: v.copy() for k, v in positions.items()}

        tetramers = self.get_tetramers(positions)
        if not tetramers:
            return new_pos

        # Select random tetramer
        tetramer = tetramers[np.random.randint(len(tetramers))]
        a_idx, b_idx, c_idx1, c_idx2 = tetramer
        particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]

        # Get coordinates and centroid
        coords = np.array([new_pos[part][idx] for part, idx in particles])
        centroid = np.mean(coords, axis=0)

        # Fixed step sizes
        trans_step = 0.25
        rot_step = 0.2

        # Translation
        displacement = np.random.normal(0.0, trans_step, 3)

        # Rotation
        rotation_axis = np.random.randn(3)
        rotation_axis /= np.linalg.norm(rotation_axis) + 1e-10
        rotation_angle = np.random.normal(0.0, rot_step)

        # Build rotation matrix
        half_angle = rotation_angle / 2.0
        qw = np.cos(half_angle)
        qx = rotation_axis[0] * np.sin(half_angle)
        qy = rotation_axis[1] * np.sin(half_angle)
        qz = rotation_axis[2] * np.sin(half_angle)
        rot_matrix = np.array([
            [1 - 2*(qy**2 + qz**2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
            [2*(qx*qy + qw*qz), 1 - 2*(qx**2 + qz**2), 2*(qy*qz - qw*qx)],
            [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx**2 + qy**2)]
        ])

        # Apply transformation
        for i, (part, idx) in enumerate(particles):
            vec = coords[i] - centroid
            rotated = rot_matrix @ vec
            final_pos = centroid + rotated + displacement
            new_pos[part][idx] = np.mod(final_pos, self.params.box_size)

        return new_pos

    def calculate_tetramer_scores_batch(self, positions, tetramers, sig, debug_logging=False):
        """Calculate scores for all tetramers."""
        if not tetramers:
            return np.array([], dtype=np.float32)
        
        # Extract indices
        a_indices = np.array([t[0] for t in tetramers], dtype=np.int32)
        b_indices = np.array([t[1] for t in tetramers], dtype=np.int32)
        c1_indices = np.array([t[2] for t in tetramers], dtype=np.int32)
        c2_indices = np.array([t[3] for t in tetramers], dtype=np.int32)
        
        # Get positions
        pos_a = positions['A'][a_indices]
        pos_b = positions['B'][b_indices]
        pos_c1 = positions['C'][c1_indices]
        pos_c2 = positions['C'][c2_indices]
        
        # Calculate distances
        ab_dists = np.sqrt(np.sum((pos_a - pos_b)**2, axis=1))
        bc1_dists = np.sqrt(np.sum((pos_b - pos_c1)**2, axis=1))
        bc2_dists = np.sqrt(np.sum((pos_b - pos_c2)**2, axis=1))
        
        # Target distances
        ab_target = self.params.pair_distances['AB']
        bc_target = self.params.pair_distances['BC']
        
        # Calculate scores
        ab_scores = ((ab_dists - ab_target)**2)/(2*sig['AB']**2) + 0.5 * np.log(2 * np.pi * sig['AB']**2)
        bc1_scores = ((bc1_dists - bc_target)**2)/(2*sig['BC']**2) + 0.5 * np.log(2 * np.pi * sig['BC']**2)
        bc2_scores = ((bc2_dists - bc_target)**2)/(2*sig['BC']**2) + 0.5 * np.log(2 * np.pi * sig['BC']**2)
        
        scores = ab_scores + bc1_scores + bc2_scores
        return scores

    def neg_log_posterior(self, positions: Dict[str, np.ndarray], tetramers: List[Tuple[int, ...]],
                         prior_penalty_from_distribution: float = 0.0, sig: Dict[str, float] = None,
                         exclusion_weight: float = 1.0, pair_weight: float = 1.0,
                         tetramer_weight: float = 1.0, debug: bool = False) -> Tuple[float, float, float, float]:
        """Calculate the negative log posterior."""
        sigma = sig if sig is not None else self.sigma
        
        if not tetramers:
            result = self.ps.calculate_score(
                positions, sigma, self.sigma_range,
                set(), prior_penalty_from_distribution
            )
            
            if isinstance(result, tuple) and len(result) >= 3:
                score, ex_score, pair_score = result[:3]
                return score, ex_score, pair_score, 0.0
            else:
                return 0.0, 0.0, 0.0, 0.0
        
        # Calculate score excluding tetramer pairs (empty set for now)
        tetramer_pairs = set()
        result = self.ps.calculate_score(
            positions, sigma, self.sigma_range,
            tetramer_pairs, prior_penalty_from_distribution
        )
        
        if isinstance(result, tuple) and len(result) >= 3:
            score, ex_score, pair_score = result[:3]
        else:
            return 0.0, 0.0, 0.0, 0.0

        # Calculate tetramer score
        scores_array = self.calculate_tetramer_scores_batch(positions, tetramers, sigma, debug_logging=debug)
        total_tet_score = scores_array.sum()
        weighted_tet_score = tetramer_weight * total_tet_score
        score += weighted_tet_score

        if debug:
            print(f"Exclusion: {ex_score:.2f}, Pair: {pair_score:.2f}, Tetramer: {weighted_tet_score:.2f}, Total: {score:.2f}")

        return score, ex_score, pair_score, weighted_tet_score
