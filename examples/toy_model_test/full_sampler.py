#----------------------------------------------------------------------
# EM_Sampler class that uses only EM density map scoring for MCMC simulations
# Self-contained with cross-correlation calculation
#----------------------------------------------------------------------
import os
import sys
import logging
import math
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import h5py
import random
from pathlib import Path
from types import SimpleNamespace

# Math and file libraries
import scipy
import scipy.ndimage
import mrcfile

from base_sampler import BaseMCSampler
from parameters import SystemParameters
import networkx as nx

# --- Backend Selection ---
BACKEND = 'cpu'
try:
    import cupy as cp
    import cupyx.scipy.ndimage
    BACKEND = 'gpu'
    logging.warning('CuPy is available, using GPU backend.')
except ImportError:
    logging.info('CuPy is unavailable, using NumPy backend for CPU.')


class FullSampler(BaseMCSampler):
    """
    Sampler that uses only EM density map scoring with self-contained CCC calculation.
    """
    def __init__(self, 
                 sampler_sequence: List[str], 
                 sequence_idx: int,
                 em_map_file: str,
                 resolution: float = 20.0,
                 base_output_dir: str = "output_analysis", 
                 positions_init=None,
                 specific_chain: int = None):
        """
        Initialize EM_Sampler with EM density map.
        
        Args:
            sampler_sequence: List of sampler names in execution order
            sequence_idx: Current position in the sampler sequence
            em_map_file: Path to the EM density map file (MRC format)
            resolution: Resolution for EM density map scoring (Angstroms)
            base_output_dir: Base directory for analysis results
            positions_init: Optional initial positions to use
            specific_chain: Specific chain to load from (if not specified, random)
        """
        super().__init__()
        
        # Initialize parameters and positions
        self.params = SystemParameters()
        self.base_output_dir = base_output_dir
        self.sampler_sequence = sampler_sequence
        self.sequence_idx = sequence_idx
        
        # Move parameters
        self.position_step = 1.0
        self.tetramer_trans_step = 0.5
        self.tetramer_rot_step = 0.3
        self.octet_trans_step = 0.25
        self.octet_rot_step = 0.2
        
        # Handle positions - either load from previous sampler or initialize
        if positions_init is None:
            if sequence_idx > 0:
                # Try to load from previous sampler
                self.positions_os = self._load_from_previous_sampler()
            else:
                # Initialize new positions
                self.positions_os = self.initialize_positions()
        else:
            self.positions_os = positions_init
            
        # Create dummy sigma values (not used but needed for compatibility)
        self.sigma = {'AA': 1.0, 'AB': 1.0, 'BC': 1.0}
        
        # Set up EM restraint
        self.em_map_file = em_map_file
        self.resolution = resolution
        
        # Load target map
        if not os.path.exists(self.em_map_file):
            raise ValueError(f"EM map file not found: {self.em_map_file}")
            
        self.target_density_map = self.parse_density(self.em_map_file)
        print(f"Successfully loaded EM map: {self.em_map_file}")

    def _load_from_previous_sampler(self) -> Dict[str, np.ndarray]:
        """Load positions from the previous sampler in the sequence."""
        if self.sequence_idx == 0:
            return self.initialize_positions()
            
        prev_sampler = self.sampler_sequence[self.sequence_idx - 1]
        prev_dir = f"{self.base_output_dir}/{prev_sampler}sampler_results_1"
        
        if not os.path.exists(prev_dir):
            print(f"Previous sampler directory not found: {prev_dir}")
            return self.initialize_positions()
            
        # Try to load the latest trajectory file
        traj_file = os.path.join(prev_dir, "trajectory_chain_1.h5")
        if os.path.exists(traj_file):
            return self._load_positions_from_h5(traj_file)
        else:
            print(f"Trajectory file not found: {traj_file}")
            return self.initialize_positions()

    def _load_positions_from_h5(self, filename: str) -> Dict[str, np.ndarray]:
        """Load the latest positions from an HDF5 trajectory file."""
        try:
            with h5py.File(filename, 'r') as f:
                if 'trajectory' in f:
                    # Get the latest state
                    state_names = list(f['trajectory'].keys())
                    if state_names:
                        latest_state = f['trajectory'][state_names[-1]]
                        positions = {}
                        if 'positions' in latest_state:
                            for key in latest_state['positions'].keys():
                                positions[key] = latest_state['positions'][key][...]
                        return positions
        except Exception as e:
            print(f"Error loading from {filename}: {e}")
        
        return self.initialize_positions()

    # =====================================================================
    # CROSS-CORRELATION CALCULATION FUNCTIONS
    # =====================================================================
    
    def resolution_to_sigma(self, resolution: float, pixel_size: float) -> float:
        """Convert resolution to sigma for Gaussian blurring."""
        rsig_ = resolution / (4 * math.sqrt(2. * math.log(2.))) / pixel_size
        return rsig_

    def parse_density(self, fname: str):
        """Parse MRC density file."""
        density = mrcfile.open(fname, permissive=True)
        try:
            assert np.allclose(density.voxel_size.x, density.voxel_size.y) and np.allclose(density.voxel_size.x, density.voxel_size.z)
        except AssertionError:
            logging.error('Non-uniform grids are not supported')
            density.close()
            sys.exit()
        return density

    def bins_from_density(self, density) -> list:
        """Generate bins from density map."""
        binsx = (np.linspace(0, density.header.nx, density.header.nx + 1) - density.header.nx/2) * density.voxel_size.x
        binsy = (np.linspace(0, density.header.ny, density.header.ny + 1) - density.header.ny/2) * density.voxel_size.y
        binsz = (np.linspace(0, density.header.nz, density.header.nz + 1) - density.header.nz/2) * density.voxel_size.z
        return (binsx, binsy, binsz)

    def calc_projection_cpu(self, coords, weights, bins, resolution: float):
        """Calculate projection using CPU."""
        img_, _ = np.histogramdd(coords, weights=weights, bins=bins)
        img_ = np.swapaxes(img_, 0, 2)
        voxel_size = bins[0][1] - bins[0][0]
        sigma = self.resolution_to_sigma(resolution, voxel_size)
        img_sm = scipy.ndimage.gaussian_filter(img_, sigma, truncate=4).astype(np.float32)
        return img_sm

    def calc_projection_gpu(self, coords, weights, bins, resolution: float):
        """Calculate projection using GPU."""
        img_, _ = cp.histogramdd(coords, weights=weights, bins=bins)
        img_ = cp.swapaxes(img_, 0, 2)
        voxel_size = bins[0][1] - bins[0][0]
        sigma = self.resolution_to_sigma(resolution, float(voxel_size))
        img_sm = cupyx.scipy.ndimage.gaussian_filter(img_, sigma, truncate=4).astype(np.float32)
        return img_sm

    def pairwise_correlation_cpu(self, A, B):
        """Calculate pairwise correlation using CPU."""
        am = A - np.mean(A)
        bm = B - np.mean(B)
        return np.sum(am * bm) / (np.sqrt(np.sum(am**2)) * np.sqrt(np.sum(bm**2)))

    def pairwise_correlation_gpu(self, A, B):
        """Calculate pairwise correlation using GPU."""
        am = A - cp.mean(A)
        bm = B - cp.mean(B)
        return cp.sum(am * bm) / (cp.sqrt(cp.sum(am**2)) * cp.sqrt(cp.sum(bm**2)))

    def compare_data_cpu(self, density, model, resolution: float) -> float:
        """Compare model to density using CPU."""
        bins = self.bins_from_density(density)
        coords = model.getCoords()
        weights = model.getMasses()
        projection = self.calc_projection_cpu(coords, weights, bins, resolution)
        ccc = self.pairwise_correlation_cpu(projection.flatten(), density.data.flatten())
        return ccc

    def compare_data_gpu(self, density, model, resolution: float) -> float:
        """Compare model to density using GPU."""
        bins = self.bins_from_density(density)
        coords = model.getCoords()
        weights = model.getMasses()
        bins_gpu = cp.asarray(bins)
        coords_gpu = cp.asarray(coords)
        weights_gpu = cp.asarray(weights)
        projection = self.calc_projection_gpu(coords_gpu, weights_gpu, bins_gpu, resolution)
        density_data_gpu = cp.asarray(density.data)
        ccc = self.pairwise_correlation_gpu(projection.flatten(), density_data_gpu.flatten())
        return float(cp.asnumpy(ccc))

    def calculate_ccc_score(self, sphere_coords, sphere_radii, target_density_map, resolution, backend='cpu'):
        """
        Calculate the CCC score for a coarse-grained model against a target map.
        """
        # Create a "mock" model that provides coordinates and masses (derived from radii)
        mock_model = SimpleNamespace()
        mock_model.getCoords = lambda: sphere_coords
        # Use the volume (radius^3) as the weight for each sphere
        mock_model.getMasses = lambda: sphere_radii**3

        # Run the comparison using the script's existing functions
        ccc = -1.0  # Default bad score
        use_gpu = (BACKEND == 'gpu' and backend == 'gpu')

        if use_gpu:
            ccc = self.compare_data_gpu(target_density_map, mock_model, resolution)
        else:
            ccc = self.compare_data_cpu(target_density_map, mock_model, resolution)

        return ccc

    # =====================================================================
    # MOVE PROPOSAL FUNCTIONS
    # =====================================================================

    def propose_position_move(self, positions: Dict[str, np.ndarray], 
                             accept_rate: float = 0.3) -> Dict[str, np.ndarray]:
        """Propose a random position move for a single particle."""
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        all_types = list(positions.keys())
        if not all_types:
            return new_pos
        
        particle_type = random.choice(all_types)
        if len(new_pos[particle_type]) == 0:
            return new_pos
            
        particle_idx = np.random.randint(len(new_pos[particle_type]))
        
        # Adjust step size based on acceptance rate
        step_size = self.position_step
        if accept_rate > 0.5:
            step_size *= 1.05
        elif accept_rate < 0.2:
            step_size *= 0.95
            
        displacement = np.random.normal(0, step_size, 3)
        new_pos[particle_type][particle_idx] += displacement
        
        return new_pos

    def get_tetramers(self, positions: Dict[str, np.ndarray]) -> List[Tuple[int, ...]]:
        """Get tetramers from current positions."""
        # Simplified tetramer identification
        tetramers = []
        n_a = len(positions['A'])
        n_b = len(positions['B'])
        n_c = len(positions['C'])
        
        for i in range(min(n_a, n_b)):
            if i * 2 + 1 < n_c:
                tetramers.append((i, i, i * 2, i * 2 + 1))
        
        return tetramers

    def propose_tetramer_move(self, positions: Dict[str, np.ndarray], 
                             accept_rate: float = 0.3) -> Dict[str, np.ndarray]:
        """Propose a tetramer move."""
        new_pos = {k: v.copy() for k, v in positions.items()}
        tetramers = self.get_tetramers(positions)
        
        if not tetramers:
            return new_pos
            
        # Select random tetramer
        tet_idx = np.random.randint(len(tetramers))
        a_idx, b_idx, c_idx1, c_idx2 = tetramers[tet_idx]
        
        # Get tetramer coordinates
        tet_coords = np.array([
            new_pos['A'][a_idx],
            new_pos['B'][b_idx], 
            new_pos['C'][c_idx1],
            new_pos['C'][c_idx2]
        ])
        
        centroid = np.mean(tet_coords, axis=0)
        
        # Apply random transformation
        displacement = np.random.normal(0, self.tetramer_trans_step, 3)
        axis = self._random_unit_vector()
        angle = np.random.normal(0, self.tetramer_rot_step)
        rotation_matrix = self._rotation_matrix(axis, angle)
        
        # Transform coordinates
        transformed_coords = (tet_coords - centroid) @ rotation_matrix.T + centroid + displacement
        
        # Update positions
        new_pos['A'][a_idx] = transformed_coords[0]
        new_pos['B'][b_idx] = transformed_coords[1]
        new_pos['C'][c_idx1] = transformed_coords[2]
        new_pos['C'][c_idx2] = transformed_coords[3]
        
        return new_pos

    def get_octets(self, positions: Dict[str, np.ndarray]) -> Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], List[Tuple[int, ...]]]:
        """Identify octets (pairs of tetramers) in the structure."""
        tetramers = self.get_tetramers(positions)

        if len(tetramers) < 2:
            return [], tetramers

        # Compute geometric centers for each tetramer
        centers = np.zeros((len(tetramers), 3), dtype=np.float64)
        for i, (a_idx, b_idx, c_idx1, c_idx2) in enumerate(tetramers):
            coords = np.vstack([
                positions['A'][a_idx],
                positions['B'][b_idx],
                positions['C'][c_idx1],
                positions['C'][c_idx2]
            ])
            centers[i] = np.mean(coords, axis=0)

        # Build a graph and find maximum weight matching
        G = nx.Graph()
        for i_t in range(len(tetramers)):
            G.add_node(i_t)
        for i_t in range(len(tetramers)):
            for j_t in range(i_t + 1, len(tetramers)):
                dist_ij = np.linalg.norm(centers[i_t] - centers[j_t])
                G.add_edge(i_t, j_t, weight=-dist_ij)  # Negative for max weight = min distance

        matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)

        octets = []
        for i_t, j_t in matching:
            i_t, j_t = sorted([i_t, j_t])
            octets.append((tetramers[i_t], tetramers[j_t]))

        return octets, tetramers

    def propose_octet_move(self, positions: Dict[str, np.ndarray], octets) -> Dict[str, np.ndarray]:
        """Propose a move for a randomly selected octet."""
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        if not octets:
            return new_pos
        
        # Select random octet
        octet_idx = np.random.randint(len(octets))
        tetramer1, tetramer2 = octets[octet_idx]
        
        # Collect all particle indices in the octet
        octet_particles = []
        for tet in [tetramer1, tetramer2]:
            octet_particles.extend([
                ('A', tet[0]), ('B', tet[1]), ('C', tet[2]), ('C', tet[3])
            ])
        
        # Group indices by particle type
        grouped_indices = {}
        for particle_type, particle_idx in octet_particles:
            grouped_indices.setdefault(particle_type, []).append(particle_idx)
        
        # Get coordinates
        octet_coords_list = []
        for particle_type, indices_list in grouped_indices.items():
            octet_coords_list.append(new_pos[particle_type][indices_list])
        
        if not octet_coords_list:
            return new_pos 
            
        all_octet_coords = np.vstack(octet_coords_list)
        centroid = np.mean(all_octet_coords, axis=0)
        
        # Generate transformation
        displacement = np.random.normal(0, self.octet_trans_step, 3)
        axis = self._random_unit_vector()
        angle = np.random.normal(0, self.octet_rot_step)
        rotation_matrix = self._rotation_matrix(axis, angle)
        
        # Apply transformation
        for particle_type, indices_list in grouped_indices.items():
            current_positions = new_pos[particle_type][indices_list]
            transformed_positions = (current_positions - centroid) @ rotation_matrix.T + centroid + displacement
            new_pos[particle_type][indices_list] = transformed_positions
            
        return new_pos

    def _random_unit_vector(self):
        """Generate a random unit vector."""
        vec = np.random.randn(3)
        vec /= np.linalg.norm(vec) + 1e-10
        return vec

    def _rotation_matrix(self, axis, theta):
        """Create a 3D rotation matrix using Rodrigues' formula."""
        axis = np.asarray(axis)
        axis = axis / np.linalg.norm(axis)
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        
        return np.array([
            [a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
            [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
            [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]
        ])

    # =====================================================================
    # SCORING AND MCMC
    # =====================================================================

    def calculate_em_score(self, positions: Dict[str, np.ndarray]) -> Tuple[float, dict]:
        """Calculate the EM density map score for the current configuration."""
        all_coords = []
        all_radii = []
        
        # Combine all particles
        if 'A' in positions and len(positions['A']) > 0:
            all_coords.append(positions['A'])
            all_radii.append(np.full(len(positions['A']), self.params.radii['A']))
            
        if 'B' in positions and len(positions['B']) > 0:
            all_coords.append(positions['B'])
            all_radii.append(np.full(len(positions['B']), self.params.radii['B']))
            
        if 'C' in positions and len(positions['C']) > 0:
            all_coords.append(positions['C'])
            all_radii.append(np.full(len(positions['C']), self.params.radii['C']))
            
        if not all_coords:
            return 0.0, {"correlation": 0.0}
            
        sphere_coords = np.vstack(all_coords)
        sphere_radii = np.concatenate(all_radii)
        
        # Calculate cross-correlation coefficient
        ccc = self.calculate_ccc_score(
            sphere_coords, sphere_radii, 
            self.target_density_map, self.resolution, 'cpu'
        )
        
        # Convert to negative log-likelihood (lower CCC = higher penalty)
        # Use transformation: score = -log(CCC + 1) where CCC is normalized to [0,1]
        normalized_ccc = (ccc + 1) / 2  # Normalize from [-1,1] to [0,1]
        score = -np.log(normalized_ccc + 1e-10)  # Add small epsilon to avoid log(0)
        
        info = {"correlation": ccc}
        return score, info

    def run_mc(self, n_steps=50000, save_freq=1000, output_dir=None):
        """Monte Carlo sampling with position, tetramer, and octet moves."""
        if output_dir is None:
            output_dir = f"{self.base_output_dir}/em_sampler_results_1/"
            
        os.makedirs(output_dir, exist_ok=True)
        trajectory_file = os.path.join(output_dir, "trajectory_chain_1.h5")
        
        # Debug log file
        em_debug_log_file = os.path.join(output_dir, "em_debug.txt")

        with open(em_debug_log_file, 'w') as f:
            f.write("Step,MoveType,Correlation,Score,Accepted,Temp\n")
        
        # Create empty trajectory file
        with h5py.File(trajectory_file, 'w') as f:
            pass

        # Initialize tracking
        best_positions = {k: v.copy() for k, v in self.positions_os.items()}
        best_score = float('inf')
        em_history = {
            'correlation': np.zeros(n_steps // save_freq + 1),
            'score': np.zeros(n_steps // save_freq + 1)
        }
        accepts = {'position': 0, 'tetramer': 0, 'octet': 0}
        attempts = {'position': 0, 'tetramer': 0, 'octet': 0}

        # Initialize state
        current_positions = {k: v.copy() for k, v in self.positions_os.items()}
        current_octets, current_tetramers = self.get_octets(current_positions)

        # Calculate initial score
        current_score, curr_em_info = self.calculate_em_score(current_positions)

        if not np.isfinite(current_score):
            print(f"FATAL: Initial score is non-finite ({current_score}). Exiting.")
            return None, None

        best_score = current_score
        em_history['correlation'][0] = curr_em_info.get('correlation', 0.0)
        em_history['score'][0] = current_score

        # Save initial state
        self.save_state_to_disk(0, current_positions, self.sigma, current_score, traj_file=trajectory_file)

        # MCMC parameters
        move_types = ['position', 'tetramer', 'octet']
        move_probs = [0.3, 0.3, 0.4]
        temp_start, temp_end = 5.0, 0.1
        temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)

        print(f"Starting MCMC sampling for {n_steps} total steps...")

        # Main MCMC loop
        for step in range(n_steps):
            temp = max(temp_start * (temp_decay ** step), temp_end)

            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1

            # Propose move
            if move_type == 'position':
                accept_rate = accepts['position'] / max(1, attempts['position'])
                proposed_positions = self.propose_position_move(current_positions, accept_rate)
            elif move_type == 'tetramer':
                accept_rate = accepts['tetramer'] / max(1, attempts['tetramer'])
                proposed_positions = self.propose_tetramer_move(current_positions, accept_rate)
            else:  # octet
                proposed_positions = self.propose_octet_move(current_positions, current_octets)

            # Calculate score
            proposed_octets, proposed_tetramers = self.get_octets(proposed_positions)
            proposed_score, prop_em_info = self.calculate_em_score(proposed_positions)

            if not np.isfinite(proposed_score) or not np.isfinite(current_score):
                with open(em_debug_log_file, 'a') as f:
                    f.write(f"{step},{move_type},NaN,NaN,REJECTED(InvalidScore),{temp:.3f}\n")
                continue

            # Metropolis criterion
            delta = proposed_score - current_score
            accept = delta < 0 or (temp > 1e-9 and np.random.random() < np.exp(-delta / temp))

            # Log step
            em_info = prop_em_info if accept else curr_em_info
            with open(em_debug_log_file, 'a') as f:
                f.write(f"{step},{move_type},"
                        f"{em_info.get('correlation', 0.0):.6f},"
                        f"{proposed_score if accept else current_score:.6f},"
                        f"{'ACCEPTED' if accept else 'REJECTED'},"
                        f"{temp:.3f}\n")

            # Update state if accepted
            if accept:
                current_positions = proposed_positions
                current_tetramers = proposed_tetramers
                current_octets = proposed_octets
                current_score = proposed_score
                curr_em_info = prop_em_info
                accepts[move_type] += 1

                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in current_positions.items()}

            # Save state periodically
            if step % save_freq == 0 and step > 0:
                save_idx = step // save_freq
                em_history['correlation'][save_idx] = curr_em_info.get('correlation', 0.0)
                em_history['score'][save_idx] = current_score
                
                self.save_state_to_disk(step, current_positions, self.sigma, current_score, traj_file=trajectory_file)
                
                acceptance_rate = sum(accepts.values()) / (step + 1)
                print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, "
                      f"AcceptRate={acceptance_rate:.2f}, EM Corr={curr_em_info.get('correlation', 0.0):.4f}")

            # Print acceptance rates
            if step % 500 == 0 and step > 0:
                pos_rate = accepts['position'] / max(1, attempts['position'])
                tet_rate = accepts['tetramer'] / max(1, attempts['tetramer'])
                oct_rate = accepts['octet'] / max(1, attempts['octet'])
                print(f"Step {step}: Acceptance rates - Pos: {pos_rate:.2f}, "
                      f"Tet: {tet_rate:.2f}, Oct: {oct_rate:.2f}")

        # Save final results
        em_history_df = pd.DataFrame(em_history)
        em_history_df.to_csv(os.path.join(output_dir, "em_history.csv"), index=False)

        print("\nSampling complete:")
        for mv_type in move_types:
            rate = accepts[mv_type] / max(1, attempts[mv_type])
            print(f"- {mv_type}: {rate:.2f} acceptance ({accepts[mv_type]}/{attempts[mv_type]})")
        print(f"- Best score: {best_score:.4f}")
        print(f"- Final EM correlation: {curr_em_info.get('correlation', 0.0):.6f}")

        return best_positions, trajectory_file

    def save_state_to_disk(self, step, positions, sigma, total_score, traj_file="trajectory.h5"):
        """Save the current state to an HDF5 file."""
        with h5py.File(traj_file, 'a') as f:
            if 'trajectory' not in f:
                traj_group = f.create_group('trajectory')
            else:
                traj_group = f['trajectory']
                
            state_name = f"state_{step:05d}"
            if state_name in traj_group:
                del traj_group[state_name]
                
            state_group = traj_group.create_group(state_name)
            
            # Save attributes
            state_group.attrs['step'] = step
            state_group.attrs['total_score'] = float(total_score)
            
            # Save sigma values
            sigma_group = state_group.create_group('sigma')
            for key, value in sigma.items():
                sigma_group.attrs[key] = float(value)
                
            # Save positions
            pos_group = state_group.create_group('positions')
            for key, value in positions.items():
                if isinstance(value, np.ndarray) and len(value) > 0:
                    pos_group.create_dataset(key, data=value)