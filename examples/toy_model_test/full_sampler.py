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
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from octet_sampler import OctetSampler
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
                 resolution: float = 50.0,
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
        
        # instantiate various samplers for move proposals
        self.pair_sampler = PairSampler(sampler_sequence, sequence_idx)
        self.tetramer_sampler = TetramerSampler(sampler_sequence, 
                                                sequence_idx)
        self.octet_sampler = OctetSampler(sampler_sequence, sequence_idx)

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
        """Calculate pairwise correlation using CPU with NaN handling."""
        # Check for empty arrays
        if len(A) == 0 or len(B) == 0:
            return 0.0
        
        # Center the arrays
        am = A - np.mean(A)
        bm = B - np.mean(B)
        
        # Calculate standard deviations
        std_a = np.sqrt(np.sum(am**2))
        std_b = np.sqrt(np.sum(bm**2))
        
        # Check for zero variance (constant arrays)
        if std_a < 1e-10 or std_b < 1e-10:
            print(f"Warning: Zero variance detected. std_a={std_a:.2e}, std_b={std_b:.2e}")
            return 0.0
        
        # Calculate correlation
        correlation = np.sum(am * bm) / (std_a * std_b)
        
        # Check for NaN result
        if not np.isfinite(correlation):
            print(f"Warning: Non-finite correlation detected: {correlation}")
            return 0.0
        
        return correlation

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
    # All of the necessary proposal functions are already implemented in various
    # files, such as base_sampler.py, tetramer_sampler.py, octet_sampler.py 
    # and can be reused here without modification.
    # =====================================================================
    # Create a new move, the full system treated as a single rigid body, and the 
    # entire system is translated or rotated as a whole.
    def propose_full_system_move(self, current_positions):
        """Propose a move that translates/rotates the entire system."""
        # Combine all positions into a single array
        all_coords = []
        for key in current_positions:
            if isinstance(current_positions[key], np.ndarray) and len(current_positions[key]) > 0:
                all_coords.append(current_positions[key])
        combined_coords = np.vstack(all_coords)

        # Calculate center of mass
        com = np.mean(combined_coords, axis=0)
        
        # Random translation
        full_trans_step = 1.0 # Translation step size
        displacement = np.random.normal(0, full_trans_step, 3)

        # Create a random vector passing through the COM
        full_rot_step = 0.5
        rand_vec = np.random.normal(size=3)
        rand_vec /= np.linalg.norm(rand_vec) + 1e-10  # Normalize and avoid division by zero
        angle = np.random.normal(0, full_rot_step)
        
        # Use Rodrigues' rotation formula
        K = np.array([[0, -rand_vec[2], rand_vec[1]],
                      [rand_vec[2], 0, -rand_vec[0]],
                      [-rand_vec[1], rand_vec[0], 0]])
        R = np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)
        
        # Apply transformation to all coordinates
        new_positions = {}
        for key in current_positions:
            if isinstance(current_positions[key], np.ndarray) and len(current_positions[key]) > 0:
                coords = current_positions[key]
                # Translate to origin
                centered_coords = coords - com
                # Rotate
                rotated_coords = centered_coords @ R.T
                # Translate back and apply displacement
                new_coords = rotated_coords + com + displacement
                new_positions[key] = new_coords
            else:
                new_positions[key] = current_positions[key]

        return new_positions
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
        
        # Instead of normalization, change to 1-ccc so that the ideal score is 0 and 
        # any deviation increase the score to be > 0.
        score = 100*(1 - ccc)
        
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
        accepts = {'position': 0, 'tetramer': 0, 'octet': 0, 'full': 0}
        attempts = {'position': 0, 'tetramer': 0, 'octet': 0, 'full': 0}

        # Initialize state
        current_positions = {k: v.copy() for k, v in self.positions_os.items()}
        current_octets, current_tetramers = self.octet_sampler.get_octets(current_positions)

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
        move_types = ['position', 'tetramer', 'octet', 'full']
        move_probs = [0.05, 0.075, 0.125, 0.75]
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
                proposed_positions = self.tetramer_sampler.propose_tetramer_move(current_positions, accept_rate)
            elif move_type == 'full':
                proposed_positions = self.propose_full_system_move(current_positions)
            else:  # octet
                proposed_positions = self.octet_sampler.propose_octet_move(current_positions, current_octets)

            # Calculate score
            proposed_octets, proposed_tetramers = self.octet_sampler.get_octets(proposed_positions)
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
