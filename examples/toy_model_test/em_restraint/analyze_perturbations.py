#!/usr/bin/env python3
"""
PerturbationAnalyzer: A class-based approach for generating and 
scoring molecular perturbations with ALL scoring functions.

Key Features:
- Clean class structure with shared state
- Minimal code redundancy 
- Well-named, focused methods
- Configurable parameters
- Comprehensive scoring (CCC+ExVol, Pair, Tetramer, Octet)
- RMSD analysis (raw vs aligned)
- Always includes unperturbed reference structure as first frame
"""

import numpy as np
import os
import sys
import h5py
from scipy.spatial.distance import cdist

# -----------------------------------------------------------------------------
# Import external scoring helpers
# -----------------------------------------------------------------------------
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from working_example import calculate_ccc_score, parse_density, create_dummy_map_from_model
from parameters import SystemParameters
from base_sampler import BaseMCSampler
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from octet_sampler import OctetSampler

class PerturbationAnalyzer:
    """
    Generates and analyzes molecular perturbations with comprehensive scoring.
    """
    
    def __init__(self, n_target=100, jitter_range=(0.1, 25.0), max_total_score=10000.0):
        """Initialize analyzer with configuration parameters."""
        # Configuration
        self.n_target = n_target
        self.jitter_range = jitter_range
        self.max_total_score = max_total_score
        self.incremental_steps = 5
        self.overlap_tolerance = 1.0
        self.max_relax_iters = 120
        self.pair_buffer = 0.25
        
        # System definition
        self._setup_molecular_system()
        
        # Analysis data
        self.ref_coords = None
        self.target_density = None
        self.results = []
        
        # Scoring components
        self.sigmas = {"AA": 1.5, "AB": 1.2, "BC": 1.0}
        self.base_sampler = BaseMCSampler()
        self.pair_sampler = PairSampler(sampler_sequence=["pair"], sequence_idx=0)
        self.tetramer_sampler = TetramerSampler(sampler_sequence=["tetramer"], sequence_idx=0)
        self.octet_sampler = OctetSampler(sampler_sequence=["octet"], sequence_idx=0)

    def _setup_molecular_system(self):
        """Define the reference molecular system."""
        self.array_A = np.array([
            [ 63.00,   0.00,   0.00], [ 44.55,  44.55,   0.00], [  0.00,  63.00,   0.00],
            [-44.55,  44.55,   0.00], [-63.00,   0.00,   0.00], [-44.55, -44.55,   0.00],
            [ -0.00, -63.00,   0.00], [ 44.55, -44.55,   0.00]
        ])
        self.array_B = np.array([
            [ 63.00,   0.00, -38.50], [ 44.55,  44.55, -38.50], [  0.00,  63.00, -38.50],
            [-44.55,  44.55, -38.50], [-63.00,   0.00, -38.50], [-44.55, -44.55, -38.50],
            [ -0.00, -63.00, -38.50], [ 44.55, -44.55, -38.50]
        ])
        self.array_C = np.array([
            [ 47.00,   0.00, -68.50], [ 79.00,   0.00, -68.50], [ 55.86,  55.86, -68.50],
            [ 33.23,  33.23, -68.50], [  0.00,  47.00, -68.50], [  0.00,  79.00, -68.50],
            [-55.86,  55.86, -68.50], [-33.23,  33.23, -68.50], [-47.00,   0.00, -68.50],
            [-79.00,   0.00, -68.50], [-55.86, -55.86, -68.50], [-33.23, -33.23, -68.50],
            [  0.00, -47.00, -68.50], [  0.00, -79.00, -68.50], [ 55.86, -55.86, -68.50],
            [ 33.23, -33.23, -68.50]
        ])
        self.n_A, self.n_B, self.n_C = len(self.array_A), len(self.array_B), len(self.array_C)
        
    def initialize_system(self):
        """Initialize reference coordinates and density map."""
        self.ref_coords = np.vstack([self.array_A, self.array_B, self.array_C])
        self._setup_density_map()
        
    def _setup_density_map(self):
        """Create or load the target density map."""
        target_map_file = "target_map.mrc"
        resolution, voxel_size, box_size = 15.0, 4.0, 800.0
        
        # Create dummy radii array for density map generation
        radii_array = np.concatenate([
            np.full(self.n_A, 24.0),  # radius A
            np.full(self.n_B, 14.0),  # radius B  
            np.full(self.n_C, 16.0)   # radius C
        ])
        
        if not os.path.exists(target_map_file):
            print("Creating target density map...")
            create_dummy_map_from_model(
                self.ref_coords, radii_array, resolution, 
                voxel_size, box_size, target_map_file
            )
        self.target_density = parse_density(target_map_file)
        
    def generate_perturbation(self, intensity):
        """Generate a perturbed structure with given jitter intensity."""
        step_sigma = intensity / np.sqrt(self.incremental_steps)
        coords = self.ref_coords.copy()
        
        for _ in range(self.incremental_steps):
            coords += np.random.normal(scale=step_sigma, size=coords.shape)
            
        return self._resolve_clashes(coords)
        
    def _resolve_clashes(self, coords):
        """Resolve overlapping particles using iterative pair separation."""
        coords = coords.copy()
        
        # Create radii array for clash detection
        radii_array = np.concatenate([
            np.full(self.n_A, 24.0),
            np.full(self.n_B, 14.0), 
            np.full(self.n_C, 16.0)
        ])
        
        for _ in range(self.max_relax_iters):
            clashes = self._find_clashes(coords, radii_array)
            if not clashes:
                break
            coords = self._separate_clashing_pairs(coords, clashes)
            
        return coords
        
    def _find_clashes(self, coords, radii_array):
        """Identify clashing particle pairs."""
        distances = cdist(coords, coords)
        np.fill_diagonal(distances, np.inf)
        min_allowed = (radii_array[:, None] + radii_array) * self.overlap_tolerance
        
        clash_mask = distances < min_allowed
        if not np.any(clash_mask):
            return None
            
        return {
            'mask': clash_mask,
            'distances': distances, 
            'min_allowed': min_allowed
        }
        
    def _separate_clashing_pairs(self, coords, clashes):
        """Separate clashing pairs by moving them apart."""
        coords = coords.copy()
        i_idx, j_idx = np.where(clashes['mask'])
        processed = set()
        
        for i, j in zip(i_idx, j_idx):
            if i >= j or (i, j) in processed:
                continue
                
            distance = clashes['distances'][i, j]
            if distance < 1e-8:
                direction = np.random.normal(size=3)
                direction /= np.linalg.norm(direction)
            else:
                direction = (coords[j] - coords[i]) / distance
                
            needed_separation = (clashes['min_allowed'][i, j] + self.pair_buffer) - distance
            if needed_separation > 0:
                shift = 0.5 * needed_separation
                coords[i] -= direction * shift
                coords[j] += direction * shift
                
            processed.add((i, j))
            
        return coords
        
    def calculate_scores(self, coords):
        """Calculate ALL scoring functions for given coordinates."""
        # Create radii array for CCC calculation
        radii_array = np.concatenate([
            np.full(self.n_A, 24.0),
            np.full(self.n_B, 14.0),
            np.full(self.n_C, 16.0)
        ])
        
        # Position dictionary for all scorers
        positions = {
            'A': coords[:self.n_A],
            'B': coords[self.n_A:self.n_A+self.n_B],
            'C': coords[self.n_A+self.n_B:]
        }
        
        # 1. CCC + Excluded Volume Score
        ccc_raw = calculate_ccc_score(
            coords, radii_array, self.target_density, 15.0, backend='cpu'
        )
        ccc_score = 500.0 * (1.0 - ccc_raw)
        
        # Get excluded volume from pair sampler
#        _, ex_score, _, _ = self.pair_sampler.calculate_score(pos=positions, sig=self.sigmas)
#        ccc_exvol_score = ccc_score + ex_score
        
        # 2. Pair Score (includes excluded volume, pair, and prior components)
        pair_total_score, ex_score, _, _ = self.pair_sampler.calculate_score(pos=positions, sig=self.sigmas)
        ccc_exvol_score = ccc_score + ex_score

        # 3. Tetramer Score
        tetramers = self.tetramer_sampler.get_tetramers(positions)
        tetramer_score, _, _, _ = self.tetramer_sampler.neg_log_posterior(
            positions=positions, tetramers=tetramers, sig=self.sigmas
        )
        
        # 4. Octet Score
        octets, tetramers = self.octet_sampler.get_octets(positions)
        octet_score, _, _, _, _ = self.octet_sampler.neg_log_posterior(
            positions=positions, tetramers=tetramers, 
            octets=octets, sig=self.sigmas
        )
        
        return {
            'ccc_exvol': ccc_exvol_score,
            'pair': pair_total_score,
            'tetramer': tetramer_score,
            'octet': octet_score,
            'ccc_only': ccc_score,
            'exvol_only': ex_score
        }
                
    def calculate_rmsd_metrics(self, coords):
        """Calculate both raw and aligned RMSD values."""
        raw_rmsd = np.sqrt(np.mean(np.sum((self.ref_coords - coords)**2, axis=1)))
        aligned_rmsd = self._calculate_aligned_rmsd(coords)
        return raw_rmsd, aligned_rmsd
        
    def _calculate_aligned_rmsd(self, coords):
        """Calculate RMSD after optimal alignment using Kabsch algorithm."""
        ref_centered = self.ref_coords - self.ref_coords.mean(0)
        coords_centered = coords - coords.mean(0)
        
        H = coords_centered.T @ ref_centered
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        if np.linalg.det(R) < 0:
            Vt[-1] *= -1
            R = Vt.T @ U.T
            
        aligned_coords = (R @ coords_centered.T).T + self.ref_coords.mean(0)
        return np.sqrt(np.mean(np.sum((self.ref_coords - aligned_coords)**2, axis=1)))
        
    def is_structure_acceptable(self, all_scores):
        """Determine if structure meets acceptance criteria based on any score."""
        # Accept if any of the main scores is below threshold
        acceptable_scores = [
            all_scores['ccc_exvol'],
            all_scores['pair'],
            all_scores['tetramer'],
            all_scores['octet']
        ]
        return any(score <= self.max_total_score for score in acceptable_scores)
        
    def create_result_record(self, coords, all_scores, rmsd_metrics, intensity, attempt_num, is_reference=False):
        """Create a standardized result record with ALL scores."""
        raw_rmsd, aligned_rmsd = rmsd_metrics
        
        record = {
            'coords_A': coords[:self.n_A],
            'coords_B': coords[self.n_A:self.n_A+self.n_B], 
            'coords_C': coords[self.n_A+self.n_B:],
            'rmsd': aligned_rmsd,
            'rmsd_aligned': aligned_rmsd,
            'rmsd_raw': raw_rmsd,
            'intensity_requested': intensity if not is_reference else 0.0,
            'attempt_index': attempt_num,
            'is_reference': is_reference
        }
        
        # Add all score types
        for score_name, score_value in all_scores.items():
            record[f'score_{score_name}'] = score_value
            
        return record
        
    def run_analysis(self):
        """Execute the complete perturbation analysis."""
        print("Perturbation Generation + Multi-Score Analysis")
        print("=" * 60)
        
        self.initialize_system()
        
        # ALWAYS start with reference structure as frame 0
        ref_scores = self.calculate_scores(self.ref_coords)
        ref_rmsd_metrics = self.calculate_rmsd_metrics(self.ref_coords)
        ref_record = self.create_result_record(
            self.ref_coords, ref_scores, ref_rmsd_metrics, 
            intensity=0.0, attempt_num=0, is_reference=True
        )
        self.results = [ref_record]
        
        self._print_reference_info(ref_scores)
        
        # Generate perturbations
        attempts = 0
        max_attempts = self.n_target * 6
        
        while len(self.results) < (self.n_target + 1) and attempts < max_attempts:  # +1 for reference
            attempts += 1
            intensity = np.random.uniform(*self.jitter_range)
            
            # Generate and score perturbation
            pert_coords = self.generate_perturbation(intensity)
            all_scores = self.calculate_scores(pert_coords)
            
            if not self.is_structure_acceptable(all_scores):
                continue
                
            rmsd_metrics = self.calculate_rmsd_metrics(pert_coords)
            record = self.create_result_record(
                pert_coords, all_scores, rmsd_metrics, intensity, attempts
            )
            self.results.append(record)
            
            if len(self.results) % 10 == 0:  # Log every 10th acceptance
                self._log_acceptance(record, len(self.results)-1, attempts)  # -1 to exclude reference
                
        self._print_final_summary(attempts)
        
    def _print_reference_info(self, ref_scores):
        """Print information about the reference system."""
        print(f"Reference Scores:")
        print(f"  CCC+ExVol: {ref_scores['ccc_exvol']:.2f}")
        print(f"  Pair:      {ref_scores['pair']:.2f}")
        print(f"  Tetramer:  {ref_scores['tetramer']:.2f}")
        print(f"  Octet:     {ref_scores['octet']:.2f}")
        print(f"Settings: overlap_tol={self.overlap_tolerance} max_total={self.max_total_score}")
        print(f"Jitter std range: {self.jitter_range} (Å) incremental_steps={self.incremental_steps}\n")
        
    def _log_acceptance(self, record, n_accepted, attempt_num):
        """Log details of an accepted perturbation."""
        print(f"Accepted {n_accepted:03d}/{self.n_target} | attempt {attempt_num:04d} | "
              f"req_int={record['intensity_requested']:5.1f} | rawRMSD={record['rmsd_raw']:6.2f} | "
              f"alignRMSD={record['rmsd_aligned']:6.2f} | "
              f"CCC+Ex={record['score_ccc_exvol']:6.1f} | Pair={record['score_pair']:6.1f} | "
              f"Tet={record['score_tetramer']:6.1f} | Oct={record['score_octet']:6.1f}")
              
    def _print_final_summary(self, total_attempts):
        """Print comprehensive analysis summary."""
        print("\nSummary")
        print("-" * 60)
        
        n_perturbations = len(self.results) - 1  # Exclude reference
        acceptance_rate = n_perturbations / max(1, total_attempts)
        print(f"Reference + {n_perturbations} perturbations = {len(self.results)} total structures")
        print(f"Attempts: {total_attempts}  Acceptance Rate: {acceptance_rate:.2%}")
        
        if len(self.results) <= 1:
            return
            
        # Extract metrics for analysis (excluding reference)
        perturbations_only = self.results[1:]  # Skip reference
        if not perturbations_only:
            return
            
        # Collect all score types
        score_types = ['ccc_exvol', 'pair', 'tetramer', 'octet']
        metrics = {
            'aligned_rmsd': np.array([r['rmsd_aligned'] for r in perturbations_only]),
            'raw_rmsd': np.array([r['rmsd_raw'] for r in perturbations_only])
        }
        
        # Add all score metrics
        for score_type in score_types:
            metrics[f'score_{score_type}'] = np.array([r[f'score_{score_type}'] for r in perturbations_only])
        
        # Print statistical summaries
        print("\nRMSD Statistics:")
        for name in ['aligned_rmsd', 'raw_rmsd']:
            values = metrics[name]
            print(f"  {name.replace('_', ' ').title():12}: mean {values.mean():.2f} ± {values.std():.2f} "
                  f"(min {values.min():.2f}, max {values.max():.2f})")
        
        print("\nScore Statistics:")
        for score_type in score_types:
            values = metrics[f'score_{score_type}']
            print(f"  {score_type.upper():12}: mean {values.mean():.2f} ± {values.std():.2f} "
                  f"(min {values.min():.2f}, max {values.max():.2f})")
                  
        # Print correlations
        if len(perturbations_only) > 1:
            print("\nCorrelations (Aligned RMSD vs Scores):")
            for score_type in score_types:
                correlation = np.corrcoef(metrics['aligned_rmsd'], metrics[f'score_{score_type}'])[0,1]
                print(f"  vs {score_type.upper():12}: {correlation:.3f}")
            
    def save_results(self, output_file="output_data_revised/perturbation_analysis_all_scores.h5"):
        """Save results to HDF5 file with ALL scores."""
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with h5py.File(output_file, 'w') as f:
            # Create groups for coordinate arrays
            groups = {
                'array_A': f.create_group('array_A'),
                'array_B': f.create_group('array_B'), 
                'array_C': f.create_group('array_C'),
                'metadata': f.create_group('metadata')
            }
            
            # Save each result
            for i, result in enumerate(self.results):
                frame_name = f"frame_{i:04d}"
                
                # Save coordinates
                groups['array_A'].create_dataset(frame_name, data=result['coords_A'])
                groups['array_B'].create_dataset(frame_name, data=result['coords_B'])
                groups['array_C'].create_dataset(frame_name, data=result['coords_C'])
                
                # Save metadata
                meta_group = groups['metadata'].create_group(frame_name)
                for key, value in result.items():
                    if not key.startswith("coords_"):
                        if isinstance(value, str):
                            meta_group.attrs[key] = value.encode('utf-8')
                        else:
                            meta_group.attrs[key] = value
                            
            # Save global attributes
            f.attrs['num_frames'] = len(self.results)
            f.attrs['description'] = "Reference + Perturbations with ALL scoring functions".encode('utf-8')
            f.attrs['overlap_tolerance'] = self.overlap_tolerance
            f.attrs['max_total_score'] = self.max_total_score
            f.attrs['has_reference_frame'] = True
            f.attrs['score_types'] = "ccc_exvol,pair,tetramer,octet,ccc_only,exvol_only".encode('utf-8')
            
        print(f"\nSaved H5: {output_file}")


def main():
    """Main execution function."""
    analyzer = PerturbationAnalyzer(
        n_target=400,
        jitter_range=(0.1, 25.0),
        max_total_score=10000.0
    )
    
    analyzer.run_analysis()
    analyzer.save_results()
    print("Done.")


if __name__ == "__main__":
    main()