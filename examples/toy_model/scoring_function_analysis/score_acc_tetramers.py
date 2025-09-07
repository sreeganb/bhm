from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
from scipy.spatial.transform import Rotation
from toy_model.parameters import SystemParameters
from toy_model.pair_sampler import PairSampler  # Assuming PairSampler is defined in pair_sampler.py
from toy_model.tetramer_sampler import TetramerSampler  # Assuming TetramerSampler is defined in tetramer_sampler.py
from toy_model.base_sampler import BaseMCSampler  # Assuming BaseMCSampler is defined in base_sampler.py

@dataclass
class PerturbedTetramerAnalysis:
    def __init__(self):
        """Initialize with system parameters and precompute constants."""
        params = SystemParameters()
        print("reached after SystemParameters")
        self.ideal_coordinates = params.ideal_coordinates
        self.radii = params.radii
        self.pair_distances = params.pair_distances
        self.bs = BaseMCSampler()  # Initialize BaseMCSampler for excluded volume calculation
        # Fixed sigma values for scoring consistency
        self.sig = {'AA': 0.5, 'AB': 0.35, 'BC': 0.25, 'CC': 0.4}
        self.sigma_range = {'AA': (0.1, 2.0), 'AB': (0.1, 2.0), 'BC': (0.1, 2.0), 'CC': (0.1, 2.0)}
        # Precompute minimum distances for overlap resolution
        self.min_dists = {
            tuple(sorted((comp1, comp2))): self.radii[comp1] + self.radii[comp2]
            for comp1 in self.radii for comp2 in self.radii
        }
        self.t_sampler = TetramerSampler(use_sigma_distribution = False, 
                                         positions_ts = self.ideal_coordinates, 
                                         sig_passed = self.sig, 
                                         sig_range_passed = self.sigma_range)  # Initialize TetramerSampler
        self.ideal_tetramers = self.t_sampler.get_tetramers(self.ideal_coordinates)
        # Initialize PairSampler for pair scoring
        self.ps = PairSampler(use_def_sig_pos = False, sig_passed = self.sig, sig_range_passed = self.sigma_range,
                              pos_passed = self.ideal_coordinates)  # Assumes default initialization is valid
        self.use_sigma_distribution = False  # Use fixed sigma
        # Precompute ideal score - score of the ideal coordinates and IDEAL tetramers
        self.ideal_score, ex_ideal, pair_ideal, tet_ideal = self.neg_log_posterior(self.ideal_coordinates,
                                                           self.ideal_tetramers, sig=self.sig)
        print(f"Ideal Score: {self.ideal_score}", f"ex_ideal: {ex_ideal}", f"pair_ideal: {pair_ideal}", f"tet_ideal: {tet_ideal}")
        print(f"Ideal Score: {self.ideal_score}")
        # from the ideal coordinates, I want to parse out A-1, B-1, C-0 and C-1 coordinates
        # pick out the coordinates of A-1, B-1, C-0, C-1 from the ideal coordinates
        self.ideal_tetramer_coords = {
        'A': np.array([self.ideal_coordinates['A'][0]]),      # First A particle
        'B': np.array([self.ideal_coordinates['B'][0]]),      # First B particle
        'C': np.array([self.ideal_coordinates['C'][0],        # First C particle
                      self.ideal_coordinates['C'][1]])        # Second C particle
        }
        print("ideal_tetramer_coords :", self.ideal_tetramer_coords)
        
        self.ideal_rmsd = 0.0
        self.ideal_accuracy = 1.0

    def perturb_ideal_coordinates(self, magnitude: float = 1.0, seed: int = None) -> Tuple[Dict[str, np.ndarray], bool]:
        if seed is not None:
            np.random.seed(seed)
        perturbed_coords = {
            comp: coords + np.random.uniform(-magnitude, magnitude, size=coords.shape)
            for comp, coords in self.ideal_coordinates.items()
        }
        return self.resolve_overlaps_vectorized(perturbed_coords)

    def resolve_overlaps_vectorized(self, coords: Dict[str, np.ndarray]) -> Tuple[Dict[str, np.ndarray], bool]:
        """Resolve overlaps between components using a vectorized approach.

        Returns:
            Tuple[Dict[str, np.ndarray], bool]: Resolved coordinates and a boolean indicating if overlaps still persist.
        """
        max_iterations = 5000
        overlap_exists = False
        for _ in range(max_iterations):
            total_displacement = {comp: np.zeros_like(coords[comp]) for comp in coords}
            overlap_exists_in_iteration = False # Reset for each iteration
            for comp1 in coords:
                for comp2 in coords:
                    if comp1 > comp2:
                        continue
                    coords1, coords2 = coords[comp1], coords[comp2]
                    same_component = comp1 == comp2
                    coords1_exp = coords1[:, np.newaxis, :]
                    coords2_exp = coords2[np.newaxis, :, :]
                    distances = np.linalg.norm(coords1_exp - coords2_exp, axis=2)
                    min_dist = self.min_dists[tuple(sorted((comp1, comp2)))]
                    if same_component:
                        np.fill_diagonal(distances, np.inf) # Avoid self-overlap for same component
                    overlap_mask = distances < min_dist
                    overlapping_indices = np.where(overlap_mask)
                    if overlapping_indices[0].size > 0:
                        overlap_exists_in_iteration = True
                        directions = (coords1_exp - coords2_exp)[overlapping_indices]
                        directions /= (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-10) # Normalize direction
                        overlap_amounts = (min_dist - distances[overlapping_indices]) / 2 * 0.12 # Displacement amount, scaled down
                        displacements = directions * overlap_amounts[:, np.newaxis]

                        # Corrected np.add.at indexing: use tuple of index arrays, directly use overlapping_indices
                        np.add.at(total_displacement[comp1], overlapping_indices[0], displacements)  
                        np.add.at(total_displacement[comp2], overlapping_indices[1], -displacements)

            if not overlap_exists_in_iteration:
                overlap_exists = False
                break # No overlaps in this iteration, resolution complete
            else:
                overlap_exists = True # Overlaps still exist
            for comp in coords:
                coords[comp] += total_displacement[comp] # Apply accumulated displacements

        return coords, overlap_exists

    def calculate_rmsd_tetramer(self, ideal_coords: np.ndarray, aligned_coords: np.ndarray) -> float:
        """
        Compute the RMSD between the ideal and aligned perturbed coordinates for a tetramer.
        
        Parameters:
            ideal_coords (np.ndarray): Ideal coordinates of shape (N, 3)
            aligned_coords (np.ndarray): Aligned perturbed coordinates of shape (N, 3)
            
        Returns:
            float: The RMSD value.
        """
        # Compute the difference for each coordinate
        diff = ideal_coords - aligned_coords
        # Compute the squared distance for each coordinate (summing over x, y, z)
        squared_distances = np.sum(diff ** 2, axis=1)
        # Compute the mean of the squared distances
        mean_squared_distance = np.mean(squared_distances)
        # Return the square root of the mean squared distance (RMSD)
        return np.sqrt(mean_squared_distance)

    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
        """Calculate the total RMSD over tetramers between perturbed and ideal coordinates.
        
        For each tetramer, the ideal and perturbed coordinates are centered,
        aligned using the Kabsch algorithm, and then their RMSD is computed.
        The final RMSD is the sum over all tetramers.
        """
        total_rmsd = 0.0
        tetramers = self.t_sampler.get_tetramers(perturbed_coords)
        
        # Loop over each tetramer defined by indices (a_idx, b_idx, c1_idx, c2_idx)
        for tetramer in tetramers:
            a_idx, b_idx, c1_idx, c2_idx = tetramer
            # Build perturbed tetramer coordinates dictionary
            perturbed_tetramer = {
                'A': np.array([perturbed_coords['A'][a_idx]]),
                'B': np.array([perturbed_coords['B'][b_idx]]),
                'C': np.array([perturbed_coords['C'][c1_idx], perturbed_coords['C'][c2_idx]])
            }
            
            # Get the ideal tetramer coordinates.
            # It is assumed that self.ideal_tetramer_coords is a dictionary with keys 'A', 'B', and 'C'
            ideal_tetramer = self.ideal_tetramer_coords
            
            # Flatten the coordinates in a consistent order by sorting the keys.
            ideal_flat = np.concatenate([ideal_tetramer[k] for k in sorted(ideal_tetramer.keys())])
            perturbed_flat = np.concatenate([perturbed_tetramer[k] for k in sorted(perturbed_tetramer.keys())])
            
            # Center both sets by subtracting their centroids
            ideal_centroid = np.mean(ideal_flat, axis=0)
            perturbed_centroid = np.mean(perturbed_flat, axis=0)
            ideal_centered = ideal_flat - ideal_centroid
            perturbed_centered = perturbed_flat - perturbed_centroid
            
            # Use the Kabsch algorithm to align perturbed tetramer to ideal tetramer
            rot = Rotation.align_vectors(perturbed_centered, ideal_centered)[0]
            aligned_perturbed = rot.apply(perturbed_centered)
            
            # Compute RMSD for this tetramer (calculate_rmsd_tetramer is assumed to be defined)
            tet_rmsd = self.calculate_rmsd_tetramer(ideal_centered, aligned_perturbed)
            total_rmsd += tet_rmsd

        return total_rmsd


    def calculate_tetramer_score(self, positions: Dict[str, np.ndarray],
                                 tetramer: Tuple[int, ...], sig: Dict[str, float]) -> float:
        """Calculate the tetramer score for a given tetramer configuration."""
        a_idx, b_idx, c1_idx, c2_idx = tetramer
        score = 0.0
        score += self.bs.pair_score_nll(positions['A'][a_idx], positions['B'][b_idx], self.pair_distances['AB'], sig['AB'])
        for c_idx in [c1_idx, c2_idx]:
            score += self.bs.pair_score_nll(positions['B'][b_idx], positions['C'][c_idx], self.pair_distances['BC'], sig['BC'])
        return score

    def neg_log_posterior(
        self,
        positions: Dict[str, np.ndarray],
        tetramers: List[Tuple[int, ...]], # Fixed tetramers from ideal coordinates passed here
        prior_penalty_from_distribution: float = 0.0,
        sig: Dict[str, float] = None,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
    ) -> Tuple[float, float, float, float]:
        """Calculate the negative log posterior score, including pair and tetramer terms."""
        tetramer_pairs = set() # Pairs involved in tetramer interactions
        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
            tetramer_pairs.add(('A', a_idx, 'B', b_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c1_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c2_idx))
            tetramer_pairs.add(('C', c1_idx, 'C', c2_idx))
        score = 0.0; ex_score = 0.0; pair_score = 0.0
        # Calculate pair score using PairSampler, excluding tetramer-internal pairs
        score, ex_score, pair_score, _ = self.ps.calculate_score(
            positions,
            sig,
            self.sigma_range,
            tetramer_pairs,
            self.use_sigma_distribution,
            prior_penalty_from_distribution
        )
        total_tet_score = 0.0
        # Sum up tetramer scores for all tetramers
        for tetramer in tetramers:
            tet_score = self.calculate_tetramer_score(positions, tetramer, sig)
            total_tet_score += tet_score
        score += tetramer_weight * total_tet_score # Add tetramer score to total score
        return score, ex_score, pair_score, tetramer_weight * total_tet_score

    def analyze_perturbations(self, magnitudes: np.ndarray, num_perturbations: int = 10, seed: int = None, max_retries: int = 5):
        """Analyze the effect of perturbations across a range of magnitudes with retry for overlap resolution."""
        results = []
        current_seed = seed
        print("Starting perturbation analysis...") # Indicate start of analysis
        for mag in magnitudes:
            for i in range(num_perturbations):
                perturbed_coords, _ = self.perturb_ideal_coordinates(magnitude=mag, seed=current_seed)
                retries = 0
                resolved_coords = None
                overlap_persists = True  # Initialize to True to enter the loop
                while overlap_persists and retries <= max_retries:
                    resolved_coords, overlap_persists = self.resolve_overlaps_vectorized(perturbed_coords)
                    if overlap_persists:
                        print(f"Warning: Overlaps persist after resolution (Retry {retries+1}/{max_retries}). Re-perturbing...")
                        current_seed += 1  # Increment seed to get a new perturbation
                        perturbed_coords, _ = self.perturb_ideal_coordinates(magnitude=mag, seed=current_seed)  # Re-perturb
                    retries += 1

                if overlap_persists:
                    print(f"Warning: Overlaps may persist after {max_retries+1} attempts for magnitude {mag:.2f}, perturbation {i}. Proceeding with potentially unresolved coordinates.")
                else:
                    print(f"Overlaps resolved after {retries-1} retries (if any) for magnitude {mag:.2f}, perturbation {i}.")


                if resolved_coords is not None: # Use resolved_coords if resolution was attempted
                    score, _, _, _ = self.neg_log_posterior(resolved_coords, self.ideal_tetramers, sig=self.sig)
                    rmsd = self.calculate_rmsd(resolved_coords)
                    accuracy = 1.0 / (1 + rmsd)
                    results.append({'Magnitude': mag, 'Score': score, 'RMSD': rmsd, 'Accuracy': accuracy, 'Perturbation_Number': i, 'Retries': retries-1 if not overlap_persists else max_retries+1})
                else: # Fallback in case of issue, though resolved_coords should always be assigned
                    print("Error: Resolved coordinates are None. This should not happen. Skipping this perturbation.")


                current_seed += 1 # Increment seed for next perturbation

        df_results = pd.DataFrame(results)
        summary_stats = df_results.groupby('Magnitude').agg(
            mean_score=('Score', 'mean'), std_score=('Score', 'std'),
            mean_rmsd=('RMSD', 'mean'), std_rmsd=('RMSD', 'std'),
            mean_accuracy=('Accuracy', 'mean'), std_accuracy=('Accuracy', 'std'),
            mean_retries=('Retries', 'mean'), std_retries=('Retries', 'std')
        ).reset_index()

        # Print summary statistics for each magnitude
        print("\nSummary Statistics:")
        for _, row in summary_stats.iterrows():
            print(f"\nMagnitude {row['Magnitude']:.2f}:")
            print(f"Mean Score: {row['mean_score']:.2f} ± {row['std_score']:.2f}")
            print(f"Mean RMSD: {row['mean_rmsd']:.4f} ± {row['std_rmsd']:.4f}")
            print(f"Mean Accuracy: {row['mean_accuracy']:.4f} ± {row['std_accuracy']:.4f}")
            print(f"Mean Retries: {row['mean_retries']:.2f} ± {row['std_retries']:.2f}")


        # Visualization using seaborn and matplotlib
        sns.set_theme(style='whitegrid')
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=df_results, x='RMSD', y='Score', hue='Magnitude', size='Magnitude', alpha=0.6)
        plt.scatter([self.ideal_rmsd], [self.ideal_score], color='red', label='Ideal', s=100) # Mark ideal point
        plt.xlabel('RMSD')
        plt.ylabel('Tetramer Score')
        plt.title('RMSD vs. Tetramer Score for Perturbed Coordinates')
        plt.legend(title='Magnitude')
        plt.tight_layout()

        # Output directory for saving results
        output_dir = 'output_analysis'
        os.makedirs(output_dir, exist_ok=True)
        # Save detailed results and summary statistics to CSV
        df_results.to_csv(os.path.join(output_dir, 'detailed_results.csv'), index=False)
        summary_stats.to_csv(os.path.join(output_dir, 'summary_stats.csv'), index=False)
        # Save the plot as a PNG image
        plt.savefig(os.path.join(output_dir, 'rmsd_vs_score.png'), dpi=300)
        plt.show()
        print("Perturbation analysis complete. Results saved to 'output_analysis/'") # Indicate end of analysis

if __name__ == "__main__":
    analyzer = PerturbedTetramerAnalysis()
    magnitudes = np.linspace(0.01, 8.0, 30)
    analyzer.analyze_perturbations(magnitudes, num_perturbations=5, seed=42, max_retries=8)
