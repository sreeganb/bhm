from dataclasses import dataclass
from typing import Dict, List
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # Moved import here
from scipy.spatial.transform import Rotation  # For Kabsch Algorithm - more efficient library

from calculate_score import ScoreCalculator  # Assuming you have your ScoreCalculator defined
from parameters import SystemParameters
import os
import pandas as pd  # Moved import here


@dataclass
class PerturbSystemParameters:
    def __init__(self):
        params = SystemParameters()
        self.ideal_coordinates = params.ideal_coordinates
        self.radii = params.radii
        self.pair_distances = params.pair_distances
        self.scalc = ScoreCalculator()  # Initialize ScoreCalculator once
        self.ideal_score, _, _, _ = self.scalc.calculate_score(self.ideal_coordinates) # Calculate ideal score once
        self.ideal_accuracy = 1.0
        self.ideal_rmsd = 0.0
        self.min_dists = {} # Pre-calculate min distances for efficiency
        for comp1 in self.radii:
            for comp2 in self.radii:
                key = tuple(sorted((comp1, comp2))) # Use sorted tuple for consistent key
                self.min_dists[key] = self.radii[comp1] + self.radii[comp2]

    def perturb_ideal_coordinates(self, magnitude=1.0, seed=None):
        """
        Introduces small random perturbations to the ideal coordinates and resolves overlaps efficiently.

        Args:
            magnitude (float): The maximum magnitude of the perturbation in each dimension.
            seed (int, optional): Random seed for reproducibility.

        Returns:
            Dict[str, np.ndarray]: A dictionary of perturbed coordinates.
        """
        if seed is not None:
            np.random.seed(seed)

        perturbed_coords = {}
        for component, coords in self.ideal_coordinates.items():
            perturbations = np.random.uniform(-magnitude, magnitude, size=coords.shape)
            perturbed_coords[component] = coords + perturbations

        return self.resolve_overlaps_vectorized(perturbed_coords) # Use vectorized overlap resolution

    def resolve_overlaps_vectorized(self, coords):
        """
        Vectorized overlap resolution, significantly faster than iterative methods.
        """
        max_iterations = 1000  # Reduced max iterations for efficiency

        for _ in range(max_iterations): # Limited iterations for vectorized approach
            total_displacement = {comp: np.zeros_like(coords[comp]) for comp in coords} # Track total displacement
            overlap_exists = False

            for comp1 in coords:
                for comp2 in coords:
                    if comp1 > comp2:
                        continue

                    coords1, coords2 = coords[comp1], coords[comp2]
                    same_component = (comp1 == comp2)

                    coords1_expanded = coords1[:, np.newaxis, :]
                    coords2_expanded = coords2[np.newaxis, :, :]
                    distances = np.linalg.norm(coords1_expanded - coords2_expanded, axis=2)

                    min_dist = self.min_dists[tuple(sorted((comp1, comp2)))]

                    if same_component:
                        np.fill_diagonal(distances, np.inf)

                    overlap_mask = distances < min_dist
                    overlapping_indices = np.where(overlap_mask)

                    if overlapping_indices[0].size > 0:
                        overlap_exists = True

                        directions = (coords1_expanded - coords2_expanded)[overlapping_indices]
                        directions /= (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-10)
                        overlap_amounts = (min_dist - distances[overlapping_indices]) / 2 * 0.1
                        displacements = directions * overlap_amounts[:, np.newaxis]

                        np.add.at(total_displacement[comp1], overlapping_indices[0], displacements)
                        np.add.at(total_displacement[comp2], overlapping_indices[1], -displacements) # Corrected to '-'


            if not overlap_exists: # Early exit if no overlaps
                break

            for comp in coords:
                coords[comp] += total_displacement[comp] # Apply accumulated displacements in one go

        if overlap_exists:
            print("Warning: Overlaps may still exist after vectorized resolution.") # Indicate if max iterations reached

        return coords


    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
        """
        Calculate RMSD using scipy's Rotation for Kabsch algorithm - more efficient and cleaner.
        """
        keys = sorted(self.ideal_coordinates.keys())
        ideal_coords_flat = np.concatenate([self.ideal_coordinates[k] for k in keys])
        perturbed_coords_flat = np.concatenate([perturbed_coords[k] for k in keys])

        if ideal_coords_flat.shape != perturbed_coords_flat.shape:
            raise ValueError("Ideal and perturbed coordinates must have the same shape.")

        ideal_centroid = np.mean(ideal_coords_flat, axis=0)
        perturbed_centroid = np.mean(perturbed_coords_flat, axis=0)
        ideal_centered = ideal_coords_flat - ideal_centroid
        perturbed_centered = perturbed_coords_flat - perturbed_centroid

        # Use scipy's Rotation for Kabsch Algorithm
        rot = Rotation.align_vectors(perturbed_centered, ideal_centered)[0] # Get rotation from scipy
        aligned_perturbed = rot.apply(perturbed_centered) # Apply rotation

        rmsd = np.sqrt(np.mean(np.sum((ideal_centered - aligned_perturbed) ** 2, axis=1))) # Efficient RMSD calculation
        return rmsd


    def analyze_perturbations(self, magnitudes, num_perturbations=10, seed=None):
        fixed_sigma = {'AA': 0.2025, 'AB': 0.5, 'BC': 0.5, 'CC': 0.2822}
        fix_sigma = True

        all_accuracies = []
        all_rmsds = []
        all_scores = []

        current_seed = seed

        results = [] # Store results for dataframe creation

        for mag in magnitudes:
            for _ in range(num_perturbations):
                perturbed_coords = self.perturb_ideal_coordinates(magnitude=mag, seed=current_seed)
                if current_seed is not None:
                    current_seed += 1

                score, _, _, _ = self.scalc.calculate_score(perturbed_coords)
                rmsd = self.calculate_rmsd(perturbed_coords)
                accuracy = 1.0 / (1 + rmsd)

                results.append({'Magnitude': mag, 'Score': score, 'RMSD': rmsd, 'Accuracy': accuracy}) # Append results

        df_results = pd.DataFrame(results) # Create DataFrame once, outside loops

        # Aggregate results for printing mean and std
        summary_stats = df_results.groupby('Magnitude').agg(
            mean_score=('Score', 'mean'), std_score=('Score', 'std'),
            mean_rmsd=('RMSD', 'mean'), std_rmsd=('RMSD', 'std'),
            mean_accuracy=('Accuracy', 'mean'), std_accuracy=('Accuracy', 'std')
        ).reset_index()


        for _, row in summary_stats.iterrows(): # Iterate over summary DataFrame
            mag = row['Magnitude']
            mean_score, std_score = row['mean_score'], row['std_score']
            mean_accuracy, std_accuracy = row['mean_accuracy'], row['std_accuracy']
            print(f"\n--- Magnitude {mag:.2f} ---")
            print(f"Mean Score: {mean_score:.2f} ± {std_score:.2f}")
            print(f"Mean Accuracy: {mean_accuracy:.4f} ± {std_accuracy:.4f}")


        # Plotting using DataFrame and seaborn
        sns.set_theme(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=1.2)
        fig_width = 10
        phi = (1 + np.sqrt(5)) / 2
        fig_height = fig_width / phi
        plt.figure(figsize=(fig_width, fig_height))

        sns.scatterplot(data=df_results, x='RMSD', y='Score', color='blue', marker='o', s=80,
                        edgecolor='black', facecolors='blue', label='Perturbed Structures')

        plt.scatter([self.ideal_rmsd], [self.ideal_score], color='lime', marker='^', s=150,
                    edgecolor='black', label='Ideal Structure')

        plt.xlabel('RMSD', fontsize=16)
        plt.ylabel('Score', fontsize=16)
        plt.title('RMSD vs. Model Score', fontsize=18)
        plt.grid(True, alpha=0.5)
        plt.legend(fontsize=14, loc='best')
        plt.tight_layout()

        output_dir = 'output_analysis/accuracy_vs_score' # Define output directory once
        os.makedirs(output_dir, exist_ok=True)
        df_results.to_csv(os.path.join(output_dir, 'rmsd_vs_score.csv'), index=False) # Save results DF
        plt.savefig(os.path.join(output_dir, 'rmsd_vs_score.png'), dpi=300)
        plt.show()


if __name__ == "__main__":
    analyze = PerturbSystemParameters()
    magnitudes = np.linspace(0.001, 5.0, 10)
    analyze.analyze_perturbations(magnitudes, num_perturbations=10, seed=52)
