from dataclasses import dataclass
from typing import Dict, List
import numpy as np
import matplotlib.pyplot as plt
from calculate_score import ScoreCalculator  # Assuming you have your ScoreCalculator defined
from parameters import SystemParameters
from visualization import visualize_3d_configuration
import os

@dataclass
class PerturbSystemParameters:
    def __init__(self):
        params = SystemParameters()
        self.ideal_coordinates = params.ideal_coordinates
        self.radii = params.radii
        self.pair_distances = params.pair_distances
        # pass the ideal structure onto the score calculator function and find out what the score is
        scalc = ScoreCalculator()
        score_1, excl_score, pair_score, prior_score = scalc.calculate_score(self.ideal_coordinates)

    def perturb_ideal_coordinates(self, magnitude=1.0, seed=None):
        """
        Introduces small random perturbations to the ideal coordinates 
        while ensuring no overlaps between particles.

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
            # Generate random perturbations within the specified magnitude
            perturbations = np.random.uniform(
                low=-magnitude, high=magnitude, size=coords.shape
            )
            perturbed_coords[component] = coords + perturbations

        # Check and resolve overlaps
        perturbed_coords = self.resolve_overlaps(perturbed_coords)
#        print("perturbed coordinates: ", perturbed_coords)

        return perturbed_coords

    def resolve_overlaps(self, coords):
            """
            Resolves overlaps between particles by iteratively moving them apart.
            Uses a more efficient approach compared to nested loops.

            Args:
                coords (Dict[str, np.ndarray]): The coordinates to check for overlaps.

            Returns:
                Dict[str, np.ndarray]: The coordinates with overlaps resolved.
            """
            overlap_exists = True
            iterations = 0
            max_iterations = 1000  # Add a maximum iteration limit

            while overlap_exists and iterations < max_iterations:
                overlap_exists = False
                iterations += 1

                for comp1, coords1 in coords.items():
                    for comp2, coords2 in coords.items():
                        # Optimization: Only check pairs where comp1 <= comp2 (lexicographically)
                        # to avoid redundant checks (A-B and B-A) and double-counting.
                        if comp1 > comp2:
                            continue

                        # Create a boolean mask to identify pairs within the same component.
                        same_component = comp1 == comp2

                        # Use broadcasting to calculate distances efficiently.
                        # Expand dimensions for broadcasting:
                        coords1_expanded = coords1[:, np.newaxis, :]  # Shape: (N1, 1, 3)
                        coords2_expanded = coords2[np.newaxis, :, :]  # Shape: (1, N2, 3)

                        distances = np.linalg.norm(coords1_expanded - coords2_expanded, axis=2)  # Shape: (N1, N2)

                        # Define minimum allowed distances based on radii.
                        if same_component:
                            min_distances = 2 * self.radii[comp1]
                            # Create a mask to exclude self-comparisons (i != j).
                            np.fill_diagonal(distances, np.inf)
                        else:
                            min_distances = self.radii[comp1] + self.radii[comp2]

                        # Find overlapping pairs.
                        overlap_mask = distances < min_distances  # Shape: (N1, N2)
                        overlapping_indices = np.where(overlap_mask)

                        if overlapping_indices[0].size > 0:  # Check if any overlaps exist
                            overlap_exists = True
                            
                            # Calculate displacement vectors for overlapping pairs.
                            directions = (coords1_expanded - coords2_expanded)[overlapping_indices]
                            directions = directions / (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-10)

                            # Calculate the amount to move each particle.
                            overlap_amounts = (min_distances - distances[overlapping_indices]) / 2 * 0.1

                            # Create displacement vectors, scaled by the overlap amount.
                            displacements = directions * overlap_amounts[:, np.newaxis]
                            
                            # Update coordinates efficiently using np.add.at.
                            np.add.at(coords1, overlapping_indices[0], displacements)
                            np.add.at(coords2, overlapping_indices[1], -displacements)

            if iterations == max_iterations:
                print("Warning: Maximum iterations reached in resolve_overlaps.")

            return coords

    import numpy as np
    from typing import Dict

    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
        """
        Calculate the RMSD between ideal and perturbed coordinates after optimal superposition
        using the quaternion-based Kabsch algorithm.

        Parameters:
            self:  The object containing the ideal_coordinates attribute (typically 'self').
                Must have an attribute `ideal_coordinates`, a dictionary mapping residue/atom keys to (N,3) numpy arrays of 3D coordinates.
            perturbed_coords: Dictionary mapping residue/atom keys to (N, 3) numpy arrays of 3D coordinates.
                            The keys should match those in self.ideal_coordinates.

        Returns:
            RMSD value as a float.
        """

        # Ensure consistent ordering by sorting the keys
        keys = sorted(self.ideal_coordinates.keys())
        ideal_coords_flat = np.concatenate([self.ideal_coordinates[k] for k in keys])
        perturbed_coords_flat = np.concatenate([perturbed_coords[k] for k in keys])

        # Validate that the two sets have the same shape
        if ideal_coords_flat.shape != perturbed_coords_flat.shape:
            raise ValueError("Ideal and perturbed coordinates must have the same shape.")

        # 1. Centering
        ideal_centroid = np.mean(ideal_coords_flat, axis=0)
        perturbed_centroid = np.mean(perturbed_coords_flat, axis=0)
        ideal_centered = ideal_coords_flat - ideal_centroid
        perturbed_centered = perturbed_coords_flat - perturbed_centroid

        # 2. Covariance Matrix
        H = ideal_centered.T @ perturbed_centered

        # 3. Construct the 4x4 Symmetric Matrix (K)
        K = np.array([
            [H[0, 0] + H[1, 1] + H[2, 2], H[1, 2] - H[2, 1], H[2, 0] - H[0, 2], H[0, 1] - H[1, 0]],
            [H[1, 2] - H[2, 1], H[0, 0] - H[1, 1] - H[2, 2], H[0, 1] + H[1, 0], H[2, 0] + H[0, 2]],
            [H[2, 0] - H[0, 2], H[0, 1] + H[1, 0], -H[0, 0] + H[1, 1] - H[2, 2], H[1, 2] + H[2, 1]],
            [H[0, 1] - H[1, 0], H[2, 0] + H[0, 2], H[1, 2] + H[2, 1], -H[0, 0] - H[1, 1] + H[2, 2]]
        ])

        # 4. Eigenvalue and Eigenvector (Quaternion)
        eigenvalues, eigenvectors = np.linalg.eigh(K)
        max_index = np.argmax(eigenvalues)
        optimal_quaternion = eigenvectors[:, max_index]

        # 5. Quaternion to Rotation Matrix
        w, x, y, z = optimal_quaternion
        R = np.array([
            [1 - 2*(y**2 + z**2), 2*(x*y - w*z), 2*(x*z + w*y)],
            [2*(x*y + w*z), 1 - 2*(x**2 + z**2), 2*(y*z - w*x)],
            [2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x**2 + y**2)]
        ])

        # 6. Rotation
        aligned_perturbed = perturbed_centered @ R

        # 7. Calculate RMSD
        rmsd = np.sqrt(np.sum((ideal_centered - aligned_perturbed) ** 2) / ideal_centered.shape[0])
        
        return rmsd

#    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
#        """
#        Calculate the RMSD between the ideal and perturbed coordinates using the Kabsch algorithm.
#        
#        Parameters:
#        perturbed_coords: Dictionary mapping keys to numpy arrays of coordinates.
#                            The keys should match those in self.ideal_coordinates.
#        
#        Returns:
#        RMSD value as a float.
#        """
#        # Ensure consistent ordering by sorting the keys
#        keys = sorted(self.ideal_coordinates.keys())
#        ideal_coords_flat = np.concatenate([self.ideal_coordinates[k] for k in keys])
#        perturbed_coords_flat = np.concatenate([perturbed_coords[k] for k in keys])
#        
#        # Validate that the two sets have the same shape
#        if ideal_coords_flat.shape != perturbed_coords_flat.shape:
#            raise ValueError("Ideal and perturbed coordinates must have the same shape.")
#        
#        # Center the coordinates
#        ideal_centroid = np.mean(ideal_coords_flat, axis=0)
#        perturbed_centroid = np.mean(perturbed_coords_flat, axis=0)
#        ideal_centered = ideal_coords_flat - ideal_centroid
#        perturbed_centered = perturbed_coords_flat - perturbed_centroid
#
#        # Compute the covariance matrix and perform SVD (Kabsch algorithm)
#        H = ideal_centered.T @ perturbed_centered
#        U, S, Vt = np.linalg.svd(H)
#        R = Vt.T @ U.T
#        # Correct for reflection if necessary
#        if np.linalg.det(R) < 0:
#            Vt[-1, :] *= -1
#            R = Vt.T @ U.T
#
#        # Apply the rotation to the perturbed coordinates
#        aligned_perturbed = perturbed_centered @ R
#
#        # Compute RMSD using the Frobenius norm
#        diff = ideal_centered - aligned_perturbed
#        rmsd = np.linalg.norm(diff, ord='fro') / np.sqrt(ideal_centered.shape[0])
#        return rmsd

    def analyze_perturbations(self, magnitudes, num_perturbations=10, seed=None):
        fixed_sigma = {'AA': 0.2025, 'AB': 0.5, 'BC': 0.5, 'CC': 0.2822}
        fix_sigma = True
        scores = ScoreCalculator()
        
        ideal_score, _, _, _ = scores.calculate_score(self.ideal_coordinates)
        ideal_accuracy = 1.0  # since ideal_rmsd is 0
        ideal_rmsd = 0.0

        all_accuracies = []
        all_rmsds = []
        all_scores = []
        
        # Cache ideal flattened coordinates outside loop for RMSD calculation
        ideal_coords_flat = np.concatenate(list(self.ideal_coordinates.values()))
        
        current_seed = seed  # only one copy is needed
        for mag in magnitudes:
            perturbation_scores = []
            perturbation_accuracies = []
            for _ in range(num_perturbations):
                perturbed_coords = self.perturb_ideal_coordinates(magnitude=mag, seed=current_seed)
                if current_seed is not None:
                    current_seed += 1

                score, excl_score, pair_score, prior_score = scores.calculate_score(perturbed_coords)
                # Compute flattened coordinates for RMSD calculation
                rmsd = self.calculate_rmsd(perturbed_coords)
                accuracy = 1.0 / (1 + rmsd)
        
                perturbation_scores.append(score)
                perturbation_accuracies.append(accuracy)
                all_scores.append(score)
                all_accuracies.append(accuracy)
                all_rmsds.append(rmsd)
        
            print(f"\n--- Magnitude {mag:.2f} ---")
            print("Scores:", perturbation_scores)
            print("Accuracies:", [f"{a:.4f}" for a in perturbation_accuracies])
            print(f"Mean Score: {np.mean(perturbation_scores):.2f} ± {np.std(perturbation_scores):.2f}")
            print(f"Mean Accuracy: {np.mean(perturbation_accuracies):.4f} ± {np.std(perturbation_accuracies):.4f}")
        
        # Plotting
        import seaborn as sns

        # Set seaborn style for publication-quality plots.
        sns.set_theme(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=1.2)

        # Define figure dimensions (width in inches; height based on golden ratio).
        fig_width = 10         # example width in inches
        phi = (1 + np.sqrt(5)) / 2  # golden ratio ≈ 1.618
        fig_height = fig_width / phi

        plt.figure(figsize=(fig_width, fig_height))

        # Plot perturbed structures.
#        plt.scatter(all_accuracies, all_scores, color='blue', marker='o', s=80,
#                    edgecolor='black', facecolors='blue', label='Perturbed Structures')
        plt.scatter(all_rmsds, all_scores, color='blue', marker='o', s=80,
                    edgecolor='black', facecolors='blue', label='Perturbed Structures')

        # Plot ideal structure.
#        plt.scatter([ideal_accuracy], [ideal_score], color='lime', marker='^', s=150,
#                    edgecolor='black', label='Ideal Structure')
        plt.scatter([ideal_rmsd], [ideal_score], color='lime', marker='^', s=150,
                    edgecolor='black', label='Ideal Structure')

#        plt.xlabel('Accuracy (1 / (1 + RMSD))', fontsize=16)
        plt.xlabel('RMSD', fontsize=16)
        plt.ylabel('Score', fontsize=16)
#        plt.title('Structural Accuracy vs. Model Score', fontsize=18)
        plt.title('RMSD vs. Model Score', fontsize=18)
#        plt.xlim(0, 30.1)
#        plt.ylim(0, 20000)
        plt.grid(True, alpha=0.5)
        plt.legend(fontsize=14, loc='best')
        plt.tight_layout()
        # check if output_analysis exists, if not create it, and then save the figure inside another 
        # folder called score_vs_accuracy
        if not os.path.exists('output_analysis'):
            os.makedirs('output_analysis')
        if not os.path.exists('output_analysis/accuracy_vs_score'):
            os.makedirs('output_analysis/accuracy_vs_score')
        # write to a csv file the rmsd and the score values
        import pandas as pd
        df = pd.DataFrame({'RMSD': all_rmsds, 'Score': all_scores})
        df.to_csv('output_analysis/accuracy_vs_score/rmsd_vs_score.csv', index=False)
        # save the figure as a png file
#        plt.savefig('output_analysis/accuracy_vs_score/accuracy_vs_score.png', dpi=300)
        plt.savefig('output_analysis/accuracy_vs_score/rmsd_vs_score.png', dpi=300)
        plt.show()

# Example usage remains the same
if __name__ == "__main__":
    params = SystemParameters()
    analyze = PerturbSystemParameters()
    magnitudes = np.linspace(0.001, 6.0, 400)
    analyze.analyze_perturbations(magnitudes, num_perturbations=1, seed=52)
