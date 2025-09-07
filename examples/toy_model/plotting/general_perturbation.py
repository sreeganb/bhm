import numpy as np
import torch
from parameters import SystemParameters
from typing import Dict
from scipy.spatial.transform import Rotation
from pair_score import ScoringSystem
# This class is used to perturb the coordinates of a system of particles

class PerturbCoords:
    def __init__(self, coords, move_size=1.0):
        self.move_size = move_size
        self.params = SystemParameters()
        self.coords = coords
        self.radii = self.params.radii
        self.ps = ScoringSystem(self.params)
        self.fix_sigma = {'AA' : 1.8, 'AB' : 1.5, 'BC': 2.0, 'CC' : 2.5}
        self.fix_sigma_range = {'AA' : [0.1, 20.0], 'AB' : [0.1, 20.0], 'BC': [0.1, 20.0], 'CC' : [0.1, 20.0]}

    def perturb(self):
        success = False
        while success == False:
            # choose how many particles to move from the set 1, 2 or 3 with probabilities 0.5, 0.3, 0.2
            n_particles = np.random.choice([1, 2, 3], p=[0.8, 0.15, 0.05])
            # randomly select n_particles from the set of particles
            # Get the particle types
            particle_types = list(self.coords.keys())
            # Randomly select n_particles from the particle types
            selected_types = np.random.choice(particle_types, n_particles, replace=False)
            
            moved_particles = []
            old_coords = {}
            
            # Randomly select particles to move
            for particle_type in selected_types:
                coords_for_type = self.coords[particle_type]
                particle_index = np.random.randint(0, len(coords_for_type))
                moved_particles.append((particle_type, particle_index))
                old_coords[(particle_type, particle_index)] = coords_for_type[particle_index].copy()
                
            # Move the selected particles
            for particle_type, particle_index in moved_particles:
                coords_for_type = self.coords[particle_type]
                direction = np.random.choice([-1, 1], size=3)
                distance = np.random.uniform(0, self.move_size)
                coords_for_type[particle_index] += direction * distance
            
            # Check for excluded volume
            overlap = False
            for p1_type, p1_index in moved_particles:
                coords_1 = self.coords[p1_type][p1_index]
                r1 = self.radii[p1_type]
                
                for other_type in self.coords.keys():
                    other_coords = self.coords[other_type]
                    for i in range(len(other_coords)):
                        r2 = self.radii[other_type]
                        dist = np.linalg.norm(coords_1 - other_coords[i])
                        
                        # Check for overlap, excluding the particle itself
                        if dist < (r1 + r2):
                            if other_type == p1_type and i == p1_index:
                                continue  # Skip self-comparison
                            
                            # this means that the particles are overlapping
                            # so go back to the original coordinates and continue the MCMC
                            # move the particle back
                            overlap = True
                            break
                    if overlap:
                        break
                if overlap:
                    break
            
            # If overlap, revert all moved particles to their original positions
            if overlap:
                for particle_type, particle_index in moved_particles:
                    self.coords[particle_type][particle_index] = old_coords[(particle_type, particle_index)]
            else:
                success = True
        return self.coords
    
    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
        keys = sorted(self.params.ideal_coordinates.keys())
        
        # Ensure all components expected in ideal_coordinates are present in perturbed_coords
        # and have the same number of atoms. This is critical.
        for k in keys:
            if k not in perturbed_coords or self.params.ideal_coordinates[k].shape != perturbed_coords[k].shape:
                raise ValueError(f"Component {k} shape mismatch or missing in perturbed_coords for RMSD.")

        ideal_coords_flat = np.concatenate([self.params.ideal_coordinates[k] for k in keys])
        perturbed_coords_flat = np.concatenate([perturbed_coords[k] for k in keys])

        if ideal_coords_flat.shape[0] == 0: # Handle empty structures
             return 0.0 if perturbed_coords_flat.shape[0] == 0 else np.inf


        ideal_centroid = np.mean(ideal_coords_flat, axis=0)
        perturbed_centroid = np.mean(perturbed_coords_flat, axis=0)
        ideal_centered = ideal_coords_flat - ideal_centroid
        perturbed_centered = perturbed_coords_flat - perturbed_centroid
        
        try:
            # Rotation.align_vectors returns rotation and RMSD
            # The old script calculated RMSD manually after alignment, this is more direct
            rot, rmsd_val = Rotation.align_vectors(perturbed_centered, ideal_centered)
            return rmsd_val 
        except Exception as e: # Catch potential errors in alignment (e.g. too few points, collinear points)
            print(f"Error during Rotation.align_vectors: {e}. Calculating non-aligned RMSD as fallback.")
            return np.sqrt(np.mean(np.sum((ideal_centered - perturbed_centered) ** 2, axis=1)))
        
    def score_calculator(self, coords):
        # Calculate the score for the perturbed coordinates
        # This is a placeholder function, replace with actual scoring logic
        total_score, exclusion_score, pairwise_score, prior_penalty = self.ps.calculate_score(coords, 
                                                                                              self.fix_sigma, 
                                                                                              self.fix_sigma_range)
        return total_score, exclusion_score, pairwise_score, prior_penalty
        
    
# Example usage
if __name__ == "__main__":
    # Initially pass the ideal coordinates for iteration 1, after that keep perturbing the structure that got 
    # returned
    n_perturbations = 20
    score_list = []
    rmsd_list = []
    coords = SystemParameters().latest_ideal()
    perturb = PerturbCoords(coords)
    for i in range(n_perturbations):
        if i == 0:
            # Use the ideal coordinates for the first iteration
            #coords = SystemParameters().ideal_coordinates
            # calculate score and rmsd for the ideal coordinates
            total_score, exclusion_score, pairwise_score, prior_penalty = perturb.score_calculator(coords)
            rmsd = perturb.calculate_rmsd(coords)
            score_list.append(total_score)
            rmsd_list.append(rmsd)
            print("ideal coordinates")
            print(f"Iteration {i+1}: Score: {total_score}, RMSD: {rmsd}")
        else:
            # Use the perturbed coordinates from the previous iteration
            coords = new_coords
        # perturb the coordinates
        new_coords = perturb.perturb()
        # calculate score and rmsd for the perturbed coordinates
        total_score, exclusion_score, pairwise_score, prior_penalty = perturb.score_calculator(new_coords)
        rmsd = perturb.calculate_rmsd(new_coords)
        score_list.append(total_score)
        rmsd_list.append(rmsd)
        print(f"Iteration {i+1}: Score: {total_score}, RMSD: {rmsd}")