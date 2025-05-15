import numpy as np
import torch
from parameters import SystemParameters

class PerturbCoords:
    def __init__(self, coords, move_size=0.1):
        self.coords = coords
        self.move_size = move_size
        self.params = SystemParameters()
        self.radii = self.params.radii

    def perturb(self):
        success = False
        while success == False:
            # Randomly select a particle type
            particle_type = np.random.choice(list(self.coords.keys()))
            # Get the coordinates for the selected particle type
            coords_for_type = self.coords[particle_type]
            # Randomly select a particle within that type to move
            particle_index = np.random.randint(0, len(coords_for_type))
            # Randomly select a direction to move in
            direction = np.random.choice([-1, 1], size=3)
            # Randomly select a distance to move
            distance = np.random.uniform(0, self.move_size)
            # Move the particle
            old_coords = coords_for_type[particle_index].copy()
            coords_for_type[particle_index] += direction * distance
            # Check for excluded volume
            overlap = False
            for other_type, other_coords in self.coords.items():
                for i in range(len(other_coords)):
                    dist = np.linalg.norm(coords_for_type[particle_index] - other_coords[i])
                    # Get radii for the particle types
                    r1 = self.radii[particle_type]
                    r2 = self.radii[other_type]
                    if dist < (r1 + r2) and (particle_type != other_type or i != particle_index):
                        # this means that the particles are overlapping
                        # so go back to the original coordinates and continue the MCMC
                        # move the particle back
                        coords_for_type[particle_index] = old_coords
                        overlap = True
                        break
                if overlap:
                    break
            if not overlap:
                success = True
        return self.coords
    
    def calculate_rmsd(self, coords):
        """
        ideal structure from parameters. this is a dictionary"""
    
# Example usage
if __name__ == "__main__":
    coords = {
        'A': np.array([
            [63.  ,   0.  ,   0.  ],
            [44.55,  44.55,   0.  ],
            [ 0.  ,  63.  ,   0.  ],
            [-44.55,  44.55,   0.  ],
            [-63.  ,   0.  ,   0.  ],
            [-44.55, -44.55,   0.  ],
            [ -0.  , -63.  ,   0.  ],
            [44.55, -44.55,   0.  ]
        ]),
        'B': np.array([
            [63.  ,   0.  , -38.5 ],
            [44.55,  44.55, -38.5 ],
            [ 0.  ,  63.  , -38.5 ],
            [-44.55,  44.55, -38.5 ],
            [-63.  ,   0.  , -38.5 ],
            [-44.55, -44.55, -38.5 ],
            [ -0.  , -63.  , -38.5 ],
            [44.55, -44.55, -38.5 ]
        ]),
        'C': np.array([
            [ 47.00,   0.00, -68.50],
            [ 79.00,   0.00, -68.50],
            [ 55.86,  55.86, -68.50],
            [ 33.23,  33.23, -68.50],
            [  0.00,  47.00, -68.50],
            [  0.00,  79.00, -68.50],
            [-55.86,  55.86, -68.50],
            [-33.23,  33.23, -68.50],
            [-47.00,   0.00, -68.50],
            [-79.00,   0.00, -68.50],
            [-55.86, -55.86, -68.50],
            [-33.23, -33.23, -68.50],
            [  0.00, -47.00, -68.50],
            [  0.00, -79.00, -68.50],
            [ 55.86, -55.86, -68.50],
            [ 33.23, -33.23, -68.50],
        ])
    }
    perturb = PerturbCoords(coords)
    new_coords = perturb.perturb()
    print(new_coords)