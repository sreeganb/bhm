from dataclasses import dataclass
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
#from generate_coordinates import find_min_distance, scale_for_zero_overlap

@dataclass
class SystemParameters:
    box_size: float = 600.0
    radii: Dict[str, float] = None
    pair_distances: Dict[str, float] = None
    component_counts: Dict[str, int] = None
    ideal_coordinates: Dict[str, np.ndarray] = None
    half_ideal_coordinates: Dict[str, np.ndarray] = None
    quarter_ideal_coordinates: Dict[str, np.ndarray] = None
    single_tetramer_ideal_coordinates: Dict[str, np.ndarray] = None

    def __post_init__(self):
        if self.radii is None:
            self.radii = {'A': 24.0, 'B': 14.0, 'C': 16.0}
        if self.pair_distances is None:
            self.pair_distances = {
                'AA': 48.22,
                'AB': 38.5,
                'BC': 34.0,  # minimum allowed C-C distance
                'CC': 32.0   # note: with the following placement BC bonds won’t be exactly this long!
            }
#        if self.radii is None:
#            self.radii = {'A': 40.0, 'B': 10.0, 'C': 16.0}
#        if self.pair_distances is None:
#            self.pair_distances = {
#                'AA': 80.0,
#                'AB': 50.1,
#                'BC': 26.1,  # minimum allowed C-C distance
#                'CC': 32.5   # note: with the following placement BC bonds won’t be exactly this long!
#            }
        if self.component_counts is None:
            self.component_counts = {'A': 8, 'B': 8, 'C': 16}
        if self.ideal_coordinates is None:
            #self.ideal_coordinates = self.new_generate_ideal_coordinates()
            self.ideal_coordinates = self.latest_ideal()
            #self.ideal_coordinates = self.torch_ideal_coords()

    def new_generate_ideal_coordinates(self, tol: float = 0.5) -> Dict[str, np.ndarray]:
        """
        Generates the ideal coordinates for 8 tetramers in an 8-fold symmetric structure.
        Each tetramer consists of:
        - A: An A particle on an octagon in the plane z=0.
        - B: A B particle directly below A with center-to-center distance AB = rA + rB + tol.
        - C: Two C particles placed below B at z = -(AB_distance + vertical_drop) 
                such that the B-to-C center separation is rB + rC + tol.
                The two C particles are initially placed along the tangent direction at B;
                then they are rotated about the vertical axis (passing through the tetramer’s center)
                by a common angle φ chosen to maximize separation between C particles of adjacent tetramers,
                ensuring that no two particles in the entire system overlap.
        
        Each tetramer’s four particles (A, B, C1, C2) will have their centers collinear
        in (x,y) so that the center-of-mass of the tetramer lies on a vertical line.
        Overlap is checked using the minimum allowed separation equal to the sum of radii + tol.
        """
        # Particle radii (example values)
        rA = self.radii['A']   # e.g. 20.0
        rB = self.radii['B']   # e.g. 16.0
        rC = self.radii['C']   # e.g. 18.0

        # Number of tetramers (placed at the vertices of an octagon)
        N = 8  
        R_A = 65.5   # radius at which A and B are placed in the xy-plane

        # Set distances based on non-overlap plus tolerance.
        AB_distance = rA + rB + tol      # e.g. 20 + 16 + 0.5 = 36.5
        BC_distance = rB + rC + tol      # e.g. 16 + 18 + 0.5 = 34.5

        # Vertical drop between B and C centers.
        vertical_drop = 25.5  # chosen value; adjust as needed.
        # Compute horizontal displacement delta such that:
        # sqrt(delta^2 + vertical_drop^2) = BC_distance
        delta = np.sqrt(BC_distance**2 - vertical_drop**2)

        # Generate angles for the octagon.
        theta = np.linspace(0, 2 * np.pi, N, endpoint=False)

        # Allocate arrays for A and B (each tetramer has one A and one B).
        array_A = np.zeros((N, 3))
        array_B = np.zeros((N, 3))
        
        # Fill in A and B positions.
        for i, angle in enumerate(theta):
            x = R_A * np.cos(angle)
            y = R_A * np.sin(angle)
            array_A[i, :] = [x, y, 0.0]
            array_B[i, :] = [x, y, -AB_distance]  # B’s center is directly below A

        # We now determine a rotation angle φ (phi) for the placement of the two C particles
        # in each tetramer. In the unrotated (default) configuration, the two C particles are
        # placed symmetrically along the tangent to the circle at B:
        #   v = [delta * (-sin(angle)), delta * cos(angle)]
        # Their preliminary positions would be:
        #   c1 = [x + v_x, y + v_y, z_C] and c2 = [x - v_x, y - v_y, z_C],
        # with z_C = -AB_distance - vertical_drop.
        # However, such placement may result in overlaps between C particles from adjacent tetramers.
        # To avoid this, we find a common rotation angle φ (about the vertical axis through each tetramer’s center)
        # so that the final positions become:
        #   new_v = [cos(φ)*v_x - sin(φ)*v_y, sin(φ)*v_x + cos(φ)*v_y]
        # and then:
        #   c1 = [x + new_v_x, y + new_v_y, z_C],
        #   c2 = [x - new_v_x, y - new_v_y, z_C].
        #
        # We choose φ from a set of candidate values (0 to π/2) that maximizes the minimum inter-tetramer
        # distance between any two C particles, with the requirement that every two C particles (from different tetramers)
        # are at least (2*rC + tol) apart.
        
        candidate_phis = np.linspace(0, np.pi/2, 50)
        best_phi = None
        best_min_dist = -np.inf
        desired_min = 2 * rC + tol  # minimum allowed center-to-center distance for C particles

        # Iterate over candidate phi values.
        for phi in candidate_phis:
            temp_array_C = np.zeros((2 * N, 3))  # temporary storage for C positions
            for i, angle in enumerate(theta):
                # Center of the tetramer (from A/B, same (x,y)); z-coordinate for C particles:
                x = R_A * np.cos(angle)
                y = R_A * np.sin(angle)
                z_C = -AB_distance - vertical_drop
                
                # Original tangent offset vector v at this tetramer.
                v = np.array([delta * (-np.sin(angle)), delta * np.cos(angle)])
                # Rotate v by candidate phi.
                new_v = np.array([
                    np.cos(phi) * v[0] - np.sin(phi) * v[1],
                    np.sin(phi) * v[0] + np.cos(phi) * v[1]
                ])
                # Compute the two C positions for tetramer i.
                c1 = np.array([x + new_v[0], y + new_v[1], z_C])
                c2 = np.array([x - new_v[0], y - new_v[1], z_C])
                temp_array_C[2 * i, :] = c1
                temp_array_C[2 * i + 1, :] = c2

            # Now, check the minimum inter-tetramer distance among C particles.
            min_dist = np.inf
            for i in range(N):
                for j in range(i + 1, N):
                    # Tetramer i: indices 2*i and 2*i+1.
                    # Tetramer j: indices 2*j and 2*j+1.
                    for a in [2 * i, 2 * i + 1]:
                        for b in [2 * j, 2 * j + 1]:
                            d = np.linalg.norm(temp_array_C[a] - temp_array_C[b])
                            if d + 1e-6 < min_dist:
                                min_dist = d
            # If this candidate phi yields a minimum distance above desired and is better than previous candidates, choose it.
            if min_dist >= desired_min and min_dist > best_min_dist:
                best_min_dist = min_dist
                best_phi = phi

        if best_phi is None:
            raise ValueError("Could not find a rotation angle phi for C particles that avoids inter-tetramer overlaps.")
        
        # With the selected best_phi, compute the final array_C.
        array_C = np.zeros((2 * N, 3))
        for i, angle in enumerate(theta):
            x = R_A * np.cos(angle)
            y = R_A * np.sin(angle)
            z_C = -AB_distance - vertical_drop
            v = np.array([delta * (-np.sin(angle)), delta * np.cos(angle)])
            new_v = np.array([
                np.cos(best_phi) * v[0] - np.sin(best_phi) * v[1],
                np.sin(best_phi) * v[0] + np.cos(best_phi) * v[1]
            ])
            c1 = [x + new_v[0], y + new_v[1], z_C]
            c2 = [x - new_v[0], y - new_v[1], z_C]
            array_C[2 * i, :] = c1
            array_C[2 * i + 1, :] = c2

        # Final check: Ensure no overlaps occur within each tetramer.
        for i, angle in enumerate(theta):
            # Tetramer particles: A_i, B_i, C1, C2.
            tetramer = np.array([array_A[i, :], array_B[i, :],
                                array_C[2 * i, :], array_C[2 * i + 1, :]])
            radii_list = [rA, rB, rC, rC]
            for j in range(4):
                for k in range(j + 1, 4):
                    dist = np.linalg.norm(tetramer[j] - tetramer[k])
                    min_allowed = radii_list[j] + radii_list[k] + tol
                    if dist + 1e-6 < min_allowed:
                        raise ValueError(f"Overlap in tetramer {i} between particles {j} and {k}: "
                                        f"distance = {dist:.2f}, expected >= {min_allowed:.2f}")

        # Optionally, one may also check globally that no two particles (from different tetramers)
        # overlap. This can be added if necessary.

        return {'A': array_A, 'B': array_B, 'C': array_C}

    def generate_ideal_coordinates(self) -> Dict[str, np.ndarray]:
        """
        Generates the ideal 8-fold symmetric coordinates.
        """
        array_A = np.array([
            [0, 100, 0],
            [-70.71, 70.71, 0],
            [-100, 0, 0],
            [-70.71, -70.71, 0],
            [0, -100, 0],
            [70.71, -70.71, 0],
            [100, 0, 0],
            [70.71, 70.71, 0]
        ])
        
        array_B = np.array([
            [0, 100, -50],
            [-70.71, 70.71, -50],
            [-100, 0, -50],
            [-70.71, -70.71, -50],
            [0, -100, -50],
            [70.71, -70.71, -50],
            [100, 0, -50],
            [70.71, 70.71, -50]
        ])
        array_C = np.array([
            [11.25, 100, -68],  # pair 1
            [-11.25, 100, -68],
            [-59.46, 70.71, -68],  # pair 2
            [-81.96, 70.71, -68],
            [-100, 11.25, -68],    # pair 3
            [-100, -11.25, -68],
            [-81.96, -70.71, -68], # pair 4
            [-59.46, -70.71, -68],
            [11.25, -100, -68],    # pair 5
            [-11.25, -100, -68],
            [81.96, -70.71, -68],  # pair 6
            [59.46, -70.71, -68],
            [100, 11.25, -68],     # pair 7
            [100, -11.25, -68],
            [81.96, 70.71, -68],   # pair 8
            [59.46, 70.71, -68]
        ])

#        array_A, array_B, array_C = scale_for_zero_overlap(array_A, array_B, array_C)
        # read in numpy arrays from files scaled_A.npy, scaled_B.npy, scaled_C.npy
#        array_A = np.load('scaled_A.npy')
#        array_B = np.load('scaled_B.npy')
#        array_C = np.load('scaled_C.npy')

        return {'A': array_A, 'B': array_B, 'C': array_C}

    def latest_ideal(self) -> Dict[str, np.ndarray]:
        array_A  = np.array([
            [ 63.  ,   0.  ,   0.  ],
            [ 44.55,  44.55,   0.  ],
            [  0.  ,  63.  ,   0.  ],
            [-44.55,  44.55,   0.  ],
            [-63.  ,   0.  ,   0.  ],
            [-44.55, -44.55,   0.  ],
            [ -0.  , -63.  ,   0.  ],
            [ 44.55, -44.55,   0.  ]
            ])
        array_B = np.array([
            [ 63.  ,   0.  , -38.5 ],
            [ 44.55,  44.55, -38.5 ],
            [  0.  ,  63.  , -38.5 ],
            [-44.55,  44.55, -38.5 ],
            [-63.  ,   0.  , -38.5 ],
            [-44.55, -44.55, -38.5 ],
            [ -0.  , -63.  , -38.5 ],
            [ 44.55, -44.55, -38.5 ]
            ]) 
        array_C = np.array([
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
        
        return {'A': array_A, 'B': array_B, 'C': array_C}
        
    def half_rotated_coordinates(self) -> Dict[str, np.ndarray]:
        """Half-rotated system coordinates with 180-degree rotation around y-axis through center of geometry"""
        original = self.latest_ideal()
        
        def rotate_180_y_about_center(coords, center):
            """Rotate coordinates 180 degrees around y-axis through a specified center"""
            # Translate to origin
            coords_centered = coords - center
            
            # 180-degree rotation around y-axis
            rotation_matrix = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
            coords_rotated = coords_centered @ rotation_matrix.T
            
            # Translate back
            return coords_rotated + center
        
        def calculate_min_translation_distance(coords1, coords2, radii1, radii2, safety_margin=19.0):
            """Calculate minimum translation distance to avoid overlaps"""
            min_required = 0.0
            
            for i, (c1, r1) in enumerate(zip(coords1, radii1)):
                for j, (c2, r2) in enumerate(zip(coords2, radii2)):
                    # Current distance
                    current_dist = np.linalg.norm(c1 - c2)  # Full 3D distance
                    # Required distance (sum of radii + safety margin)
                    required_dist = r1 + r2 + safety_margin
                    # Additional translation needed
                    additional_needed = max(0, required_dist - current_dist)
                    min_required = max(min_required, additional_needed)
            
            return min_required
        
        # Calculate center of geometry for the entire system
        all_coords = np.vstack([original['A'], original['B'], original['C']])
        center_of_geometry = np.mean(all_coords, axis=0)
        print(f"Center of geometry: {center_of_geometry}")
        
        # Split coordinates - first half vs second half
        half1_A, half2_A = original['A'][:4], original['A'][4:]  # A particles 0-3 vs 4-7
        half1_B, half2_B = original['B'][:4], original['B'][4:]  # B particles 0-3 vs 4-7
        half1_C, half2_C = original['C'][:8], original['C'][8:]  # C particles 0-7 vs 8-15
        
        # Rotate the ENTIRE second half (A4B4C8 system) by 180 degrees around y-axis through center of geometry
        half2_A_rotated = rotate_180_y_about_center(half2_A, center_of_geometry)
        half2_B_rotated = rotate_180_y_about_center(half2_B, center_of_geometry)
        half2_C_rotated = rotate_180_y_about_center(half2_C, center_of_geometry)
        
        print(f"Rotated second half (A4B4C8 system) by 180 degrees around y-axis through center of geometry")
        
        # Combine all coordinates and radii for overlap checking
        half1_coords = np.vstack([half1_A, half1_B, half1_C])
        half2_coords = np.vstack([half2_A_rotated, half2_B_rotated, half2_C_rotated])
        
        half1_radii = [self.radii['A']] * 4 + [self.radii['B']] * 4 + [self.radii['C']] * 8
        half2_radii = [self.radii['A']] * 4 + [self.radii['B']] * 4 + [self.radii['C']] * 8
        
        # Calculate minimum translation distance to avoid overlaps
        min_translation = calculate_min_translation_distance(
            half1_coords, half2_coords, half1_radii, half2_radii
        )
        
        # Add extra safety margin
        translation_distance = min_translation + 15.0
        
        # For 180-degree rotation around y-axis, translate in x-direction to separate the halves
        translation_vector = np.array([translation_distance, 0.0, 0.0])
        
        # Apply translation to the rotated second half
        half2_A_translated = half2_A_rotated + translation_vector
        half2_B_translated = half2_B_rotated + translation_vector
        half2_C_translated = half2_C_rotated + translation_vector
        
        print(f"Applied translation of {translation_distance:.2f} units in x-direction to avoid overlaps")
        
        return {
            'A': np.vstack([half1_A, half2_A_translated]),
            'B': np.vstack([half1_B, half2_B_translated]),
            'C': np.vstack([half1_C, half2_C_translated])
        }
    
    def single_tetramer_coordinates(self) -> Dict[str, np.ndarray]:
        """Single tetramer (ABCC) coordinates"""
        original = self.latest_ideal()
        return {
            'A': original['A'][0:1],
            'B': original['B'][0:1], 
            'C': original['C'][0:2]
        }
    
    def octamer_coordinates(self) -> Dict[str, np.ndarray]:
        """Octamer coordinates (first two tetramers)"""
        original = self.latest_ideal()
        return {
            'A': original['A'][0:2],
            'B': original['B'][0:2],
            'C': original['C'][0:4]
        }