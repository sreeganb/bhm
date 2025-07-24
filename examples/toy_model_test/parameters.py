from dataclasses import dataclass
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
#from generate_coordinates import find_min_distance, scale_for_zero_overlap

@dataclass
class SystemParameters:
    box_size: float = 800.0
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
                'BC': 34.0,   # minimum allowed C-C distance
                #'CC': 32.0   # note: with the following placement BC bonds won’t be exactly this long!
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
            self.component_counts = {'A': 8, 'B': 8, 'C': 14}
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
        
                
#        array_A = np.array([
#            [ 104.63,    0.00,    0.00],
#            [  73.97,   73.97,    0.00],
#            [   0.00,  104.63,    0.00],
#            [ -73.97,   73.97,    0.00],
#            [-104.63,    0.00,    0.00],
#            [ -73.97,  -73.97,    0.00],
#            [   0.00, -104.63,    0.00],
#            [  73.97,  -73.97,    0.00],
#        ])
#
#        array_B = np.array([
#            [ 104.63,    0.00,  -50.20],
#            [  73.97,   73.97,  -50.20],
#            [   0.00,  104.63,  -50.20],
#            [ -73.97,   73.97,  -50.20],
#            [-104.63,    0.00,  -50.20],
#            [ -73.97,  -73.97,  -50.20],
#            [   0.00, -104.63,  -50.20],
#            [  73.97,  -73.97,  -50.20],
#        ])
#
#        array_C = np.array([
#            [ 115.99,   11.39,  -73.20],
#            [  93.28,  -11.39,  -73.20],
#            [  85.28,   85.28,  -73.20],
#            [  62.58,   62.58,  -73.20],
#            [  11.39,  115.99,  -73.20],
#            [ -11.39,   93.28,  -73.20],
#            [ -62.58,   85.28,  -73.20],
#            [ -85.28,   62.58,  -73.20],
#            [-115.99,   11.39,  -73.20],
#            [ -93.28,  -11.39,  -73.20],
#            [ -85.28,  -85.28,  -73.20],
#            [ -62.58,  -62.58,  -73.20],
#            [ -11.39, -115.99,  -73.20],
#            [  11.39,  -93.28,  -73.20],
#            [  62.58,  -85.28,  -73.20],
#            [  85.28,  -62.58,  -73.20],
#        ])
#
        
        return {'A': array_A, 'B': array_B, 'C': array_C}

