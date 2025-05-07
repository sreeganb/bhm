from dataclasses import dataclass, field
from typing import Dict, Union
import torch
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

@dataclass
class SystemParametersTorch:
    """
    Defines an 8-fold symmetric tetramer system using PyTorch for coordinate
    calculations and Matplotlib for a true-to-scale 3D sphere visualization.
    """
    box_size: float = 500.0
    radii: Dict[str, float] = field(default_factory=lambda: {'A': 20.0, 'B': 16.0, 'C': 18.0})
    pair_distances: Dict[str, float] = field(default_factory=lambda: {
        'AA': 40.5,
        'AB': 36.5,
        'BC': 34.5,
        'CC': 36.5
    })
    component_counts: Dict[str, int] = field(default_factory=lambda: {'A': 8, 'B': 8, 'C': 16})
    octagon_radius: float = 80.0
    tetramer_vertical_drop: float = 20.0
    overlap_tolerance: float = 0.5
    device: Union[str, torch.device] = field(default_factory=lambda: torch.device('cuda' if torch.cuda.is_available() else 'cpu'))
    ideal_coordinates: Dict[str, torch.Tensor] = field(init=False)

    def __post_init__(self):
        self.device = torch.device(self.device)
        self.ideal_coordinates = self._generate_ideal_coordinates()

    def _generate_ideal_coordinates(self) -> Dict[str, torch.Tensor]:
        # Load radii as tensors
        rA, rB, rC = [torch.tensor(self.radii[k], device=self.device, dtype=torch.float32) for k in ['A', 'B', 'C']]
        N = self.component_counts['A']
        assert self.component_counts['B'] == N and self.component_counts['C'] == 2 * N, \
            "Component counts must satisfy B=A and C=2A."

        # Compute distances
        AB = rA + rB + self.overlap_tolerance
        BC = rB + rC + self.overlap_tolerance
        drop = torch.tensor(self.tetramer_vertical_drop, device=self.device, dtype=torch.float32)
        assert drop < BC, f"Vertical drop {drop} must be less than BC {BC}."

        delta = torch.sqrt(BC**2 - drop**2)
        angles = torch.arange(N, device=self.device, dtype=torch.float32) * (2 * torch.pi / N)

        # Positions for A and B
        xy = torch.stack((torch.cos(angles), torch.sin(angles)), dim=1) * self.octagon_radius
        A_pos = torch.cat((xy, torch.zeros(N, 1, device=self.device)), dim=1)
        B_pos = torch.cat((xy, -AB * torch.ones(N, 1, device=self.device)), dim=1)

        # Prepare for C placement optimization
        phis = torch.linspace(0, torch.pi/2, 50, device=self.device)
        base_tangent = delta * torch.stack((-torch.sin(angles), torch.cos(angles)), dim=1)
        zC = -AB - drop
        desired_min = 2 * rC + self.overlap_tolerance

        best_phi = None
        best_min_dist = -1.0

        # Search best phi
        for phi in phis:
            rot = torch.tensor([[torch.cos(phi), -torch.sin(phi)], [torch.sin(phi), torch.cos(phi)]], device=self.device)
            v = (rot @ base_tangent.T).T
            C1 = torch.cat((xy + v, zC.expand(N, 1)), dim=1)
            C2 = torch.cat((xy - v, zC.expand(N, 1)), dim=1)

            interleave = torch.empty((2 * N, 3), device=self.device)
            interleave[0::2] = C1
            interleave[1::2] = C2

            # Mask self-distances
            d = torch.cdist(interleave, interleave)
            mask = ~torch.eye(2*N, device=self.device, dtype=torch.bool)
            min_dist = d[mask].min()

            if min_dist >= desired_min and min_dist > best_min_dist:
                best_min_dist = min_dist
                best_phi = phi

        assert best_phi is not None, "No valid phi found; consider adjusting parameters."

        # Final C coordinates
        rot = torch.tensor([[torch.cos(best_phi), -torch.sin(best_phi)], [torch.sin(best_phi), torch.cos(best_phi)]], device=self.device)
        v = (rot @ base_tangent.T).T
        C1 = torch.cat((xy + v, zC.expand(N, 1)), dim=1)
        C2 = torch.cat((xy - v, zC.expand(N, 1)), dim=1)
        C_pos = torch.empty((2 * N, 3), device=self.device)
        C_pos[0::2] = C1
        C_pos[1::2] = C2

        return {'A': A_pos, 'B': B_pos, 'C': C_pos}

    def plot_3d_spheres(self, resolution: int = 20):
        """
        Renders each particle as a true-to-scale sphere in 3D,
        using distinct uniform colors per species.
        """
        coords = {k: v.cpu().numpy() for k, v in self.ideal_coordinates.items()}
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect([1, 1, 1])

        # Sphere mesh parameters
        u = np.linspace(0, 2 * np.pi, resolution)
        v = np.linspace(0, np.pi, resolution)

        # Define uniform colors per label
        color_map = {'A': 'red', 'B': 'green', 'C': 'blue'}

        for label, arr in coords.items():
            # Determine color by first character ('C' covers C1/C2)
            col = color_map[label[0]] if label[0] in color_map else 'gray'
            r = self.radii[label[0]] if label[0] in self.radii else self.radii['C']
            for (x0, y0, z0) in arr:
                X = x0 + r * np.outer(np.cos(u), np.sin(v))
                Y = y0 + r * np.outer(np.sin(u), np.sin(v))
                Z = z0 + r * np.outer(np.ones_like(u), np.cos(v))
                ax.plot_surface(X, Y, Z, color=col, alpha=0.6, linewidth=0)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('True-scale 3D Spheres of Tetramer System')
        plt.show()


# Example Usage:
if __name__ == '__main__':
    try:
        # --- Default Parameters ---
        print("--- Default Parameters ---")
        system_params = SystemParametersTorch()
        print(f"Generated A coordinates shape: {system_params.ideal_coordinates['A'].shape}, "
              f"device: {system_params.ideal_coordinates['A'].device}")
        print(f"Generated B coordinates shape: {system_params.ideal_coordinates['B'].shape}")
        print(f"Generated C coordinates shape: {system_params.ideal_coordinates['C'].shape}")
        print(f"Example A coordinates:\n{system_params.ideal_coordinates['A'][:2]}")
        print(f"Example B coordinates:\n{system_params.ideal_coordinates['B'][:2]}")
        print(f"Example C coordinates:\n{system_params.ideal_coordinates['C'][:4]}")  # first two C‑pairs

        # Show the 3D plot (true‑scale spheres)
        system_params.plot_3d_spheres(resolution=30)

        # --- Testing Potential Vertical Drop Issue (Expect AssertionError) ---
        print("\n--- Testing Bad Vertical Drop (Expect AssertionError) ---")
        try:
            SystemParametersTorch(tetramer_vertical_drop=40.0)
        except AssertionError as e:
            print(f"Caught expected assertion: {e}")

        # --- Testing Potential C–C Overlap (Expect AssertionError) ---
        print("\n--- Testing Small Octagon (Expect AssertionError) ---")
        try:
            SystemParametersTorch(octagon_radius=50.0)
        except AssertionError as e:
            print(f"Caught expected assertion: {e}")

        # --- Testing Intra‑tetramer Overlap (tiny drop) ---
        print("\n--- Testing Tiny Vertical Drop (Expect AssertionError) ---")
        try:
            SystemParametersTorch(tetramer_vertical_drop=1.0)
        except AssertionError as e:
            print(f"Caught expected assertion: {e}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        import traceback; traceback.print_exc()
