import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from typing import Dict
import os

sns.set_style("whitegrid")

def plot_3d(positions: Dict[str, np.ndarray], 
                   title: str = "3D Molecular System",
                   save_name: str = "plot"):
    """
    Create amazing 3D visualization and save as PNG/PDF.
    
    Args:
        positions: Dict with 'A', 'B', 'C' keys and coordinate arrays
        title: Plot title
        save_name: Filename (without extension)
    """
    # Setup
    radii = {'A': 24, 'B': 14, 'C': 16}
    colors = {'A': '#e74c3c', 'B': '#3498db', 'C': '#2ecc71'}
    
    # Create figure
    fig = plt.figure(figsize=(12, 10), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    
    # Style
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.set_alpha(0.1)
        pane.set_edgecolor('lightgray')
    
    # Draw spheres
    for particle_type, coords in positions.items():
        if len(coords) == 0:
            continue
        radius = radii[particle_type]
        color = colors[particle_type]
        
        for x, y, z in coords:
            u = np.linspace(0, 2*np.pi, 30)
            v = np.linspace(0, np.pi, 30)
            sphere_x = x + radius * np.outer(np.cos(u), np.sin(v))
            sphere_y = y + radius * np.outer(np.sin(u), np.sin(v))
            sphere_z = z + radius * np.outer(np.ones(len(u)), np.cos(v))
            
            ax.plot_surface(sphere_x, sphere_y, sphere_z, 
                          color=color, alpha=0.8, linewidth=0, 
                          antialiased=True, shade=True)
    
    # Formatting
    ax.set_xlabel('X', fontsize=14, fontweight='bold')
    ax.set_ylabel('Y', fontsize=14, fontweight='bold')
    ax.set_zlabel('Z', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
    
    # Legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=colors[pt], markersize=15, 
                                 label=f'{pt} (r={radii[pt]})') 
                      for pt in positions.keys() if len(positions[pt]) > 0]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=12)
    
    # Set limits and view
    all_coords = np.vstack([coords for coords in positions.values() if len(coords) > 0])
    if len(all_coords) > 0:
        padding = max(radii.values()) * 1.5
        for i, dim in enumerate(['x', 'y', 'z']):
            center = np.mean(all_coords[:, i])
            range_val = np.ptp(all_coords[:, i]) / 2 + padding
            getattr(ax, f'set_{dim}lim')(center - range_val, center + range_val)
    
    ax.view_init(elev=25, azim=45)
    plt.tight_layout()
    
    # Save
    os.makedirs(os.path.dirname(save_name) if os.path.dirname(save_name) else '.', exist_ok=True)
    plt.savefig(f"{save_name}.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{save_name}.pdf", bbox_inches='tight', facecolor='white')
    print(f"Saved: {save_name}.png and {save_name}.pdf")
    
    plt.show()
    return fig

# Usage
if __name__ == "__main__":
    from parameters import SystemParameters
    params = SystemParameters()
    plot_3d(params.latest_ideal(), "Original System", "plots/original")