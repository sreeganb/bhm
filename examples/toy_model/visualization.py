import numpy as np
from scipy.spatial.distance import cdist
from typing import Dict, List, Tuple
from parameters import SystemParameters
import pickle
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import for 3D plotting


def draw_sphere(ax, x0, y0, z0, radius, color, alpha=0.5, resolution=30):
    """
    Draw a 3D sphere at (x0, y0, z0) with the given radius using parametric surfaces.

    Args:
        ax (Axes3D): A 3D matplotlib axis to plot on.
        x0, y0, z0 (float): Center of the sphere.
        radius (float): Radius of the sphere.
        color (str): Matplotlib-compatible color.
        alpha (float): Transparency of the sphere surface.
        resolution (int): Controls the resolution of the parametric grid.
    """
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    u, v = np.meshgrid(u, v)
    x = x0 + radius * np.cos(u) * np.sin(v)
    y = y0 + radius * np.sin(u) * np.sin(v)
    z = z0 + radius * np.cos(v)
    ax.plot_surface(x, y, z, color=color, alpha=alpha, shade=True)

def visualize_3d_configuration(positions, radii, title="System Configuration"):
    """
    Plots the system's components as true scaled spheres in 3D based on their radii.

    Args:
        positions (Dict[str, np.ndarray]): Positions of each component type.
        radii (Dict[str, float]): Radii of each component type.
        title (str): Title of the plot.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    colors = {'A': 'red', 'B': 'blue', 'C': 'green'}

    for type_name, coords_array in positions.items():
        c = colors.get(type_name, 'gray')
        r = radii[type_name]
        for coords in coords_array:
            draw_sphere(ax, coords[0], coords[1], coords[2], r, c)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    plt.show()
    
def visualize_trajectory(trajectory: List[Dict], box_size: float):
    """
    Animates the trajectory of the Monte Carlo simulation.

    Args:
        trajectory (List[Dict]): The trajectory data.
        box_size (float): The size of the simulation box.
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Initialize empty scatter plots for each component type
    scatter_plots = {}
    for type_name in trajectory[0]['positions'].keys():
        scatter_plots[type_name] = ax.scatter([], [], [], label=type_name)

    # Set plot labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Monte Carlo Trajectory')
    ax.set_xlim(0, box_size)
    ax.set_ylim(0, box_size)
    ax.set_zlim(0, box_size)
    ax.legend()

    def update(frame):
        """Updates the scatter plot for each frame."""
        state = trajectory[frame]
        for type_name, scatter_plot in scatter_plots.items():
            scatter_plot._offsets3d = (
                state['positions'][type_name][:, 0],
                state['positions'][type_name][:, 1],
                state['positions'][type_name][:, 2]
            )
        return scatter_plots.values()  # Return updated plots

    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(trajectory), interval=100, blit=False)
    #plt.show()


from typing import List, Dict
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def visualize_trajectory_plotly(trajectory: List[Dict], box_size: float, radii: Dict[str, float]):
    """
    Visualizes the trajectory of the Monte Carlo simulation using Plotly,
    with varying radii for different component types and distinct colors.

    Args:
        trajectory (List[Dict]): The trajectory data.
        box_size (float): The size of the simulation box.
        radii (Dict[str, float]): Radii of each component type.
    """

    # Define colors for each component type
    colors = {
        'A': 'red',
        'B': 'blue',
        'C': 'green'
    }

    fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

    # Create scatter plot traces for each component type
    for type_name in trajectory[0]['positions'].keys():
        x_coords = [state['positions'][type_name][:, 0] for state in trajectory]
        y_coords = [state['positions'][type_name][:, 1] for state in trajectory]
        z_coords = [state['positions'][type_name][:, 2] for state in trajectory]

        # Flatten the coordinates for all time steps
        x_flat = [x[i] for x in x_coords for i in range(len(x_coords[0]))]
        y_flat = [y[i] for y in y_coords for i in range(len(y_coords[0]))]
        z_flat = [z[i] for z in z_coords for i in range(len(z_coords[0]))]

        # Create a scatter3d trace for the current component type
        fig.add_trace(
            go.Scatter3d(
                x=x_flat,
                y=y_flat,
                z=z_flat,
                mode='markers',
                marker=dict(
                    size=radii[type_name] * 1.0,  # Scale radii for visualization
                    color=colors.get(type_name, 'gray')  # Use default gray if type_name not in colors
                ),
                name=type_name
            )
        )

    # Set layout options for the plot
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[0, box_size]),
            yaxis=dict(range=[0, box_size]),
            zaxis=dict(range=[0, box_size]),
            aspectmode='cube'  # Maintain equal aspect ratio
        ),
        updatemenus=[
            dict(
                type="buttons",
                buttons=[
                    dict(
                        label="Play",
                        method="animate",
                        args=[None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True}]  # Adjust duration for speed
                    ),
                    dict(
                        label="Pause",
                        method="animate",
                        args=[[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate", "transition": {"duration": 0}}]
                    ),
                    dict(
                        label="Step Forward",
                        method="animate",
                        args=[None, {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}]
                    ),
                    dict(
                        label="Jump to End",
                        method="animate",
                        args=[None, {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "to": -1}]  # -1 refers to the last frame
                    )

                ]
            )
        ]
    )

    # Create frames for the animation
    frames = [go.Frame(data=[
        go.Scatter3d(
            x=state['positions'][type_name][:, 0],
            y=state['positions'][type_name][:, 1],
            z=state['positions'][type_name][:, 2],
            mode='markers',
            marker=dict(
                size=radii[type_name] * 0.9,  # Scale radii for visualization
                color=colors.get(type_name, 'gray')  # Use default gray if type_name not in colors
            ),
            name=type_name
        ) for type_name in state['positions']
    ]) for state in trajectory]
    fig.frames = frames

    #fig.show()


def visualize_distances(positions, radii):
    """
    Given final positions (a dict of {type_name: np.ndarray}),
    compute and plot histograms of distances for each allowed pair type.

    Args:
        positions (Dict[str, np.ndarray]): Dictionary of component types to their coordinates.
        radii (Dict[str, float]): Radii (not directly used here, but kept for compatibility).

    Returns:
        None
    """
    allowed_pairs = ["AA", "AB", "BC", "CC"]

    # Flatten all positions and track types
    all_coords = []
    all_types = []
    for tname, coords_array in positions.items():
        for coords in coords_array:
            all_coords.append(coords)
            all_types.append(tname)

    # Calculate distances for each pair type
    pair_distances = {pair_type: [] for pair_type in allowed_pairs}
    for i in range(len(all_coords) - 1):
        for j in range(i + 1, len(all_coords)):
            type_i = all_types[i]
            type_j = all_types[j]
            pair_type = ''.join(sorted([type_i, type_j]))
            if pair_type in allowed_pairs:
                dist = np.linalg.norm(all_coords[i] - all_coords[j])
                pair_distances[pair_type].append(dist)

    # Plot histograms for each allowed pair type
    plt.figure(figsize=(10, 6))
    for pair_type, distances in pair_distances.items():
        plt.hist(distances, bins=20, alpha=0.5, label=pair_type)

    plt.title('Pairwise Distances (Filtered by Allowed Pairs)')
    plt.xlabel('Distance')
    plt.ylabel('Count')
    plt.legend()
    #plt.show()
    
def count_satisfied_restraints(positions, sigma=5.0):
    """
    Takes a single set of coordinates (positions) and computes how many distance
    restraints are satisfied. 'positions' is expected to be a dict of {type_name: np.ndarray}.

    Args:
        positions (Dict[str, np.ndarray]): Mapping from component types to their coordinates.
        pair_distances (Dict[str, float]): Ideal distances keyed by pair types (e.g. "AB").
        sigma (float): Threshold for satisfying the distance restraint.

    Returns:
        int: The number of satisfied distance restraints.
    """
    params = SystemParameters()
    pair_distances = params.pair_distances
    all_positions = []
    all_types = []
    for tname, coords_array in positions.items():
        for coords in coords_array:
            all_positions.append(coords)
            all_types.append(tname)

    num_satisfied = 0
    # Check all pairs
    for i in range(len(all_positions) - 1):
        for j in range(i + 1, len(all_positions)):
            type1 = all_types[i]
            type2 = all_types[j]
            pair_type = ''.join(sorted([type1, type2]))
            if pair_type in pair_distances:
                ideal_distance = pair_distances[pair_type]
                distance = np.linalg.norm(all_positions[i] - all_positions[j])
                #print(f"Distance between {type1} and {type2}: {distance}")
                if abs(distance - ideal_distance) <= sigma:
                    print(f"distance between {type1} and {type2} satisfied: ", distance)
                    num_satisfied += 1
    print(f"Satisfied restraints: {num_satisfied}")

    return num_satisfied
