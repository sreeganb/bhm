import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import List, Dict, Optional, Union, Tuple
import h5py
import os
import argparse
from tqdm.auto import tqdm
import pandas as pd
import plotly.express as px

def load_trajectory_from_hdf5(filename: str, max_frames: Optional[int] = None, step: int = 1) -> List[Dict]:
    """
    Loads MCMC trajectory data from an HDF5 file with options for downsampling.
    
    Args:
        filename: Path to the HDF5 file
        max_frames: Maximum number of frames to load (None for all)
        step: Sample every nth frame for efficiency
    
    Returns:
        List of state dictionaries
    """
    trajectory = []
    
    with h5py.File(filename, 'r') as f:
        traj_grp = f['trajectory']
        state_names = sorted(traj_grp.keys(), key=lambda x: int(x.split('_')[1]))
        
        # Apply downsampling if specified
        if max_frames:
            if len(state_names) > max_frames:
                step = max(1, len(state_names) // max_frames)
        
        state_names = state_names[::step]
        
        print(f"Loading {len(state_names)} states from trajectory...")
        for state_name in tqdm(state_names, desc="Loading frames"):
            state_grp = traj_grp[state_name]
            state = {
                "step": state_grp.attrs["step"],
                "total_score": state_grp.attrs["total_score"],
                "prior_score": state_grp.attrs["prior_score"],
                "pair_score": state_grp.attrs["pair_score"],
                "exvol_score": state_grp.attrs["exvol_score"],
                "sigma": {},
                "positions": {},
                "types": {},
                "bead_numbers": {}
            }
            
            # Load sigma values
            sigma_grp = state_grp['sigma']
            for key in sigma_grp.attrs:
                state["sigma"][key] = sigma_grp.attrs[key]
            
            # Load positions
            pos_grp = state_grp['positions']
            for type_name in pos_grp:
                state["positions"][type_name] = pos_grp[type_name][:]
            
            # Load metadata
            try:
                types_keys_dataset = state_grp['types_keys'][:]
                types_vals_dataset = state_grp['types_vals'][:]
                bead_keys_dataset = state_grp['bead_keys'][:]
                bead_vals_dataset = state_grp['bead_vals'][:]

                state["types"] = {types_keys_dataset[i].decode('utf-8'): types_vals_dataset[i].decode('utf-8') 
                                for i in range(len(types_keys_dataset))}
                state["bead_numbers"] = {int(bead_keys_dataset[i]): int(bead_vals_dataset[i]) 
                                        for i in range(len(bead_keys_dataset))}
            except (KeyError, ValueError) as e:
                print(f"Warning: Error loading metadata for state {state_name}: {e}")
                # Create default values if data is missing
                state["types"] = {}
                state["bead_numbers"] = {}
            
            trajectory.append(state)
    
    return trajectory

def extract_trajectory_statistics(trajectory: List[Dict]) -> pd.DataFrame:
    """Extract statistical data from the trajectory for plotting."""
    stats = []
    for state in trajectory:
        stats.append({
            'step': state['step'],
            'total_score': state['total_score'],
            'prior_score': state['prior_score'],
            'pair_score': state['pair_score'],
            'exvol_score': state['exvol_score'],
        })
    
    return pd.DataFrame(stats)

def visualize_trajectory_plotly(
    trajectory: List[Dict], 
    box_size: float, 
    particle_config: Dict[str, Dict],
    animation_speed: int = 200,
    show_statistics: bool = True,
    display_bonds: bool = False,
    bond_distance_cutoff: float = 2.0,
    frame_range: Optional[Tuple[int, int]] = None
):
    """
    Visualizes the trajectory using Plotly with enhanced options.

    Args:
        trajectory: List of state dictionaries
        box_size: Size of the simulation box
        particle_config: Configuration for particle rendering, format:
            {particle_type: {'radius': float, 'color': str, 'opacity': float}}
        animation_speed: Animation speed in ms per frame
        show_statistics: Whether to show statistics panel
        display_bonds: Whether to display bonds between particles
        bond_distance_cutoff: Maximum distance for drawing bonds
        frame_range: Optional tuple (start, end) to display specific frame range
    """
    if frame_range:
        start, end = frame_range
        trajectory = trajectory[start:end]

    # Calculate a scaling factor to make particle sizes proportional to box
    # This converts physical units to reasonable display sizes
    scaling_factor = box_size / 1000  # Adjust divisor to taste
    
    # Create initial frame data - this ensures particles show correctly on page load
    frame_data = []
    initial_state = trajectory[0]
    
    # Add particles for each type to initial view
    for type_name in initial_state['positions'].keys():
        config = particle_config.get(type_name, {'radius': 0.5, 'color': 'gray', 'opacity': 0.9})
        
        # Scale marker size to preserve relative proportions with box size
        scaled_size = config['radius'] / scaling_factor
        
        frame_data.append(go.Scatter3d(
            x=initial_state['positions'][type_name][:, 0],
            y=initial_state['positions'][type_name][:, 1],
            z=initial_state['positions'][type_name][:, 2],
            mode='markers',
            marker=dict(
                size=scaled_size,  # Apply scaling
                color=config['color'],
                opacity=config['opacity'],
                symbol='circle',
                sizemode='diameter',  # Confirm we're treating size as diameter
            ),
            name=f"{type_name} (N={len(initial_state['positions'][type_name])})"
        ))
        
    # Initialize figure with subplots if statistics are displayed
    if show_statistics:
        fig = make_subplots(
            rows=1, cols=2,
            specs=[[{'type': 'scatter3d'}, {'type': 'scatter'}]],
            column_widths=[0.7, 0.3],
            subplot_titles=("Particle Simulation", "Statistics")
        )
    else:
        fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

    # Define default animation duration
    frame_duration = animation_speed
    
    # Create initial frame data - this ensures particles show correctly on page load
    frame_data = []
    initial_state = trajectory[0]
    
    # Add particles for each type to initial view
    for type_name in initial_state['positions'].keys():
        config = particle_config.get(type_name, {'radius': 0.5, 'color': 'gray', 'opacity': 0.9})
        
        frame_data.append(go.Scatter3d(
            x=initial_state['positions'][type_name][:, 0],
            y=initial_state['positions'][type_name][:, 1],
            z=initial_state['positions'][type_name][:, 2],
            mode='markers',
            marker=dict(
                size=config['radius'],
                color=config['color'],
                opacity=config['opacity'],
                symbol='circle',
                sizemode='diameter',
            ),
            name=f"{type_name} (N={len(initial_state['positions'][type_name])})"
        ))
    
    # Add initial frame data to figure
    for trace in frame_data:
        fig.add_trace(trace, row=1, col=1)
    
    # Create bonds for the initial frame if requested
    if display_bonds:
        bonds = compute_bonds(initial_state['positions'], bond_distance_cutoff, particle_config)
        for bond_trace in bonds:
            fig.add_trace(bond_trace, row=1, col=1)

    # Create animation frames
    print("Creating animation frames...")
    frames = []
    for i, state in enumerate(tqdm(trajectory)):
        frame_data = []
        
        # Add particles
        for type_name in state['positions'].keys():
            config = particle_config.get(type_name, {'radius': 0.5, 'color': 'gray', 'opacity': 0.9})
            
            # Scale marker size here too
            scaled_size = config['radius'] / scaling_factor
            
            frame_data.append(go.Scatter3d(
                x=state['positions'][type_name][:, 0],
                y=state['positions'][type_name][:, 1],
                z=state['positions'][type_name][:, 2],
                mode='markers',
                marker=dict(
                    size=scaled_size,  # Apply scaling
                    color=config['color'],
                    opacity=config['opacity'],
                    symbol='circle',
                    sizemode='diameter',
                ),
                name=f"{type_name} (N={len(state['positions'][type_name])})"
            ))
            
    # Add statistics plot if enabled
    if show_statistics:
        stats_df = extract_trajectory_statistics(trajectory)
        
        fig.add_trace(
            go.Scatter(
                x=stats_df['step'], 
                y=stats_df['total_score'],
                mode='lines+markers',
                name='Total Score',
                line=dict(width=2, color='black')
            ),
            row=1, col=2
        )
        
        fig.add_trace(
            go.Scatter(
                x=stats_df['step'], 
                y=stats_df['pair_score'],
                mode='lines',
                name='Pair Score',
                line=dict(width=1, color='blue', dash='dot')
            ),
            row=1, col=2
        )
        
        fig.add_trace(
            go.Scatter(
                x=stats_df['step'], 
                y=stats_df['exvol_score'],
                mode='lines',
                name='Exclusion Score',
                line=dict(width=1, color='red', dash='dot')
            ),
            row=1, col=2
        )
        
        # Add vertical indicator for current frame
        fig.add_trace(
            go.Scatter(
                x=[initial_state['step'], initial_state['step']],
                y=[stats_df['total_score'].min(), stats_df['total_score'].max()],
                mode='lines',
                line=dict(color='gray', width=2, dash='dash'),
                name='Current Frame',
                showlegend=False
            ),
            row=1, col=2
        )
    
    # Create animation frames
    print("Creating animation frames...")
    frames = []
    for i, state in enumerate(tqdm(trajectory)):
        frame_data = []
        
        # Add particles
        for type_name in state['positions'].keys():
            config = particle_config.get(type_name, {'radius': 0.5, 'color': 'gray', 'opacity': 0.9})
            
            frame_data.append(go.Scatter3d(
                x=state['positions'][type_name][:, 0],
                y=state['positions'][type_name][:, 1],
                z=state['positions'][type_name][:, 2],
                mode='markers',
                marker=dict(
                    size=config['radius'],
                    color=config['color'],
                    opacity=config['opacity'],
                    symbol='circle',
                    sizemode='diameter',
                ),
                name=f"{type_name} (N={len(state['positions'][type_name])})"
            ))
        
        # Add bonds if requested
        if display_bonds:
            bonds = compute_bonds(state['positions'], bond_distance_cutoff, particle_config)
            for bond_trace in bonds:
                frame_data.append(bond_trace)
        
        # Add statistics updates
        if show_statistics:
            # Frame-specific stats traces
            step = state['step']
            
            # Position indicator for current frame
            frame_data.append(
                go.Scatter(
                    x=[step, step],
                    y=[stats_df['total_score'].min(), stats_df['total_score'].max()],
                    mode='lines',
                    line=dict(color='gray', width=2, dash='dash'),
                    name='Current Frame',
                    showlegend=False
                )
            )
        
        frames.append(go.Frame(data=frame_data, name=str(i)))
    
    fig.frames = frames
    
    # Update layout with improved controls
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[0, box_size], title="X"),
            yaxis=dict(range=[0, box_size], title="Y"),
            zaxis=dict(range=[0, box_size], title="Z"),
            aspectmode='cube',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.1,
                y=0,
                buttons=[
                    dict(label="▶️ Play",
                         method="animate",
                         args=[None, {"frame": {"duration": frame_duration, "redraw": True}, "fromcurrent": True}]),
                    dict(label="⏸️ Pause",
                         method="animate",
                         args=[[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}]),
                    dict(label="⏮️ Start",
                         method="animate",
                         args=[["0"], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}]),
                    dict(label="⏭️ End",
                         method="animate",
                         args=[[str(len(frames)-1)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}]),
                ]
            )
        ],
        sliders=[dict(
            active=0,
            steps=[dict(
                label=f"Step {trajectory[i]['step']}",
                method="animate",
                args=[[str(i)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}]
            ) for i in range(0, len(trajectory), max(1, len(trajectory)//20))],  # Show fewer labels for readability
            currentvalue={"prefix": "Frame: ", "visible": True},
            pad={"t": 50}
        )],
        title=f"Particle Simulation Trajectory",
        hovermode='closest',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
    )
    
    # Optimize for large datasets
    if len(trajectory) > 50:
        print("Large trajectory detected, optimizing for performance...")
        fig.update_layout(uirevision='true')  # Preserve UI state on updates
    
    # Show the figure
    fig.show(config={"scrollZoom": True})
    
    return fig

def compute_bonds(positions: Dict[str, np.ndarray], cutoff: float, particle_config: Dict) -> List:
    """Compute bonds between particles based on distance cutoff."""
    bond_traces = []
    
    # Create all pairs of particle types
    types = list(positions.keys())
    
    for i, type1 in enumerate(types):
        for type2 in types[i:]:  # Only upper triangle to avoid duplicates
            pos1 = positions[type1]
            pos2 = positions[type2]
            
            # Compute pairwise distances
            if type1 == type2:
                # Same type - need to avoid self-connections
                for i in range(len(pos1)):
                    for j in range(i+1, len(pos1)):
                        dist = np.linalg.norm(pos1[i] - pos1[j])
                        if dist <= cutoff:
                            bond_traces.append(
                                go.Scatter3d(
                                    x=[pos1[i,0], pos1[j,0]],
                                    y=[pos1[i,1], pos1[j,1]],
                                    z=[pos1[i,2], pos1[j,2]],
                                    mode='lines',
                                    line=dict(
                                        color='rgba(100,100,100,0.7)',
                                        width=2
                                    ),
                                    showlegend=False
                                )
                            )
            else:
                # Different types
                for i in range(len(pos1)):
                    for j in range(len(pos2)):
                        dist = np.linalg.norm(pos1[i] - pos2[j])
                        if dist <= cutoff:
                            bond_traces.append(
                                go.Scatter3d(
                                    x=[pos1[i,0], pos2[j,0]],
                                    y=[pos1[i,1], pos2[j,1]],
                                    z=[pos1[i,2], pos2[j,2]],
                                    mode='lines',
                                    line=dict(
                                        color='rgba(100,100,100,0.7)', 
                                        width=2
                                    ),
                                    showlegend=False
                                )
                            )
    
    return bond_traces

def create_trajectory_movie(
    h5_file_path: str,
    max_frames: Optional[int] = None,
    frame_step: int = 1,
    animation_speed: int = 200,
    show_statistics: bool = True,
    display_bonds: bool = False,
    bond_distance_cutoff: float = 2.0,
    custom_colors: Dict[str, str] = None,
    custom_radii: Dict[str, float] = None,
    custom_opacity: Dict[str, float] = None,
    frame_range: Optional[Tuple[int, int]] = None
):
    """
    Create an interactive trajectory visualization with customizable parameters.
    
    Args:
        h5_file_path: Path to HDF5 trajectory file
        max_frames: Maximum frames to load (None for all)
        frame_step: Sample every nth frame
        animation_speed: Animation speed (ms per frame)
        show_statistics: Whether to show statistics panel
        display_bonds: Whether to display bonds between nearby particles
        bond_distance_cutoff: Maximum distance for bond display
        custom_colors: Optional dict mapping particle types to colors
        custom_radii: Optional dict mapping particle types to radii
        custom_opacity: Optional dict mapping particle types to opacity values
        frame_range: Optional tuple (start, end) to display specific frame range
    """
    if not os.path.exists(h5_file_path):
        raise FileNotFoundError(f"HDF5 file not found: {h5_file_path}")
    
    # Load trajectory with potential downsampling
    trajectory = load_trajectory_from_hdf5(h5_file_path, max_frames, frame_step)
    
    if not trajectory:
        print("Trajectory data is empty.")
        return
    
    # Get all particle types
    first_state = trajectory[0]
    all_types = list(first_state['positions'].keys())
    
    # Calculate maximum coordinate for box size
    max_coord = 0.0
    for state in trajectory:
        for type_name in all_types:
            if len(state['positions'][type_name]) > 0:
                max_coord = max(max_coord, np.max(state['positions'][type_name]))
    
    # Set up default particle configuration
    default_colors = {'A': 'red', 'B': 'blue', 'C': 'green'}
    default_radii = {'A': 16.0, 'B': 9.333, 'C': 10.666}
    default_opacity = {'A': 0.9, 'B': 0.9, 'C': 0.9}
    
    # Apply custom settings if provided
    if custom_colors:
        default_colors.update(custom_colors)
    if custom_radii:
        default_radii.update(custom_radii)
    if custom_opacity:
        default_opacity.update(custom_opacity)
    
    # Build particle configuration dictionary
    particle_config = {}
    for type_name in all_types:
        particle_config[type_name] = {
            'radius': default_radii.get(type_name, 3.0),
            'color': default_colors.get(type_name, px.colors.qualitative.Plotly[hash(type_name) % len(px.colors.qualitative.Plotly)]),
            'opacity': default_opacity.get(type_name, 0.8)
        }
    
    # Create and show the visualization
    fig = visualize_trajectory_plotly(
        trajectory, 
        max_coord * 1.2,  # Add 10% padding to box size
        particle_config,
        animation_speed,
        show_statistics,
        display_bonds,
        bond_distance_cutoff,
        frame_range
    )
    
    return fig

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a 3D visualization of a molecular trajectory.')
    parser.add_argument('trajectory_file', type=str, help='Path to the HDF5 trajectory file')
    parser.add_argument('--max-frames', type=int, default=None, help='Maximum number of frames to load')
    parser.add_argument('--step', type=int, default=1, help='Sample every nth frame')
    parser.add_argument('--speed', type=int, default=200, help='Animation speed (ms per frame)')
    parser.add_argument('--no-stats', action='store_false', dest='show_stats', help='Hide statistics panel')
    parser.add_argument('--show-bonds', action='store_true', help='Display bonds between nearby particles')
    parser.add_argument('--bond-cutoff', type=float, default=2.0, help='Maximum distance for bond display')
    
    args = parser.parse_args()
    
    # Example custom settings (you can modify these or implement command-line options)
    custom_colors = {'A': 'crimson', 'B': 'royalblue', 'C': 'forestgreen'}
    custom_radii = {'A': 16.0, 'B': 9.333, 'C': 10.666}
    
    create_trajectory_movie(
        args.trajectory_file,
        max_frames=args.max_frames,
        frame_step=args.step,
        animation_speed=args.speed,
        show_statistics=args.show_stats,
        display_bonds=args.show_bonds,
        bond_distance_cutoff=args.bond_cutoff,
        custom_colors=custom_colors,
        custom_radii=custom_radii
    )
