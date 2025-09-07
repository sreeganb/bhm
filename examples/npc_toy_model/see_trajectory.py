import numpy as np
import plotly.graph_objects as go
from typing import List, Dict, Optional
import h5py
import os
import argparse
from tqdm.auto import tqdm

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
        if max_frames and len(state_names) > max_frames:
            step = max(1, len(state_names) // max_frames)
        
        state_names = state_names[::step]
        
        print(f"Loading {len(state_names)} states from trajectory...")
        for state_name in tqdm(state_names, desc="Loading frames"):
            state_grp = traj_grp[state_name]
            state = {
                "step": state_grp.attrs.get("step", 0),
                "positions": {}
            }
            
            # Load positions
            pos_grp = state_grp['positions']
            for type_name in pos_grp:
                state["positions"][type_name] = pos_grp[type_name][:]
            
            trajectory.append(state)
    
    return trajectory

def visualize_trajectory_plotly(trajectory: List[Dict], box_size: float, radii: Dict[str, float] = None):
    """
    Visualizes the trajectory using Plotly.
    
    Args:
        trajectory: List of state dictionaries
        box_size: Size of the simulation box
        radii: Dictionary mapping particle types to radii
    """
    if radii is None:
        radii = {'A': 1.0, 'B': 0.6, 'C': 0.8, 'D': 0.5, 'E': 1.0}
    
    # Define colors for each component type
    colors = {
        'A': 'red',
        'B': 'blue',
        'C': 'green',
        'D': 'purple',
        'E': 'orange'
    }
    
    # Create the 3D plot
    fig = go.Figure()
    
    # Add initial state traces
    initial_state = trajectory[0]
    for type_name in initial_state['positions'].keys():
        radius = radii.get(type_name, 0.5)
        color = colors.get(type_name, 'gray')
        
        fig.add_trace(go.Scatter3d(
            x=initial_state['positions'][type_name][:, 0],
            y=initial_state['positions'][type_name][:, 1],
            z=initial_state['positions'][type_name][:, 2],
            mode='markers',
            marker=dict(
                size=radius * 10,  # Scale up for visibility
                color=color,
                opacity=0.8,
                symbol='circle'
            ),
            name=f"{type_name} (n={len(initial_state['positions'][type_name])})"
        ))
    
    # Create animation frames
    frames = []
    for i, state in enumerate(trajectory):
        frame_data = []
        for type_name in state['positions'].keys():
            radius = radii.get(type_name, 0.5)
            color = colors.get(type_name, 'gray')
            
            frame_data.append(go.Scatter3d(
                x=state['positions'][type_name][:, 0],
                y=state['positions'][type_name][:, 1],
                z=state['positions'][type_name][:, 2],
                mode='markers',
                marker=dict(
                    size=radius * 10,  # Scale up for visibility
                    color=color,
                    opacity=0.8,
                    symbol='circle'
                ),
                name=f"{type_name} (n={len(state['positions'][type_name])})"
            ))
        frames.append(go.Frame(data=frame_data, name=str(i)))
    
    fig.frames = frames
    
    # Configure animation buttons
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[0, box_size], title="X"),
            yaxis=dict(range=[0, box_size], title="Y"),
            zaxis=dict(range=[0, box_size], title="Z"),
            aspectmode='cube'
        ),
        updatemenus=[
            dict(
                type="buttons",
                buttons=[
                    dict(label="▶️ Play",
                         method="animate",
                         args=[None, {"frame": {"duration": 200, "redraw": True}, "fromcurrent": True}]),
                    dict(label="⏸️ Pause",
                         method="animate",
                         args=[[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}]),
                    dict(label="⏮️ Start",
                         method="animate", 
                         args=[["0"], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}]),
                    dict(label="⏭️ End",
                         method="animate",
                         args=[[str(len(frames)-1)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}])
                ]
            )
        ],
        title="Particle Simulation"
    )
    
    return fig

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize a molecular trajectory')
    parser.add_argument('trajectory_file', type=str, help='Path to the HDF5 trajectory file')
    parser.add_argument('--max-frames', type=int, default=None, help='Maximum number of frames to load')
    parser.add_argument('--step', type=int, default=1, help='Sample every nth frame')
    args = parser.parse_args()
    
    # Load trajectory data
    trajectory = load_trajectory_from_hdf5(args.trajectory_file, args.max_frames, args.step)
    
    if not trajectory:
        print("No trajectory data found.")
        exit(1)
    
    # Determine box size from data
    max_coord = 0
    for state in trajectory:
        for pos in state['positions'].values():
            if len(pos) > 0:
                max_pos = np.max(pos)
                if max_pos > max_coord:
                    max_coord = max_pos
    
    box_size = max_coord * 1.1
    
    # Define radii for different particle types
    radii = {'A': 2.4, 'B': 1.4, 'C': 1.6, 'D': 1.0, 'E': 2.0}
    
    # Visualize the trajectory
    fig = visualize_trajectory_plotly(trajectory, box_size, radii)
    fig.show()