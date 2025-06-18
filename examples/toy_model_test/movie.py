import numpy as np
import pyvista as pv
from typing import List, Dict, Optional, Tuple
import h5py
import os
import argparse
from tqdm.auto import tqdm
import pandas as pd
import time
import threading

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
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"HDF5 file not found: {filename}")
    
    with h5py.File(filename, 'r') as f:
        if 'trajectory' not in f:
            raise KeyError(f"No 'trajectory' group found in {filename}")
            
        traj_grp = f['trajectory']
        
        # Debug: Print all available keys
        all_keys = list(traj_grp.keys())
        print(f"Found {len(all_keys)} states in trajectory file")
        if len(all_keys) > 0:
            print(f"First few keys: {all_keys[:5]}")
            print(f"Last few keys: {all_keys[-5:]}")
        
        # More robust sorting function
        def extract_step_number(state_name):
            try:
                # Handle both "state_XXXXX" and other formats
                if '_' in state_name:
                    return int(state_name.split('_')[-1])  # Use last part after split
                else:
                    # Try to extract number from the string
                    import re
                    numbers = re.findall(r'\d+', state_name)
                    return int(numbers[0]) if numbers else 0
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not extract step number from '{state_name}': {e}")
                return 0
        
        # Sort state names properly
        state_names = sorted(all_keys, key=extract_step_number)
        
        # Debug: Verify sorting worked correctly
        if len(state_names) > 1:
            first_step = extract_step_number(state_names[0])
            last_step = extract_step_number(state_names[-1])
            print(f"Step range: {first_step} to {last_step}")
            
            # Check for gaps or issues
            if last_step < first_step:
                print("WARNING: Last step is smaller than first step - sorting may have failed!")
        
        # Apply downsampling if specified
        if max_frames:
            if len(state_names) > max_frames:
                step = max(1, len(state_names) // max_frames)
        
        state_names = state_names[::step]
        
        print(f"Loading {len(state_names)} states from trajectory...")
        
        for i, state_name in enumerate(tqdm(state_names, desc="Loading frames")):
            try:
                state_grp = traj_grp[state_name]
                
                # More robust attribute reading with error handling
                state = {
                    "step": state_grp.attrs.get("step", 0),
                    "total_score": state_grp.attrs.get("total_score", 0.0),
                    "prior_score": state_grp.attrs.get("prior_score", 0.0),
                    "pair_score": state_grp.attrs.get("pair_score", 0.0),
                    "exvol_score": state_grp.attrs.get("exvol_score", 0.0),
                    "tet_score": state_grp.attrs.get("tet_score", 0.0),      # Added missing field
                    "oct_score": state_grp.attrs.get("oct_score", 0.0),      # Added missing field
                    "sigma": {},
                    "positions": {},
                    "types": {},
                    "bead_numbers": {}
                }
                
                # Debug: Print step info for first and last frames
                if i == 0 or i == len(state_names) - 1:
                    print(f"Frame {i} ({state_name}): step={state['step']}, score={state['total_score']:.4f}")
                
                # Load sigma values with error handling
                if 'sigma' in state_grp:
                    sigma_grp = state_grp['sigma']
                    for key in sigma_grp.attrs:
                        state["sigma"][key] = float(sigma_grp.attrs[key])
                else:
                    print(f"Warning: No sigma group in state {state_name}")
                
                # Load positions with error handling and validation
                if 'positions' in state_grp:
                    pos_grp = state_grp['positions']
                    for type_name in pos_grp:
                        positions = pos_grp[type_name][:]
                        state["positions"][type_name] = positions
                        
                        # Debug: Print position info for first frame
                        if i == 0:
                            print(f"  {type_name}: {len(positions)} particles")
                            if len(positions) > 0:
                                center = np.mean(positions, axis=0)
                                print(f"    Center: {center}")
                                print(f"    Range: {np.min(positions, axis=0)} to {np.max(positions, axis=0)}")
                else:
                    print(f"Warning: No positions group in state {state_name}")
                
                # Load metadata with better error handling
                try:
                    if 'types_keys' in state_grp and 'types_vals' in state_grp:
                        types_keys_dataset = state_grp['types_keys'][:]
                        types_vals_dataset = state_grp['types_vals'][:]
                        
                        if len(types_keys_dataset) == len(types_vals_dataset):
                            state["types"] = {}
                            for j in range(len(types_keys_dataset)):
                                k = types_keys_dataset[j]
                                v = types_vals_dataset[j]
                                if isinstance(k, bytes):
                                    k = k.decode('utf-8', errors='ignore')
                                if isinstance(v, bytes):
                                    v = v.decode('utf-8', errors='ignore')
                                state["types"][k] = v
                    
                    if 'bead_keys' in state_grp and 'bead_vals' in state_grp:
                        bead_keys_dataset = state_grp['bead_keys'][:]
                        bead_vals_dataset = state_grp['bead_vals'][:]
                        
                        if len(bead_keys_dataset) == len(bead_vals_dataset):
                            state["bead_numbers"] = {}
                            for j in range(len(bead_keys_dataset)):
                                bkey = bead_keys_dataset[j]
                                bval = bead_vals_dataset[j]
                                state["bead_numbers"][int(bkey)] = int(bval)
                                
                except (KeyError, ValueError, TypeError) as e:
                    print(f"Warning: Error loading metadata for state {state_name}: {e}")
                    state["types"] = {}
                    state["bead_numbers"] = {}
                
                trajectory.append(state)
                
            except Exception as e:
                print(f"Error loading state {state_name}: {e}")
                continue  # Skip this frame and continue
    
    print(f"Successfully loaded {len(trajectory)} states from {filename}")
    
    # Final validation
    if len(trajectory) > 1:
        first_step = trajectory[0]['step']
        last_step = trajectory[-1]['step']
        print(f"Trajectory spans steps {first_step} to {last_step}")
        
        # Check for step consistency
        steps = [state['step'] for state in trajectory]
        if steps != sorted(steps):
            print("WARNING: Trajectory steps are not in ascending order!")
            print(f"First 10 steps: {steps[:10]}")
    
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

class TrajectoryVisualizer:
    """PyVista-based trajectory visualizer with enhanced controls."""
    
    def __init__(self, trajectory, box_size, particle_config, show_statistics=True):
        self.trajectory = trajectory
        self.box_size = box_size
        self.particle_config = particle_config
        self.show_statistics = show_statistics
        self.current_frame = 0
        self.is_playing = False
        self.animation_speed = 200  # ms
        self.auto_loop = True  # Loop animation automatically
        
        # Initialize plotter with enhanced settings
        self.plotter = pv.Plotter(
            window_size=[1400, 800],
            off_screen=False,
            notebook=False
        )
        self.plotter.set_background('white')
        
        # Store actors for each particle type
        self.particle_actors = {}
        self.box_actor = None
        self.slider = None
        
        # Statistics data
        if self.show_statistics:
            self.stats_df = extract_trajectory_statistics(trajectory)
        
        try:
            self.setup_scene()
            self.create_ui_controls()
            self.add_slider_widget()
        except Exception as e:
            print(f"Error setting up scene: {e}")
            raise
        
    def setup_scene(self):
        """Set up the 3D scene with initial frame."""
        # Create simulation box
        self.create_simulation_box()
        
        # Add particles for initial frame
        self.update_particles(self.trajectory[0])
        
        # Set up camera and lighting
        self.setup_camera_and_lighting()
        
        # Add coordinate axes
        self.plotter.add_axes(
            xlabel='X (Å)', 
            ylabel='Y (Å)', 
            zlabel='Z (Å)',
            line_width=3,
            labels_off=False
        )
        
        # Add frame information text
        self.frame_text = self.plotter.add_text(
            f"Frame: {self.current_frame}/{len(self.trajectory)-1}\n"
            f"Step: {self.trajectory[self.current_frame]['step']}\n"
            f"Score: {self.trajectory[self.current_frame]['total_score']:.2f}",
            position='upper_left',
            font_size=12,
            color='black'
        )
            
    def create_simulation_box(self):
        """Create a wireframe box showing simulation boundaries."""
        # Create box wireframe
        box = pv.Box(bounds=[0, self.box_size, 0, self.box_size, 0, self.box_size])
        box_wireframe = box.extract_all_edges()
        
        self.box_actor = self.plotter.add_mesh(
            box_wireframe, 
            color='gray', 
            line_width=2, 
            opacity=0.3,
            name='simulation_box'
        )
        
    def add_slider_widget(self):
        """Add a slider widget for frame control."""
        def slider_callback(value):
            """Callback function for slider interaction."""
            self.is_playing = False  # Stop auto-play when manually using slider
            new_frame = int(value)
            if new_frame != self.current_frame:
                self.current_frame = new_frame
                self.update_frame()
        
        # Add slider widget
        self.slider = self.plotter.add_slider_widget(
            callback=slider_callback,
            rng=[0, len(self.trajectory) - 1],
            value=0,
            title="Frame",
            pointa=(0.1, 0.1),
            pointb=(0.9, 0.1),
            style='modern'
        )
        
    def update_particles(self, state):
        """Update particle positions for the current state."""
        # Fix 1: Create a copy of the dictionary values to avoid iteration issues
        actors_to_remove = list(self.particle_actors.values())
        for actor in actors_to_remove:
            try:
                self.plotter.remove_actor(actor)
            except:
                pass  # Ignore if actor doesn't exist
        self.particle_actors.clear()
        
        # Add particles for each type
        for type_name, positions in state['positions'].items():
            if len(positions) == 0:
                continue
                
            config = self.particle_config.get(type_name, {
                'radius': 5.0, 
                'color': 'gray', 
                'opacity': 0.8
            })
            
            # More robust sphere creation with proper error handling
            try:
                # Create points mesh
                points = pv.PolyData(positions)
                
                # Create sphere geometry with reduced complexity for better performance
                sphere_glyph = pv.Sphere(
                    radius=config['radius'],
                    phi_resolution=12,  # Reduced from default 20
                    theta_resolution=12  # Reduced from default 20
                )
                
                # Create spheres at all positions with better error handling
                spheres = points.glyph(
                    geom=sphere_glyph, 
                    scale=False,
                    orient=False
                )
                
                # Ensure proper normals for OpenGL rendering
                if not spheres.point_normals.size:
                    spheres = spheres.compute_normals(
                        cell_normals=True,
                        point_normals=True,
                        split_vertices=False,
                        flip_normals=False
                    )
                
                # Add to scene with safer material properties
                actor = self.plotter.add_mesh(
                    spheres,
                    color=config['color'],
                    opacity=config['opacity'],
                    specular=0.3,
                    specular_power=10,
                    ambient=0.4,
                    diffuse=0.6,
                    smooth_shading=True,
                    name=f'particles_{type_name}',
                    reset_camera=False,
                    render=False
                )
                
                self.particle_actors[type_name] = actor
                
            except Exception as e:
                print(f"Error creating spheres for {type_name}: {e}")
                # Fallback to simple points if sphere creation fails
                try:
                    actor = self.plotter.add_points(
                        positions,
                        color=config['color'],
                        point_size=config['radius'],
                        name=f'particles_{type_name}_fallback',
                        reset_camera=False,
                        render=False
                    )
                    self.particle_actors[type_name] = actor
                except Exception as e2:
                    print(f"Fallback also failed for {type_name}: {e2}")
            
    def setup_camera_and_lighting(self):
        """Set up camera position and lighting for optimal visualization."""
        # Set camera position
        self.plotter.camera_position = [
            (self.box_size * 1.5, self.box_size * 1.5, self.box_size * 1.5),
            (self.box_size / 2, self.box_size / 2, self.box_size / 2),
            (0, 0, 1)
        ]
        
        # Add additional lighting for better sphere rendering
        self.plotter.add_light(pv.Light(
            position=(self.box_size, self.box_size, self.box_size * 2),
            focal_point=(self.box_size/2, self.box_size/2, self.box_size/2),
            color='white',
            intensity=0.8
        ))
        
        self.plotter.add_light(pv.Light(
            position=(0, 0, self.box_size * 2),
            focal_point=(self.box_size/2, self.box_size/2, self.box_size/2),
            color='white',
            intensity=0.3
        ))
        
    def create_ui_controls(self):
        """Create UI controls for animation playback."""
        def play_pause_callback():
            self.is_playing = not self.is_playing
            if self.is_playing:
                print("▶️ Starting auto-play animation")
                self.start_animation()
            else:
                print("⏸️ Pausing animation")
            
        def next_frame_callback():
            self.is_playing = False
            self.next_frame()
            
        def prev_frame_callback():
            self.is_playing = False
            self.previous_frame()
            
        def reset_callback():
            self.is_playing = False
            self.current_frame = 0
            self.update_frame()
            self.update_slider()
            
        def end_callback():
            self.is_playing = False
            self.current_frame = len(self.trajectory) - 1
            self.update_frame()
            self.update_slider()
            
        def speed_up_callback():
            self.animation_speed = max(50, self.animation_speed - 50)
            print(f"🐰 Animation speed increased: {self.animation_speed}ms per frame")
            
        def speed_down_callback():
            self.animation_speed = min(1000, self.animation_speed + 50)
            print(f"🐌 Animation speed decreased: {self.animation_speed}ms per frame")
            
        def toggle_loop_callback():
            self.auto_loop = not self.auto_loop
            print(f"🔄 Auto-loop: {'ON' if self.auto_loop else 'OFF'}")
        
        # Add keyboard controls
        self.plotter.add_key_event('space', play_pause_callback)
        self.plotter.add_key_event('Right', next_frame_callback)
        self.plotter.add_key_event('Left', prev_frame_callback)
        self.plotter.add_key_event('r', reset_callback)
        self.plotter.add_key_event('e', end_callback)
        self.plotter.add_key_event('Up', speed_up_callback)
        self.plotter.add_key_event('Down', speed_down_callback)
        self.plotter.add_key_event('l', toggle_loop_callback)
        
        # Enhanced instructions
        instructions = (
            "🎮 CONTROLS:\n"
            "SPACE: ▶️/⏸️ Play/Pause auto-run\n"
            "← →: Previous/Next frame\n"
            "R: ⏮️ Reset to start\n"
            "E: ⏭️ Jump to end\n"
            "↑ ↓: 🐰🐌 Speed up/down\n"
            "L: 🔄 Toggle auto-loop\n"
            "🖱️ Mouse: Rotate/Zoom\n"
            "📊 Slider: Scrub frames"
        )
        
        self.plotter.add_text(
            instructions,
            position='lower_left',
            font_size=9,
            color='dark_blue'
        )
        
        # Add status text
        self.status_text = self.plotter.add_text(
            "Press SPACE to start auto-play",
            position='lower_right',
            font_size=10,
            color='green'
        )
        
    def update_slider(self):
        """Update slider position to match current frame."""
        if self.slider is not None:
            try:
                # Update slider value without triggering callback
                self.slider.GetRepresentation().SetValue(self.current_frame)
                self.plotter.render()
            except:
                pass  # Ignore slider update errors
        
    def update_frame(self):
        """Update the visualization to show the current frame."""
        if 0 <= self.current_frame < len(self.trajectory):
            state = self.trajectory[self.current_frame]
            
            # Update particles with error handling
            try:
                self.update_particles(state)
            except Exception as e:
                print(f"Error updating particles for frame {self.current_frame}: {e}")
                return
            
            # Update frame information
            frame_info = (
                f"Frame: {self.current_frame}/{len(self.trajectory)-1}\n"
                f"Step: {state['step']}\n"
                f"Total Score: {state['total_score']:.2f}\n"
                f"Pair Score: {state['pair_score']:.2f}\n"
                f"ExVol Score: {state['exvol_score']:.2f}"
            )
            
            # Update status
            status_info = (
                f"{'▶️ PLAYING' if self.is_playing else '⏸️ PAUSED'} | "
                f"Speed: {self.animation_speed}ms | "
                f"Loop: {'ON' if self.auto_loop else 'OFF'}"
            )
            
            # Safer text update
            try:
                self.plotter.remove_actor(self.frame_text)
                self.plotter.remove_actor(self.status_text)
            except:
                pass
                
            self.frame_text = self.plotter.add_text(
                frame_info,
                position='upper_left',
                font_size=12,
                color='black'
            )
            
            self.status_text = self.plotter.add_text(
                status_info,
                position='lower_right',
                font_size=10,
                color='green' if self.is_playing else 'orange'
            )
            
            # Update slider position
            self.update_slider()
            
            # Force render update with error handling
            try:
                self.plotter.render()
            except Exception as e:
                print(f"Render error: {e}")
            
    def next_frame(self):
        """Go to next frame with auto-loop support."""
        if self.current_frame < len(self.trajectory) - 1:
            self.current_frame += 1
        elif self.auto_loop and self.is_playing:
            self.current_frame = 0  # Loop back to start
        else:
            self.is_playing = False  # Stop at end if not looping
            
        self.update_frame()
            
    def previous_frame(self):
        """Go to previous frame."""
        if self.current_frame > 0:
            self.current_frame -= 1
            self.update_frame()
            
    def start_animation(self):
        """Start the enhanced animation loop."""
        def animation_loop():
            try:
                while self.is_playing:
                    time.sleep(self.animation_speed / 1000.0)
                    if self.is_playing:  # Check again in case it was paused
                        try:
                            self.next_frame()
                            # If we're at the end and not looping, stop
                            if (self.current_frame == len(self.trajectory) - 1 and 
                                not self.auto_loop):
                                self.is_playing = False
                                print("📍 Animation completed - reached end")
                                break
                        except Exception as e:
                            print(f"Error in animation frame {self.current_frame}: {e}")
                            self.is_playing = False
                            break
            except Exception as e:
                print(f"Animation loop error: {e}")
            finally:
                if self.is_playing:  # Update status if we stopped due to reaching end
                    self.is_playing = False
                    self.update_frame()  # Update status display
            
        # Run animation in separate thread to avoid blocking
        animation_thread = threading.Thread(target=animation_loop)
        animation_thread.daemon = True
        animation_thread.start()
        
    def show(self):
        """Display the visualization with enhanced instructions."""
        print("\n🎬 PyVista Enhanced Trajectory Visualizer")
        print("=" * 50)
        print("🎮 CONTROLS:")
        print("  SPACEBAR: ▶️/⏸️ Play/Pause auto-run animation")
        print("  ← →: Previous/Next frame manually")
        print("  R: ⏮️ Reset to start")
        print("  E: ⏭️ Jump to end")
        print("  ↑ ↓: 🐰🐌 Speed up/down animation")
        print("  L: 🔄 Toggle auto-loop (restart at end)")
        print("  🖱️ Mouse: Rotate/Zoom view")
        print("  📊 Slider: Click and drag to scrub through frames")
        print("=" * 50)
        print("💡 TIP: Press SPACE to start auto-play, then use slider")
        print("    to jump to specific frames while pausing auto-play!")
        print("=" * 50)
        
        self.plotter.show()

def visualize_trajectory_pyvista(
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
    Visualizes the trajectory using PyVista with enhanced controls.
    """
    if frame_range:
        start, end = frame_range
        trajectory = trajectory[start:end]
    
    print(f"Creating enhanced PyVista visualization with {len(trajectory)} frames...")
    
    # Create and show the visualizer
    visualizer = TrajectoryVisualizer(
        trajectory, 
        box_size, 
        particle_config, 
        show_statistics
    )
    
    visualizer.animation_speed = animation_speed
    visualizer.show()
    
    return visualizer

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
    Create an enhanced interactive trajectory visualization using PyVista.
    """
    if not os.path.exists(h5_file_path):
        raise FileNotFoundError(f"HDF5 file not found: {h5_file_path}")
    
    # Load trajectory with potential downsampling
    trajectory = load_trajectory_from_hdf5(h5_file_path, max_frames, frame_step)
    
    if not trajectory:
        print("Trajectory data is empty.")
        return
    
    # Get all particle types from first state
    first_state = trajectory[0]
    all_types = list(first_state['positions'].keys())
    
    # Calculate box size more intelligently
    all_coords = []
    for state in trajectory[:min(10, len(trajectory))]:
        for type_name in all_types:
            if len(state['positions'][type_name]) > 0:
                all_coords.extend(state['positions'][type_name].flatten())
    
    if all_coords:
        min_coord = np.min(all_coords)
        max_coord = np.max(all_coords)
        
        # Adjust coordinates to start from 0 if needed
        if min_coord < 0:
            print(f"Shifting coordinates by {-min_coord} to ensure positive values")
            for state in trajectory:
                for type_name in all_types:
                    state['positions'][type_name] += abs(min_coord)
            box_size = max_coord - min_coord + 100
        else:
            box_size = max_coord + 100
    else:
        box_size = 1000
    
    print(f"Using box size: {box_size:.1f} Å")
    
    # Set up particle configuration
    default_colors = {
        'A': 'red', 
        'B': 'blue', 
        'C': 'green', 
        'D': 'orange', 
        'E': 'purple'
    }
    default_radii = {'A': 24.0, 'B': 14.0, 'C': 16.0}
    default_opacity = {'A': 0.8, 'B': 0.8, 'C': 0.8}
    
    # Apply custom settings if provided
    if custom_colors:
        default_colors.update(custom_colors)
    if custom_radii:
        default_radii.update(custom_radii)
    if custom_opacity:
        default_opacity.update(custom_opacity)
    
    # Build particle configuration dictionary
    particle_config = {}
    color_options = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow']
    
    for i, type_name in enumerate(all_types):
        particle_config[type_name] = {
            'radius': default_radii.get(type_name, 5.0),
            'color': default_colors.get(type_name, color_options[i % len(color_options)]),
            'opacity': default_opacity.get(type_name, 0.8)
        }
    
    print("Particle configuration:")
    for type_name, config in particle_config.items():
        print(f"  {type_name}: radius={config['radius']:.1f}Å, color={config['color']}")
    
    # Create and show the visualization
    visualizer = visualize_trajectory_pyvista(
        trajectory, 
        box_size,
        particle_config,
        animation_speed,
        show_statistics,
        display_bonds,
        bond_distance_cutoff,
        frame_range
    )
    
    return visualizer

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create an enhanced 3D PyVista visualization of a molecular trajectory.')
    parser.add_argument('trajectory_file', type=str, help='Path to the HDF5 trajectory file')
    parser.add_argument('--max-frames', type=int, default=None, help='Maximum number of frames to load')
    parser.add_argument('--step', type=int, default=1, help='Sample every nth frame')
    parser.add_argument('--speed', type=int, default=200, help='Animation speed (ms per frame)')
    parser.add_argument('--no-stats', action='store_false', dest='show_stats', help='Hide statistics panel')
    parser.add_argument('--show-bonds', action='store_true', help='Ignored (bonds removed)')
    parser.add_argument('--bond-cutoff', type=float, default=2.0, help='Ignored (bonds removed)')
    
    args = parser.parse_args()
    
    # Example custom settings
    custom_colors = {'A': 'red', 'B': 'blue', 'C': 'green'}
    custom_radii = {'A': 24.0, 'B': 14.0, 'C': 10.666}
    
    try:
        visualizer = create_trajectory_movie(
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
        print("🎉 Enhanced visualization completed successfully!")
        
    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()