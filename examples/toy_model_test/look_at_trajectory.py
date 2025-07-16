import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import h5py
import os
import sys
from typing import List, Dict, Optional
from tqdm.auto import tqdm
import re
# Imports for 3D spheres
from mpl_toolkits.mplot3d import art3d
from matplotlib.colors import LightSource

# Import SystemParameters from parameters.py
from parameters import SystemParameters

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
                    "tet_score": state_grp.attrs.get("tet_score", 0.0),
                    "oct_score": state_grp.attrs.get("oct_score", 0.0),
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

class MatplotlibTrajectoryVisualizer:
    def __init__(self, trajectory, params):
        self.trajectory = trajectory
        self.params = params
        self.box_size = params.box_size
        self.particle_radii = params.radii
        self.current_frame = 0
        
        # Set up the figure and 3D axis with improved resolution
        self.fig = plt.figure(figsize=(12, 10), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # Initialize collections for spheres
        self.sphere_collections = {}
        
        # Define particle colors with metallic appearance
        self.particle_colors = {
            'A': '#FF3B30',  # Bright red
            'B': '#147EFB',  # Bright blue
            'C': '#28CD41',  # Bright green
            'D': '#FF9500',  # Bright orange
            'E': '#AF52DE',  # Bright purple
            'F': '#A2845E'   # Rich brown
        }
        
        # For metallic appearance
        self.light_source = LightSource(azdeg=315, altdeg=45)
        
        # Quality of sphere rendering - higher is better quality but slower
        # Adjust based on performance needs
        self.sphere_resolution = 15
        
        # Reference to the legend for updating
        self.legend_elements = []
        
        self.setup_scene()
        
    def setup_scene(self):
        """Set up the 3D scene with improved initial view."""
        # Get first frame to determine proper view
        first_state = self.trajectory[0]
        
        # Calculate center and extent of particles
        all_positions = []
        for type_name, pos in first_state['positions'].items():
            if len(pos) > 0:
                all_positions.append(pos)
        
        if all_positions:
            all_positions = np.vstack(all_positions)
            center = np.mean(all_positions, axis=0)
            max_radius = max(self.particle_radii.values()) if self.particle_radii else 10.0
            
            # Focus on particles with closer zoom
            # Calculate range based on particle distribution
            pos_range = np.max(all_positions, axis=0) - np.min(all_positions, axis=0)
            max_range = max(pos_range) * 0.65  # Use 65% of the max range for better focus
            
            # Set limits with better focus on particles
            self.ax.set_xlim(center[0] - max_range, center[0] + max_range)
            self.ax.set_ylim(center[1] - max_range, center[1] + max_range)
            self.ax.set_zlim(center[2] - max_range, center[2] + max_range)
        else:
            # Fallback to box center if no particles
            center = self.box_size / 2
            self.ax.set_xlim(0, self.box_size)
            self.ax.set_ylim(0, self.box_size)
            self.ax.set_zlim(0, self.box_size)
        
        # Improve axis appearance
        self.ax.set_xlabel('X (Å)', fontweight='bold', fontsize=12)
        self.ax.set_ylabel('Y (Å)', fontweight='bold', fontsize=12)
        self.ax.set_zlabel('Z (Å)', fontweight='bold', fontsize=12)
        
        # Set optimal viewing angle
        self.ax.view_init(elev=35, azim=45)
        
        # Add title
        self.title = self.ax.set_title('Molecular Trajectory Visualization', 
                                      fontsize=14, fontweight='bold')
        
        # Create wireframe box
        self.create_box_wireframe()
        
        # Enable better 3D depth perception
        self.ax.set_facecolor('#F5F5F5')  # Light gray background
        self.fig.patch.set_facecolor('#E8E8E8')  # Slightly darker for figure
        
        # Configure lighting for better 3D effect
        self.ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
        
    def create_sphere(self, center, radius, color):
        """Create a 3D sphere with metallic appearance."""
        # Create sphere coordinates
        u = np.linspace(0, 2 * np.pi, self.sphere_resolution)
        v = np.linspace(0, np.pi, self.sphere_resolution)
        x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
        
        # Create metallic appearance
        # Create gradient for metallic effect
        rgb = np.array(plt.cm.colors.to_rgb(color))
        
        # Lighten color for highlights
        highlight = rgb + 0.4
        highlight = np.clip(highlight, 0, 1)
        
        # Create color array with gradient for metallic effect
        facecolors = np.zeros((u.size, v.size, 3))
        for i in range(u.size):
            # Create highlight effect based on angle
            gradient_factor = 0.7 + 0.3 * np.sin(u[i] + np.pi/4) * np.sin(v)
            for j in range(v.size):
                # Interpolate between base color and highlight
                blend = (1 - gradient_factor[j]) * rgb + gradient_factor[j] * highlight
                facecolors[i, j] = blend
        
        # Plot the sphere surface with the metallic gradient
        return self.ax.plot_surface(
            x, y, z,
            facecolors=facecolors,
            rstride=1, cstride=1,
            linewidth=0, antialiased=True,
            alpha=0.95  # Almost no transparency
        )
        
    def create_box_wireframe(self):
        """Create simulation box wireframe."""
        # Define box edges
        edges = [
            [[0, self.box_size], [0, 0], [0, 0]],
            [[0, self.box_size], [self.box_size, self.box_size], [0, 0]],
            [[0, self.box_size], [0, 0], [self.box_size, self.box_size]],
            [[0, self.box_size], [self.box_size, self.box_size], [self.box_size, self.box_size]],
            [[0, 0], [0, self.box_size], [0, 0]],
            [[self.box_size, self.box_size], [0, self.box_size], [0, 0]],
            [[0, 0], [0, self.box_size], [self.box_size, self.box_size]],
            [[self.box_size, self.box_size], [0, self.box_size], [self.box_size, self.box_size]],
            [[0, 0], [0, 0], [0, self.box_size]],
            [[self.box_size, self.box_size], [0, 0], [0, self.box_size]],
            [[0, 0], [self.box_size, self.box_size], [0, self.box_size]],
            [[self.box_size, self.box_size], [self.box_size, self.box_size], [0, self.box_size]]
        ]
        
        for edge in edges:
            self.ax.plot3D(*edge, 'k-', alpha=0.4, linewidth=1)
    
    def update_frame(self, frame_num):
        """Update visualization for given frame with 3D spheres."""
        self.current_frame = frame_num
        state = self.trajectory[frame_num]
        
        # Clear previous spheres
        for collection in self.sphere_collections.values():
            for sphere in collection:
                sphere.remove()
        self.sphere_collections.clear()
        
        # Clear previous legend
        for elem in self.legend_elements:
            if elem in self.ax.collections:
                elem.remove()
        self.legend_elements = []
        
        # For tracking legend entries
        legend_handles = []
        legend_labels = []
        
        # Plot particles for each type
        for type_name, positions in state['positions'].items():
            if len(positions) == 0:
                continue
            
            # Get particle radius from SystemParameters
            radius = self.particle_radii.get(type_name, 10.0)
            
            # Get particle color from color map or default to gray
            color = self.particle_colors.get(type_name, 'gray')
            
            # Limit number of spheres for performance
            max_spheres = 50  # Adjust based on your computer's performance
            plot_positions = positions
            if len(positions) > max_spheres:
                # Sample uniformly for consistent appearance
                indices = np.linspace(0, len(positions)-1, max_spheres, dtype=int)
                plot_positions = positions[indices]
                print(f"Note: Rendering {max_spheres} of {len(positions)} particles for type {type_name}")
            
            # Create collection for this particle type
            self.sphere_collections[type_name] = []
            
            # Add a dummy scatter point for legend (hidden)
            legend_point = self.ax.scatter(
                [-1000], [-1000], [-1000],  # Place off-screen
                c=color, 
                s=100,
                alpha=1.0
            )
            legend_handles.append(legend_point)
            legend_labels.append(f"{type_name}: r={radius}Å ({len(positions)} particles)")
            self.legend_elements.append(legend_point)
            
            # Create spheres for particles
            for pos in plot_positions:
                sphere = self.create_sphere(pos, radius, color)
                self.sphere_collections[type_name].append(sphere)
        
        # Update title with frame info
        self.title.set_text(
            f'Frame: {frame_num}/{len(self.trajectory)-1} | '
            f'Step: {state["step"]} | '
            f'Score: {state["total_score"]:.2f}'
        )
        
        # Update legend
        self.ax.legend(legend_handles, legend_labels, loc='upper right', framealpha=0.7)
        
        # Return list of artist objects for animation
        return self.legend_elements
    
    def create_animation(self, interval=200, save_path=None):
        """Create animated visualization."""
        anim = FuncAnimation(
            self.fig, 
            self.update_frame,
            frames=len(self.trajectory),
            interval=interval,
            blit=False,  # Cannot use blit with 3D plot_surface objects
            repeat=True
        )
        
        if save_path:
            print(f"Saving animation to {save_path}...")
            anim.save(save_path, writer='ffmpeg', fps=5, 
                     dpi=100, extra_args=['-vcodec', 'libx264'])
        
        return anim
    
    def show_interactive(self):
        """Show interactive plot with slider."""
        from matplotlib.widgets import Slider, Button
        
        # Adjust subplot to make room for slider
        plt.subplots_adjust(bottom=0.25)
        
        # Add slider
        ax_slider = plt.axes([0.1, 0.1, 0.65, 0.03])
        slider = Slider(
            ax_slider, 'Frame', 0, len(self.trajectory)-1, 
            valinit=0, valfmt='%d'
        )
        
        def update_slider(val):
            frame = int(slider.val)
            self.update_frame(frame)
            self.fig.canvas.draw_idle()  # Less CPU-intensive redraw
        
        slider.on_changed(update_slider)
        
        # Add play/pause buttons
        ax_play = plt.axes([0.8, 0.1, 0.05, 0.04])
        ax_pause = plt.axes([0.86, 0.1, 0.05, 0.04])
        
        btn_play = Button(ax_play, '▶️ Play')
        btn_pause = Button(ax_pause, '⏸️ Pause')
        
        self.animation = None
        
        def play_animation(event):
            if self.animation is None:
                self.animation = self.create_animation(interval=200)
            plt.show()
        
        def pause_animation(event):
            if self.animation:
                self.animation.event_source.stop()
        
        btn_play.on_clicked(play_animation)
        btn_pause.on_clicked(pause_animation)
        
        # Initialize with first frame
        self.update_frame(0)
        plt.show()

def main():
    """Main function to run the trajectory visualizer."""
    
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("Usage: python trajectory_visualizer.py <trajectory.h5> [max_frames] [step]")
        print("Example: python trajectory_visualizer.py output.h5 1000 2")
        sys.exit(1)
    
    trajectory_file = sys.argv[1]
    max_frames = int(sys.argv[2]) if len(sys.argv) > 2 else None
    step = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    
    print(f"🎬 Matplotlib Trajectory Visualizer (3D Spheres Edition)")
    print(f"=====================================================")
    print(f"📁 Loading trajectory: {trajectory_file}")
    if max_frames:
        print(f"🎯 Max frames: {max_frames}")
    if step > 1:
        print(f"⏭️  Sampling every {step} frames")
    
    # Load trajectory
    try:
        trajectory = load_trajectory_from_hdf5(trajectory_file, max_frames, step)
    except Exception as e:
        print(f"❌ Error loading trajectory: {e}")
        sys.exit(1)
    
    if len(trajectory) == 0:
        print("❌ No trajectory data loaded!")
        sys.exit(1)
    
    # Load system parameters directly from parameters.py
    print(f"🔧 Loading system parameters from parameters.py")
    params = SystemParameters()
    
    print(f"📏 Box size: {params.box_size} Å")
    print(f"⚛️ Particle configuration:")
    for type_name, radius in params.radii.items():
        print(f"  {type_name}: radius={radius}Å")
    
    # Create and run visualizer
    print(f"\n🎨 Creating visualization with {len(trajectory)} frames...")
    print(f"🔍 Using 3D spheres with metallic appearance")
    
    try:
        visualizer = MatplotlibTrajectoryVisualizer(trajectory, params)
        
        print(f"\n🎮 CONTROLS:")
        print(f"  🖱️  Mouse: Rotate/Zoom view")
        print(f"  📊 Slider: Drag to scrub through frames")
        print(f"  ▶️  Play button: Start auto-play animation")
        print(f"  ⏸️  Pause button: Stop auto-play")
        print(f"\n💡 NOTE: For better performance, large numbers of particles are subsampled")
        print(f"\n🎉 Starting interactive visualization...")
        
        visualizer.show_interactive()
        
    except Exception as e:
        print(f"❌ Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()