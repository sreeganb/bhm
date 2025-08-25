#!/usr/bin/env python3
"""
Simple 2-panel visualization: 3D structure + RMSD vs Score plot with play controls.
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button
import numpy as np
import h5py
import sys
import os

class SimpleVisualization:
    def __init__(self, h5_file):
        self.h5_file = h5_file
        self.frames_data = []
        self.current_frame = 0
        
        # Particle properties
        self.radii = {'A': 24.0, 'B': 14.0, 'C': 16.0}
        self.colors = {'A': '#FF3B30', 'B': '#147EFB', 'C': '#28CD41'}
        
        self.load_data()
        
        # Set up plots
        self.fig = None
        self.ax_3d = None
        self.ax_plot = None
        self.scatter_objects = {}
        self.animation = None
        
    def load_data(self):
        """Load data from H5 file."""
        print(f"Loading data from {self.h5_file}...")
        
        with h5py.File(self.h5_file, 'r') as f:
            num_frames = f.attrs.get('num_frames', 0)
            print(f"Found {num_frames} frames")
            
            for frame_idx in range(num_frames):
                frame_name = f'frame_{frame_idx:04d}'
                
                # Load coordinates
                coords_A = f[f'array_A/{frame_name}'][:]
                coords_B = f[f'array_B/{frame_name}'][:]
                coords_C = f[f'array_C/{frame_name}'][:]
                
                # Load metadata
                frame_meta = f[f'metadata/{frame_name}']
                
                # Get scores
                if 'total_score' in frame_meta.attrs:
                    total_score = frame_meta.attrs['total_score']
                else:
                    total_score = frame_meta.attrs.get('score', 0.0)
                
                rmsd = frame_meta.attrs['rmsd']
                
                frame_data = {
                    'coordinates': {'A': coords_A, 'B': coords_B, 'C': coords_C},
                    'total_score': float(total_score),
                    'rmsd': float(rmsd)
                }
                
                self.frames_data.append(frame_data)
        
        print(f"Loaded {len(self.frames_data)} frames")
    
    def setup_plots(self):
        """Set up the two-panel layout."""
        self.fig = plt.figure(figsize=(16, 8))
        
        # 3D structure plot (left panel)
        self.ax_3d = self.fig.add_subplot(1, 2, 1, projection='3d')
        
        # RMSD vs Score plot (right panel)
        self.ax_plot = self.fig.add_subplot(1, 2, 2)
        
        self.setup_3d_scene()
        self.setup_rmsd_plot()
        self.setup_controls()
    
    def setup_3d_scene(self):
        """Set up 3D scene."""
        # Calculate bounds
        all_coords = []
        for frame_data in self.frames_data:
            for coords in frame_data['coordinates'].values():
                all_coords.extend(coords)
        
        all_coords = np.array(all_coords)
        margin = 20.0
        
        x_min, x_max = all_coords[:, 0].min() - margin, all_coords[:, 0].max() + margin
        y_min, y_max = all_coords[:, 1].min() - margin, all_coords[:, 1].max() + margin
        z_min, z_max = all_coords[:, 2].min() - margin, all_coords[:, 2].max() + margin
        
        self.ax_3d.set_xlim(x_min, x_max)
        self.ax_3d.set_ylim(y_min, y_max)
        self.ax_3d.set_zlim(z_min, z_max)
        
        self.ax_3d.set_xlabel('X (Å)', fontsize=12)
        self.ax_3d.set_ylabel('Y (Å)', fontsize=12)
        self.ax_3d.set_zlabel('Z (Å)', fontsize=12)
        self.ax_3d.set_title('3D Structure', fontsize=14, fontweight='bold')
        self.ax_3d.view_init(elev=20, azim=45)
    
    def setup_rmsd_plot(self):
        """Set up RMSD vs Score plot."""
        # Extract all data for plotting
        rmsds = [frame['rmsd'] for frame in self.frames_data]
        scores = [frame['total_score'] for frame in self.frames_data]
        
        # Plot all points
        self.ax_plot.scatter(rmsds, scores, alpha=0.6, s=50, color='lightblue', 
                           label='All Structures', zorder=1)
        
        # Add ideal structure (RMSD=0, lowest score)
        ideal_score = min(scores) if scores else 0
        self.ax_plot.scatter([0], [ideal_score], color='red', s=200, marker='*', 
                           label='Ideal Structure', edgecolor='black', linewidth=2, zorder=3)
        
        # Current frame marker (will be updated)
        self.current_marker = self.ax_plot.scatter([0], [0], color='orange', s=150, 
                                                  marker='o', label='Current Frame', 
                                                  edgecolor='black', linewidth=2, zorder=2)
        
        self.ax_plot.set_xlabel('RMSD (Å)', fontsize=12)
        self.ax_plot.set_ylabel('Total Score', fontsize=12)
        self.ax_plot.set_title('RMSD vs Total Score', fontsize=14, fontweight='bold')
        self.ax_plot.grid(True, alpha=0.3)
        self.ax_plot.legend(fontsize=10)
        
        # Add correlation
        if len(rmsds) > 1:
            correlation = np.corrcoef(rmsds, scores)[0, 1]
            self.ax_plot.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                            transform=self.ax_plot.transAxes, fontsize=11,
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    def setup_controls(self):
        """Set up play controls."""
        plt.subplots_adjust(bottom=0.15)
        
        # Frame slider
        ax_slider = plt.axes([0.1, 0.05, 0.65, 0.03])
        self.frame_slider = Slider(ax_slider, 'Frame', 0, len(self.frames_data)-1, 
                                  valinit=0, valfmt='%d')
        self.frame_slider.on_changed(self.update_frame_from_slider)
        
        # Play/Pause buttons
        ax_play = plt.axes([0.8, 0.05, 0.05, 0.04])
        ax_pause = plt.axes([0.86, 0.05, 0.05, 0.04])
        
        self.btn_play = Button(ax_play, '▶️')
        self.btn_pause = Button(ax_pause, '⏸️')
        
        self.btn_play.on_clicked(self.play_animation)
        self.btn_pause.on_clicked(self.pause_animation)
    
    def update_frame(self, frame_idx):
        """Update both plots for given frame."""
        if frame_idx >= len(self.frames_data):
            return
        
        self.current_frame = frame_idx
        frame_data = self.frames_data[frame_idx]
        
        # Clear previous 3D scatter plots
        for scatter in self.scatter_objects.values():
            scatter.remove()
        self.scatter_objects.clear()
        
        # Plot 3D structure
        for type_name in ['A', 'B', 'C']:
            coords = frame_data['coordinates'][type_name]
            if len(coords) > 0:
                scatter = self.ax_3d.scatter(
                    coords[:, 0], coords[:, 1], coords[:, 2],
                    color=self.colors[type_name],
                    s=self.radii[type_name] * 8,
                    alpha=0.8,
                    label=f'{type_name}'
                )
                self.scatter_objects[type_name] = scatter
        
        # Update 3D title
        self.ax_3d.set_title(
            f'Frame {frame_idx}/{len(self.frames_data)-1} | '
            f'RMSD: {frame_data["rmsd"]:.2f}Å | '
            f'Score: {frame_data["total_score"]:.1f}',
            fontsize=12, fontweight='bold'
        )
        
        # Update current marker in RMSD plot
        self.current_marker.set_offsets([[frame_data['rmsd'], frame_data['total_score']]])
        
        # Update slider
        if abs(self.frame_slider.val - frame_idx) > 0.1:
            self.frame_slider.set_val(frame_idx)
        
        self.fig.canvas.draw_idle()
    
    def update_frame_from_slider(self, val):
        """Update frame from slider."""
        frame_idx = int(self.frame_slider.val)
        if frame_idx != self.current_frame:
            self.update_frame(frame_idx)
    
    def play_animation(self, event):
        """Start animation."""
        if self.animation is None or not hasattr(self.animation, 'event_source') or not self.animation.event_source:
            self.animation = FuncAnimation(
                self.fig, 
                lambda frame: self.update_frame(frame % len(self.frames_data)),
                frames=len(self.frames_data),
                interval=500,  # 0.5 seconds per frame
                repeat=True
            )
    
    def pause_animation(self, event):
        """Pause animation."""
        if self.animation and hasattr(self.animation, 'event_source') and self.animation.event_source:
            self.animation.event_source.stop()
            self.animation = None
    
    def show(self):
        """Show the visualization."""
        self.setup_plots()
        self.update_frame(0)
        
        print("🎮 CONTROLS:")
        print("  🖱️  Mouse: Rotate/Zoom 3D view")
        print("  📊 Slider: Scrub through frames")
        print("  ▶️  Play: Start animation")
        print("  ⏸️  Pause: Stop animation")
        print("  🔶 Orange dot: Current frame position")
        
        plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python simple_visualization.py <perturbation_analysis.h5>")
        sys.exit(1)
    
    h5_file = sys.argv[1]
    
    if not os.path.exists(h5_file):
        print(f"Error: File {h5_file} not found!")
        sys.exit(1)
    
    print("🎬 Simple Structure Visualization")
    print("================================")
    
    visualizer = SimpleVisualization(h5_file)
    visualizer.show()

if __name__ == "__main__":
    main()