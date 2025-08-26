#!/usr/bin/env python3
"""
Multi-Score Visualization: 3D structure + 4-panel RMSD vs Score plots with play controls.
Compatible with ALL scoring functions format (CCC+ExVol, Pair, Tetramer, Octet).
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button
import numpy as np
import h5py
import sys
import os

class MultiScoreVisualization:
    def __init__(self, h5_file):
        self.h5_file = h5_file
        self.frames_data = []
        self.current_frame = 0
        
        # Particle properties
        self.radii = {'A': 24.0, 'B': 14.0, 'C': 16.0}
        self.colors = {'A': '#FF3B30', 'B': '#147EFB', 'C': '#28CD41'}
        
        # Score types and their display properties
        self.score_types = {
            'ccc_exvol': {'name': 'CCC + ExVol', 'color': '#FF6B6B'},
            'pair': {'name': 'Pair Score', 'color': '#4ECDC4'},
            'tetramer': {'name': 'Tetramer Score', 'color': '#45B7D1'},
            'octet': {'name': 'Octet Score', 'color': '#96CEB4'}
        }
        
        self.load_data()
        
        # Set up plots
        self.fig = None
        self.ax_3d = None
        self.ax_plots = {}  # Dictionary to hold score plot axes
        self.scatter_objects = {}
        self.current_markers = {}  # Dictionary to hold current frame markers for each plot
        self.animation = None
        
    def load_data(self):
        """Load data from H5 file with ALL scoring functions format."""
        print(f"Loading data from {self.h5_file}...")
        
        with h5py.File(self.h5_file, 'r') as f:
            num_frames = f.attrs.get('num_frames', 0)
            has_reference = f.attrs.get('has_reference_frame', False)
            
            # Handle both string and bytes for score_types
            score_types_attr = f.attrs.get('score_types', '')
            if isinstance(score_types_attr, bytes):
                score_types_str = score_types_attr.decode('utf-8')
            else:
                score_types_str = str(score_types_attr)
            
            print(f"Found {num_frames} frames (includes reference: {has_reference})")
            print(f"Available score types: {score_types_str}")
            
            for frame_idx in range(num_frames):
                frame_name = f'frame_{frame_idx:04d}'
                
                # Load coordinates
                coords_A = f[f'array_A/{frame_name}'][:]
                coords_B = f[f'array_B/{frame_name}'][:]
                coords_C = f[f'array_C/{frame_name}'][:]
                
                # Load metadata
                frame_meta = f[f'metadata/{frame_name}']
                
                # Get all scores with new naming convention
                scores = {}
                for score_type in self.score_types.keys():
                    score_key = f'score_{score_type}'
                    scores[score_type] = float(frame_meta.attrs.get(score_key, 0.0))
                
                # Also get individual CCC and ExVol if available
                scores['ccc_only'] = float(frame_meta.attrs.get('score_ccc_only', 0.0))
                scores['exvol_only'] = float(frame_meta.attrs.get('score_exvol_only', 0.0))
                
                # Get RMSD
                rmsd = frame_meta.attrs.get('rmsd', 0.0)
                rmsd_aligned = frame_meta.attrs.get('rmsd_aligned', rmsd)
                rmsd_raw = frame_meta.attrs.get('rmsd_raw', rmsd)
                
                # Check if this is the reference structure
                is_reference = frame_meta.attrs.get('is_reference', False)
                intensity = frame_meta.attrs.get('intensity_requested', 0.0)
                
                frame_data = {
                    'coordinates': {'A': coords_A, 'B': coords_B, 'C': coords_C},
                    'scores': scores,
                    'rmsd': float(rmsd_aligned),  # Use aligned RMSD for plotting
                    'rmsd_aligned': float(rmsd_aligned),
                    'rmsd_raw': float(rmsd_raw),
                    'is_reference': bool(is_reference),
                    'intensity': float(intensity)
                }
                
                self.frames_data.append(frame_data)
        
        print(f"Loaded {len(self.frames_data)} frames")
        
        # Print reference structure info
        if self.frames_data and self.frames_data[0]['is_reference']:
            ref = self.frames_data[0]
            print(f"Reference structure (frame 0) scores:")
            for score_type, score_info in self.score_types.items():
                print(f"  {score_info['name']}: {ref['scores'][score_type]:.2f}")
    
    def setup_plots(self):
        """Set up the multi-panel layout: 3D structure + 4 score plots."""
        self.fig = plt.figure(figsize=(20, 12))
        
        # 3D structure plot (top left, spans 2 rows)
        self.ax_3d = plt.subplot2grid((3, 4), (0, 0), rowspan=2, projection='3d')
        
        # Score plots (2x2 grid on the right)
        plot_positions = [(0, 2), (0, 3), (1, 2), (1, 3)]
        
        for i, (score_type, score_info) in enumerate(self.score_types.items()):
            row, col = plot_positions[i]
            self.ax_plots[score_type] = plt.subplot2grid((3, 4), (row, col))
        
        self.setup_3d_scene()
        self.setup_score_plots()
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
    
    def setup_score_plots(self):
        """Set up RMSD vs Score plots for all scoring functions."""
        # Separate reference and perturbations
        reference_frames = [f for f in self.frames_data if f['is_reference']]
        perturbation_frames = [f for f in self.frames_data if not f['is_reference']]
        
        for score_type, score_info in self.score_types.items():
            ax = self.ax_plots[score_type]
            
            # Plot perturbations
            if perturbation_frames:
                pert_rmsds = [f['rmsd'] for f in perturbation_frames]
                pert_scores = [f['scores'][score_type] for f in perturbation_frames]
                ax.scatter(pert_rmsds, pert_scores, alpha=0.6, s=30, 
                          color='lightgray', label=f'Perturbations ({len(perturbation_frames)})', zorder=1)
            
            # Plot reference structure(s)
            if reference_frames:
                ref_rmsds = [f['rmsd'] for f in reference_frames]
                ref_scores = [f['scores'][score_type] for f in reference_frames]
                ax.scatter(ref_rmsds, ref_scores, color='red', s=150, marker='*', 
                          label='Reference', edgecolor='black', linewidth=1.5, zorder=3)
            
            # Current frame marker (will be updated)
            self.current_markers[score_type] = ax.scatter([0], [0], color=score_info['color'], s=100, 
                                                         marker='o', label='Current', 
                                                         edgecolor='black', linewidth=1.5, zorder=2)

            ax.set_xlabel('RMSD (Å)', fontsize=10)
            ax.set_ylabel(f'{score_info["name"]}', fontsize=10)
            ax.set_title(f'RMSD vs {score_info["name"]}', fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            # Add correlation (only for perturbations)
            if len(perturbation_frames) > 1:
                pert_rmsds = [f['rmsd'] for f in perturbation_frames]
                pert_scores = [f['scores'][score_type] for f in perturbation_frames]
                correlation = np.corrcoef(pert_rmsds, pert_scores)[0, 1]
                ax.text(0.05, 0.95, f'r = {correlation:.3f}', 
                       transform=ax.transAxes, fontsize=9,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Add score range info
            if self.frames_data:
                all_scores = [f['scores'][score_type] for f in self.frames_data]
                score_range = [min(all_scores), max(all_scores)]
                ax.text(0.05, 0.85, f'Range: {score_range[0]:.1f} - {score_range[1]:.1f}', 
                       transform=ax.transAxes, fontsize=8,
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
        
        # Add legend to the first plot only
        first_ax = list(self.ax_plots.values())[0]
        first_ax.legend(fontsize=8, loc='upper right')
    
    def setup_controls(self):
        """Set up play controls."""
        plt.subplots_adjust(bottom=0.15)
        
        # Frame slider (spans most of the bottom)
        ax_slider = plt.axes([0.1, 0.05, 0.6, 0.03])
        self.frame_slider = Slider(ax_slider, 'Frame', 0, len(self.frames_data)-1, 
                                  valinit=0, valfmt='%d')
        self.frame_slider.on_changed(self.update_frame_from_slider)
        
        # Play/Pause buttons
        ax_play = plt.axes([0.75, 0.05, 0.04, 0.04])
        ax_pause = plt.axes([0.8, 0.05, 0.04, 0.04])
        
        self.btn_play = Button(ax_play, '▶️')
        self.btn_pause = Button(ax_pause, '⏸️')
        
        self.btn_play.on_clicked(self.play_animation)
        self.btn_pause.on_clicked(self.pause_animation)
    
    def update_frame(self, frame_idx):
        """Update all plots for given frame."""
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
                    s=self.radii[type_name] * 6,
                    alpha=0.8,
                    label=f'{type_name}'
                )
                self.scatter_objects[type_name] = scatter
        
        # Create detailed title with all score breakdown
        frame_type = "REF" if frame_data['is_reference'] else "PERT"
        title_parts = [
            f'Frame {frame_idx}/{len(self.frames_data)-1} ({frame_type})',
            f'RMSD: {frame_data["rmsd"]:.2f}Å'
        ]
        
        # Add key scores to title
        for score_type, score_info in self.score_types.items():
            score_val = frame_data['scores'][score_type]
            title_parts.append(f'{score_type.upper()}: {score_val:.1f}')
        
        if not frame_data['is_reference']:
            title_parts.append(f'Int: {frame_data["intensity"]:.1f}')
        
        # Split title into multiple lines if too long
        title_line1 = ' | '.join(title_parts[:3])
        title_line2 = ' | '.join(title_parts[3:])
        full_title = title_line1 + '\n' + title_line2 if len(title_parts) > 3 else title_line1
        
        self.ax_3d.set_title(full_title, fontsize=9, fontweight='bold')
        
        # Update current markers in all score plots
        for score_type in self.score_types.keys():
            marker = self.current_markers[score_type]
            score_val = frame_data['scores'][score_type]
            marker.set_offsets([[frame_data['rmsd'], score_val]])
        
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
                interval=800,  # 0.8 seconds per frame (slower for 4 plots)
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
        
        print("🎮 MULTI-SCORE ANALYSIS CONTROLS:")
        print("  🖱️  Mouse: Rotate/Zoom 3D view")
        print("  📊 Slider: Scrub through frames")
        print("  ▶️  Play: Start animation") 
        print("  ⏸️  Pause: Stop animation")
        print("  🔶 Colored dots: Current frame position in each score plot")
        print("  ⭐ Red star: Reference structure")
        print("  ⚪ Gray dots: Perturbations")
        print(f"  📈 Analyzing {len(self.score_types)} scoring functions")
        
        plt.tight_layout()
        plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python simplified_visualization.py <perturbation_analysis_all_scores.h5>")
        sys.exit(1)
    
    h5_file = sys.argv[1]
    
    if not os.path.exists(h5_file):
        print(f"Error: File {h5_file} not found!")
        sys.exit(1)
    
    print("🎬 Multi-Score Structure Visualization")
    print("=====================================")
    
    visualizer = MultiScoreVisualization(h5_file)
    visualizer.show()

if __name__ == "__main__":
    main()