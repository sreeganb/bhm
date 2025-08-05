#------------------------------------------------------------------
# Script to compute the number of restraints that are satisfied
# for a given set of coordinates.
#-------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd 
import os
import matplotlib.pyplot as plt
import seaborn as sns
#-------------------------------------------------------------------
# Just add the parent directory to Python's path - that's it!
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

import parameters
import pair_sampler
import tetramer_sampler
import octet_sampler

from parameters import SystemParameters
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from octet_sampler import OctetSampler
#-----------------------------------------------------
# Instantiate the parameters object
#-----------------------------------------------------
params = SystemParameters()
# get the ideal coordinates
ideal_coords = params.latest_ideal()
print(f"Ideal coordinates: {ideal_coords}")
#-----------------------------------------------------
def analyze_restraints(coords, tolerance=0.1):
    """
    Analyze distances and restraint satisfaction with modern visualizations.
    
    Args:
        coords: Dictionary with coordinates for A, B, C particles
        tolerance: Fractional tolerance to consider a restraint satisfied (±10%)
    """
    # Create output directory
    output_dir = "output_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    # Get reference distances
    params = SystemParameters()
    ideal_aa = params.pair_distances['AA']
    ideal_ab = params.pair_distances['AB']
    ideal_bc = params.pair_distances['BC']
    
    # Extract coordinates
    a_coords = coords['A']
    b_coords = coords['B']
    c_coords = coords['C']
    
    # Calculate distances
    aa_dists = np.array([np.linalg.norm(a_coords[i] - a_coords[j]) 
                for i in range(len(a_coords)) 
                for j in range(i+1, len(a_coords))])
    
    ab_dists = np.array([np.linalg.norm(a_coords[i] - b_coords[j]) 
                for i in range(len(a_coords)) 
                for j in range(len(b_coords))])
    
    bc_dists = np.array([np.linalg.norm(b_coords[i] - c_coords[j]) 
                for i in range(len(b_coords)) 
                for j in range(len(c_coords))])
    
    # Calculate restraint satisfaction
    aa_satisfied = np.sum(np.abs(aa_dists - ideal_aa) <=  ideal_aa)
    ab_satisfied = np.sum(np.abs(ab_dists - ideal_ab) <=  ideal_ab)
    bc_satisfied = np.sum(np.abs(bc_dists - ideal_bc) <=  ideal_bc)
    
    # Create statistics dictionary
    stats_dict = {
        'Pair_Type': ['A-A', 'A-B', 'B-C'],
        'Ideal_Distance': [ideal_aa, ideal_ab, ideal_bc],
        'Count': [len(aa_dists), len(ab_dists), len(bc_dists)],
        'Mean': [np.mean(aa_dists), np.mean(ab_dists), np.mean(bc_dists)],
        'Std_Dev': [np.std(aa_dists), np.std(ab_dists), np.std(bc_dists)],
        'Min': [np.min(aa_dists), np.min(ab_dists), np.min(bc_dists)],
        'Max': [np.max(aa_dists), np.max(ab_dists), np.max(bc_dists)],
        'Satisfied_Count': [aa_satisfied, ab_satisfied, bc_satisfied],
        'Satisfied_Percent': [
            100 * aa_satisfied / len(aa_dists), 
            100 * ab_satisfied / len(ab_dists), 
            100 * bc_satisfied / len(bc_dists)
        ]
    }
    
    # Create and save statistics DataFrame
    stats_df = pd.DataFrame(stats_dict)
    csv_path = os.path.join(output_dir, 'restraint_statistics.csv')
    stats_df.to_csv(csv_path, index=False)
    print(f"Statistics saved to {csv_path}")
    
    # Set seaborn style for modern plots
    sns.set_theme(style="whitegrid")
    
    # Create a single figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Data to plot and plot settings
    all_dists = [aa_dists, ab_dists, bc_dists]
    titles = ['A-A Distances', 'A-B Distances', 'B-C Distances']
    colors = ['dodgerblue', 'forestgreen', 'darkorange']
    ideal_vals = [ideal_aa, ideal_ab, ideal_bc]
    sat_pcts = stats_df['Satisfied_Percent'].values
    
    # Create plots with a simple loop
    for i, (dist, title, color, ideal, sat_pct) in enumerate(zip(all_dists, titles, colors, ideal_vals, sat_pcts)):
        # Modern histogram with KDE overlay
        # Modern histogram with KDE overlay - with more bins
        sns.histplot(dist, kde=True, color=color, alpha=0.7, ax=axes[i], bins=30)  # Try different values (25-40)
        
        # Add ideal line
        axes[i].axvline(ideal, color='red', linestyle='--', linewidth=2, 
                       label=f'Ideal: {ideal:.2f} Å')
        
        # Add stats as text
        stats_text = (f"Mean: {stats_dict['Mean'][i]:.2f} Å\n"
                      f"Std: {stats_dict['Std_Dev'][i]:.2f} Å\n"
                      f"Satisfied: {sat_pct:.1f}%")
        
        axes[i].text(0.05, 0.95, stats_text, transform=axes[i].transAxes,
                    va='top', bbox=dict(facecolor='white', alpha=0.7))
        
        # Customize plot
        axes[i].set_title(title, fontsize=12, fontweight='bold')
        axes[i].set_xlabel('Distance (Å)')
        if i == 0:  # Only add y-label to first plot
            axes[i].set_ylabel('Frequency')
        else:
            axes[i].set_ylabel('')
        axes[i].legend()
    
    # Add overall title
    fig.suptitle(f'Distance Distributions and Restraint Satisfaction (±{tolerance*100}% tolerance)', 
                fontsize=14, fontweight='bold', y=1.05)
    
    plt.tight_layout()
    
    # Save high-quality PDF
    pdf_path = os.path.join(output_dir, 'distance_analysis.pdf')
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    print(f"Plot saved to {pdf_path}")
    
    # Return the statistics for further analysis
    return stats_df

# Run the analysis
stats = analyze_restraints(ideal_coords)
print("\nRestraint Satisfaction Summary:")
print(f"A-A: {stats['Satisfied_Percent'][0]:.1f}% satisfied")
print(f"A-B: {stats['Satisfied_Percent'][1]:.1f}% satisfied") 
print(f"B-C: {stats['Satisfied_Percent'][2]:.1f}% satisfied")