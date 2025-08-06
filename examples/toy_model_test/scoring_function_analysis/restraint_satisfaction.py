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
import h5py
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
def load_trajectory_from_hdf5(filename: str) -> list:
    """
    Loads MCMC trajectory data from an HDF5 file, matching the save_state_to_disk format.
    Safely handles missing keys and empty datasets, then sorts the states by step number.
    """
    all_states = []
    if not os.path.isfile(filename):
        print(f"File not found: {filename}")
        return all_states

    with h5py.File(filename, 'r') as f:
        # If the 'trajectory' group doesn't exist, bail out
        if 'trajectory' not in f:
            print(f"No 'trajectory' group in file: {filename}")
            return all_states
        
        traj_grp = f['trajectory']
        for state_name in traj_grp:
            state_grp = traj_grp[state_name]
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
            
            # Safely read sigma values
            if 'sigma' in state_grp:
                sigma_grp = state_grp['sigma']
                for key in sigma_grp.attrs:
                    state["sigma"][key] = float(sigma_grp.attrs[key])

            # Safely read positions
            if 'positions' in state_grp:
                pos_grp = state_grp['positions']
                for type_name in pos_grp:
                    state["positions"][type_name] = pos_grp[type_name][()]

            # Safely read types
            if 'types_keys' in state_grp and 'types_vals' in state_grp:
                types_keys_dataset = state_grp['types_keys'][()]
                types_vals_dataset = state_grp['types_vals'][()]
                if len(types_keys_dataset) == len(types_vals_dataset):
                    for i in range(len(types_keys_dataset)):
                        k = types_keys_dataset[i]
                        v = types_vals_dataset[i]
                        if isinstance(k, bytes):
                            k = k.decode('utf-8', errors='ignore')
                        if isinstance(v, bytes):
                            v = v.decode('utf-8', errors='ignore')
                        state["types"][k] = v
            
            # Safely read bead_numbers
            if 'bead_keys' in state_grp and 'bead_vals' in state_grp:
                bead_keys_dataset = state_grp['bead_keys'][()]
                bead_vals_dataset = state_grp['bead_vals'][()]
                if len(bead_keys_dataset) == len(bead_vals_dataset):
                    for i in range(len(bead_keys_dataset)):
                        bkey = bead_keys_dataset[i]
                        bval = bead_vals_dataset[i]
                        state["bead_numbers"][int(bkey)] = int(bval)

            all_states.append(state)

    # Sort the collected states by their step number
    all_states.sort(key=lambda s: s["step"])
    print(f"Loaded {len(all_states)} states from {filename}")

    return all_states
#------------------------------------------------------
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

def analyze_trajectory_restraints(filename: str, burnin_fraction=0.2, tolerance=0.1):
    """
    Load trajectory data and compute restraint satisfaction for each frame after burnin.
    
    Args:
        filename: Path to HDF5 trajectory file
        burnin_fraction: Fraction of trajectory to discard as burnin (default 0.2 = 20%)
        tolerance: Fractional tolerance to consider a restraint satisfied (±10%)
    """
    # Create output directory
    output_dir = "output_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load all states from trajectory
    all_states = load_trajectory_from_hdf5(filename)
    
    if not all_states:
        print("No states found in trajectory file")
        return None
    
    # Calculate burnin cutoff
    total_states = len(all_states)
    burnin_cutoff = int(burnin_fraction * total_states)
    
    # Keep only post-burnin states
    post_burnin_states = all_states[burnin_cutoff:]
    
    print(f"Total states: {total_states}")
    print(f"Burnin fraction: {burnin_fraction} ({burnin_cutoff} states)")
    print(f"Analyzing {len(post_burnin_states)} post-burnin states...")
    
    if not post_burnin_states:
        print("No states remaining after burnin - try reducing burnin_fraction")
        return None
    
    # Get reference distances
    params = SystemParameters()
    ideal_aa = params.pair_distances['AA']
    ideal_ab = params.pair_distances['AB']
    ideal_bc = params.pair_distances['BC']
    
    # Storage for results
    results = {
        'step': [],
        'frame_number': [],
        'aa_satisfied_percent': [],
        'ab_satisfied_percent': [],
        'bc_satisfied_percent': [],
        'total_satisfied_percent': [],
        'aa_mean_distance': [],
        'ab_mean_distance': [],
        'bc_mean_distance': [],
        'aa_std_distance': [],
        'ab_std_distance': [],
        'bc_std_distance': []
    }
    
    # Process each post-burnin state
    for frame_idx, state in enumerate(post_burnin_states):
        step = state['step']
        coords = state['positions']
        
        print(f"Processing frame {frame_idx+1}/{len(post_burnin_states)} (step {step})...")
        
        # Check if we have required data
        if not coords:
            print(f"  Skipping step {step} - missing coordinates")
            continue
            
        if not all(ptype in coords for ptype in ['A', 'B', 'C']):
            print(f"  Skipping step {step} - missing particle types")
            continue
        
        try:
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
            aa_satisfied = np.sum(np.abs(aa_dists - ideal_aa) <= tolerance * ideal_aa)
            ab_satisfied = np.sum(np.abs(ab_dists - ideal_ab) <= tolerance * ideal_ab)
            bc_satisfied = np.sum(np.abs(bc_dists - ideal_bc) <= tolerance * ideal_bc)
            
            # Calculate percentages
            aa_satisfied_pct = 100 * aa_satisfied / len(aa_dists)
            ab_satisfied_pct = 100 * ab_satisfied / len(ab_dists)
            bc_satisfied_pct = 100 * bc_satisfied / len(bc_dists)
            
            # Total satisfied percentage (weighted by number of restraints)
            total_restraints = len(aa_dists) + len(ab_dists) + len(bc_dists)
            total_satisfied = aa_satisfied + ab_satisfied + bc_satisfied
            total_satisfied_pct = 100 * total_satisfied / total_restraints
            
            # Store results
            results['step'].append(step)
            results['frame_number'].append(frame_idx + 1)
            results['aa_satisfied_percent'].append(aa_satisfied_pct)
            results['ab_satisfied_percent'].append(ab_satisfied_pct)
            results['bc_satisfied_percent'].append(bc_satisfied_pct)
            results['total_satisfied_percent'].append(total_satisfied_pct)
            results['aa_mean_distance'].append(np.mean(aa_dists))
            results['ab_mean_distance'].append(np.mean(ab_dists))
            results['bc_mean_distance'].append(np.mean(bc_dists))
            results['aa_std_distance'].append(np.std(aa_dists))
            results['ab_std_distance'].append(np.std(ab_dists))
            results['bc_std_distance'].append(np.std(bc_dists))
            
            print(f"  Frame {frame_idx+1}: AA={aa_satisfied_pct:.1f}%, AB={ab_satisfied_pct:.1f}%, BC={bc_satisfied_pct:.1f}%, Total={total_satisfied_pct:.1f}%")
            
        except Exception as e:
            print(f"  Error processing step {step}: {e}")
            continue
    
    print(f"\nSuccessfully processed {len(results['step'])} post-burnin frames")
    
    # Convert to DataFrame for easy analysis
    import pandas as pd
    df = pd.DataFrame(results)
    
    # Save results with burnin info in filename
    output_file = os.path.join(output_dir, f"trajectory_restraints_burnin{int(burnin_fraction*100)}.csv")
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
    # Create comprehensive plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Restraint satisfaction vs frame
    axes[0, 0].plot(df['frame_number'], df['aa_satisfied_percent'], 'o-', label='A-A', alpha=0.7)
    axes[0, 0].plot(df['frame_number'], df['ab_satisfied_percent'], 's-', label='A-B', alpha=0.7)
    axes[0, 0].plot(df['frame_number'], df['bc_satisfied_percent'], '^-', label='B-C', alpha=0.7)
    axes[0, 0].plot(df['frame_number'], df['total_satisfied_percent'], 'D-', label='Total', alpha=0.9, linewidth=2)
    axes[0, 0].set_xlabel('Frame Number')
    axes[0, 0].set_ylabel('Restraints Satisfied (%)')
    axes[0, 0].set_title(f'Restraint Satisfaction vs Frame (±{tolerance*100}% tolerance)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Mean distances vs frame
    axes[0, 1].plot(df['frame_number'], df['aa_mean_distance'], 'o-', label=f'A-A (ideal: {ideal_aa:.1f})', alpha=0.7)
    axes[0, 1].plot(df['frame_number'], df['ab_mean_distance'], 's-', label=f'A-B (ideal: {ideal_ab:.1f})', alpha=0.7)
    axes[0, 1].plot(df['frame_number'], df['bc_mean_distance'], '^-', label=f'B-C (ideal: {ideal_bc:.1f})', alpha=0.7)
    axes[0, 1].axhline(ideal_aa, color='blue', linestyle='--', alpha=0.5)
    axes[0, 1].axhline(ideal_ab, color='orange', linestyle='--', alpha=0.5)
    axes[0, 1].axhline(ideal_bc, color='green', linestyle='--', alpha=0.5)
    axes[0, 1].set_xlabel('Frame Number')
    axes[0, 1].set_ylabel('Mean Distance (Å)')
    axes[0, 1].set_title('Mean Distances vs Frame')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Standard deviation of distances
    axes[0, 2].plot(df['frame_number'], df['aa_std_distance'], 'o-', label='A-A', alpha=0.7)
    axes[0, 2].plot(df['frame_number'], df['ab_std_distance'], 's-', label='A-B', alpha=0.7)
    axes[0, 2].plot(df['frame_number'], df['bc_std_distance'], '^-', label='B-C', alpha=0.7)
    axes[0, 2].set_xlabel('Frame Number')
    axes[0, 2].set_ylabel('Std Dev of Distance (Å)')
    axes[0, 2].set_title('Distance Variability vs Frame')
    axes[0, 2].legend()
    axes[0, 2].grid(True, alpha=0.3)
    
    # Plot 4: Distribution of total satisfaction
    axes[1, 0].hist(df['total_satisfied_percent'], bins=20, alpha=0.7, color='purple')
    axes[1, 0].axvline(df['total_satisfied_percent'].mean(), color='red', linestyle='--', 
                      label=f'Mean: {df["total_satisfied_percent"].mean():.1f}%')
    axes[1, 0].set_xlabel('Total Restraints Satisfied (%)')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Total Restraint Satisfaction')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 5: Running average of total satisfaction
    window_size = max(1, len(df) // 10)
    running_avg = df['total_satisfied_percent'].rolling(window=window_size).mean()
    axes[1, 1].plot(df['frame_number'], df['total_satisfied_percent'], 'o-', alpha=0.3, label='Individual')
    axes[1, 1].plot(df['frame_number'], running_avg, '-', linewidth=2, label=f'Running avg (window={window_size})')
    axes[1, 1].set_xlabel('Frame Number')
    axes[1, 1].set_ylabel('Total Restraints Satisfied (%)')
    axes[1, 1].set_title('Convergence of Restraint Satisfaction')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    # Plot 6: Correlation matrix of satisfaction percentages
    corr_data = df[['aa_satisfied_percent', 'ab_satisfied_percent', 'bc_satisfied_percent']].corr()
    im = axes[1, 2].imshow(corr_data, cmap='coolwarm', vmin=-1, vmax=1)
    axes[1, 2].set_xticks(range(3))
    axes[1, 2].set_yticks(range(3))
    axes[1, 2].set_xticklabels(['A-A', 'A-B', 'B-C'])
    axes[1, 2].set_yticklabels(['A-A', 'A-B', 'B-C'])
    axes[1, 2].set_title('Restraint Satisfaction Correlation')
    
    # Add correlation values to the plot
    for i in range(3):
        for j in range(3):
            axes[1, 2].text(j, i, f'{corr_data.iloc[i, j]:.2f}', 
                           ha='center', va='center', fontweight='bold')
    
    plt.colorbar(im, ax=axes[1, 2])
    
    # Overall title
    fig.suptitle(f'Trajectory Restraint Analysis (burnin={burnin_fraction}, tolerance=±{tolerance*100}%)', 
                fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    
    # Save high-quality plot
    plot_file = os.path.join(output_dir, f'trajectory_restraints_burnin{int(burnin_fraction*100)}.pdf')
    plt.savefig(plot_file, format='pdf', dpi=300, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    
    # Print summary statistics
    print(f"\n=== Restraint Satisfaction Summary (burnin={burnin_fraction}) ===")
    print(f"A-A Satisfaction: mean={df['aa_satisfied_percent'].mean():.1f}%, std={df['aa_satisfied_percent'].std():.1f}%")
    print(f"A-B Satisfaction: mean={df['ab_satisfied_percent'].mean():.1f}%, std={df['ab_satisfied_percent'].std():.1f}%")
    print(f"B-C Satisfaction: mean={df['bc_satisfied_percent'].mean():.1f}%, std={df['bc_satisfied_percent'].std():.1f}%")
    print(f"Total Satisfaction: mean={df['total_satisfied_percent'].mean():.1f}%, std={df['total_satisfied_percent'].std():.1f}%")
    
    plt.show()
    
    return df

# Usage example
if __name__ == "__main__":
    # Specify your trajectory file
    trajectory_file = "../output_analysis/tetramersampler_results_1/trajectory_chain_8.h5"  # Update this path
    
    # Analyze restraint satisfaction over trajectory
    restraint_results = analyze_trajectory_restraints(
        trajectory_file, 
        burnin_fraction=0.5, 
        tolerance=0.01
    )
    
    if restraint_results is not None:
        print(f"\nProcessed {len(restraint_results)} frames successfully!")
    
    # Run the analysis
    stats = analyze_restraints(ideal_coords)
    print("\nRestraint Satisfaction Summary:")
    print(f"A-A: {stats['Satisfied_Percent'][0]:.1f}% satisfied")
    print(f"A-B: {stats['Satisfied_Percent'][1]:.1f}% satisfied") 
    print(f"B-C: {stats['Satisfied_Percent'][2]:.1f}% satisfied")