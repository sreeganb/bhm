#!/usr/bin/env python3
"""
Simple restraint satisfaction analysis - compact version.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from parameters import SystemParameters

def load_trajectory(filename):
    """Load trajectory from HDF5 file."""
    states = []
    with h5py.File(filename, 'r') as f:
        for state_name in f['trajectory']:
            state_grp = f['trajectory'][state_name]
            positions = {}
            for ptype in state_grp['positions']:
                positions[ptype] = state_grp['positions'][ptype][()]
            states.append({
                'step': state_grp.attrs['step'],
                'positions': positions
            })
    return sorted(states, key=lambda x: x['step'])

def check_restraints(coords, thresholds=None):
    """
    Check if restraints are satisfied.
    Returns: dict with counts of satisfied restraints.
    """
    if thresholds is None:
        params = SystemParameters()
        thresholds = params.pair_distances  # {'AA': 48.22, 'AB': 38.5, 'BC': 34.0}
        tolerance = 5.0  # Tolerance for distance checks
        thresholds = {k: v + tolerance for k, v in thresholds.items()}
    
    # Calculate distances
    aa_dists = [np.linalg.norm(coords['A'][i] - coords['A'][j]) 
                for i in range(len(coords['A'])) 
                for j in range(i+1, len(coords['A']))]
    
    ab_dists = [np.linalg.norm(coords['A'][i] - coords['B'][j]) 
                for i in range(len(coords['A'])) 
                for j in range(len(coords['B']))]
    
    bc_dists = [np.linalg.norm(coords['B'][i] - coords['C'][j]) 
                for i in range(len(coords['B'])) 
                for j in range(len(coords['C']))]
    #print("AA distances: ", aa_dists)
    # Count satisfied restraints (distance < threshold)
    return {
        'aa_satisfied': sum(1 for d in aa_dists if d < thresholds['AA']),
        'ab_satisfied': sum(1 for d in ab_dists if d < thresholds['AB']),
        'bc_satisfied': sum(1 for d in bc_dists if d < thresholds['BC']),
        'aa_total': len(aa_dists),
        'ab_total': len(ab_dists),
        'bc_total': len(bc_dists)
    }

def analyze_single_structure(coords):
    """Analyze restraint satisfaction for a single structure."""
    results = check_restraints(coords)
    total_satisfied = results['aa_satisfied'] + results['ab_satisfied'] + results['bc_satisfied']
    total_restraints = results['aa_total'] + results['ab_total'] + results['bc_total']
    
    print(f"Restraint Satisfaction:")
    print(f"  A-A: {results['aa_satisfied']}/{results['aa_total']}")
    print(f"  A-B: {results['ab_satisfied']}/{results['ab_total']}")
    print(f"  B-C: {results['bc_satisfied']}/{results['bc_total']}")
    print(f"  Total: {total_satisfied}/{total_restraints} ({100*total_satisfied/total_restraints:.1f}%)")
    
    return results

def analyze_trajectory(filename, burnin=0.2):
    """Analyze restraint satisfaction over trajectory."""
    states = load_trajectory(filename)
    start_idx = int(burnin * len(states))
    post_burnin = states[start_idx:]
    
    print(f"Analyzing {len(post_burnin)} frames (after {burnin*100}% burnin)")
    
    # Collect data
    data = []
    for i, state in enumerate(post_burnin):
        results = check_restraints(state['positions'])
        total_satisfied = results['aa_satisfied'] + results['ab_satisfied'] + results['bc_satisfied']
        total_restraints = results['aa_total'] + results['ab_total'] + results['bc_total']
        
        data.append({
            'frame': i+1,
            'step': state['step'],
            'total_satisfied': total_satisfied,
            'total_restraints': total_restraints,
            'aa_satisfied': results['aa_satisfied'],
            'ab_satisfied': results['ab_satisfied'],
            'bc_satisfied': results['bc_satisfied']
        })
    
    df = pd.DataFrame(data)
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Satisfied restraints vs frame
    #ax1.plot(df['frame'], df['total_satisfied'], 'o-', label='Total satisfied')
    ax1.plot(df['frame'], df['aa_satisfied'], 'o-', label='A-A satisfied')
    ax1.plot(df['frame'], df['ab_satisfied'], 'o-', label='A-B satisfied')
    ax1.plot(df['frame'], df['bc_satisfied'], 'o-', label='B-C satisfied')
    #ax1.axhline(df['total_restraints'].iloc[0], color='red', linestyle='--', label='Maximum possible')
    ax1.set_xlabel('Frame')
    ax1.set_ylabel('Satisfied Restraints')
    ax1.set_title('Restraint Satisfaction vs Frame')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Distribution
    ax2.hist(df['total_satisfied'], bins=20, alpha=0.7)
    ax2.axvline(df['total_satisfied'].mean(), color='red', linestyle='--', 
                label=f'Mean: {df["total_satisfied"].mean():.1f}')
    ax2.set_xlabel('Total Satisfied Restraints')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('restraint_analysis.png', dpi=150)
    plt.show()
    
    # Save data
    df.to_csv('restraint_data.csv', index=False)
    
    # Summary
    mean_satisfied = df['total_satisfied'].mean()
    total_possible = df['total_restraints'].iloc[0]
    print(f"\nSummary:")
    print(f"  Mean satisfied: {mean_satisfied:.1f}/{total_possible} ({100*mean_satisfied/total_possible:.1f}%)")
    print(f"  Data saved to: restraint_data.csv")
    print(f"  Plot saved to: restraint_analysis.png")
    
    return df

if __name__ == "__main__":
    # Example 1: Analyze ideal coordinates
    params = SystemParameters()
    ideal_coords = params.latest_ideal()
    print("=== Ideal Structure Analysis ===")
    analyze_single_structure(ideal_coords)
    
    # Example 2: Analyze trajectory
    trajectory_file = "../output_analysis/tetramersampler_results_1/trajectory_chain_7.h5"
    if os.path.exists(trajectory_file):
        print("\n=== Trajectory Analysis ===")
        df = analyze_trajectory(trajectory_file, burnin=0.2)
    else:
        print(f"\nTrajectory file not found: {trajectory_file}")