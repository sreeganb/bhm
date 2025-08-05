#!/usr/bin/env python3
"""
Simple scoring script - no fancy imports, just add the path and import directly.
"""

import sys
import os
import numpy as np

# Just add the parent directory to Python's path - that's it!
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

# Now import normally - Python will find the files
import parameters
import pair_sampler
import tetramer_sampler
import octet_sampler
import h5py

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
# create instance of the pair sampler
#-----------------------------------------------------
def pair_sampler_score(coords, sigmas, sigma_ranges):
    sampler_sequence = ["pair_sampler"]
    sampler_idx = 0
    sampler = PairSampler(sampler_sequence, sampler_idx)
    tot_score, ex_score, pair_score, prior_score = sampler.calculate_score(pos = coords, sig = sigmas, sig_range=sigma_ranges)
    return tot_score, ex_score, pair_score, prior_score

def tetramer_sampler_score(coords, sigmas, sigma_ranges):
    sampler_sequence = ["tetramer_sampler"]
    sampler_idx = 0
    sampler = TetramerSampler(sampler_sequence, sampler_idx, positions_ts=coords)
    tetramers = sampler.get_tetramers(positions =coords)
    tot_score, ex_score, pair_score, tet_score = sampler.neg_log_posterior(positions = coords,
                                                                             tetramers = tetramers, 
                                                                           sig = sigmas)
    return tot_score, ex_score, pair_score, tet_score

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


def analyze_trajectory_scores(filename: str, burnin_fraction=0.2):
    """
    Load trajectory data and compute scores for each state after burnin period.
    
    Args:
        filename: Path to HDF5 trajectory file
        burnin_fraction: Fraction of trajectory to discard as burnin (default 0.2 = 20%)
    """
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
    
    # Storage for results
    results = {
        'step': [],
        'total_score': [],
        'pair_total': [],
        'exclusion_score': [],
        'pairwise_score': [],
        'prior_score': [],
        'tetramer_score': [],
        'sigma_AA': [],
        'sigma_AB': [],
        'sigma_BC': []
    }
    
    # Define sigma ranges (you can adjust these)
    sigma_ranges = {
        'AA': (0.1, 2.2),
        'AB': (0.1, 2.3),
        'BC': (0.2, 2.4),
    }
    
    # Process each post-burnin state
    for i, state in enumerate(post_burnin_states):
        step = state['step']
        coords = state['positions']
        sigmas = state['sigma']
        
        print(f"Processing step {step} (post-burnin {i+1}/{len(post_burnin_states)})...")
        
        # Check if we have required data
        if not coords or not sigmas:
            print(f"  Skipping step {step} - missing coordinates or sigmas")
            continue
            
        if not all(ptype in coords for ptype in ['A', 'B', 'C']):
            print(f"  Skipping step {step} - missing particle types")
            continue
            
        if not all(stype in sigmas for stype in ['AA', 'AB', 'BC']):
            print(f"  Skipping step {step} - missing sigma values")
            continue
        
        try:
            # Compute pair scores
            tot_pair_score, ex_score, pair_score, prior_score = pair_sampler_score(
                coords, sigmas, sigma_ranges
            )
            
#            # Compute tetramer score
#            _, _, _, tet_score = tetramer_sampler_score(
#                coords, sigmas, sigma_ranges
#            )
            
#            # Total score
#            total_score = tot_pair_score + tet_score
            
            # Store results
            results['step'].append(step)
#            results['total_score'].append(total_score)
            results['pair_total'].append(tot_pair_score)
            results['exclusion_score'].append(ex_score)
            results['pairwise_score'].append(pair_score)
            results['prior_score'].append(prior_score)
            results['tetramer_score'].append(tet_score)
            results['sigma_AA'].append(sigmas['AA'])
            results['sigma_AB'].append(sigmas['AB'])
            results['sigma_BC'].append(sigmas['BC'])
            
#            print(f"  Step {step}: Total={total_score:.3f}, Pair={tot_pair_score:.3f}, Tet={tet_score:.3f}")
            
        except Exception as e:
            print(f"  Error processing step {step}: {e}")
            continue
    
    print(f"\nSuccessfully processed {len(results['step'])} post-burnin states")
    
    # Convert to DataFrame for easy analysis
    import pandas as pd
    df = pd.DataFrame(results)
    
    # Save results with burnin info in filename
    output_file = f"trajectory_scores_burnin{int(burnin_fraction*100)}.csv"
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
    # Print summary statistics
    print(f"\n=== Post-Burnin Summary Statistics (burnin={burnin_fraction}) ===")
    print(f"Total Score: mean={df['total_score'].mean():.3f}, std={df['total_score'].std():.3f}")
    print(f"Pair Score: mean={df['pair_total'].mean():.3f}, std={df['pair_total'].std():.3f}")
    print(f"Tetramer Score: mean={df['tetramer_score'].mean():.3f}, std={df['tetramer_score'].std():.3f}")
    print(f"Sigma AA: mean={df['sigma_AA'].mean():.3f}, std={df['sigma_AA'].std():.3f}")
    print(f"Sigma AB: mean={df['sigma_AB'].mean():.3f}, std={df['sigma_AB'].std():.3f}")
    print(f"Sigma BC: mean={df['sigma_BC'].mean():.3f}, std={df['sigma_BC'].std():.3f}")
    
    return df

# Usage example with different burnin fractions
if __name__ == "__main__":
    # Specify your trajectory file
    trajectory_file = "../output_analysis/pairsampler_results_1/trajectory_chain_2.h5"  # Update this path
    
    # Try different burnin fractions
    burnin_options = [0.1, 0.2, 0.3]  # 10%, 20%, 30%
    
    for burnin in burnin_options:
        print(f"\n{'='*60}")
        print(f"ANALYZING WITH BURNIN FRACTION: {burnin}")
        print(f"{'='*60}")
        
        # Analyze the trajectory with specific burnin
        results_df = analyze_trajectory_scores(trajectory_file, burnin_fraction=burnin)
        
        if results_df is not None:
            # Create plots with burnin info
            import matplotlib.pyplot as plt
            
            plt.figure(figsize=(15, 10))
            
            # Plot 1: Total score vs step
            plt.subplot(2, 3, 1)
            plt.plot(results_df['step'], results_df['total_score'])
            plt.xlabel('Step')
            plt.ylabel('Total Score')
            plt.title(f'Total Score vs Step (burnin={burnin})')
            
            # Plot 2: Score components
            plt.subplot(2, 3, 2)
            plt.plot(results_df['step'], results_df['pair_total'], label='Pair')
            plt.plot(results_df['step'], results_df['tetramer_score'], label='Tetramer')
            plt.xlabel('Step')
            plt.ylabel('Score')
            plt.title('Score Components')
            plt.legend()
            
            # Plot 3: Sigma values
            plt.subplot(2, 3, 3)
            plt.plot(results_df['step'], results_df['sigma_AA'], label='AA')
            plt.plot(results_df['step'], results_df['sigma_AB'], label='AB')
            plt.plot(results_df['step'], results_df['sigma_BC'], label='BC')
            plt.xlabel('Step')
            plt.ylabel('Sigma')
            plt.title('Sigma Values')
            plt.legend()
            
            # Plot 4: Score distribution
            plt.subplot(2, 3, 4)
            plt.hist(results_df['total_score'], bins=20, alpha=0.7)
            plt.xlabel('Total Score')
            plt.ylabel('Frequency')
            plt.title('Score Distribution')
            
            # Plot 5: Running average (to check convergence)
            plt.subplot(2, 3, 5)
            running_avg = results_df['total_score'].rolling(window=max(1, len(results_df)//10)).mean()
            plt.plot(results_df['step'], running_avg)
            plt.xlabel('Step')
            plt.ylabel('Running Average Score')
            plt.title('Convergence Check')
            
            # Plot 6: Autocorrelation of total score
            plt.subplot(2, 3, 6)
            from scipy.stats import pearsonr
            lags = range(1, min(50, len(results_df)//4))
            autocorr = [pearsonr(results_df['total_score'][:-lag], 
                               results_df['total_score'][lag:])[0] for lag in lags]
            plt.plot(lags, autocorr)
            plt.xlabel('Lag')
            plt.ylabel('Autocorrelation')
            plt.title('Score Autocorrelation')
            plt.axhline(y=0, color='r', linestyle='--', alpha=0.5)
            
            plt.tight_layout()
            plot_filename = f'trajectory_analysis_burnin{int(burnin*100)}.png'
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            print(f"Plots saved as {plot_filename}")
            plt.show()

    # Create a dictionary of sigma values AA, AB, BC
    sigmas = {
        'AA': 0.9,
        'AB': 0.7,
        'BC': 0.1,
    }
    # this is how I want the ranges to be defined sigma_ranges: [str, Tuple[float, float]]]
    sigma_ranges = {
        'AA': (0.1, 2.2),
        'AB': (0.1, 2.3),
        'BC': (0.2, 2.4),
    }
    tot_score, ex_score, pair_score, prior_score = pair_sampler_score(ideal_coords, sigmas, sigma_ranges)

    print(f"Total Score: {tot_score}, Excluded Volume Score: {ex_score}, Pair Score: {pair_score}, Prior Score: {prior_score}")

    tot_score, ex_score, pair_score, tet_score = tetramer_sampler_score(ideal_coords, sigmas, sigma_ranges)
    print(f"Total Score: {tot_score}, Excluded Volume Score: {ex_score}, Pair Score: {pair_score}, Tetramer Score: {tet_score}")
