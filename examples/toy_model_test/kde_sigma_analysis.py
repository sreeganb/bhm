#============================================================
# Analysis script to find well-converged trajectories and plot
# KDE distributions of sigma values for each sampler type
#============================================================
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import os
import pandas as pd
import h5py
from scipy.stats import kstest
from scipy import signal

#--------------------------------------------------------------
# Function to read the sampler sequence from a JSON file
#--------------------------------------------------------------
def read_sampler_sequence(file_path):
    import json
    with open(file_path, 'r') as f:
        data = json.load(f)
    print(f"Sampler sequence: {data['sequence']}")
    print(f"Number of chains: {data['n_chains']}")
    return data['sequence'], data['n_chains']

#--------------------------------------------------------------
# Extract sigma time series from HDF5 file
#--------------------------------------------------------------
def load_sigma_timeseries(filename: str) -> pd.DataFrame:
    """
    Loads only the sigma parameters (AA, AB, BC) from an HDF5 trajectory file.
    Returns a pandas DataFrame with columns for step and sigma values.
    """
    steps = []
    sigma_AA = []
    sigma_AB = []
    sigma_BC = []
    
    if not os.path.isfile(filename):
        print(f"File not found: {filename}")
        return pd.DataFrame()

    with h5py.File(filename, 'r') as f:
        if 'trajectory' not in f:
            print(f"No 'trajectory' group in file: {filename}")
            return pd.DataFrame()
        
        traj_grp = f['trajectory']
        for state_name in traj_grp:
            state_grp = traj_grp[state_name]
            step = state_grp.attrs.get("step", 0)
            
            if 'sigma' in state_grp:
                sigma_grp = state_grp['sigma']
                steps.append(step)
                sigma_AA.append(sigma_grp.attrs.get('AA', float('nan')))
                sigma_AB.append(sigma_grp.attrs.get('AB', float('nan')))
                sigma_BC.append(sigma_grp.attrs.get('BC', float('nan')))
    
    df = pd.DataFrame({
        'step': steps,
        'sigma_AA': sigma_AA,
        'sigma_AB': sigma_AB,
        'sigma_BC': sigma_BC
    })
    
    df = df.sort_values('step').reset_index(drop=True)
    return df

#--------------------------------------------------------------
# Diagnostic function to understand the data
#--------------------------------------------------------------
def diagnose_trajectory_data(df: pd.DataFrame, chain_id: str, sampler_label: str):
    """
    Diagnose what's wrong with the trajectory data.
    """
    print(f"\n--- Diagnosing {sampler_label} Chain {chain_id} ---")
    print(f"Total samples: {len(df)}")
    print(f"Steps range: {df['step'].min()} to {df['step'].max()}")
    
    for pair_type in ['AA', 'AB', 'BC']:
        sigma_col = f'sigma_{pair_type}'
        if sigma_col in df.columns:
            data = df[sigma_col].dropna()
            if len(data) > 0:
                print(f"{pair_type}: {len(data)} samples, mean={data.mean():.3f}, std={data.std():.3f}, range=({data.min():.3f}, {data.max():.3f})")
                
                # Check for constant values
                if data.std() < 1e-6:
                    print(f"  WARNING: {pair_type} appears constant!")
                
                # Check for extreme values
                if data.max() > 100 or data.min() < 0:
                    print(f"  WARNING: {pair_type} has extreme values!")
                    
            else:
                print(f"{pair_type}: No valid data")
        else:
            print(f"{pair_type}: Column missing")

#--------------------------------------------------------------
# Relaxed convergence assessment
#--------------------------------------------------------------
def assess_convergence_relaxed(df: pd.DataFrame, pair_type: str, min_samples: int = 50) -> dict:
    """
    More relaxed convergence assessment for challenging data.
    """
    sigma_col = f'sigma_{pair_type}'
    if sigma_col not in df.columns or df[sigma_col].isna().all():
        return {'converged': False, 'reason': 'No data', 'score': 0}
    
    values = df[sigma_col].dropna().values
    if len(values) < min_samples:
        return {'converged': False, 'reason': f'Insufficient samples ({len(values)} < {min_samples})', 'score': 0}
    
    # Check for constant values
    if np.std(values) < 1e-6:
        return {'converged': False, 'reason': 'Constant values', 'score': 0}
    
    # Apply burn-in (remove first 30% for relaxed criteria)
    burn_in = max(10, int(len(values) * 0.3))  # At least 10 samples burn-in
    post_burnin = values[burn_in:]
    
    if len(post_burnin) < 20:
        return {'converged': False, 'reason': 'Insufficient post-burnin samples', 'score': 0}
    
    metrics = {}
    score_components = []
    
    # 1. Basic stationarity (less strict)
    if len(post_burnin) >= 20:
        mid_point = len(post_burnin) // 2
        first_half = post_burnin[:mid_point]
        second_half = post_burnin[mid_point:]
        
        try:
            ks_stat, ks_pvalue = kstest(first_half, second_half)
            metrics['ks_pvalue'] = ks_pvalue
            score_components.append(min(ks_pvalue * 2, 1.0))  # Boost low p-values
        except:
            metrics['ks_pvalue'] = 0.01
            score_components.append(0.1)
    
    # 2. Simple mixing check (standard deviation should be reasonable)
    std_ratio = np.std(post_burnin) / (np.mean(np.abs(post_burnin)) + 1e-6)
    metrics['std_ratio'] = std_ratio
    score_components.append(min(std_ratio * 2, 1.0))
    
    # 3. Range check (values should vary)
    value_range = np.max(post_burnin) - np.min(post_burnin)
    relative_range = value_range / (np.mean(np.abs(post_burnin)) + 1e-6)
    metrics['relative_range'] = relative_range
    score_components.append(min(relative_range, 1.0))
    
    # 4. Simple trend check
    if len(post_burnin) > 10:
        x = np.arange(len(post_burnin))
        try:
            slope, _, _, p_trend, _ = stats.linregress(x, post_burnin)
            metrics['trend_pvalue'] = p_trend
            metrics['abs_slope'] = abs(slope)
            score_components.append(min(p_trend * 2, 1.0))
        except:
            metrics['trend_pvalue'] = 0.5
            score_components.append(0.5)
    
    # Relaxed convergence criteria
    converged = (
        len(post_burnin) >= 20 and
        metrics.get('ks_pvalue', 0) > 0.01 and  # Much more relaxed
        std_ratio > 0.01 and  # Must have some variation
        relative_range > 0.01  # Must have some range
    )
    
    metrics['converged'] = converged
    metrics['score'] = np.mean(score_components) if score_components else 0
    metrics['n_samples'] = len(values)
    metrics['n_post_burnin'] = len(post_burnin)
    
    return metrics

#--------------------------------------------------------------
# Find best trajectory with relaxed criteria
#--------------------------------------------------------------
def find_best_trajectory_relaxed(folder_name: str, n_chains: int) -> tuple:
    """
    Find the best trajectory with relaxed convergence criteria.
    """
    best_chain = None
    best_df = None
    best_score = -1
    all_metrics = {}
    
    print(f"\nAnalyzing {n_chains} chains in {folder_name}")
    
    for chain in range(n_chains):
        filename = f"trajectory_chain_{chain+1}.h5"
        df = load_sigma_timeseries(filename)
        
        if df.empty:
            print(f"Chain {chain+1}: No data")
            continue
        
        # Diagnose the data
        diagnose_trajectory_data(df, str(chain+1), folder_name)
        
        # Assess convergence for all pair types
        chain_metrics = {}
        total_score = 0
        valid_pairs = 0
        
        for pair_type in ['AA', 'AB', 'BC']:
            metrics = assess_convergence_relaxed(df, pair_type)
            chain_metrics[pair_type] = metrics
            if metrics['score'] > 0:
                valid_pairs += 1
                total_score += metrics['score']
        
        # Overall score: average of valid pairs
        if valid_pairs > 0:
            overall_score = total_score / valid_pairs
            chain_metrics['overall_score'] = overall_score
            chain_metrics['valid_pairs'] = valid_pairs
            
            print(f"Chain {chain+1}: Score = {overall_score:.3f}, Valid pairs = {valid_pairs}/3")
            
            if overall_score > best_score:
                best_score = overall_score
                best_chain = chain + 1
                best_df = df.copy()
        else:
            print(f"Chain {chain+1}: No valid pairs")
        
        all_metrics[f'chain_{chain+1}'] = chain_metrics
    
    return best_chain, best_df, all_metrics

#--------------------------------------------------------------
# Plot trace plots for diagnosis
#--------------------------------------------------------------
def plot_trace_diagnostics(selected_trajectories):
    """
    Plot trace plots to visually assess the selected trajectories.
    """
    n_samplers = len(selected_trajectories)
    if n_samplers == 0:
        return
    
    fig, axes = plt.subplots(3, n_samplers, figsize=(5*n_samplers, 12))
    if n_samplers == 1:
        axes = axes.reshape(-1, 1)
    
    for j, (sampler_label, traj_info) in enumerate(selected_trajectories.items()):
        df = traj_info['df']
        chain_id = traj_info['chain']
        
        # Apply burn-in for plotting
        burn_in = max(10, int(len(df) * 0.3))
        df_plot = df.iloc[burn_in:]
        
        for i, pair_type in enumerate(['AA', 'AB', 'BC']):
            ax = axes[i, j]
            sigma_col = f'sigma_{pair_type}'
            
            if sigma_col in df_plot.columns:
                data = df_plot[sigma_col].dropna()
                if len(data) > 0:
                    ax.plot(data.index, data.values, alpha=0.7)
                    ax.set_title(f'{sampler_label} Chain {chain_id}\nSigma {pair_type}')
                    ax.set_ylabel(f'Sigma {pair_type}')
                    if i == 2:  # Last row
                        ax.set_xlabel('Sample Index')
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
            else:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    plt.savefig('sigma_plots/trace_diagnostics.pdf', dpi=300, bbox_inches='tight')
    plt.show()

#--------------------------------------------------------------
# Modified KDE fitting with better error handling
#--------------------------------------------------------------
def fit_kde_to_sigma_robust(df: pd.DataFrame, bandwidth: float = 0.2):
    """
    Robust KDE fitting with better error handling.
    """
    kdes = {}
    burn_in_fraction = 0.3  # Reduced burn-in
    burn_in_index = max(10, int(len(df) * burn_in_fraction))
    df_post_burnin = df.iloc[burn_in_index:]
    
    for pair_type in ['AA', 'AB', 'BC']:
        sigma_col = f'sigma_{pair_type}'
        if sigma_col not in df_post_burnin.columns:
            continue
        
        data = df_post_burnin[sigma_col].dropna()
        if len(data) > 5 and np.std(data) > 1e-6:  # More lenient requirements
            try:
                # Use smaller bandwidth for potentially noisy data
                kde = stats.gaussian_kde(data, bw_method=bandwidth)
                kdes[pair_type] = kde
                print(f"  KDE fitted for {pair_type}: {len(data)} samples")
            except Exception as e:
                print(f"  Failed to fit KDE for {pair_type}: {e}")
    
    return kdes

#--------------------------------------------------------------
# Main analysis function with better fallbacks
#--------------------------------------------------------------
def analyze_converged_trajectories(sequence, n_chains):
    """
    Find best trajectory for each sampler with relaxed criteria.
    """
    # Create output directory
    if not os.path.exists("sigma_plots"):
        os.makedirs("sigma_plots")
    
    # Track sampler instances
    sampler_counts = {'pair': 0, 'tetramer': 0, 'octet': 0}
    selected_trajectories = {}
    all_convergence_info = {}
    
    # Process each sampler in sequence
    for i, sampler in enumerate(sequence):
        sampler_counts[sampler] += 1
        sampler_id = sampler_counts[sampler]
        sampler_label = f"{sampler}_{sampler_id}"
        
        folder_name = f"./output_analysis/{sampler}sampler_results_{sampler_id}/"
        print(f"\n{'='*50}")
        print(f"=== Analyzing {sampler_label} ===")
        print(f"Folder: {folder_name}")
        
        if not os.path.exists(folder_name):
            print(f"Folder {folder_name} does not exist. Skipping.")
            continue
        
        # Change to sampler directory
        parent_dir = os.getcwd()
        os.chdir(folder_name)
        
        # Find best trajectory with relaxed criteria
        best_chain, best_df, metrics = find_best_trajectory_relaxed(".", n_chains)
        
        # Return to parent directory
        os.chdir(parent_dir)
        
        if best_chain is not None:
            print(f"\n✓ Selected chain {best_chain} for {sampler_label}")
            conv_info = metrics[f'chain_{best_chain}']
            print(f"Overall score: {conv_info['overall_score']:.3f}")
            print(f"Valid pairs: {conv_info['valid_pairs']}/3")
            
            for pair_type in ['AA', 'AB', 'BC']:
                pair_metrics = conv_info[pair_type]
                status = "✓" if pair_metrics['converged'] else "~"
                print(f"  {pair_type}: {status} (score: {pair_metrics['score']:.3f})")
            
            selected_trajectories[sampler_label] = {
                'df': best_df,
                'chain': best_chain,
                'folder': folder_name
            }
            all_convergence_info[sampler_label] = metrics
        else:
            print(f"\n✗ No usable trajectories found for {sampler_label}")
    
    # Generate outputs if we have any trajectories
    if selected_trajectories:
        print(f"\n{'='*50}")
        print(f"Found {len(selected_trajectories)} usable trajectories")
        
        # Plot trace diagnostics
        plot_trace_diagnostics(selected_trajectories)
        
        # Plot KDEs for selected trajectories
        plot_selected_kdes_robust(selected_trajectories, all_convergence_info)
        save_convergence_report(all_convergence_info)
    else:
        print("\n" + "="*50)
        print("No usable trajectories found!")
        print("This might indicate:")
        print("1. MCMC chains are not mixing well")
        print("2. Insufficient sampling")
        print("3. Model/data issues")
        print("4. Try running longer chains or different initialization")

#--------------------------------------------------------------
# Robust plotting function
#--------------------------------------------------------------
def plot_selected_kdes_robust(selected_trajectories, convergence_info):
    """
    Plot KDEs for selected trajectories with robust error handling.
    """
    # Fit KDEs for all selected trajectories
    all_kdes = {'AA': {}, 'AB': {}, 'BC': {}}
    
    print("\nFitting KDEs...")
    for sampler_label, traj_info in selected_trajectories.items():
        print(f"Processing {sampler_label}...")
        kdes = fit_kde_to_sigma_robust(traj_info['df'])
        for pair_type, kde in kdes.items():
            all_kdes[pair_type][sampler_label] = kde
    
    # Check if we have any KDEs to plot
    total_kdes = sum(len(pair_kdes) for pair_kdes in all_kdes.values())
    if total_kdes == 0:
        print("No KDEs could be fitted. Check your data quality.")
        return
    
    # Create combined plot
    plt.figure(figsize=(15, 12))
    colors = plt.cm.Set1(np.linspace(0, 1, max(len(selected_trajectories), 3)))
    
    for i, pair_type in enumerate(['AA', 'AB', 'BC']):
        plt.subplot(3, 1, i+1)
        
        if len(all_kdes[pair_type]) > 0:
            x = np.linspace(0, 20, 1000)
            for j, (sampler_label, kde) in enumerate(all_kdes[pair_type].items()):
                color = colors[j % len(colors)]
                chain_id = selected_trajectories[sampler_label]['chain']
                label = f"{sampler_label} (chain {chain_id})"
                
                plt.plot(x, kde(x), label=label, color=color, linewidth=2)
                plt.fill_between(x, kde(x), alpha=0.2, color=color)
            
            plt.title(f'Sigma {pair_type} KDE - Available Trajectories', fontsize=14)
            plt.xlabel('Sigma Value')
            plt.ylabel('Density')
            plt.legend()
            plt.grid(True, alpha=0.3)
        else:
            plt.text(0.5, 0.5, f'No KDE available for Sigma {pair_type}', 
                    ha='center', va='center', transform=plt.gca().transAxes)
            plt.title(f'Sigma {pair_type} - No Data Available', fontsize=14)
    
    plt.tight_layout()
    plt.savefig('sigma_plots/available_trajectories_kde.pdf', dpi=300, bbox_inches='tight')
    print(f"Saved KDE plot to sigma_plots/available_trajectories_kde.pdf")
    plt.show()

def save_convergence_report(convergence_info):
    """
    Save detailed convergence report.
    """
    with open('sigma_plots/convergence_report.txt', 'w') as f:
        f.write("=== MCMC Convergence Analysis Report (Relaxed Criteria) ===\n\n")
        
        for sampler_label, metrics in convergence_info.items():
            f.write(f"\n{sampler_label}:\n")
            f.write("-" * 40 + "\n")
            
            for chain_id, chain_metrics in metrics.items():
                if chain_id.startswith('chain_'):
                    f.write(f"\n{chain_id}:\n")
                    f.write(f"  Overall Score: {chain_metrics.get('overall_score', 'N/A'):.4f}\n")
                    f.write(f"  Valid Pairs: {chain_metrics.get('valid_pairs', 0)}/3\n")
                    
                    for pair_type in ['AA', 'AB', 'BC']:
                        pair_metrics = chain_metrics[pair_type]
                        f.write(f"  {pair_type}: {'✓' if pair_metrics['converged'] else '~'}\n")
                        f.write(f"    Score: {pair_metrics['score']:.4f}\n")
                        f.write(f"    Samples: {pair_metrics.get('n_samples', 'N/A')}\n")
                        if 'ks_pvalue' in pair_metrics:
                            f.write(f"    KS p-value: {pair_metrics['ks_pvalue']:.4f}\n")
    
    print("Saved detailed convergence report to sigma_plots/convergence_report.txt")

if __name__ == "__main__":
    sequence, n_chains = read_sampler_sequence('sampler_sequence.json')
    analyze_converged_trajectories(sequence, n_chains)