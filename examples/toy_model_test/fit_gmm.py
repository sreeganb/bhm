from email import parser
import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # Fix Qt platform plugin warning
import numpy as np
import matplotlib.pyplot as plt
import h5py
import arviz as az
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from sklearn.mixture import GaussianMixture
import json
import argparse
import sys
from diagnostics import AdvancedConvergenceDiagnostics

# Replace your load_trajectory_from_hdf5 function with this corrected version:

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
                "step": int(state_grp.attrs.get("step", 0)),
                "total_score": float(state_grp.attrs.get("total_score", 0.0)),
                "prior_score": float(state_grp.attrs.get("prior_score", 0.0)),
                "pair_score": float(state_grp.attrs.get("pair_score", 0.0)),
                "exvol_score": float(state_grp.attrs.get("exvol_score", 0.0)),
                "tet_score": float(state_grp.attrs.get("tet_score", 0.0)),
                "oct_score": float(state_grp.attrs.get("oct_score", 0.0)),
                "sigma": {},
                "positions": {},
                "types": {},
                "bead_numbers": {}
            }
            
            # FIX: Safely read sigma values - handle HDF5 attributes properly
            if 'sigma' in state_grp:
                sigma_grp = state_grp['sigma']
                for key in sigma_grp.attrs:
                    # Get the raw value and ensure it's a Python float
                    raw_value = sigma_grp.attrs[key]
                    if hasattr(raw_value, 'item'):  # numpy scalar
                        state["sigma"][key] = float(raw_value.item())
                    elif isinstance(raw_value, (int, float, np.integer, np.floating)):
                        state["sigma"][key] = float(raw_value)
                    else:
                        print(f"Warning: unexpected sigma value type {type(raw_value)} for {key}")
                        state["sigma"][key] = 0.0

            # Safely read positions
            if 'positions' in state_grp:
                pos_grp = state_grp['positions']
                for type_name in pos_grp:
                    # Read the dataset and convert to numpy array
                    dataset = pos_grp[type_name]
                    state["positions"][type_name] = np.array(dataset[:])

            # Safely read types
            if 'types_keys' in state_grp and 'types_vals' in state_grp:
                try:
                    types_keys_data = state_grp['types_keys'][:]
                    types_vals_data = state_grp['types_vals'][:]
                    if len(types_keys_data) == len(types_vals_data):
                        for i in range(len(types_keys_data)):
                            k = types_keys_data[i]
                            v = types_vals_data[i]
                            if isinstance(k, bytes):
                                k = k.decode('utf-8', errors='ignore')
                            if isinstance(v, bytes):
                                v = v.decode('utf-8', errors='ignore')
                            state["types"][str(k)] = str(v)
                except Exception as e:
                    print(f"Warning: failed to read types: {e}")
            
            # Safely read bead_numbers
            if 'bead_keys' in state_grp and 'bead_vals' in state_grp:
                try:
                    bead_keys_data = state_grp['bead_keys'][:]
                    bead_vals_data = state_grp['bead_vals'][:]
                    if len(bead_keys_data) == len(bead_vals_data):
                        for i in range(len(bead_keys_data)):
                            bkey = bead_keys_data[i]
                            bval = bead_vals_data[i]
                            # Convert to Python int
                            if hasattr(bkey, 'item'):
                                bkey = int(bkey.item())
                            else:
                                bkey = int(bkey)
                            if hasattr(bval, 'item'):
                                bval = int(bval.item())
                            else:
                                bval = int(bval)
                            state["bead_numbers"][bkey] = bval
                except Exception as e:
                    print(f"Warning: failed to read bead_numbers: {e}")

            all_states.append(state)

    # Sort the collected states by their step number
    all_states.sort(key=lambda s: s["step"])
    print(f"Loaded {len(all_states)} states from {filename}")

    return all_states

# Add this import at the top with your other imports
from statsmodels.tsa.stattools import acf

# Add this function after load_trajectory_from_hdf5 function
def calculate_autocorrelation(data, nlags=50):
    """Calculate autocorrelation function for time series data."""
    try:
        # Use statsmodels acf function to calculate autocorrelation
        acf_values = acf(data, nlags=min(nlags, len(data)//3), fft=True)
        return acf_values
    except Exception as e:
        print(f"Error calculating autocorrelation: {e}")
        return None

def calculate_effective_sample_size(data):
    """Estimate effective sample size using autocorrelation."""
    try:
        # Remove any inf/nan values
        data = np.array(data)
        data = data[np.isfinite(data)]
        
        if len(data) < 10:
            return None
            
        acf_values = acf(data, nlags=min(100, len(data)//3), fft=True)
        # Find first negative autocorrelation or stop at acf_values < 0.05
        tau = 1  # Default correlation time
        for i, acf_val in enumerate(acf_values[1:], 1):  # Skip lag 0
            if acf_val < 0.05:
                tau = i
                break
        # ESS = N / (1 + 2 * sum(ACF))
        ess = len(data) / (1 + 2 * np.sum(acf_values[1:tau]))
        return max(1, ess)  # Ensure ESS is at least 1
    except Exception as e:
        print(f"Error calculating effective sample size: {e}")
        return None

def fit_gmm_simple(data, sigma_type: str, max_components: int = 4, min_components: int = 1,
                   criterion: str = "bic", random_state: int = 42):
    """
    Fit a 1D Gaussian Mixture Model to data using BIC (default) or AIC to 
    pick component count.

    Args:
        data: 1D iterable of samples
        sigma_type: label (used only for logging / filenames)
        max_components: upper bound on number of mixture components
        min_components: lower bound (>=1)
        criterion: 'bic' or 'aic'
        random_state: RNG seed for reproducibility

    Returns:
        best_gmm (sklearn.mixture.GaussianMixture) or None if not enough data
    """
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]             # remove NaN/inf
    if data.ndim != 1:
        data = data.ravel()

    if len(data) < 5:
        print(f"[GMM] {sigma_type}: insufficient data ({len(data)} samples)")
        return None

    X = data.reshape(-1, 1)

    best_score = np.inf
    best_gmm = None
    chosen_k = None

    max_k = min(max_components, len(data))  # cannot exceed number of samples
    for k in range(min_components, max_k + 1):
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type='diag',
                random_state=random_state,
                max_iter=200,
                reg_covar=1e-8
            ).fit(X)
            if criterion == "aic":
                score = gmm.aic(X)
            else:
                score = gmm.bic(X)  # default BIC
            if score < best_score:
                best_score = score
                best_gmm = gmm
                chosen_k = k
        except Exception as e:
            print(f"[GMM] {sigma_type}: failed for k={k}: {e}")

    if best_gmm is None:
        print(f"[GMM] {sigma_type}: no valid model fit")
        return None

    print(f"[GMM] {sigma_type}: selected k={chosen_k} ({criterion.upper()}={best_score:.2f})")
    return best_gmm


def save_gmm_simple(gmm, sigma_type: str, chain_id: str, output_dir: str, data):
    """
    Save minimal GMM parameters + basic data stats to JSON.

    Args:
        gmm: fitted GaussianMixture
        sigma_type: label
        chain_id: identifier for chain / level
        output_dir: directory to write JSON
        data: original 1D data array
    """
    if gmm is None:
        print(f"[GMM] {sigma_type}_{chain_id}: nothing to save")
        return False

    data = np.asarray(data)
    params = {
        "sigma_type": sigma_type,
        "chain_id": chain_id,
        "n_components": int(gmm.n_components),
        "means": gmm.means_.ravel().tolist(),
        "variances": gmm.covariances_.ravel().tolist(),  # diag variances
        "weights": gmm.weights_.ravel().tolist(),
        "bic": float(gmm.bic(data.reshape(-1, 1))),
        "aic": float(gmm.aic(data.reshape(-1, 1))),
        "n_iter": int(gmm.n_iter_),
        "converged": bool(getattr(gmm, "converged_", True)),
        "data_n": int(len(data)),
        "data_mean": float(np.mean(data)),
        "data_std": float(np.std(data)),
        "data_min": float(np.min(data)),
        "data_max": float(np.max(data))
    }

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"gmm_fit_{sigma_type}_{chain_id}.json")
    try:
        with open(out_path, "w") as f:
            json.dump(params, f, indent=2)
        return True
    except Exception as e:
        print(f"[GMM] {sigma_type}_{chain_id}: save failed: {e}")
        return False


# Backward-compatible aliases (so existing calls still work)
fit_gmm_robust = fit_gmm_simple
fit_gmm = fit_gmm_simple
save_gmm_with_validation = save_gmm_simple

def get_sampler_folder_name(sampler_type: str, sampler_position: int) -> str:
    """
    Generate the folder name for a specific sampler at a given position in the sequence.
    
    Args:
        sampler_type (str): Type of sampler ('pair', 'tetramer', 'octet', 'full')
        sampler_position (int): Position in the sequence (1-based)
    
    Returns:
        str: Folder name (e.g., 'pairsampler_results_1')
    """
    return f"{sampler_type}sampler_results_{sampler_position}"

def count_sampler_occurrences(sampler_sequence: list, target_sampler: str, position: int) -> int:
    """
    Count how many times a sampler appears up to a specific position in the sequence.
    
    Args:
        sampler_sequence (list): List of sampler names
        target_sampler (str): The sampler we're counting
        position (int): Position in sequence (0-based)
    
    Returns:
        int: The occurrence number (1-based)
    """
    count = 0
    for i in range(position + 1):
        if sampler_sequence[i] == target_sampler:
            count += 1
    return count

def plot_combined_gmm(all_data: dict, all_gmms: dict, sigma_type: str, sampler_name: str, output_dir: str, pdf: PdfPages):
    """Plots combined GMM fits for all chains of a single sigma component."""

    plt.figure(figsize=(12, 8))
    sns.set(style="darkgrid")
    palette = sns.color_palette("husl", len(all_data))  # Colors for each chain

    for i, (chain_id, data) in enumerate(all_data.items()):
        gmm = all_gmms[chain_id]
        x = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)
        logprob = gmm.score_samples(x)
        pdf_ind = np.exp(logprob)

        sns.histplot(data, bins=30, kde=False, stat='density', label=f'Chain {chain_id}', color=palette[i], alpha=0.4)
        plt.plot(x, pdf_ind, '-', label=f'GMM Fit - Chain {chain_id}', color=palette[i])

    plt.title(f'Combined GMM Fit for {sigma_type} ({sampler_name})', fontsize=18)
    plt.xlabel('Sigma Value', fontsize=16)
    plt.ylabel('Density', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.tight_layout()
    pdf.savefig()  # Save to PDF
    plt.savefig(os.path.join(output_dir, f"gmm_combined_plot_{sigma_type}_{sampler_name}.png"))
    plt.close()

def analyze_mcmc_data(output_folder: str, sampler_type: str, sampler_position: int, 
                     burnin=0.4, do_trace_plots=True, do_gmm_fits=True, 
                     do_diagnostics=True):
    """Analyzes MCMC data with improved burnin handling."""
    
    # Setup paths and validation
    sampler_folder_name = get_sampler_folder_name(sampler_type, sampler_position)
    sampler_output_dir = os.path.join(output_folder, sampler_folder_name)
    
    if not os.path.exists(sampler_output_dir):
        print(f"Directory not found: {sampler_output_dir}")
        return
    
    trajectory_files = [f for f in os.listdir(sampler_output_dir) 
                        if f.startswith("trajectory_") and f.endswith(".h5")]
    if not trajectory_files:
        print(f"No trajectory files found in: {sampler_output_dir}")
        return

    sampler_display_name = f"{sampler_type.capitalize()}Sampler_{sampler_position}"
    analysis_output_dir = sampler_output_dir
    
    # Load all trajectory data first (without burn-in applied)
    all_sigma_histories = defaultdict(lambda: defaultdict(list))
    all_scores = defaultdict(lambda: defaultdict(list)) if do_trace_plots else None
    
    for traj_file in trajectory_files:
        level_name = traj_file.replace("trajectory_", "").replace(".h5", "")
        traj_path = os.path.join(sampler_output_dir, traj_file)
        trajectory_data = load_trajectory_from_hdf5(traj_path)
        
        for state in trajectory_data:
            # Defensive sigma handling
            if 'sigma' in state and isinstance(state['sigma'], dict):
                for sigma_type, sigma_value in state['sigma'].items():
                    try:
                        if hasattr(sigma_value, 'item'):
                            sigma_value = float(sigma_value.item())
                        else:
                            sigma_value = float(sigma_value)
                        all_sigma_histories[level_name][sigma_type].append(sigma_value)
                    except (ValueError, TypeError) as e:
                        print(f"Warning: Could not convert sigma value to float: {e}")
                        continue
            
            # Defensive score handling
            if do_trace_plots:
                try:
                    all_scores['total_score'][level_name].append(float(state.get('total_score', 0)))
                    all_scores['prior_score'][level_name].append(float(state.get('prior_score', 0)))
                    all_scores['pair_score'][level_name].append(float(state.get('pair_score', 0)))
                    all_scores['exvol_score'][level_name].append(float(state.get('exvol_score', 0)))
                    all_scores['tet_score'][level_name].append(float(state.get('tet_score', 0)))      
                    all_scores['oct_score'][level_name].append(float(state.get('oct_score', 0)))
                except (ValueError, TypeError) as e:
                    print(f"Warning: Could not convert score values to float: {e}")
                    continue
    
    # Determine burn-in fraction
    burnin_fraction = burnin
    if isinstance(burnin, str) and burnin.lower() == "auto":
        diag_chains = {}
        for chain_id, sig_dict in all_sigma_histories.items():
            for sigma_type, series in sig_dict.items():
                diag_chains.setdefault(sigma_type, []).append(np.array(series))
        
        if diag_chains:
            diag = AdvancedConvergenceDiagnostics(diag_chains)
            burnin_fraction = diag.auto_burnin_fraction()
            print(f"[Diagnostics] Auto burn-in fraction selected: {burnin_fraction:.3f}")
        else:
            burnin_fraction = 0.4  # fallback
    
    # Apply burn-in to all data
    def apply_burnin(data_dict, fraction):
        result = defaultdict(lambda: defaultdict(list))
        for chain_id, sig_dict in data_dict.items():
            for sigma_type, series in sig_dict.items():
                n = len(series)
                cut = int(fraction * n)
                result[chain_id][sigma_type] = series[cut:]
        return result
    
    all_sigma_histories = apply_burnin(all_sigma_histories, burnin_fraction)
    if do_trace_plots:
        all_scores = apply_burnin(all_scores, burnin_fraction)
    
    # Generate reports
    if do_trace_plots:
        _generate_trace_plots_and_diagnostics(all_sigma_histories, all_scores, 
                                            sampler_display_name, analysis_output_dir)
    
    if do_diagnostics:
        _generate_convergence_diagnostics(all_sigma_histories, sampler_display_name, 
                                        analysis_output_dir, burnin_fraction)
    
    if do_gmm_fits:
        _generate_gmm_fits(all_sigma_histories, sampler_display_name, analysis_output_dir)


def _generate_trace_plots_and_diagnostics(all_sigma_histories, all_scores, sampler_display_name, analysis_output_dir):
    """Generate trace plots, R-hat, and autocorrelation analysis."""
    pdf_filename = os.path.join(analysis_output_dir, f"{sampler_display_name}_report.pdf")
    
    with PdfPages(pdf_filename) as pdf:
        sns.set(style="darkgrid")
        palette = sns.color_palette("husl", 4)
        
        # Sigma trace plots
        for chain_idx, chain_id in enumerate(all_sigma_histories):
            plt.figure(figsize=(12, 8))
            for i, sigma_type in enumerate(["AA", "AB", "BC"]):
                if sigma_type in all_sigma_histories[chain_id]:
                    plt.plot(all_sigma_histories[chain_id][sigma_type], 
                            label=sigma_type, color=palette[i])
            
            plt.title(f'Sigma Trace Plot - Chain {chain_idx} ({sampler_display_name})', fontsize=16)
            plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)
            plt.ylabel('Sigma Value', fontsize=14)
            plt.legend(fontsize=12)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        
        # Score plots
        if all_scores:
            for score_type in all_scores:
                plt.figure(figsize=(12, 8))
                for chain_id in all_scores[score_type]:
                    plt.plot(all_scores[score_type][chain_id], label=chain_id, alpha=0.7)
                plt.title(f'{score_type} ({sampler_display_name})', fontsize=16)
                plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)
                plt.ylabel('Score Value', fontsize=14)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
        
        # R-hat analysis
        _add_rhat_analysis(pdf, all_sigma_histories, sampler_display_name, analysis_output_dir)
        
        # Autocorrelation analysis
        _add_autocorr_analysis(pdf, all_sigma_histories, all_scores, sampler_display_name, analysis_output_dir)


def _add_rhat_analysis(pdf, all_sigma_histories, sampler_display_name, analysis_output_dir):
    """Add R-hat analysis to PDF and save text file."""
    print(f"\nR-hat Statistics for Sigma Components ({sampler_display_name}):")
    
    combined_sigma_data = defaultdict(list)
    for chain_data in all_sigma_histories.values():
        for sigma_type, values in chain_data.items():
            combined_sigma_data[sigma_type].append(values)
    
    rhat_values = {}
    rhat_filename = os.path.join(analysis_output_dir, f"rhat_statistics_{sampler_display_name}.txt")
    
    with open(rhat_filename, 'w') as f:
        f.write(f"R-hat Statistics for Sigma Components ({sampler_display_name}):\n")
        
        for sigma_type, histories in combined_sigma_data.items():
            if len(histories) < 2:
                message = f"  {sigma_type}: R-hat not calculated - insufficient chains"
                print(message)
                f.write(message + "\n")
                rhat_values[sigma_type] = None
            else:
                try:
                    # Use equal length chains
                    min_len = min(len(chain) for chain in histories)
                    truncated = [chain[:min_len] for chain in histories]
                    r_hat = float(az.rhat(np.array(truncated)))
                    rhat_values[sigma_type] = r_hat
                    message = f"  {sigma_type}: {r_hat:.3f} using {len(truncated)} chains (length={min_len})"
                    print(message)
                    f.write(message + "\n")
                except Exception as e:
                    message = f"  {sigma_type}: R-hat calculation failed - {str(e)}"
                    print(message)
                    f.write(message + "\n")
                    rhat_values[sigma_type] = None
    
    # Add R-hat table to PDF
    fig, ax = plt.subplots(figsize=(6, 2 + 0.3*len(rhat_values)))
    fig.suptitle(f"R-hat Values ({sampler_display_name})", fontsize=14)
    ax.axis('off')
    
    rows = [[sigma_type, f"{r_hat:.3f}" if r_hat is not None else "N/A"] 
            for sigma_type, r_hat in rhat_values.items()]
    
    table = ax.table(cellText=rows, colLabels=["Sigma Type", "R-hat"], loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    pdf.savefig(fig)
    plt.close(fig)


def _add_autocorr_analysis(pdf, all_sigma_histories, all_scores, sampler_display_name, analysis_output_dir):
    """Add autocorrelation analysis to PDF."""
    print("\nCalculating autocorrelation and effective sample sizes...")
    
    # Get representative trajectory length for max_lag calculation
    sample_chain = next(iter(all_sigma_histories.values()))
    sample_series = next(iter(sample_chain.values()))
    max_lag = min(50, len(sample_series) // 10)
    
    ess_values = {}
    combined_sigma_data = defaultdict(list)
    for chain_data in all_sigma_histories.values():
        for sigma_type, values in chain_data.items():
            combined_sigma_data[sigma_type].append(values)
    
    # Autocorr plots for sigma parameters
    for sigma_type in combined_sigma_data:
        plt.figure(figsize=(12, 8))
        for chain_idx, chain_values in enumerate(combined_sigma_data[sigma_type]):
            if len(chain_values) > 10:
                autocorr = calculate_autocorrelation(chain_values, nlags=max_lag)
                if autocorr is not None:
                    plt.plot(autocorr, label=f'Chain {chain_idx}', alpha=0.7)
                    
                    ess = calculate_effective_sample_size(chain_values)
                    if ess is not None:
                        ess_values.setdefault(sigma_type, []).append(ess)
        
        plt.axhline(y=0, color='r', linestyle='--')
        plt.title(f'Autocorrelation for {sigma_type} ({sampler_display_name})', fontsize=16)
        plt.xlabel('Lag', fontsize=14)
        plt.ylabel('Autocorrelation', fontsize=14)
        plt.legend(fontsize=12)
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    
    # ESS table
    if ess_values:
        fig, ax = plt.subplots(figsize=(8, 2 + 0.3*len(ess_values)))
        fig.suptitle(f"Effective Sample Size ({sampler_display_name})", fontsize=14)
        ax.axis('off')
        
        rows = []
        for sigma_type, ess_list in ess_values.items():
            if ess_list:
                rows.append([sigma_type, f"{np.mean(ess_list):.1f}", 
                           f"{np.min(ess_list):.1f}", f"{np.max(ess_list):.1f}", str(len(ess_list))])
        
        table = ax.table(cellText=rows, 
                        colLabels=["Sigma Type", "Mean ESS", "Min ESS", "Max ESS", "N Chains"], 
                        loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        pdf.savefig(fig)
        plt.close(fig)
        
        # Save ESS to file
        ess_filename = os.path.join(analysis_output_dir, f"ess_statistics_{sampler_display_name}.txt")
        with open(ess_filename, 'w') as f:
            f.write(f"Effective Sample Size for Sigma Components ({sampler_display_name}):\n")
            for sigma_type, ess_list in ess_values.items():
                if ess_list:
                    f.write(f"{sigma_type}: Mean={np.mean(ess_list):.1f}, "
                           f"Min={np.min(ess_list):.1f}, Max={np.max(ess_list):.1f}, "
                           f"Chains={len(ess_list)}\n")
        print(f"ESS statistics saved to: {ess_filename}")

def _generate_convergence_diagnostics(all_sigma_histories, sampler_display_name, analysis_output_dir, burnin_fraction):
    """Generate advanced convergence diagnostics."""
    diag_chains = {}
    for chain_id, sig_dict in all_sigma_histories.items():
        for sigma_type, series in sig_dict.items():
            # Ensure series is a proper numpy array of floats
            clean_series = []
            for val in series:
                try:
                    if hasattr(val, 'item'):
                        clean_series.append(float(val.item()))
                    else:
                        clean_series.append(float(val))
                except (ValueError, TypeError) as e:
                    print(f"Warning: Could not convert value {val} to float: {e}")
                    continue
            
            if clean_series:  # Only add if we have valid data
                diag_chains.setdefault(sigma_type, []).append(np.array(clean_series))
    
    if diag_chains:
        try:
            diag = AdvancedConvergenceDiagnostics(diag_chains)
            
            # Generate PDF and CSV reports with error handling
            try:
                pdf_path = diag.generate_pdf_report(analysis_output_dir, sampler_display_name, burnin_fraction)
                print(f"[Diagnostics] PDF report saved: {pdf_path}")
            except Exception as e:
                print(f"[Diagnostics] PDF generation failed: {e}")
            
            try:
                csv_path = diag.save_diagnostics_csv(analysis_output_dir, sampler_display_name, burnin_fraction)
                print(f"[Diagnostics] CSV report saved: {csv_path}")
            except Exception as e:
                print(f"[Diagnostics] CSV generation failed: {e}")
                
        except Exception as e:
            print(f"[Diagnostics] Failed to create diagnostics: {e}")
    else:
        print("[Diagnostics] No valid chain data for diagnostics")


def _generate_gmm_fits(all_sigma_histories, sampler_display_name, analysis_output_dir):
    """Generate GMM fits and plots."""
    pdf_filename = os.path.join(analysis_output_dir, f"{sampler_display_name}_combined_gmm_plots.pdf")
    
    with PdfPages(pdf_filename) as pdf:
        for sigma_type in ["AA", "AB", "BC"]:
            print(f"\nFitting GMM for {sigma_type}...")
            all_data_for_type = {}
            all_gmms_for_type = {}

            for chain_id, chain_data in all_sigma_histories.items():
                if sigma_type in chain_data:
                    data = np.array(chain_data[sigma_type])
                    print(f"  Chain {chain_id}: {len(data)} samples")
                    
                    gmm = fit_gmm_robust(data, sigma_type, max_components=4)
                    
                    if gmm is not None:
                        all_data_for_type[chain_id] = data
                        all_gmms_for_type[chain_id] = gmm
                        save_gmm_with_validation(gmm, sigma_type, chain_id, analysis_output_dir, data)

            if all_data_for_type:
                plot_combined_gmm(all_data_for_type, all_gmms_for_type, 
                                sigma_type, sampler_display_name, analysis_output_dir, pdf)
                print(f"Combined GMM plot for {sigma_type} saved.")

def analyze_sampler_in_sequence(sampler_sequence: list, sampler_type: str, sampler_position: int, 
                               output_folder: str = "output_analysis", burnin: float = 0.4, 
                               do_trace_plots: bool = True, do_gmm_fits: bool = True, 
                               do_diagnostics: bool = True):
    """
    Analyze a specific sampler at a given position in a sequence.
    
    Args:
        sampler_sequence (list): Complete sequence of samplers (e.g., ['pair', 'tetramer', 'pair', 'octet'])
        sampler_type (str): Type of sampler to analyze ('pair', 'tetramer', 'octet')
        sampler_position (int): 0-based position in the sequence
        output_folder (str): Base directory for output files
        burnin (float): Fraction of initial samples to discard
        do_trace_plots (bool): Whether to generate trace plots
        do_gmm_fits (bool): Whether to perform GMM fitting
    
    Returns:
        bool: True if analysis was successful, False otherwise
    """
    # Validate inputs
    if sampler_position >= len(sampler_sequence):
        print(f"Error: Position {sampler_position} is out of range for sequence length {len(sampler_sequence)}")
        return False
    
    if sampler_sequence[sampler_position] != sampler_type:
        print(f"Error: Expected '{sampler_type}' at position {sampler_position}, but found '{sampler_sequence[sampler_position]}'")
        return False
    
    # Count occurrences to determine the numbered folder
    occurrence_number = count_sampler_occurrences(sampler_sequence, sampler_type, sampler_position)
    
    print(f"\nAnalyzing {sampler_type} sampler at position {sampler_position} (occurrence #{occurrence_number}) in sequence {sampler_sequence}")
    print(f"Loading trajectories from folder: {get_sampler_folder_name(sampler_type, occurrence_number)}")
    
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Run the analysis
    try:
        analyze_mcmc_data(output_folder, sampler_type, occurrence_number, burnin, 
                          do_trace_plots, do_gmm_fits, do_diagnostics)
        print(f"Analysis for {sampler_type} sampler (occurrence #{occurrence_number}) complete.")
        return True
    except Exception as e:
        print(f"Error analyzing {sampler_type} sampler: {str(e)}")
        return False

def main():
    """Process command line arguments and run analysis."""
    parser = argparse.ArgumentParser(description='Analyze MCMC results and fit GMM models for sampler sequences.')
    
    parser.add_argument('--sequence', type=str, required=True,
                        help='Comma-separated sampler sequence (e.g., "pair,tetramer,pair,octet, full")')

    parser.add_argument('--sampler', type=str, required=True, choices=['pair', 'tetramer', 'octet', 'full'],
                        help='Type of sampler to analyze (pair, tetramer, octet, or full)')

    parser.add_argument('--position', type=int, required=True,
                        help='0-based position of the sampler in the sequence')
    
    parser.add_argument('--output', '-o', type=str, default='output_analysis',
                        help='Output folder for analysis results (default: output_analysis)')
    
    parser.add_argument('--no-traces', action='store_true',
                        help='Disable trace plots')
    
    parser.add_argument('--no-gmm', action='store_true',
                        help='Disable GMM fitting')

    parser.add_argument('--no-diagnostics', action='store_true',
                        help='Disable convergence diagnostics')
    parser.add_argument('--burnin', '-b', type=str, default='0.4',
                       help='Burn-in fraction (e.g. 0.4) or "auto" for automatic detection')
    
    args = parser.parse_args()
    
    # Fix the burnin handling - parse ONCE and correctly
    if args.burnin and args.burnin.lower() == "auto":
        burnin = "auto"  # Keep as string literal
    else:
        try:
            burnin = float(args.burnin)
            if not (0.0 <= burnin < 1.0):
                print(f"Warning: burnin {burnin} not in [0,1), using 0.4")
                burnin = 0.4
        except (ValueError, TypeError):
            print(f"Warning: invalid burnin value '{args.burnin}', using 0.4")
            burnin = 0.4
                
    # Parse the sequence
    try:
        sampler_sequence = [s.strip().lower() for s in args.sequence.split(',')]
        # Validate sequence
        valid_samplers = ['pair', 'tetramer', 'octet', 'full']
        for sampler in sampler_sequence:
            if sampler not in valid_samplers:
                print(f"Error: Invalid sampler '{sampler}' in sequence. Must be one of {valid_samplers}")
                return
    except Exception as e:
        print(f"Error parsing sequence: {e}")
        return
    

            
    # Run analysis with specified parameters
    success = analyze_sampler_in_sequence(
        sampler_sequence=sampler_sequence,
        sampler_type=args.sampler,
        sampler_position=args.position,
        output_folder=args.output, 
        burnin=burnin,
        do_trace_plots=not args.no_traces,
        do_gmm_fits=not args.no_gmm,
        do_diagnostics=not args.no_diagnostics
    )
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()