import os
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

def analyze_mcmc_data(output_folder: str, sampler_type: str, sampler_position: int, burnin: float = 0.4, 
                     do_trace_plots: bool = True, do_gmm_fits: bool = True, 
                     do_diagnostics: bool = True):
    """
    Analyzes MCMC data, including trace plots, R-hat statistics, and GMM fitting.

    Args:
        output_folder (str): Base directory for output files.
        sampler_type (str): Type of sampler ('pair', 'tetramer', 'octet').
        sampler_position (int): Position of this sampler in the sequence (1-based).
        burnin (float): Fraction of initial samples to discard.
        do_trace_plots (bool): Whether to generate trace plots and score plots.
        do_gmm_fits (bool): Whether to perform GMM fitting and plotting.
    """

    if not 0.0 <= burnin < 1.0:
        raise ValueError("burnin must be between 0.0 and 1.0")

    # Generate folder name based on sampler type and position
    sampler_folder_name = get_sampler_folder_name(sampler_type, sampler_position)
    sampler_output_dir = os.path.join(output_folder, sampler_folder_name)
    
    if not os.path.exists(sampler_output_dir):
        print(f"Directory not found: {sampler_output_dir}")
        return

    # Create analysis output directory within the sampler folder
    #analysis_output_dir = os.path.join(sampler_output_dir, "analysis")
    analysis_output_dir = sampler_output_dir
    if not os.path.exists(analysis_output_dir):
        os.makedirs(analysis_output_dir)

    trajectory_files = [f for f in os.listdir(sampler_output_dir) if f.startswith("trajectory_") and f.endswith(".h5")]
    if not trajectory_files:
        print(f"No trajectory files found in: {sampler_output_dir}")
        return

    # Generate descriptive sampler name for plots
    sampler_display_name = f"{sampler_type.capitalize()}Sampler_{sampler_position}"

    all_sigma_histories = defaultdict(lambda: defaultdict(list))
    all_scores = defaultdict(lambda: defaultdict(list))  # Only used if do_trace_plots

    for traj_file in trajectory_files:
        level_name = traj_file.replace("trajectory_", "").replace(".h5", "")
        traj_path = os.path.join(sampler_output_dir, traj_file)
        trajectory_data = load_trajectory_from_hdf5(traj_path)

        num_samples = len(trajectory_data)
        burnin_samples = int(num_samples * burnin)
        trajectory_data = trajectory_data[burnin_samples:]

        for state in trajectory_data:
            for sigma_type, sigma_value in state['sigma'].items():
                all_sigma_histories[level_name][sigma_type].append(sigma_value)
            if do_trace_plots:
                all_scores['total_score'][level_name].append(state['total_score'])
                all_scores['prior_score'][level_name].append(state['prior_score'])
                all_scores['pair_score'][level_name].append(state['pair_score'])
                all_scores['exvol_score'][level_name].append(state['exvol_score'])
                all_scores['tet_score'][level_name].append(state['tet_score'])      
                all_scores['oct_score'][level_name].append(state['oct_score'])      
    
    # AUTO burn-in adjustment
    if isinstance(burnin, str) and burnin.lower() == "auto":
        # Build chains dict for diagnostics (per sigma param)
        diag_chains = {}
        for chain_id, sig_dict in all_sigma_histories.items():
            for sigma_type, series in sig_dict.items():
                diag_chains.setdefault(sigma_type, []).append(np.array(series))
        diag = AdvancedConvergenceDiagnostics(diag_chains)
        auto_frac = diag.auto_burnin_fraction()
        print(f"[Diagnostics] Auto burn-in fraction selected: {auto_frac:.3f}")
        burnin = auto_frac
        # Re-apply burn-in to histories
        new_hist = {}
        for chain_id, sig_dict in all_sigma_histories.items():
            new_hist[chain_id] = {}
            for sigma_type, series in sig_dict.items():
                n = len(series)
                cut = int(burnin * n)
                new_hist[chain_id][sigma_type] = series[cut:]
        all_sigma_histories = new_hist
        
    # --- Trace Plots and R-hat (if requested) ---
    if do_trace_plots:
        pdf_filename_trace = os.path.join(analysis_output_dir, f"{sampler_display_name}_report.pdf")
        with PdfPages(pdf_filename_trace) as pdf:
            sns.set(style="darkgrid")
            palette = sns.color_palette("husl", 4)

            # --- Sigma trace plots ---
            for chain_idx, chain_id in enumerate(all_sigma_histories):
                plt.figure(figsize=(12, 8))
                #for i, sigma_type in enumerate(["AA", "AB", "BC", "CC"]):
                for i, sigma_type in enumerate(["AA", "AB", "BC"]):
                    if sigma_type in all_sigma_histories[chain_id]:
                        plt.plot(all_sigma_histories[chain_id][sigma_type], label=sigma_type, color=palette[i])
                plt.title(f'Combined Sigma Trace Plot - Chain {chain_idx} ({sampler_display_name})', fontsize=16)
                plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)
                plt.ylabel('Sigma Value', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.savefig(os.path.join(analysis_output_dir, f'trace_sigma_combined_chain{chain_idx}_{sampler_display_name}.png'))
                plt.close()

            # --- Score plots ---
            for score_type in all_scores:
                plt.figure(figsize=(12, 8))
                for chain_id in all_scores[score_type]:
                    plt.plot(all_scores[score_type][chain_id], label=chain_id, alpha=0.7)
                plt.title(f'Score vs. Frame: {score_type} ({sampler_display_name})', fontsize=16)
                plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)
                plt.ylabel('Score Value', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.savefig(os.path.join(analysis_output_dir, f'score_vs_frame_{score_type}_{sampler_display_name}.png'))
                plt.close()

            print(f"\nR-hat Statistics for Sigma Components ({sampler_display_name}):")
            rhat_values = {}
            rhat_filename = os.path.join(analysis_output_dir, f"rhat_statistics_{sampler_display_name}.txt")

            combined_sigma_data = defaultdict(list)
            for chain_data in all_sigma_histories.values():
                for sigma_type, values in chain_data.items():
                    combined_sigma_data[sigma_type].append(values)
            
            # Replace the R-hat calculation section with this version:
            with open(rhat_filename, 'w') as f:
                f.write(f"R-hat Statistics for Sigma Components ({sampler_display_name}):\n")
                rhat_values = {}
                for sigma_type, histories in combined_sigma_data.items():
                    chain_lengths = [len(chain) for chain in histories]
                    max_length = max(chain_lengths)
                    # Filter out trajectories that are shorter than the maximum length.
                    filtered_histories = [chain for chain in histories if len(chain) == max_length]
                    dropped = len(histories) - len(filtered_histories)
                    
                    if len(filtered_histories) < 2:
                        message = (f"  {sigma_type}: R-hat not calculated - only "
                                   f"{len(filtered_histories)} chain(s) available after filtering "
                                   f"(dropped {dropped} chain(s) out of {len(histories)})")
                        print(message)
                        f.write(message + "\n")
                        rhat_values[sigma_type] = None
                    else:
                        try:
                            # All filtered chains have the same length (== max_length)
                            r_hat = az.rhat(np.array(filtered_histories))
                            rhat_values[sigma_type] = r_hat
                            message = (f"  {sigma_type}: {r_hat:.3f} calculated using "
                                       f"{len(filtered_histories)} chains (max_length = {max_length}, "
                                       f"dropped {dropped} chain(s) out of {len(histories)})")
                            print(message)
                            f.write(message + "\n")
                        except Exception as e:
                            message = f"  {sigma_type}: R-hat calculation failed - {str(e)}"
                            print(message)
                            f.write(message + "\n")
                            rhat_values[sigma_type] = None

            # --- Add R-hat table to the same PDF ---
            fig, ax = plt.subplots(figsize=(6, 2 + 0.3*len(rhat_values)))  # Adjust figure size to fit the table
            fig.suptitle(f"R-hat Values for Sigma Components ({sampler_display_name})", fontsize=14)
            ax.axis('off')

            # Build table data
            header = ["Sigma Type", "R-hat"]
            rows = []
            for sigma_type, r_hat in rhat_values.items():
                # Format r_hat as a float; handle the None case with 'N/A'
                val_str = f"{r_hat:.3f}" if r_hat is not None else "N/A"
                rows.append([sigma_type, val_str])

            # Create table in the axes
            table = ax.table(
                cellText=rows,
                colLabels=header,
                loc='center'
            )
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1, 1.5)

            pdf.savefig(fig)
            plt.close(fig)
            
            # --- Add autocorrelation plots to PDF ---
            print("\nCalculating autocorrelation and effective sample sizes...")
            ess_values = {}
            
            # Analyze each sigma parameter across chains
            for sigma_type in combined_sigma_data:
                # Process each chain separately
                plt.figure(figsize=(12, 8))
                max_lag = min(50, int(len(trajectory_data) / 10))  # Use reasonable max lag
                
                for chain_idx, chain_values in enumerate(combined_sigma_data[sigma_type]):
                    if len(chain_values) > 10:  # Need sufficient data points
                        autocorr = calculate_autocorrelation(chain_values, nlags=max_lag)
                        if autocorr is not None:
                            plt.plot(autocorr, label=f'Chain {chain_idx}', alpha=0.7)
                            
                            # Calculate ESS
                            ess = calculate_effective_sample_size(chain_values)
                            if ess is not None:
                                if sigma_type not in ess_values:
                                    ess_values[sigma_type] = []
                                ess_values[sigma_type].append(ess)
                
                plt.axhline(y=0, color='r', linestyle='--')
                plt.title(f'Autocorrelation for {sigma_type} ({sampler_display_name})', fontsize=16)
                plt.xlabel('Lag', fontsize=14)
                plt.ylabel('Autocorrelation', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.savefig(os.path.join(analysis_output_dir, f'autocorr_{sigma_type}_{sampler_display_name}.png'))
                plt.close()
                
            # --- Add ESS table to PDF ---
            if ess_values:
                fig, ax = plt.subplots(figsize=(8, 2 + 0.3*len(ess_values)))
                fig.suptitle(f"Effective Sample Size (ESS) for Sigma Components ({sampler_display_name})", fontsize=14)
                ax.axis('off')
                
                # Build table data
                header = ["Sigma Type", "Mean ESS", "Min ESS", "Max ESS", "N Chains"]
                rows = []
                
                for sigma_type, ess_list in ess_values.items():
                    if ess_list:
                        rows.append([
                            sigma_type,
                            f"{np.mean(ess_list):.1f}",
                            f"{np.min(ess_list):.1f}",
                            f"{np.max(ess_list):.1f}",
                            f"{len(ess_list)}"
                        ])
                    else:
                        rows.append([sigma_type, "N/A", "N/A", "N/A", "0"])
                
                # Create table in the axes
                table = ax.table(
                    cellText=rows,
                    colLabels=header,
                    loc='center'
                )
                table.auto_set_font_size(False)
                table.set_fontsize(10)
                table.scale(1, 1.5)
                
                pdf.savefig(fig)
                plt.close(fig)
                
                # Save ESS values to file
                ess_filename = os.path.join(analysis_output_dir, f"ess_statistics_{sampler_display_name}.txt")
                with open(ess_filename, 'w') as f:
                    f.write(f"Effective Sample Size (ESS) for Sigma Components ({sampler_display_name}):\n")
                    for sigma_type, ess_list in ess_values.items():
                        if ess_list:
                            f.write(f"{sigma_type}: Mean={np.mean(ess_list):.1f}, Min={np.min(ess_list):.1f}, " 
                                   f"Max={np.max(ess_list):.1f}, Chains={len(ess_list)}\n")
                        else:
                            f.write(f"{sigma_type}: No valid ESS values\n")
                print(f"ESS statistics saved to: {ess_filename}")
            
            # --- Also add autocorrelation for score components ---
            if all_scores:
                for score_type in all_scores:
                    plt.figure(figsize=(12, 8))
                    max_lag = min(50, int(len(trajectory_data) / 10))
                    
                    for chain_id, score_values in all_scores[score_type].items():
                        if len(score_values) > 10:
                            autocorr = calculate_autocorrelation(score_values, nlags=max_lag)
                            if autocorr is not None:
                                plt.plot(autocorr, label=f'Chain {chain_id}', alpha=0.7)
                    
                    plt.axhline(y=0, color='r', linestyle='--')
                    plt.title(f'Autocorrelation for {score_type} ({sampler_display_name})', fontsize=16)
                    plt.xlabel('Lag', fontsize=14)
                    plt.ylabel('Autocorrelation', fontsize=14)
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.legend(fontsize=12)
                    plt.tight_layout()
                    pdf.savefig()
                    plt.savefig(os.path.join(analysis_output_dir, f'autocorr_{score_type}_{sampler_display_name}.png'))
                    plt.close()
                    
    # Run convergence diagnostics (post-burn-in) if requested
    if do_diagnostics:
        diag_chains2 = {}
        for chain_id, sig_dict in all_sigma_histories.items():
            for sigma_type, series in sig_dict.items():
                diag_chains2.setdefault(sigma_type, []).append(np.array(series))
        diag2 = AdvancedConvergenceDiagnostics(diag_chains2)
        df_diag = diag2.compute_summary_dataframe(burnin_fraction=0.0)
        diag_csv = os.path.join(analysis_output_dir, f"{sampler_display_name}_diagnostics.csv")
        df_diag.to_csv(diag_csv, index=False)
        print(f"[Diagnostics] Saved summary to {diag_csv}")
            
    # --- GMM Fits and Plots (if requested) ---
    if do_gmm_fits:
        pdf_filename_gmm = os.path.join(analysis_output_dir, f"{sampler_display_name}_combined_gmm_plots.pdf")
        with PdfPages(pdf_filename_gmm) as pdf:
            sigma_types = ["AA", "AB", "BC"]
            
            for sigma_type in sigma_types:
                print(f"\nFitting GMM for {sigma_type}...")
                all_data_for_type = {}
                all_gmms_for_type = {}

                for chain_id, chain_data in all_sigma_histories.items():
                    if sigma_type in chain_data:
                        data = np.array(chain_data[sigma_type])
                        print(f"  Chain {chain_id}: {len(data)} samples")
                        
                        # Use robust GMM fitting
                        gmm = fit_gmm_robust(data, sigma_type, max_components=4)  # Reduced max components
                        
                        if gmm is not None:
                            all_data_for_type[chain_id] = data
                            all_gmms_for_type[chain_id] = gmm
                            
                            # Save with validation
                            save_gmm_with_validation(gmm, sigma_type, chain_id, analysis_output_dir, data)
                        else:
                            print(f"  Failed to fit GMM for {sigma_type} chain {chain_id}")

                if all_data_for_type:
                    plot_combined_gmm(all_data_for_type, all_gmms_for_type, sigma_type, sampler_display_name, analysis_output_dir, pdf)
                    print(f"Combined GMM plot for {sigma_type} saved.")
                else:
                    print(f"No valid GMMs for {sigma_type} - skipping plot")

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

    parser.add_argument('--burnin', '-b', type=str, default="auto",
                        help='Burn-in fraction (e.g. 0.4) or "auto" for automatic detection')

    parser.add_argument('--no-diagnostics', action='store_true',
                        help='Disable convergence diagnostics')
    
    # Parse arguments
    args = parser.parse_args()
    
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
    
# After parsing:
    # Parse burnin numeric or auto
    try:
        if args.burnin.lower() != "auto":
            burnin_val = float(args.burnin)
        else:
            burnin_val = "auto"
    except:
        print("Invalid --burnin value; use float or 'auto'")
        return    
    
    # Run analysis with specified parameters
    success = analyze_sampler_in_sequence(
        sampler_sequence=sampler_sequence,
        sampler_type=args.sampler,
        sampler_position=args.position,
        output_folder=args.output, 
        burnin=burnin_val,
        do_trace_plots=not args.no_traces,
        do_gmm_fits=not args.no_gmm,
        do_diagnostics=not args.no_diagnostics
    )
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()