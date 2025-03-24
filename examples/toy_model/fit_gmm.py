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

def fit_gmm(data: np.ndarray, max_components: int = 5):
    """Fits a Gaussian Mixture Model to the data."""
    data = data.reshape(-1, 1)
    best_gmm = None
    best_bic = np.inf

    for n_components in range(1, max_components + 1):
        gmm = GaussianMixture(n_components=n_components, covariance_type='full', random_state=42)
        gmm.fit(data)
        bic = gmm.bic(data)
        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm

    return best_gmm

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

def analyze_mcmc_data(output_folder: str, sampler_name: str, burnin: float = 0.3, do_trace_plots: bool = True, do_gmm_fits: bool = True):
    """
    Analyzes MCMC data, including trace plots, R-hat statistics, and GMM fitting.

    Args:
        output_folder (str): Base directory for output files.
        sampler_name (str): Name of the sampler (e.g., "PairSampler").
        burnin (float): Fraction of initial samples to discard.
        do_trace_plots (bool): Whether to generate trace plots and score plots.
        do_gmm_fits (bool):  Whether to perform GMM fitting and plotting.
    """

    if not 0.0 <= burnin < 1.0:
        raise ValueError("burnin must be between 0.0 and 1.0")

    sampler_output_dir = os.path.join(output_folder, f"{sampler_name.lower()}_results")
    if not os.path.exists(sampler_output_dir):
        os.makedirs(sampler_output_dir)

    trajectory_files = [f for f in os.listdir(sampler_output_dir) if f.startswith("trajectory_") and f.endswith(".h5")]
    if not trajectory_files:
        print(f"No trajectory files found in: {sampler_output_dir}")
        return

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

    # --- Trace Plots and R-hat (if requested) ---
    if do_trace_plots:
        pdf_filename_trace = os.path.join(sampler_output_dir, f"{sampler_name}_report.pdf")
        with PdfPages(pdf_filename_trace) as pdf:
            sns.set(style="darkgrid")
            palette = sns.color_palette("husl", 4)

            # --- Sigma trace plots ---
            for chain_idx, chain_id in enumerate(all_sigma_histories):
                plt.figure(figsize=(12, 8))
                for i, sigma_type in enumerate(["AA", "AB", "BC", "CC"]):
                    if sigma_type in all_sigma_histories[chain_id]:
                        plt.plot(all_sigma_histories[chain_id][sigma_type], label=sigma_type, color=palette[i])
                plt.title(f'Combined Sigma Trace Plot - Chain {chain_idx} ({sampler_name})', fontsize=16)
                plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)
                plt.ylabel('Sigma Value', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.savefig(os.path.join(sampler_output_dir, f'trace_sigma_combined_chain{chain_idx}_{sampler_name}.png'))
                plt.close()

            # --- Score plots ---
            for score_type in all_scores:
                plt.figure(figsize=(12, 8))
                for chain_id in all_scores[score_type]:
                    plt.plot(all_scores[score_type][chain_id], label=chain_id, alpha=0.7)
                plt.title(f'Score vs. Frame: {score_type} ({sampler_name})', fontsize=16)
                plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)
                plt.ylabel('Score Value', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.savefig(os.path.join(sampler_output_dir, f'score_vs_frame_{score_type}_{sampler_name}.png'))
                plt.close()

            print(f"\nR-hat Statistics for Sigma Components ({sampler_name}):")
            rhat_values = {}
            rhat_filename = os.path.join(sampler_output_dir, f"rhat_statistics_{sampler_name}.txt")

            combined_sigma_data = defaultdict(list)
            for chain_data in all_sigma_histories.values():
                for sigma_type, values in chain_data.items():
                    combined_sigma_data[sigma_type].append(values)
            
            # Replace the R-hat calculation section with this version:
            with open(rhat_filename, 'w') as f:
                f.write(f"R-hat Statistics for Sigma Components ({sampler_name}):\n")
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
            fig.suptitle(f"R-hat Values for Sigma Components ({sampler_name})", fontsize=14)
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
                plt.title(f'Autocorrelation for {sigma_type} ({sampler_name})', fontsize=16)
                plt.xlabel('Lag', fontsize=14)
                plt.ylabel('Autocorrelation', fontsize=14)
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.legend(fontsize=12)
                plt.tight_layout()
                pdf.savefig()
                plt.savefig(os.path.join(sampler_output_dir, f'autocorr_{sigma_type}_{sampler_name}.png'))
                plt.close()
                
            # --- Add ESS table to PDF ---
            if ess_values:
                fig, ax = plt.subplots(figsize=(8, 2 + 0.3*len(ess_values)))
                fig.suptitle(f"Effective Sample Size (ESS) for Sigma Components ({sampler_name})", fontsize=14)
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
                ess_filename = os.path.join(sampler_output_dir, f"ess_statistics_{sampler_name}.txt")
                with open(ess_filename, 'w') as f:
                    f.write(f"Effective Sample Size (ESS) for Sigma Components ({sampler_name}):\n")
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
                    plt.title(f'Autocorrelation for {score_type} ({sampler_name})', fontsize=16)
                    plt.xlabel('Lag', fontsize=14)
                    plt.ylabel('Autocorrelation', fontsize=14)
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.legend(fontsize=12)
                    plt.tight_layout()
                    pdf.savefig()
                    plt.savefig(os.path.join(sampler_output_dir, f'autocorr_{score_type}_{sampler_name}.png'))
                    plt.close()

    # --- GMM Fits and Plots (if requested) ---
    if do_gmm_fits:
        pdf_filename_gmm = os.path.join(sampler_output_dir, f"{sampler_name}_combined_gmm_plots.pdf")
        with PdfPages(pdf_filename_gmm) as pdf:
            sigma_types = ["AA", "AB", "BC", "CC"]
            for sigma_type in sigma_types:
                all_data_for_type = {}
                all_gmms_for_type = {}

                for chain_id, chain_data in all_sigma_histories.items():
                    if sigma_type in chain_data:
                        data = np.array(chain_data[sigma_type])
                        all_data_for_type[chain_id] = data
                        gmm = fit_gmm(data)
                        all_gmms_for_type[chain_id] = gmm

                        gmm_params = {
                            "n_components": gmm.n_components,
                            "means": gmm.means_.flatten().tolist(),
                            "covariances": gmm.covariances_.flatten().tolist(),
                            "weights": gmm.weights_.tolist(),
                        }
                        json_filename = os.path.join(sampler_output_dir, f"gmm_fit_{sigma_type}_{chain_id}.json")
                        with open(json_filename, "w") as f:
                            json.dump(gmm_params, f, indent=4)

                if all_data_for_type:
                    plot_combined_gmm(all_data_for_type, all_gmms_for_type, sigma_type, sampler_name, sampler_output_dir, pdf)
                    print(f"Combined GMM plot for {sigma_type} saved.")

def analyze_sampler(sampler_key, output_folder="output_analysis", burnin=0.4, 
                   do_trace_plots=True, do_gmm_fits=True):
    """
    Analyze a specific sampler's MCMC results.
    
    Args:
        sampler_key (str): Name of the sampler ('pair_sampler', 'tetramer_sampler', or 'octet_sampler')
        output_folder (str): Base directory for output files
        burnin (float): Fraction of initial samples to discard
        do_trace_plots (bool): Whether to generate trace plots
        do_gmm_fits (bool): Whether to perform GMM fitting
    
    Returns:
        bool: True if analysis was successful, False otherwise
    """
    valid_samplers = ["pair_sampler", "tetramer_sampler", "octet_sampler"]
    
    if sampler_key not in valid_samplers:
        print(f"Error: Invalid sampler '{sampler_key}'. Must be one of {valid_samplers}")
        return False
    
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    # Convert sampler_key to class name (e.g., "pair_sampler" -> "PairSampler")
    sampler_name = sampler_key.replace("_sampler", "").capitalize() + "Sampler"
    
    print(f"\nAnalyzing data for {sampler_name} with burn-in fraction: {burnin}...")
    
    # Run the analysis
    try:
        analyze_mcmc_data(output_folder, sampler_name, burnin, do_trace_plots, do_gmm_fits)
        print(f"Analysis for {sampler_name} complete.")
        return True
    except Exception as e:
        print(f"Error analyzing {sampler_name}: {str(e)}")
        return False

def main():
    """Process command line arguments and run analysis."""
    parser = argparse.ArgumentParser(description='Analyze MCMC results and fit GMM models.')
    
    parser.add_argument('sampler', type=str, choices=['pair_sampler', 'tetramer_sampler', 'octet_sampler'],
                        help='Sampler to analyze (pair_sampler, tetramer_sampler, or octet_sampler)')
    
    parser.add_argument('--output', '-o', type=str, default='output_analysis',
                        help='Output folder for analysis results (default: output_analysis)')
    
    parser.add_argument('--burnin', '-b', type=float, default=0.4,
                        help='Burn-in fraction (default: 0.4)')
    
    parser.add_argument('--no-traces', action='store_true',
                        help='Disable trace plots')
    
    parser.add_argument('--no-gmm', action='store_true',
                        help='Disable GMM fitting')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run analysis with specified parameters
    analyze_sampler(
        args.sampler, 
        output_folder=args.output, 
        burnin=args.burnin,
        do_trace_plots=not args.no_traces,
        do_gmm_fits=not args.no_gmm
    )

if __name__ == "__main__":
    main()
