import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import arviz as az
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

def load_trajectory_from_hdf5(filename: str) -> list:
    """Loads MCMC trajectory data from an HDF5 file."""
    trajectory = []
    with h5py.File(filename, 'r') as f:
        traj_grp = f['trajectory']
        for state_name in traj_grp:
            state_grp = traj_grp[state_name]
            state = {
                "step": state_grp.attrs["step"],
                "total_score": state_grp.attrs["total_score"],
                "prior_score": state_grp.attrs["prior_score"],
                "pair_score": state_grp.attrs["pair_score"],
                "exvol_score": state_grp.attrs["exvol_score"],
                "sigma": {},
                "positions": {},
                "types": {},
                "bead_numbers": {}
            }
            sigma_grp = state_grp['sigma']
            for key in sigma_grp.attrs:
                state["sigma"][key] = sigma_grp.attrs[key]
            pos_grp = state_grp['positions']
            for type_name in pos_grp:
                state["positions"][type_name] = pos_grp[type_name][:]
            types_keys_dataset = state_grp['types_keys'][:]
            types_vals_dataset = state_grp['types_vals'][:]
            bead_keys_dataset = state_grp['bead_keys'][:]
            bead_vals_dataset = state_grp['bead_vals'][:]

            state["types"] = {types_keys_dataset[i].decode('utf-8'): types_vals_dataset[i].decode('utf-8') for i in range(len(types_keys_dataset))}
            state["bead_numbers"] = {int(bead_keys_dataset[i]): int(bead_vals_dataset[i]) for i in range(len(bead_keys_dataset))}
            trajectory.append(state)
    return trajectory


def analyze_sampler_results(output_folder: str, sampler_name: str, burnin: float = 0.0):
    """Analyzes MCMC results, generates plots, and compiles a PDF report.

    Args:
        output_folder (str): Base directory for output files.
        sampler_name (str): Name of the sampler (e.g., "PairSampler").
        burnin (float): Fraction of initial samples to discard as burn-in (0.0 to 1.0).
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
    all_scores = defaultdict(lambda: defaultdict(list))

    for traj_file in trajectory_files:
        level_name = traj_file.replace("trajectory_", "").replace(".h5", "")
        traj_path = os.path.join(sampler_output_dir, traj_file)
        trajectory_data = load_trajectory_from_hdf5(traj_path)

        # Apply burn-in
        num_samples = len(trajectory_data)
        burnin_samples = int(num_samples * burnin)
        trajectory_data = trajectory_data[burnin_samples:]

        for state in trajectory_data:
            for sigma_type, sigma_value in state['sigma'].items():
                all_sigma_histories[level_name][sigma_type].append(sigma_value)
            all_scores['total_score'][level_name].append(state['total_score'])
            all_scores['prior_score'][level_name].append(state['prior_score'])
            all_scores['pair_score'][level_name].append(state['pair_score'])
            all_scores['exvol_score'][level_name].append(state['exvol_score'])


    # --- Generate Combined Sigma Trace Plots (one per chain) ---
    pdf_filename = os.path.join(sampler_output_dir, f"{sampler_name}_report.pdf")
    with PdfPages(pdf_filename) as pdf:
        sns.set(style="darkgrid")
        palette = sns.color_palette("husl", 4)

        for chain_idx, chain_id in enumerate(all_sigma_histories):
            plt.figure(figsize=(12, 8))
            for i, sigma_type in enumerate(["AA", "AB", "BC", "CC"]):
                if sigma_type in all_sigma_histories[chain_id]:
                    plt.plot(all_sigma_histories[chain_id][sigma_type], label=sigma_type, color=palette[i])
            plt.title(f'Combined Sigma Trace Plot - Chain {chain_idx} ({sampler_name})', fontsize=16)
            plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14)  # Updated x-axis label
            plt.ylabel('Sigma Value', fontsize=14)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.legend(fontsize=12)
            plt.tight_layout()
            pdf.savefig()
            plt.savefig(os.path.join(sampler_output_dir, f'trace_sigma_combined_chain{chain_idx}_{sampler_name}.png'))
            plt.close()

        # --- Generate Score vs. Frame Plots (keep as separate plots) ---
        for score_type in all_scores:
            plt.figure(figsize=(12, 8))
            for chain_id in all_scores[score_type]:
                plt.plot(all_scores[score_type][chain_id], label=chain_id, alpha=0.7)
            plt.title(f'Score vs. Frame: {score_type} ({sampler_name})', fontsize=16)
            plt.xlabel('MCMC Step (Post Burn-in)', fontsize=14) # Updated x-axis label
            plt.ylabel('Score Value', fontsize=14)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.legend(fontsize=12)
            plt.tight_layout()
            pdf.savefig()
            plt.savefig(os.path.join(sampler_output_dir, f'score_vs_frame_{score_type}_{sampler_name}.png'))
            plt.close()

        # --- Calculate and Print/Save R-hat Statistic ---
        print(f"\nR-hat Statistics for Sigma Components ({sampler_name}):")
        rhat_values = {}
        rhat_filename = os.path.join(sampler_output_dir, f"rhat_statistics_{sampler_name}.txt")

        combined_sigma_data = defaultdict(list)
        for chain_data in all_sigma_histories.values():
            for sigma_type, values in chain_data.items():
                combined_sigma_data[sigma_type].append(values)

        with open(rhat_filename, 'w') as f:
            f.write(f"R-hat Statistics for Sigma Components ({sampler_name}):\n")
            for sigma_type, histories in combined_sigma_data.items():
                if len(histories) >= 2:
                    r_hat = az.rhat(np.array(histories))
                    rhat_values[sigma_type] = r_hat
                    print(f"  {sigma_type}: {r_hat:.3f}")
                    f.write(f"  {sigma_type}: {r_hat:.3f}\n")
                else:
                    print(f"  {sigma_type}: R-hat not calculated (only {len(histories)} chain(s)).")
                    f.write(f"  {sigma_type}: R-hat not calculated (only {len(histories)} chain(s)).\n")
                    rhat_values[sigma_type] = None
    print(f"R-hat statistics saved to: {rhat_filename}")


def main():
    output_folder = "output_analysis"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    sampler_sequence = ["pair_sampler"]
    burnin = 0.3  # Set the burn-in fraction here

    for sampler_key in sampler_sequence:
        sampler_name = sampler_key.replace("_sampler", "").capitalize() + "Sampler"
        print(f"\nAnalyzing results for {sampler_name} with burn-in fraction: {burnin}...")
        analyze_sampler_results(output_folder, sampler_name, burnin)
        print(f"Analysis for {sampler_name} complete.  Report saved to '{output_folder}/{sampler_name.lower()}_results/{sampler_name}_report.pdf'")

if __name__ == "__main__":
    main()