import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import argparse
import json
import sys
from collections import defaultdict
from typing import Dict, List

import arviz as az
import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.mixture import GaussianMixture


def load_trajectory_from_hdf5(filename: str) -> List[dict]:
    """Load states written by save_state_to_disk."""
    states = []
    if not os.path.isfile(filename):
        print(f"[load] Missing file: {filename}")
        return states

    with h5py.File(filename, "r") as fh:
        if "trajectory" not in fh:
            print(f"[load] No 'trajectory' group in {filename}")
            return states

        for name, grp in fh["trajectory"].items():
            state = {
                "step": int(grp.attrs.get("step", 0)),
                "sigma": {},
            }

            if "sigma" in grp:
                for key, val in grp["sigma"].attrs.items():
                    if hasattr(val, "item"):
                        state["sigma"][key] = float(val.item())
                    elif isinstance(val, (int, float, np.integer, np.floating)):
                        state["sigma"][key] = float(val)
            states.append(state)

    states.sort(key=lambda s: s["step"])
    return states


def fit_gmm_simple(
    data: np.ndarray,
    sigma_type: str,
    max_components: int = 4,
    min_components: int = 1,
    random_state: int = 42,
) -> GaussianMixture | None:
    """Fit a 1D GaussianMixture using BIC selection."""
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if data.ndim != 1:
        data = data.ravel()

    if len(data) < 5:
        print(f"[GMM] {sigma_type}: need >=5 samples (have {len(data)})")
        return None

    X = data.reshape(-1, 1)
    best_model = None
    best_bic = np.inf
    max_k = min(max_components, len(data))

    for k in range(min_components, max_k + 1):
        try:
            model = GaussianMixture(
                n_components=k,
                covariance_type="diag",
                random_state=random_state,
                max_iter=200,
                reg_covar=1e-8,
            ).fit(X)
            bic = model.bic(X)
            if bic < best_bic:
                best_bic = bic
                best_model = model
        except Exception as exc:
            print(f"[GMM] {sigma_type}: fit failed for k={k}: {exc}")

    if best_model is None:
        print(f"[GMM] {sigma_type}: no valid mixture fit")
        return None

    print(f"[GMM] {sigma_type}: selected k={best_model.n_components}, BIC={best_bic:.2f}")
    return best_model


def save_gmm_simple(
    gmm: GaussianMixture | None,
    sigma_type: str,
    chain_id: str,
    output_dir: str,
    data: np.ndarray,
) -> bool:
    """Persist fitted GMM parameters to JSON."""
    if gmm is None:
        return False

    os.makedirs(output_dir, exist_ok=True)
    payload = {
        "sigma_type": sigma_type,
        "chain_id": chain_id,
        "n_components": int(gmm.n_components),
        "means": gmm.means_.ravel().tolist(),
        "variances": gmm.covariances_.ravel().tolist(),
        "weights": gmm.weights_.ravel().tolist(),
        "bic": float(gmm.bic(data.reshape(-1, 1))),
        "aic": float(gmm.aic(data.reshape(-1, 1))),
        "n_iter": int(gmm.n_iter_),
        "converged": bool(getattr(gmm, "converged_", True)),
        "n_samples": int(len(data)),
        "data_mean": float(np.mean(data)),
        "data_std": float(np.std(data)),
        "data_min": float(np.min(data)),
        "data_max": float(np.max(data)),
    }

    outfile = os.path.join(output_dir, f"gmm_fit_{sigma_type}_{chain_id}.json")
    with open(outfile, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return True


# Backward-compatible aliases
fit_gmm_robust = fit_gmm_simple
fit_gmm = fit_gmm_simple
save_gmm_with_validation = save_gmm_simple


def get_sampler_folder_name(sampler_type: str, sampler_position: int) -> str:
    return f"{sampler_type}sampler_results_{sampler_position}"


def count_sampler_occurrences(sequence: List[str], target: str, position: int) -> int:
    return sum(1 for i in range(position + 1) if sequence[i] == target)


def apply_burnin_with_steps(series: List[float], steps: List[int], fraction: float) -> tuple[List[float], List[int]]:
    """Apply burn-in and return both values and corresponding step numbers."""
    if not series:
        return series, steps
    cut = int(len(series) * fraction)
    return series[cut:], steps[cut:]


def compute_rhat_table(all_sigma_histories: Dict[str, Dict[str, List[float]]]) -> Dict[str, float | None]:
    """Return R-hat per sigma type (requires >=2 chains)."""
    sigma_by_type = defaultdict(list)
    for chain_data in all_sigma_histories.values():
        for sigma_type, values in chain_data.items():
            sigma_by_type[sigma_type].append(np.asarray(values, dtype=float))

    rhat = {}
    for sigma_type, chains in sigma_by_type.items():
        chains = [ch[np.isfinite(ch)] for ch in chains if len(ch) > 10]
        if len(chains) < 2:
            rhat[sigma_type] = None
            continue
        min_len = min(len(ch) for ch in chains)
        if min_len < 10:
            rhat[sigma_type] = None
            continue
        truncated = np.array([ch[:min_len] for ch in chains])
        try:
            rhat[sigma_type] = float(az.rhat(truncated))
        except Exception as exc:
            print(f"[Rhat] {sigma_type}: failed ({exc})")
            rhat[sigma_type] = None
    return rhat


def add_rhat_page(pdf: PdfPages, rhat: Dict[str, float | None], sampler_label: str, out_dir: str, burnin_frac: float) -> None:
    """Write R-hat table to PDF and text file (Page 1)."""
    txt_path = os.path.join(out_dir, f"rhat_{sampler_label}.txt")
    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write(f"R-hat values for {sampler_label} (burn-in: {burnin_frac*100:.0f}%)\n")
        for k, val in rhat.items():
            line = f"{k}: {val:.3f}" if val is not None else f"{k}: N/A"
            fh.write(line + "\n")
    print(f"[Rhat] Summary written to {txt_path}")

    fig, ax = plt.subplots(figsize=(6, 2.5 + 0.4 * max(len(rhat), 1)))
    fig.suptitle(f"R-hat Convergence Diagnostics ({sampler_label})\nBurn-in: {burnin_frac*100:.0f}%", 
                 fontsize=14, fontweight='bold')
    ax.axis("off")

    rows = [[sigma_type, f"{val:.3f}" if val is not None else "N/A"] for sigma_type, val in rhat.items()]
    table = ax.table(cellText=rows, colLabels=["Sigma Type", "R-hat"], loc="center", cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.6)
    
    # Color code R-hat values
    for i, (_, val) in enumerate(rhat.items(), start=1):
        if val is not None:
            if val < 1.01:
                color = '#90EE90'  # Light green
            elif val < 1.05:
                color = '#FFFFE0'  # Light yellow
            else:
                color = '#FFB6C1'  # Light red
            table[(i, 1)].set_facecolor(color)

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def add_individual_chain_trace_page(
    pdf: PdfPages,
    chain_id: str,
    chain_data: Dict[str, List[float]],
    chain_steps: Dict[str, List[int]],
    sampler_label: str,
    burnin_step: int
) -> None:
    """Create one page with trace plots for all sigma types in a single chain."""
    sns.set(style="whitegrid")
    
    sigma_types = ["AA", "AB", "BC"]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 9))
    fig.suptitle(f"Trace Plots - Chain {chain_id} ({sampler_label})\nPost burn-in (step ≥ {burnin_step})", 
                 fontsize=14, fontweight='bold')
    
    for ax, sigma_type, color in zip(axes, sigma_types, colors):
        values = chain_data.get(sigma_type, [])
        steps = chain_steps.get(sigma_type, [])
        if values and steps:
            ax.plot(steps, values, color=color, linewidth=1.2, alpha=0.8)
            ax.axvline(x=burnin_step, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Burn-in cutoff')
            ax.set_ylabel(f"σ_{{{sigma_type}}}", fontsize=11, fontweight='bold')
            ax.set_xlabel("MCMC Step", fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_title(f"Sigma {sigma_type}", fontsize=11)
            ax.legend(loc='best', fontsize=9)
        else:
            ax.text(0.5, 0.5, f"No data for {sigma_type}", 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_ylabel(f"σ_{{{sigma_type}}}", fontsize=11)
    
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def add_combined_trace_page(
    pdf: PdfPages,
    sigma_type: str,
    all_chain_data: Dict[str, List[float]],
    all_chain_steps: Dict[str, List[int]],
    sampler_label: str,
    burnin_step: int
) -> None:
    """Create one page showing trace plots from all chains for a single sigma type."""
    sns.set(style="whitegrid")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle(f"Combined Trace Plot - Sigma {sigma_type} ({sampler_label})\nPost burn-in (step ≥ {burnin_step})", 
                 fontsize=14, fontweight='bold')
    
    palette = sns.color_palette("husl", len(all_chain_data))
    
    for (chain_id, values), steps, color in zip(
        sorted(all_chain_data.items()),
        [all_chain_steps[cid] for cid in sorted(all_chain_data.keys())],
        palette
    ):
        if values and steps:
            ax.plot(steps, values, color=color, linewidth=1.0, alpha=0.7, label=f"Chain {chain_id}")
    
    ax.axvline(x=burnin_step, color='red', linestyle='--', linewidth=2.0, alpha=0.8, label='Burn-in cutoff')
    ax.set_xlabel("MCMC Step", fontsize=12, fontweight='bold')
    ax.set_ylabel(f"σ_{{{sigma_type}}} Value", fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def add_combined_distribution_page(
    pdf: PdfPages,
    sigma_type: str,
    all_chain_data: Dict[str, np.ndarray],
    all_chain_models: Dict[str, GaussianMixture],
    sampler_label: str
) -> None:
    """Create one page showing distributions from all chains for a single sigma type."""
    sns.set(style="whitegrid")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle(f"Distribution Comparison - Sigma {sigma_type} ({sampler_label})", 
                 fontsize=14, fontweight='bold')
    
    palette = sns.color_palette("husl", len(all_chain_data))
    
    for (chain_id, data), color in zip(sorted(all_chain_data.items()), palette):
        # Plot histogram
        sns.histplot(
            data,
            bins=30,
            stat="density",
            alpha=0.25,
            label=f"Chain {chain_id}",
            color=color,
            ax=ax
        )
        
        # Overlay GMM curve if available
        gmm = all_chain_models.get(chain_id)
        if gmm is not None:
            x_range = np.linspace(data.min(), data.max(), 400).reshape(-1, 1)
            pdf_vals = np.exp(gmm.score_samples(x_range))
            ax.plot(x_range, pdf_vals, color=color, linewidth=2.5, alpha=0.9)
    
    ax.set_xlabel(f"σ_{{{sigma_type}}} Value", fontsize=12, fontweight='bold')
    ax.set_ylabel("Density", fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def analyze_mcmc_data(
    output_folder: str,
    sampler_type: str,
    sampler_position: int,
    burnin: float = 0.4,
    do_trace_plots: bool = True,
    do_gmm_fits: bool = True,
    do_diagnostics: bool = True,
) -> None:
    """
    Generate comprehensive PDF report with structured pages:
    - Page 1: R-hat table
    - Pages 2-N: Individual chain trace plots (one page per chain, showing step numbers)
    - Pages N+1, N+2, N+3: Combined trace plots for AA, AB, BC (all chains together)
    - Pages N+4, N+5, N+6: Combined distributions for AA, AB, BC
    """
    sampler_folder = get_sampler_folder_name(sampler_type, sampler_position)
    sampler_dir = os.path.join(output_folder, sampler_folder)
    if not os.path.isdir(sampler_dir):
        print(f"[analyze] Missing directory: {sampler_dir}")
        return

    trajectory_files = sorted(
        f for f in os.listdir(sampler_dir) if f.startswith("trajectory_") and f.endswith(".h5")
    )
    if not trajectory_files:
        print(f"[analyze] No trajectory files in {sampler_dir}")
        return

    # Load all sigma histories WITH step numbers
    sigma_histories = defaultdict(lambda: defaultdict(list))
    sigma_steps = defaultdict(lambda: defaultdict(list))
    
    for filename in trajectory_files:
        chain_id = filename.replace("trajectory_chain_", "").replace(".h5", "")
        states = load_trajectory_from_hdf5(os.path.join(sampler_dir, filename))
        for state in states:
            step = state["step"]
            for sigma_type, value in state["sigma"].items():
                sigma_histories[chain_id][sigma_type].append(float(value))
                sigma_steps[chain_id][sigma_type].append(int(step))

    if isinstance(burnin, str):
        print("[analyze] Automatic burn-in detection disabled; using 0.4")
        burnin = 0.4

    # Calculate burn-in step (use first chain's steps as reference)
    first_chain = list(sigma_histories.keys())[0]
    first_sigma = list(sigma_histories[first_chain].keys())[0]
    total_samples = len(sigma_histories[first_chain][first_sigma])
    burnin_idx = int(total_samples * burnin)
    burnin_step = sigma_steps[first_chain][first_sigma][burnin_idx] if burnin_idx < total_samples else 0
    
    print(f"[Burn-in] Removing first {burnin*100:.0f}% of samples (step < {burnin_step})")

    # Apply burn-in while keeping step information
    for chain_id in sigma_histories:
        for sigma_type in list(sigma_histories[chain_id].keys()):
            values, steps = apply_burnin_with_steps(
                sigma_histories[chain_id][sigma_type],
                sigma_steps[chain_id][sigma_type],
                float(burnin)
            )
            sigma_histories[chain_id][sigma_type] = values
            sigma_steps[chain_id][sigma_type] = steps

    sampler_label = f"{sampler_type.capitalize()}Sampler_{sampler_position}"
    pdf_path = os.path.join(sampler_dir, f"{sampler_label}_report.pdf")

    with PdfPages(pdf_path) as pdf:
        # PAGE 1: R-hat table
        rhat = compute_rhat_table(sigma_histories)
        add_rhat_page(pdf, rhat, sampler_label, sampler_dir, burnin)
        print(f"[Report] Page 1: R-hat table")

        # PAGES 2-N: Individual chain traces (one page per chain, with step numbers)
        for chain_id in sorted(sigma_histories.keys()):
            add_individual_chain_trace_page(
                pdf, chain_id, 
                sigma_histories[chain_id], 
                sigma_steps[chain_id],
                sampler_label,
                burnin_step
            )
            print(f"[Report] Added trace page for Chain {chain_id}")

        # PAGES N+1, N+2, N+3: Combined trace plots (all chains together)
        sigma_types = ["AA", "AB", "BC"]
        for sigma_type in sigma_types:
            chain_data = {}
            chain_step_data = {}
            for chain_id in sigma_histories.keys():
                if sigma_type in sigma_histories[chain_id]:
                    chain_data[chain_id] = sigma_histories[chain_id][sigma_type]
                    chain_step_data[chain_id] = sigma_steps[chain_id][sigma_type]
            
            if chain_data:
                add_combined_trace_page(pdf, sigma_type, chain_data, chain_step_data, sampler_label, burnin_step)
                print(f"[Report] Added combined trace plot for Sigma {sigma_type}")

        # Fit GMMs and prepare combined distribution data
        if do_gmm_fits:
            for sigma_type in sigma_types:
                chain_data = {}
                chain_models = {}
                
                for chain_id, chain_vals in sigma_histories.items():
                    values = np.asarray(chain_vals.get(sigma_type, []), dtype=float)
                    values = values[np.isfinite(values)]
                    if len(values) == 0:
                        continue
                    
                    chain_data[chain_id] = values
                    model = fit_gmm_simple(values, sigma_type)
                    chain_models[chain_id] = model
                    save_gmm_simple(model, sigma_type, chain_id, sampler_dir, values)
                
                # PAGES N+4, N+5, N+6: Combined distributions
                if chain_data:
                    add_combined_distribution_page(pdf, sigma_type, chain_data, chain_models, sampler_label)
                    print(f"[Report] Added combined distribution page for Sigma {sigma_type}")

    print(f"[Report] Complete PDF saved to {pdf_path}")


def analyze_sampler_in_sequence(
    sampler_sequence: List[str],
    sampler_type: str,
    sampler_position: int,
    output_folder: str = "output_analysis",
    burnin: float = 0.4,
    do_trace_plots: bool = True,
    do_gmm_fits: bool = True,
    do_diagnostics: bool = True,
) -> bool:
    """Thin wrapper preserving previous API."""
    if sampler_position >= len(sampler_sequence):
        print("[sequence] Position outside sequence")
        return False
    if sampler_sequence[sampler_position] != sampler_type:
        print("[sequence] Sampler type mismatch")
        return False

    occurrence = count_sampler_occurrences(sampler_sequence, sampler_type, sampler_position)
    os.makedirs(output_folder, exist_ok=True)
    analyze_mcmc_data(
        output_folder,
        sampler_type,
        occurrence,
        burnin,
        do_trace_plots,
        do_gmm_fits,
        do_diagnostics,
    )
    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze sigma trajectories: R-hat + GMM fits."
    )
    parser.add_argument(
        "--sequence",
        type=str,
        required=True,
        help='Comma-separated sampler sequence, e.g. "pair,tetramer,pair"',
    )
    parser.add_argument(
        "--sampler",
        type=str,
        required=True,
        choices=["pair", "tetramer", "octet", "full"],
        help="Sampler type to analyze",
    )
    parser.add_argument(
        "--position",
        type=int,
        required=True,
        help="0-based index of sampler occurrence in the sequence",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="output_analysis",
        help="Directory containing sampler outputs (default: output_analysis)",
    )
    parser.add_argument(
        "--burnin",
        "-b",
        type=str,
        default="0.4",
        help='Burn-in fraction (0–1).',
    )
    parser.add_argument(
        "--no-gmm",
        action="store_true",
        help="Skip GMM fitting (kept for compatibility)",
    )

    args = parser.parse_args()

    try:
        burnin = float(args.burnin)
        if not (0.0 <= burnin < 1.0):
            print("[cli] Burn-in outside [0,1); using 0.4")
            burnin = 0.4
    except ValueError:
        print("[cli] Invalid burn-in; using 0.4")
        burnin = 0.4

    sampler_sequence = [token.strip().lower() for token in args.sequence.split(",")]
    valid = {"pair", "tetramer", "octet", "full"}
    if any(token not in valid for token in sampler_sequence):
        print("[cli] Invalid sampler in sequence")
        sys.exit(1)

    ok = analyze_sampler_in_sequence(
        sampler_sequence=sampler_sequence,
        sampler_type=args.sampler,
        sampler_position=args.position,
        output_folder=args.output,
        burnin=burnin,
        do_trace_plots=True,
        do_gmm_fits=not args.no_gmm,
        do_diagnostics=True,
    )
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()