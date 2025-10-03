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


def apply_burnin(series: List[float], fraction: float) -> List[float]:
    if not series:
        return series
    cut = int(len(series) * fraction)
    return series[cut:]


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


def add_rhat_page(pdf: PdfPages, rhat: Dict[str, float | None], sampler_label: str, out_dir: str) -> None:
    """Write R-hat table to PDF and text file."""
    txt_path = os.path.join(out_dir, f"rhat_{sampler_label}.txt")
    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write(f"R-hat values for {sampler_label}\n")
        for k, val in rhat.items():
            line = f"{k}: {val:.3f}" if val is not None else f"{k}: N/A"
            fh.write(line + "\n")
    print(f"[Rhat] Summary written to {txt_path}")

    fig, ax = plt.subplots(figsize=(5, 1.5 + 0.3 * max(len(rhat), 1)))
    fig.suptitle(f"R-hat ({sampler_label})", fontsize=14)
    ax.axis("off")

    rows = [[sigma_type, f"{val:.3f}" if val is not None else "N/A"] for sigma_type, val in rhat.items()]
    table = ax.table(cellText=rows, colLabels=["Sigma Type", "R-hat"], loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.4)

    pdf.savefig(fig)
    plt.close(fig)


def plot_combined_gmm(
    pdf: PdfPages,
    sigma_type: str,
    sampler_label: str,
    chain_data: Dict[str, np.ndarray],
    chain_models: Dict[str, GaussianMixture],
    output_dir: str,
) -> None:
    """Histogram + mixture curves for a sigma component."""
    plt.figure(figsize=(10, 7))
    sns.set(style="darkgrid")
    palette = sns.color_palette("husl", len(chain_data))

    for (chain_id, data), color in zip(chain_data.items(), palette):
        sns.histplot(
            data,
            bins=30,
            stat="density",
            alpha=0.35,
            label=f"Chain {chain_id}",
            color=color,
        )

        gmm = chain_models.get(chain_id)
        if gmm is not None:
            x = np.linspace(data.min(), data.max(), 400).reshape(-1, 1)
            pdf_vals = np.exp(gmm.score_samples(x))
            plt.plot(x, pdf_vals, color=color, lw=2)

    plt.title(f"{sigma_type} — {sampler_label}")
    plt.xlabel("Sigma")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()

    pdf.savefig()
    out_png = os.path.join(output_dir, f"gmm_{sigma_type}_{sampler_label}.png")
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[Plot] {sigma_type} saved to {out_png}")

def plot_sigma_traces(
    pdf: PdfPages,
    sigma_histories: Dict[str, Dict[str, List[float]]],
    sampler_label: str,
) -> None:
    """Trace plots for all sigma types."""
    sns.set(style="darkgrid")
    all_sigma_types = set()
    for chain_data in sigma_histories.values():
        all_sigma_types.update(chain_data.keys())
    all_sigma_types = sorted(all_sigma_types)

    n_types = len(all_sigma_types)
    n_cols = 2
    n_rows = (n_types + 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 3 * n_rows), squeeze=False)
    fig.suptitle(f"Sigma Traces ({sampler_label})", fontsize=16)

    for idx, sigma_type in enumerate(all_sigma_types):
        ax = axes[idx // n_cols][idx % n_cols]
        for chain_id, chain_data in sigma_histories.items():
            values = chain_data.get(sigma_type, [])
            if not values:
                continue
            ax.plot(values, label=f"Chain {chain_id}", alpha=0.7)
        ax.set_title(sigma_type)
        ax.set_xlabel("Sample Index")
        ax.set_ylabel("Sigma Value")
        ax.legend()

    for j in range(idx + 1, n_rows * n_cols):
        fig.delaxes(axes[j // n_cols][j % n_cols])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
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
    """Collect sigma trajectories, compute R-hat, fit GMMs, and save plots."""
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

    sigma_histories = defaultdict(lambda: defaultdict(list))
    for filename in trajectory_files:
        chain_id = filename.removeprefix("trajectory_").removesuffix(".h5")
        states = load_trajectory_from_hdf5(os.path.join(sampler_dir, filename))
        for state in states:
            for sigma_type, value in state["sigma"].items():
                sigma_histories[chain_id][sigma_type].append(float(value))

    if isinstance(burnin, str):
        print("[analyze] Automatic burn-in detection disabled; using 0.4")
        burnin = 0.4

    for chain_id in sigma_histories:
        for sigma_type in list(sigma_histories[chain_id].keys()):
            sigma_histories[chain_id][sigma_type] = apply_burnin(
                sigma_histories[chain_id][sigma_type], float(burnin)
            )

    sampler_label = f"{sampler_type.capitalize()}Sampler_{sampler_position}"
    pdf_path = os.path.join(sampler_dir, f"{sampler_label}_summary.pdf")

    with PdfPages(pdf_path) as pdf:
        rhat = compute_rhat_table(sigma_histories)
        add_rhat_page(pdf, rhat, sampler_label, sampler_dir)

        if do_gmm_fits:
            for sigma_type in ["AA", "AB", "BC"]:
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

                if chain_data:
                    plot_combined_gmm(pdf, sigma_type, sampler_label, chain_data, chain_models, sampler_dir)
                    plot_sigma_traces(pdf, sigma_histories, sampler_label)

    print(f"[analyze] Report saved to {pdf_path}")


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
        help='Burn-in fraction (0–1). "auto" falls back to 0.4.',
    )
    parser.add_argument(
        "--no-gmm",
        action="store_true",
        help="Skip GMM fitting (kept for compatibility)",
    )

    args = parser.parse_args()

    if args.burnin.lower() == "auto":
        burnin = "auto"
    else:
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
        do_trace_plots=False,
        do_gmm_fits=not args.no_gmm,
        do_diagnostics=False,
    )
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()