import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List

from scipy.spatial.transform import Rotation

# Use your OctetScorer instead of TetramerScorer
from parameters import SystemParameters  # Using the placeholder
from octet_score import OctetScorer     # Now referencing the OctetScorer class
from tetramer_score import TetramerScorer  # Assuming this is the original class

@dataclass
class PerturbSystemParameters:
    max_overlap_iterations: int = 2000
    overlap_correction_scale: float = 0.5
    log_high_score_details: bool = True
    debug_overlap_failure: bool = False

    ideal_coordinates: Dict[str, np.ndarray] = field(default_factory=dict)
    radii: Dict[str, float] = field(default_factory=dict)
    pair_distances: Dict = field(default_factory=dict)
    min_dists: Dict[Tuple[str, str], float] = field(default_factory=dict)
    fixed_sigma: Dict[str, float] = field(default_factory=lambda: {'AA': 1.5, 'AB': 0.9, 'BC': 0.8, 'CC': 1.0})

    # Replace TetramerScorer with OctetScorer
    tscorer: OctetScorer = field(default_factory=OctetScorer)
    tets_ideal: Dict = field(default_factory=dict)
    octets_ideal: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = field(default_factory=list)
    ideal_score: float = 0.0
    ideal_rmsd: float = 0.0
    ideal_accuracy: float = 1.0

    def __post_init__(self):
        params = SystemParameters()
        self.ideal_coordinates = params.ideal_coordinates
        self.radii = params.radii
        self.pair_distances = params.pair_distances
        self.ts = TetramerScorer()

        for comp1 in self.radii:
            for comp2 in self.radii:
                key = tuple(sorted((comp1, comp2)))
                self.min_dists[key] = self.radii[comp1] + self.radii[comp2]
        
        # Use the OctetScorer's default tetramer generator, then form octets
        #self.tets_ideal = self.tscorer._get_default_tetramers(self.ideal_coordinates)
        #self.tets_ideal = self.ts.get_tetramers(self.ideal_coordinates)
        #self.octets_ideal, _ = self.tscorer.get_octets(self.ideal_coordinates)

        score_tuple = self.tscorer.neg_log_posterior(
            positions=self.ideal_coordinates,
            sig=self.fixed_sigma,
            exclusion_weight=1.0,
            pair_weight=1.0,
            octet_weight=1.0
        )
        #self.tscorer.calculate_pairwise_distances(self.ideal_coordinates)
        self.ideal_score = score_tuple[0]
        print(f"Initialized Ideal Structure: Score = {self.ideal_score:.4f}")

    def _check_final_overlaps(self, current_coords: Dict[str, np.ndarray], print_details: bool = False) -> bool:
        # ...existing code...
        final_overlap_exists = False
        min_detected_dist_overall = np.inf
        overlapping_pair_info = ""

        for comp1 in current_coords:
            for comp2 in current_coords:
                coords1, coords2 = current_coords[comp1], current_coords[comp2]
                min_dist_allowed = self.min_dists[tuple(sorted((comp1, comp2)))]

                if comp1 == comp2:
                    if coords1.shape[0] < 2: continue
                    diff = coords1[:, np.newaxis, :] - coords1[np.newaxis, :, :]
                    distances = np.linalg.norm(diff, axis=2)
                    np.fill_diagonal(distances, np.inf)
                    if np.any(distances < min_dist_allowed):
                        final_overlap_exists = True
                        if print_details:
                            min_val_here = np.min(distances)
                            if min_val_here < min_detected_dist_overall:
                                min_detected_dist_overall = min_val_here
                                p_indices = np.unravel_index(np.argmin(distances), distances.shape)
                                overlapping_pair_info = f"Same comp {comp1}: p{p_indices[0]}-p{p_indices[1]}, dist {min_val_here:.3f} vs allowed {min_dist_allowed:.3f}"
                        if not print_details:
                            break
                elif comp1 < comp2:
                    coords1_expanded = coords1[:, np.newaxis, :]
                    coords2_expanded = coords2[np.newaxis, :, :]
                    distances = np.linalg.norm(coords1_expanded - coords2_expanded, axis=2)
                    if np.any(distances < min_dist_allowed):
                        final_overlap_exists = True
                        if print_details:
                            min_val_here = np.min(distances)
                            if min_val_here < min_detected_dist_overall:
                                min_detected_dist_overall = min_val_here
                                p_indices = np.unravel_index(np.argmin(distances), distances.shape)
                                overlapping_pair_info = f"Diff comps {comp1}-{comp2}: p{p_indices[0]}-p{p_indices[1]}, dist {min_val_here:.3f} vs allowed {min_dist_allowed:.3f}"
                        if not print_details:
                            break
            if final_overlap_exists and not print_details:
                break

        if final_overlap_exists and print_details:
            print(f"    Overlap Debug: Smallest distance violation found involving: {overlapping_pair_info}")
        return final_overlap_exists

    def resolve_overlaps_vectorized(self, coords: Dict[str, np.ndarray]) -> Tuple[Dict[str, np.ndarray], bool]:
        # ...existing code...
        current_coords = {comp: np.copy(arr) for comp, arr in coords.items()}

        for iteration in range(self.max_overlap_iterations):
            total_displacement = {comp: np.zeros_like(current_coords[comp]) for comp in current_coords}
            overlap_found_in_iteration = False

            for comp1 in current_coords:
                for comp2 in current_coords:
                    coords1, coords2 = current_coords[comp1], current_coords[comp2]
                    if comp1 == comp2:
                        if coords1.shape[0] < 2: continue
                        diff = coords1[:, np.newaxis, :] - coords1[np.newaxis, :, :]
                        distances = np.linalg.norm(diff, axis=2)
                        np.fill_diagonal(distances, np.inf)
                        min_dist_allowed = self.min_dists[tuple(sorted((comp1, comp1)))]
                        overlap_mask = distances < min_dist_allowed
                        overlapping_indices = np.where(overlap_mask)
                        if overlapping_indices[0].size > 0:
                            overlap_found_in_iteration = True
                            idx1, idx2 = overlapping_indices
                            valid_pairs = idx1 < idx2
                            idx1, idx2 = idx1[valid_pairs], idx2[valid_pairs]
                            if idx1.size > 0:
                                directions = coords1[idx1] - coords1[idx2]
                                norms = np.linalg.norm(directions, axis=1, keepdims=True)
                                directions /= (norms + 1e-10)
                                overlap_amounts = (min_dist_allowed - distances[idx1, idx2])
                                displacement_magnitudes = overlap_amounts * self.overlap_correction_scale / 2.0
                                displacements = directions * displacement_magnitudes[:, np.newaxis]
                                np.add.at(total_displacement[comp1], idx1, displacements)
                                np.add.at(total_displacement[comp1], idx2, -displacements)
                    elif comp1 < comp2:
                        coords1_expanded = coords1[:, np.newaxis, :]
                        coords2_expanded = coords2[np.newaxis, :, :]
                        distances = np.linalg.norm(coords1_expanded - coords2_expanded, axis=2)
                        min_dist_allowed = self.min_dists[tuple(sorted((comp1, comp2)))]
                        overlap_mask = distances < min_dist_allowed
                        overlapping_indices = np.where(overlap_mask)
                        if overlapping_indices[0].size > 0:
                            overlap_found_in_iteration = True
                            idx1, idx2 = overlapping_indices
                            directions = current_coords[comp1][idx1] - current_coords[comp2][idx2]
                            norms = np.linalg.norm(directions, axis=1, keepdims=True)
                            directions /= (norms + 1e-10)
                            overlap_amounts = (min_dist_allowed - distances[idx1, idx2])
                            displacement_magnitudes = overlap_amounts * self.overlap_correction_scale / 2.0
                            displacements = directions * displacement_magnitudes[:, np.newaxis]
                            np.add.at(total_displacement[comp1], idx1, displacements)
                            np.add.at(total_displacement[comp2], idx2, -displacements)

            if not overlap_found_in_iteration:
                return current_coords, True

            for comp in current_coords:
                current_coords[comp] += total_displacement[comp]

        # Final check
        if self._check_final_overlaps(current_coords, print_details=self.debug_overlap_failure):
            if self.debug_overlap_failure:
                print(f"    Warning: Overlaps persist after {self.max_overlap_iterations} iterations (vectorized).")
            return current_coords, False
        else:
            return current_coords, True

    def perturb_ideal_coordinates(
        self,
        magnitude: float = 1.0,
        seed: int = None,
        max_attempts: int = 20
    ) -> Optional[Dict[str, np.ndarray]]:
        # ...existing code...
        if seed is not None:
            np.random.seed(seed)

        num_mcmc_steps = 50
        step_scale = magnitude * 0.2
        for attempt in range(max_attempts):
            current_coords = {
                component: np.array(coords) for component, coords in self.ideal_coordinates.items()
            }
            for _ in range(num_mcmc_steps):
                for component in current_coords:
                    coords_array = current_coords[component]
                    idx_to_move = np.random.randint(0, coords_array.shape[0])
                    trial_coords = coords_array.copy()

                    move = np.random.normal(loc=0.0, scale=step_scale, size=(3,))
                    trial_coords[idx_to_move] += move

                    backup = coords_array[idx_to_move].copy()
                    coords_array[idx_to_move] = trial_coords[idx_to_move]

                    if self._check_final_overlaps(current_coords, print_details=False):
                        coords_array[idx_to_move] = backup

            if not self._check_final_overlaps(current_coords, print_details=False):
                return current_coords

        if self.debug_overlap_failure:
            print("Warning: MCMC approach yielded overlap after multiple attempts.")
        return None

    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
        # ...existing code...
        keys = sorted(self.ideal_coordinates.keys())
        try:
            ideal_coords_flat = np.concatenate([self.ideal_coordinates[k] for k in keys])
            perturbed_coords_flat = np.concatenate([perturbed_coords[k] for k in keys])
        except ValueError as e:
            print(f"Error concatenating coordinates for RMSD calculation: {e}")
            ideal_lengths = {k: len(self.ideal_coordinates[k]) for k in keys}
            perturbed_lengths = {k: len(perturbed_coords.get(k, [])) for k in keys}
            print(f"Ideal lengths: {ideal_lengths}, Perturbed lengths: {perturbed_lengths}")
            return np.inf

        if ideal_coords_flat.shape[0] == 0 or perturbed_coords_flat.shape[0] == 0:
            return np.inf
        if ideal_coords_flat.shape != perturbed_coords_flat.shape:
            raise ValueError("Mismatch in coordinate structure for RMSD calculation.")

        ideal_centroid = np.mean(ideal_coords_flat, axis=0)
        perturbed_centroid = np.mean(perturbed_coords_flat, axis=0)
        ideal_centered = ideal_coords_flat - ideal_centroid
        perturbed_centered = perturbed_coords_flat - perturbed_centroid
        try:
            rot, rmsd_val = Rotation.align_vectors(perturbed_centered, ideal_centered)
            return rmsd_val
        except Exception:
            return np.sqrt(np.mean(np.sum((ideal_centered - perturbed_centered) ** 2, axis=1)))

    def analyze_perturbations(self, magnitudes: List[float], num_perturbations: int = 10,
                              base_seed: Optional[int] = None,
                              score_threshold: Optional[float] = None,
                              exclusion_w: float = 1.0, pair_w: float = 1.0, octet_w: float = 1.0,
                              plot_means_only: bool = False,
                              add_error_bars_to_mean_plot: bool = True
                              ):
        # ...existing code...
        sns.set_theme(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=1.1)
        results = []

        best_coords_for_mag = {}
        best_score_for_mag = {}

        current_global_seed = base_seed
        total_skipped_overlap_failure = 0
        total_skipped_score_threshold = 0

        for mag_idx, mag in enumerate(magnitudes):
            num_successful_perturbations_for_mag = 0
            best_score_for_mag[mag] = float('inf')
            best_coords_for_mag[mag] = None

            print(f"\nProcessing Magnitude: {mag:.3f} ({mag_idx+1}/{len(magnitudes)})")
            for i in range(num_perturbations):
                current_perturb_seed = (current_global_seed + i) if base_seed is not None else None

                if self.debug_overlap_failure:
                    print(f"  Attempting Pert {i+1}/{num_perturbations}, Mag {mag:.3f}, Seed {current_perturb_seed}")

                perturbed_coords = self.perturb_ideal_coordinates(magnitude=mag, seed=current_perturb_seed)
                if perturbed_coords is None:
                    print(f"  Skipped Pert {i+1} (Mag {mag:.3f}, Seed {current_perturb_seed}): Overlap resolution failed.")
                    total_skipped_overlap_failure += 1
                    continue

                num_successful_perturbations_for_mag += 1
                # Default tetramers, then octets from those
                #tets_perturbed = self.tscorer._get_default_tetramers(perturbed_coords)
                #tets_perturbed = self.ts.get_tetramers(perturbed_coords)
                #octets_perturbed, _ = self.tscorer.get_octets(perturbed_coords)

                score_vals = self.tscorer.neg_log_posterior(
                    positions=perturbed_coords,
                    sig=self.fixed_sigma,
                    exclusion_weight=exclusion_w,
                    pair_weight=pair_w,
                    octet_weight=octet_w
                )
                score = score_vals[0]

                if score_threshold is not None and score > score_threshold:
                    print(f"  Skipped Pert {i+1} (Mag {mag:.3f}, Seed {current_perturb_seed}): Score {score:.2f} > threshold {score_threshold}.")
                    if self.log_high_score_details:
                        print(f"    Details: Excl={score_vals[1]:.2f}, Pair={score_vals[2]:.2f}, Octet={score_vals[3]:.2f}")
                    total_skipped_score_threshold += 1
                    continue

                rmsd = self.calculate_rmsd(perturbed_coords)
                accuracy = 1.0 / (1.0 + rmsd) if rmsd != np.inf else 0

                results.append({
                    'Magnitude': mag,
                    'Score': score,
                    'RMSD': rmsd,
                    'Accuracy': accuracy,
                    'ExclusionScore': score_vals[1],
                    'PairScore': score_vals[2],
                    'OctetScore': score_vals[3]
                })

                if score < best_score_for_mag[mag]:
                    best_score_for_mag[mag] = score
                    best_coords_for_mag[mag] = perturbed_coords

            if base_seed is not None:
                current_global_seed += num_perturbations

            print(f"  Completed {num_successful_perturbations_for_mag}/{num_perturbations} valid perturbations for magnitude {mag:.3f}.")

        if not results:
            print("No valid results generated. Check parameters or overlap resolution.")
            if total_skipped_overlap_failure > 0:
                print(f"Total configurations skipped due to overlap resolution failure: {total_skipped_overlap_failure}")
            if total_skipped_score_threshold > 0:
                print(f"Total configurations skipped due to exceeding score threshold: {total_skipped_score_threshold}")
            return None

        df_results = pd.DataFrame(results)

        print(f"\n--- Analysis Summary ---")
        if total_skipped_overlap_failure > 0:
            print(f"Total configurations skipped due to overlap resolution failure: {total_skipped_overlap_failure}")
        if total_skipped_score_threshold > 0:
            print(f"Total configurations skipped due to exceeding score threshold: {total_skipped_score_threshold}")
        total_skipped_configs = total_skipped_overlap_failure + total_skipped_score_threshold
        if total_skipped_configs > 0:
            print(f"Total configurations skipped overall: {total_skipped_configs}")

        summary_stats = df_results.groupby('Magnitude').agg(
            mean_score=('Score', 'mean'), std_score=('Score', 'std'),
            mean_rmsd=('RMSD', 'mean'), std_rmsd=('RMSD', 'std'),
            mean_accuracy=('Accuracy', 'mean'), std_accuracy=('Accuracy', 'std'),
            count=('Score', 'count')
        ).reset_index()

        print("\nSummary Statistics per Magnitude (based on non-skipped configurations):")
        print(summary_stats)

        fig_width = 12
        phi = (1 + np.sqrt(5)) / 2
        fig_height = fig_width / phi
        plt.figure(figsize=(fig_width, fig_height))

        plt.rcParams.update({
            'font.size': 14,
            'font.weight': 'bold',
            'axes.titlesize': 22,
            'axes.labelsize': 18,
            'axes.titleweight': 'bold',
            'axes.labelweight': 'bold',
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'legend.fontsize': 14,
            'legend.title_fontsize': 16
        })

        plot_title = 'RMSD vs. Score for Octet Sampler'
        plot_filename_suffix = 'detailed'

        if plot_means_only:
            plot_title = 'Mean RMSD vs. Mean Score (Error Bars: Std Dev)'
            plot_filename_suffix = 'means'
            if summary_stats.empty:
                print("No summary statistics to plot for means. Skipping plot.")
                return df_results

            plt.errorbar(summary_stats['mean_rmsd'], summary_stats['mean_score'],
                         xerr=summary_stats['std_rmsd'] if add_error_bars_to_mean_plot else None,
                         yerr=summary_stats['std_score'] if add_error_bars_to_mean_plot else None,
                         fmt='-o', capsize=6, elinewidth=2, markeredgewidth=2,
                         label='Mean Values (per Magnitude)', color='dodgerblue',
                         markerfacecolor='skyblue', markersize=12)
        else:
            plt.scatter(
                df_results['RMSD'], df_results['Score'],
                color='blue', s=30, alpha=0.7,
                label='Perturbed Structures'
            )

        plt.scatter(
            [self.ideal_rmsd], [self.ideal_score],
            color='red', marker='*', s=400, edgecolor='black', linewidth=2,
            label='Ground Truth', zorder=5
        )

        plt.xlabel(r'RMSD ($\AA$)', fontsize=20, fontweight='bold')
        plt.ylabel('Score (Negative Log Posterior)', fontsize=20, fontweight='bold')
        plt.title(plot_title, fontsize=22, fontweight='bold', pad=15)
        plt.grid(True, linestyle='-', alpha=0.3, linewidth=1.5)

        plt.tick_params(axis='both', which='major', labelsize=16, width=2, length=6)

        legend = plt.legend(fontsize=16, loc='best', framealpha=0.9, edgecolor='black')
        plt.tight_layout()

        output_dir = 'output_analysis/accuracy_vs_score'
        os.makedirs(output_dir, exist_ok=True)

        df_results.to_csv(os.path.join(output_dir, 'perturbation_analysis_results_raw.csv'), index=False)
        if plot_means_only and not summary_stats.empty:
            summary_stats.to_csv(os.path.join(output_dir, 'perturbation_analysis_summary_stats.csv'), index=False)

        plot_filename = os.path.join(output_dir, f'rmsd_vs_score_{plot_filename_suffix}_octet.png')
        plot_filename_pdf = os.path.join(output_dir, f'rmsd_vs_score_{plot_filename_suffix}_octet.svg')
        plt.savefig(plot_filename, dpi=300)
        plt.savefig(plot_filename_pdf)

        print(f"\nPlot saved to {plot_filename} and {plot_filename_pdf}")
        print(f"Raw results saved to {os.path.join(output_dir, 'perturbation_analysis_results_raw.csv')}")
        if plot_means_only and not summary_stats.empty:
            print(f"Summary statistics saved to {os.path.join(output_dir, 'perturbation_analysis_summary_stats.csv')}")
        plt.show()

        zoomed_df = df_results[df_results['RMSD'].between(0, 2.0)]
        if zoomed_df.empty:
            print("No results with RMSD between 0 and 2.0, skipping zoomed plot.")
        else:
            zoom_width = 8
            zoom_height = zoom_width / phi
            plt.figure(figsize=(zoom_width, zoom_height))
            plt.rcParams.update({
                'font.size': 16,
                'font.weight': 'bold',
                'axes.titlesize': 24,
                'axes.labelsize': 22,
                'axes.titleweight': 'bold',
                'axes.labelweight': 'bold',
                'xtick.labelsize': 18,
                'ytick.labelsize': 18,
                'legend.fontsize': 16,
                'legend.title_fontsize': 18
            })

            plt.scatter(
                zoomed_df['RMSD'], zoomed_df['Score'],
                color='blue', s=30, alpha=0.8,
                label='Perturbed Structures'
            )
            if 0 <= self.ideal_rmsd <= 2.0:
                plt.scatter(
                    [self.ideal_rmsd], [self.ideal_score],
                    color='red', marker='*', s=500, edgecolor='black', linewidth=2.5,
                    label='Ground Truth', zorder=10
                )

            plt.xlabel(r'RMSD ($\AA$)', fontsize=22, fontweight='bold')
            plt.ylabel('Score (Negative Log Posterior)', fontsize=22, fontweight='bold')
            plt.title("Zoomed RMSD ∈ [0, 2.0]", fontsize=24, fontweight='bold', pad=15)
            plt.grid(True, linestyle='-', alpha=0.3, linewidth=1.5)
            plt.tick_params(axis='both', which='major', labelsize=18, width=2.5, length=8)
            zoom_legend = plt.legend(fontsize=16, loc='best', framealpha=0.95, edgecolor='black', borderpad=1)
            plt.tight_layout()

            zoomed_plot_path_png = os.path.join(output_dir, "zoomed_octet_rmsd_vs_score.png")
            zoomed_plot_path_svg = os.path.join(output_dir, "zoomed_octet_rmsd_vs_score.svg")
            plt.savefig(zoomed_plot_path_png, dpi=300)
            plt.savefig(zoomed_plot_path_svg)
            print(f"\nZoomed plot saved to {zoomed_plot_path_png} and {zoomed_plot_path_svg}")
            plt.show()

        self.save_lowest_scoring_structures(df_results, best_coords_for_mag)
        return df_results

    def save_lowest_scoring_structures(self,
                                       df_results: pd.DataFrame,
                                       best_coords_for_mag: Dict[float, Dict[str, np.ndarray]],
                                       h5_filename: str = "lowest_scoring_structures.h5",
                                       csv_filename: str = "lowest_scoring_summary.csv"):
        # ...existing code...
        import h5py
        os.makedirs('output_analysis/lowest_structures', exist_ok=True)
        grouped = df_results.groupby('Magnitude', as_index=False)
        lowest_rows = grouped.apply(lambda g: g.loc[g['Score'].idxmin()]).reset_index(drop=True)

        full_h5_path = os.path.join('output_analysis/lowest_structures', h5_filename)
        with h5py.File(full_h5_path, 'w') as f:
            for row in lowest_rows.itertuples():
                mag = getattr(row, 'Magnitude')
                coords_dict = best_coords_for_mag.get(mag, None)
                if coords_dict is None:
                    continue
                group_name = f"magnitude_{mag:.3f}"
                grp = f.create_group(group_name)
                for comp_name, arr in coords_dict.items():
                    grp.create_dataset(comp_name, data=arr, compression="gzip", compression_opts=9)

        csv_path = os.path.join('output_analysis/lowest_structures', csv_filename)
        lowest_rows[['Magnitude', 'Score', 'RMSD']].to_csv(csv_path, index=False)

        print(f"\nSaved lowest-scoring structures per magnitude to HDF5: {full_h5_path}")
        print(f"Saved lowest-scoring summary CSV to: {csv_path}")


if __name__ == "__main__":
    magnitudes_to_test = np.concatenate([
        np.linspace(0.05, 0.5, 20),
        np.linspace(0.6, 2.5, 20),
        np.linspace(2.6, 10.0, 20)
    ])
    num_perturbations_per_mag = 50

    max_iters_overlap = 5000
    overlap_scale = 0.5
    debug_overlaps = False

    score_cap = 1000

    plot_only_mean_values = False
    show_error_bars_on_mean_plot = True

    analyzer = PerturbSystemParameters(
        max_overlap_iterations=max_iters_overlap,
        overlap_correction_scale=overlap_scale,
        log_high_score_details=True,
        debug_overlap_failure=debug_overlaps
    )

    print(f"Ideal Structure Score: {analyzer.ideal_score:.4f}")
    print(f"Overlap resolution: max_iterations={analyzer.max_overlap_iterations}, correction_scale={analyzer.overlap_correction_scale}")
    if score_cap is None:
        print("Score thresholding is DISABLED. All successfully resolved structures will be included.")
    else:
        print(f"Score threshold is ENABLED. Structures with score > {score_cap} will be skipped.")

    results_df = analyzer.analyze_perturbations(
        magnitudes=magnitudes_to_test,
        num_perturbations=num_perturbations_per_mag,
        base_seed=42,
        score_threshold=score_cap,
        exclusion_w=1.0,
        pair_w=1.0,
        octet_w=1.0,
        plot_means_only=plot_only_mean_values,
        add_error_bars_to_mean_plot=show_error_bars_on_mean_plot
    )

    if results_df is not None and not results_df.empty:
        print("\n--- First 5 rows of raw results ---")
        print(results_df.head())
    elif results_df is not None and results_df.empty:
        print("Analysis ran, but all configurations were skipped. Review parameters and logs.")
    else:
        print("Analysis did not produce any results DataFrame.")