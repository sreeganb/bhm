import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List

from scipy.spatial.transform import Rotation

# Assuming parameters.py and pair_score.py are in the same directory or accessible
from parameters import SystemParameters
from pair_score import ScoringSystem # Using ScoringSystem as per your first script
from visualization import plot_3d


@dataclass
class PerturbSystemParameters:
    # Configurable parameters from the newer version
    max_overlap_iterations: int = 2000
    overlap_correction_scale: float = 0.5 # This was 0.1 in the original displacement calc, now configurable
    log_high_score_details: bool = True
    debug_overlap_failure: bool = False

    # Fields to be initialized
    ideal_coordinates: Dict[str, np.ndarray] = field(default_factory=dict)
    radii: Dict[str, float] = field(default_factory=dict)
    pair_distances_from_params: Dict = field(default_factory=dict) # To store params.pair_distances
    min_dists: Dict[Tuple[str, str], float] = field(default_factory=dict)
    
    # Scoring related attributes from your first script's __init__
    scalc: Optional[ScoringSystem] = None # Will be initialized in __post_init__
    #fixed_sigma_for_scoring: Dict[str, float] = field(default_factory=lambda: {'AA': 0.2025, 'AB': 0.5, 'BC': 0.5, 'CC': 0.2822})
    fixed_sigma_for_scoring: Dict[str, float] = field(default_factory=lambda: {'AA': 3.0, 'AB': 2.0, 'BC': 1.0, 'CC': 1.5})
    sig_range_for_scoring: Dict[str, Tuple[float, float]] = field(default_factory=lambda: {"AA": (0.01, 10.0), "AB": (0.01, 10.0), "BC": (0.01, 10.0), "CC": (0.01, 10.0)})

    ideal_score: float = 0.0
    ideal_rmsd: float = 0.0
    ideal_accuracy: float = 1.0

    def __post_init__(self):
        params = SystemParameters()
        self.ideal_coordinates = params.ideal_coordinates
        self.radii = params.radii
        self.pair_distances_from_params = params.pair_distances # Store original pair distances

        # Initialize ScoringSystem (self.scalc) as in your original script
        self.scalc = ScoringSystem(system_params=params, exclusion_weight=1.0, pair_weight=1.0, prior_type='jeffreys')

        # Pre-calculate min distances for overlap resolution
        for comp1 in self.radii:
            for comp2 in self.radii:
                key = tuple(sorted((comp1, comp2)))
                self.min_dists[key] = self.radii[comp1] + self.radii[comp2]
        
        # Calculate ideal score using self.scalc and class-defined sigmas
        # The method name is calculate_score in your ScoringSystem
        score_tuple = self.scalc.calculate_score(
            self.ideal_coordinates, 
            self.fixed_sigma_for_scoring, # Use the class attribute
            self.sig_range_for_scoring    # Use the class attribute
        )
        self.ideal_score = score_tuple[0]
        self.ideal_accuracy = 1.0 # By definition
        self.ideal_rmsd = 0.0     # By definition
        print(f"Initialized Ideal Structure: Score = {self.ideal_score:.4f}")
        # plot the ideal structure
        #plot_3d(self.ideal_coordinates, "Ideal Structure", "output_analysis/ideal_structure")
        
        # Calculate the score for the half rotated system
        half_rotated_coords = params.half_rotated_coordinates()
        half_score_tuple = self.scalc.calculate_score(
            half_rotated_coords,
            self.fixed_sigma_for_scoring,
            self.sig_range_for_scoring
        )
        # plot the half rotated structure
        plot_3d(half_rotated_coords, "Half Rotated Structure", "output_analysis/half_rotated_structure")
        
        # score for a tetramer * 8 to get the full system score
        tet_coords = params.single_tetramer_coordinates()
        tet_score_tuple = self.scalc.calculate_score(
            tet_coords,
            self.fixed_sigma_for_scoring,
            self.sig_range_for_scoring
        )
        print(f"Tetramer Structure: Score = {tet_score_tuple[0]:.4f}")
        plot_3d(tet_coords, "Tetramer Structure", "output_analysis/tetramer_structure")
        
        oct_coords = params.octamer_coordinates()
        oct_score_tuple = self.scalc.calculate_score(
            oct_coords,
            self.fixed_sigma_for_scoring,
            self.sig_range_for_scoring
        )
        print(f"Octamer Structure: Score = {oct_score_tuple[0]:.4f}")
        plot_3d(oct_coords, "Octamer Structure", "output_analysis/octamer_structure")

        # write this score and ideal score to file
        with open('output_analysis/ideal_and_half_scores.txt', 'w') as f:
            f.write(f"Ideal Structure Score: {self.ideal_score:.4f}\n")
            f.write(f"Half Rotated Structure Score: {half_score_tuple[0]:.4f}\n")
            f.write(f"Half rotated structure excluded volume score: {half_score_tuple[1]:.4f}\n")
            f.write(f"Half rotated structure pair score: {half_score_tuple[2]:.4f}\n")
            f.write(f"Tetramer Structure Score: {tet_score_tuple[0]*8:.4f}\n")
            f.write(f"Tetramer Structure Excluded Volume Score: {tet_score_tuple[1]*8:.4f}\n")
            f.write(f"Tetramer Structure Pair Score: {tet_score_tuple[2]*8:.4f}\n")
            f.write(f"Octamer Structure Score: {oct_score_tuple[0]*4:.4f}\n")
            f.write(f"Octamer Structure Excluded Volume Score: {oct_score_tuple[1]*4:.4f}\n")
            f.write(f"Octamer Structure Pair Score: {oct_score_tuple[2]*4:.4f}\n")
        print(f"Half Rotated Structure: Score = {half_score_tuple[0]:.4f}")

    def _check_final_overlaps(self, current_coords: Dict[str, np.ndarray], print_details: bool = False) -> bool:
        """Helper to check for overlaps and optionally print details."""
        final_overlap_exists = False
        min_detected_dist_overall = np.inf
        overlapping_pair_info = ""

        for comp1 in current_coords:
            for comp2 in current_coords:
                coords1, coords2 = current_coords[comp1], current_coords[comp2]
                min_dist_allowed = self.min_dists[tuple(sorted((comp1, comp2)))]

                if comp1 == comp2:
                    if coords1.shape[0] < 2:
                        continue
                    diff = coords1[:, np.newaxis, :] - coords1[np.newaxis, :, :]
                    distances = np.linalg.norm(diff, axis=2)
                    np.fill_diagonal(distances, np.inf)
                    if np.any(distances < min_dist_allowed):
                        final_overlap_exists = True
                        if print_details:
                            min_val_here = np.min(distances[distances < min_dist_allowed]) \
                                if np.any(distances < min_dist_allowed) else np.inf
                            if min_val_here < min_detected_dist_overall:
                                min_detected_dist_overall = min_val_here
                                p_indices = np.unravel_index(np.argmin(distances), distances.shape)
                                overlapping_pair_info = (
                                    f"Same comp {comp1}: p{p_indices[0]}-p{p_indices[1]}, "
                                    f"dist {min_val_here:.3f} vs allowed {min_dist_allowed:.3f}"
                                )
                        if not print_details:
                            break
                elif comp1 < comp2:
                    coords1_expanded = coords1[:, np.newaxis, :]
                    coords2_expanded = coords2[np.newaxis, :, :]
                    distances = np.linalg.norm(coords1_expanded - coords2_expanded, axis=2)
                    if np.any(distances < min_dist_allowed):
                        final_overlap_exists = True
                        if print_details:
                            min_val_here = np.min(distances[distances < min_dist_allowed]) \
                                if np.any(distances < min_dist_allowed) else np.inf
                            if min_val_here < min_detected_dist_overall:
                                min_detected_dist_overall = min_val_here
                                p_indices = np.unravel_index(np.argmin(distances), distances.shape)
                                overlapping_pair_info = (
                                    f"Diff comps {comp1}-{comp2}: p{p_indices[0]}-p{p_indices[1]}, "
                                    f"dist {min_val_here:.3f} vs allowed {min_dist_allowed:.3f}"
                                )
                        if not print_details:
                            break
            if final_overlap_exists and not print_details:
                break

        if final_overlap_exists and print_details:
            print(f"    Overlap Debug: Smallest distance violation found involving: {overlapping_pair_info}")
        return final_overlap_exists

    def perturb_ideal_coordinates(
        self,
        magnitude: float = 1.0,
        seed: int = None,
        max_attempts: int = 20
    ) -> Optional[Dict[str, np.ndarray]]:
        """
        MCMC-style perturbation: for each attempt, we perform a short MCMC chain
        moving each particle randomly and accepting if no overlap is introduced.
        This avoids large collisions often seen in purely random jumps.
        """
        if seed is not None:
            np.random.seed(seed)

        num_mcmc_steps = 50  # Increase to refine final structure
        step_scale = magnitude * 0.2  # Smaller scale per step
        for attempt in range(max_attempts):
            current_coords = {
                component: np.array(coords) for component, coords in self.ideal_coordinates.items()
            }
            for _ in range(num_mcmc_steps):
                for component in current_coords:
                    coords_array = current_coords[component]
                    idx_to_move = np.random.randint(0, coords_array.shape[0])
                    trial_coords = coords_array.copy()

                    # Change size=(3,) instead of (1,3)
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

# ...remaining code below...
    def calculate_rmsd(self, perturbed_coords: Dict[str, np.ndarray]) -> float:
        keys = sorted(self.ideal_coordinates.keys())
        
        # Ensure all components expected in ideal_coordinates are present in perturbed_coords
        # and have the same number of atoms. This is critical.
        for k in keys:
            if k not in perturbed_coords or self.ideal_coordinates[k].shape != perturbed_coords[k].shape:
                raise ValueError(f"Component {k} shape mismatch or missing in perturbed_coords for RMSD.")

        ideal_coords_flat = np.concatenate([self.ideal_coordinates[k] for k in keys])
        perturbed_coords_flat = np.concatenate([perturbed_coords[k] for k in keys])

        if ideal_coords_flat.shape[0] == 0: # Handle empty structures
             return 0.0 if perturbed_coords_flat.shape[0] == 0 else np.inf


        ideal_centroid = np.mean(ideal_coords_flat, axis=0)
        perturbed_centroid = np.mean(perturbed_coords_flat, axis=0)
        ideal_centered = ideal_coords_flat - ideal_centroid
        perturbed_centered = perturbed_coords_flat - perturbed_centroid
        
        try:
            # Rotation.align_vectors returns rotation and RMSD
            # The old script calculated RMSD manually after alignment, this is more direct
            rot, rmsd_val = Rotation.align_vectors(perturbed_centered, ideal_centered)
            return rmsd_val 
        except Exception as e: # Catch potential errors in alignment (e.g. too few points, collinear points)
            print(f"Error during Rotation.align_vectors: {e}. Calculating non-aligned RMSD as fallback.")
            return np.sqrt(np.mean(np.sum((ideal_centered - perturbed_centered) ** 2, axis=1)))


    def analyze_perturbations(self, magnitudes: List['float'], num_perturbations: int = 10,
                            base_seed: Optional['int'] = None,
                            score_threshold: Optional['float'] = None,
                            exclusion_w_scorer: float = 1.0,
                            pair_w_scorer: float = 1.0,
                            plot_means_only: bool = False,
                            add_error_bars_to_mean_plot: bool = True,
                            max_rmsd_to_plot: float = 15.0  # <-- New parameter to limit RMSD in plots
                            ):
        sns.set_theme(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=1.1)
        results = []

        best_coords_for_mag = {}
        best_score_for_mag = {}

        current_global_seed = base_seed
        total_skipped_overlap_failure = 0
        total_skipped_score_threshold = 0

        for mag_idx, mag in enumerate(magnitudes):
            num_successful_perturbations_for_mag = 0
            best_score_for_mag[mag] = float('inf')  # Fixed: use mag as key, not 'mag' string
            best_coords_for_mag[mag] = None         # Fixed: use mag as key, not 'mag' string

            print(f"\nProcessing Magnitude: {mag:.3f} ({mag_idx+1}/{len(magnitudes)})")
            for i in range(num_perturbations):
                current_perturb_seed = (current_global_seed + i) if base_seed is not None else None

                if self.debug_overlap_failure:
                    print(f"  Attempting Pert {i+1}/{num_perturbations}, Mag {mag:.3f}, Seed {current_perturb_seed}")

                perturbed_coords = self.perturb_ideal_coordinates(magnitude=mag, seed=current_perturb_seed)

                if perturbed_coords is None:
                    if not self.debug_overlap_failure:
                        print(f"  Skipped Pert {i+1} (Mag {mag:.3f}, Seed {current_perturb_seed}): Overlap resolution failed.")
                    else:
                        print(f"  Skipped Pert {i+1} (Mag {mag:.3f}, Seed {current_perturb_seed}) due to previously noted overlap resolution failure.")
                    total_skipped_overlap_failure += 1
                    continue

                num_successful_perturbations_for_mag += 1

                # Calculate score using the perturbed coordinates
                score_tuple = self.scalc.calculate_score(
                    perturbed_coords,
                    self.fixed_sigma_for_scoring,
                    self.sig_range_for_scoring
                )
                score = score_tuple[0]
                
                # If score_tuple has additional components, access them by position
                exclusion_val = score_tuple[1] if len(score_tuple) > 1 else np.nan
                pair_val = score_tuple[2] if len(score_tuple) > 2 else np.nan
                other_val = score_tuple[3] if len(score_tuple) > 3 else np.nan

                if score_threshold is not None and score > score_threshold:
                    print(f"  Skipped Pert {i+1} (Mag {mag:.3f}, Seed {current_perturb_seed}): Score {score:.2f} > threshold {score_threshold}.")
                    if self.log_high_score_details:
                        print(f"    Details: Excl={exclusion_val:.2f}, Pair={pair_val:.2f}, Other={other_val:.2f}")
                    total_skipped_score_threshold += 1
                    continue

                rmsd = self.calculate_rmsd(perturbed_coords)
                accuracy = 1.0 / (1.0 + rmsd) if rmsd != np.inf else 0.0

                results.append({
                    'Magnitude': mag, 'Score': score, 'RMSD': rmsd, 'Accuracy': accuracy,
                    'ExclusionScore': exclusion_val, 'PairScore': pair_val, 'OtherScore': other_val
                })

                # Update best structure if this is the lowest score so far for this magnitude
                if score < best_score_for_mag[mag]:  # Fixed: use mag as key, not 'mag' string
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

        # Filter to plot only RMSD <= max_rmsd_to_plot
        df_results_for_plot = df_results[df_results['RMSD'] <= max_rmsd_to_plot]

        summary_stats = df_results_for_plot.groupby('Magnitude').agg(
            mean_score=('Score', 'mean'), std_score=('Score', 'std'),
            mean_rmsd=('RMSD', 'mean'), std_rmsd=('RMSD', 'std'),
            mean_accuracy=('Accuracy', 'mean'), std_accuracy=('Accuracy', 'std'),
            count=('Score', 'count')
        ).reset_index()

        print("\nSummary Statistics per Magnitude (based on non-skipped configurations):")
        print(summary_stats)

        # Square dimensions for main plot
        fig_width = 10
        fig_height = 10
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

        # Removed plot title
        plot_filename_suffix = 'detailed'

        if plot_means_only:
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
                df_results_for_plot['RMSD'], df_results_for_plot['Score'],
                color='blue', s=30, alpha=0.7,
                label='Perturbed Structures'
            )

        # Ground truth star
        if self.ideal_rmsd <= max_rmsd_to_plot:
            plt.scatter(
                [self.ideal_rmsd], [self.ideal_score],
                color='red', marker='*', s=400, edgecolor='black', linewidth=2,
                label='Ground Truth', zorder=5
            )

        plt.xlabel('RMSD [$\\AA$]', fontsize=20, fontweight='bold')
        plt.ylabel('Score [Negative Log Posterior]', fontsize=20, fontweight='bold')
        plt.grid(True, linestyle='-', alpha=0.3, linewidth=1.5)

        plt.tick_params(axis='both', which='major', labelsize=16, width=2, length=6)

        legend = plt.legend(fontsize=16, loc='best', framealpha=0.9, edgecolor='black')

        plt.tight_layout()

        output_dir = 'output_analysis/accuracy_vs_score'
        os.makedirs(output_dir, exist_ok=True)

        df_results.to_csv(os.path.join(output_dir, 'perturbation_analysis_results_raw.csv'), index=False)
        if plot_means_only and not summary_stats.empty:
            summary_stats.to_csv(os.path.join(output_dir, 'perturbation_analysis_summary_stats.csv'), index=False)

        plot_filename = os.path.join(output_dir, f'rmsd_vs_score_{plot_filename_suffix}_pair.png')
        plot_filename_pdf = os.path.join(output_dir, f'rmsd_vs_score_{plot_filename_suffix}_pair.svg')

        plt.savefig(plot_filename, dpi=300)
        plt.savefig(plot_filename_pdf)
        print(f"\nPlot saved to {plot_filename} and {plot_filename_pdf}")
        print(f"Raw results saved to {os.path.join(output_dir, 'perturbation_analysis_results_raw.csv')}")
        if plot_means_only and not summary_stats.empty:
            print(f"Summary statistics saved to {os.path.join(output_dir, 'perturbation_analysis_summary_stats.csv')}")

        plt.show()

        # Create a second "zoomed in" plot for RMSD in [0, 2.0]
        zoomed_df = df_results_for_plot[df_results_for_plot['RMSD'].between(0, 2.0)]
        if zoomed_df.empty:
            print("No results with RMSD between 0 and 2.0, skipping zoomed plot.")
        else:
            # Square dimensions for zoomed plot
            zoom_width = 8
            zoom_height = 8
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

            # Removed plot title
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

            plt.xlabel('RMSD [$\\AA$]', fontsize=22, fontweight='bold')
            plt.ylabel('Score [Negative Log Posterior]', fontsize=22, fontweight='bold')
            plt.grid(True, linestyle='-', alpha=0.3, linewidth=1.5)

            plt.tick_params(axis='both', which='major', labelsize=18, width=2.5, length=8)

            zoom_legend = plt.legend(fontsize=16, loc='best', framealpha=0.95, edgecolor='black', borderpad=1)

            plt.tight_layout()

            zoomed_plot_path_png = os.path.join(output_dir, "zoomed_pair_rmsd_vs_score.png")
            zoomed_plot_path_svg = os.path.join(output_dir, "zoomed_pair_rmsd_vs_score.svg")
            plt.savefig(zoomed_plot_path_png, dpi=300)
            plt.savefig(zoomed_plot_path_svg)
            print(f"\nZoomed plot saved to {zoomed_plot_path_png} and {zoomed_plot_path_svg}")
            plt.show()

        self.save_lowest_scoring_structures(df_results, best_coords_for_mag)
        return df_results

    def save_lowest_scoring_structures(self,
                                    df_results: pd.DataFrame,
                                    best_coords_for_mag: Dict['float', Dict['str', np.ndarray]],
                                    h5_filename: str = "lowest_scoring_structures_pair.h5",
                                    csv_filename: str = "lowest_scoring_summary_pair.csv"):
        """
        Save the lowest-score structure per magnitude to HDF5,
        and also output a CSV summarizing Magnitude, Score, RMSD.
        Each magnitude-group's minimal entry from df_results is used.
        """
        import h5py

        # 1) Find the row with the lowest Score for each magnitude in df_results
        # Use idxmin() on the Score column directly - much simpler approach
        idx_min = df_results.groupby('Magnitude')['Score'].idxmin()
        lowest_rows = df_results.loc[idx_min]

        # 2) Write each best coordinate set to one HDF5 file, keyed by "mag_XXXX"
        os.makedirs('output_analysis/lowest_structures', exist_ok=True)
        full_h5_path = os.path.join('output_analysis/lowest_structures', h5_filename)

        with h5py.File(full_h5_path, 'w') as f:
            for row in lowest_rows.itertuples():
                mag = getattr(row, 'Magnitude')
                coords_dict = best_coords_for_mag.get(mag, None)
                if coords_dict is None:
                    continue
                group_name = f"magnitude_{mag:.3f}"
                grp = f.create_group(group_name)
                # Save each component as a dataset
                for comp_name, arr in coords_dict.items():
                    grp.create_dataset(comp_name, data=arr, compression="gzip", compression_opts=9)

        # 3) Save a CSV with the minimal Score entries
        #    containing Magnitude, Score, RMSD, etc.
        csv_path = os.path.join('output_analysis/lowest_structures', csv_filename)
        lowest_rows[['Magnitude', 'Score', 'RMSD']].to_csv(csv_path, index=False)

        print(f"\nSaved lowest-scoring structures per magnitude to HDF5: {full_h5_path}")
        print(f"Saved lowest-scoring summary CSV to: {csv_path}")


if __name__ == "__main__":
    # --- Configuration for the analysis (similar to newer script) ---
    magnitudes_to_test = np.concatenate([
        np.linspace(0.05, 0.5, 30), # Start with smaller magnitudes
        np.linspace(0.6, 2.5, 30),
        np.linspace(2.6, 300.0, 30)  # Max magnitude from your old script's linspace
    ])
    # num_perturbations_per_mag from your old script was 20, new script example was 50. Using 20.
    num_perturbations_per_mag = 30 

    # --- Parameters for Overlap Resolution and Scoring ---
    max_iters_overlap = 5000    # Increased from original 1000, tunable
    # Original displacement factor was effectively 0.1 * 0.5 = 0.05.
    # Here, overlap_correction_scale is directly the factor (e.g., 0.5 means move half the overlap distance)
    # A value of 0.1 here would match the old script's displacement magnitude if self.overlap_correction_scale=0.1
    # For the overlap_amounts / 2 * 0.1 line, this means overlap_correction_scale = 0.1.
    # Let's use a slightly larger default for better resolution.
    overlap_scale = 0.2        # Tunable. 0.1 would be closer to old script's fixed factor.
    debug_overlaps = False     # Set to True for detailed overlap failure messages

    # Score Threshold: Set to None to disable filtering by score
    score_cap = None # Example: 75000, or None to disable

    # --- Plotting Choice ---
    plot_only_mean_values = False # Set to True to plot only mean values 
    show_error_bars_on_mean_plot = True


    analyzer = PerturbSystemParameters(
        max_overlap_iterations=max_iters_overlap,
        overlap_correction_scale=overlap_scale,
        log_high_score_details=True,
        debug_overlap_failure=debug_overlaps
    )

    print(f"Overlap resolution: max_iterations={analyzer.max_overlap_iterations}, correction_scale={analyzer.overlap_correction_scale}")
    if score_cap is None:
        print("Score thresholding is DISABLED. All successfully resolved structures will be included.")
    else:
        print(f"Score threshold is ENABLED. Structures with score > {score_cap} will be skipped.")

    # Call analyze_perturbations. Note: exclusion_w_scorer and pair_w_scorer are not used
    # unless you modify ScoringSystem to accept these dynamically or re-instantiate.
    # The weights are currently set during ScoringSystem initialization.
    results_df = analyzer.analyze_perturbations(
        magnitudes=magnitudes_to_test,
        num_perturbations=num_perturbations_per_mag,
        base_seed=52, # From your original script
        score_threshold=score_cap,
        # These weights would be for the scorer if it's designed to take them here.
        # exclusion_w_scorer=1.0, 
        # pair_w_scorer=1.0,
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
