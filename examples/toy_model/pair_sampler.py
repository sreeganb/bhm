#----------------------------------------------------------------------
# PairSampler class which will use the base sampler functions 
# and run the MCMC simulations for the pairwise sampling and 
# track sigma values
#----------------------------------------------------------------------
from typing import List, Dict, Tuple
import os
import numpy as np
import pandas as pd
from base_sampler import BaseMCSampler, Priors  # Import BaseMCSampler and related classes
from parameters import SystemParameters
from visualization import visualize_3d_configuration
from scipy.spatial.distance import cdist
#----------------------------------------------------------------------

class PairSampler(BaseMCSampler):
    """
    Sampler for pair-level interactions, inheriting from BaseMCSampler.
    Includes run_mc method.
    """
    def __init__(self, use_sigma_distribution=False,
                 ex_weight: float = 1.0, pair_weight: float = 1.0, 
                 use_def_sig_pos : bool = True, 
                 sig_passed: Dict[str, float] = None, 
                 sig_range_passed: Dict[str, Tuple[float, float]] = None,
                 pos_passed: Dict[str, np.ndarray] = None):
        super().__init__()  # Call BaseMCSampler constructor
        self.use_sigma_distribution = use_sigma_distribution
        self.exclusion_weight = ex_weight
        self.pair_weight = pair_weight
        self.params = SystemParameters()  # Initialize system parameters
        if not self.use_sigma_distribution:
            self.priors = Priors("jeffreys")  # Initialize priors

        if use_def_sig_pos:
            # Initialize sigma, sigma_range, and priors.
            self.sigma, self.sigma_range = self.initialize_sigma()
            self.positions_ps = self.initialize_positions()  # Initialize positions
            visualize_3d_configuration(self.positions_ps, self.params.radii)  # Visualize initial configuration
            print("PairSampler initialized with default sigma values.")
        else:
            self.sigma = sig_passed
            self.sigma_range = sig_range_passed
            self.positions_ps = pos_passed
            print("PairSampler initialized with passed sigma values.")        

    def calculate_score(
        self,
        pos: Dict[str, np.ndarray],
        sig: Dict[str, float],
        sig_range: Dict[str, Tuple[float, float]] = None,
        excluded_pairs=None,
        use_sigma_distribution=False,
        prior_penalty_from_distribution=0.0,
        step: int = 0
    ) -> Tuple[float, float, float, float]:
        """Calculate the log posterior for the pair-level interactions."""
        # 1) Excluded volume contribution
        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos, step, 
                                                                           log_file="excluded_vol_log.csv", 
                                                                           debug = True)
        
        # 2) Pairwise negative log-likelihood
        pairwise_score = 0.0
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
        
        # Collect logged pairs for printing below
        chosen_pairs_log = []
        
        for type1, type2 in pair_types:
            pair_key = f"{type1}{type2}"
            if pair_key in self.params.pair_distances:
                target_dist = self.params.pair_distances[pair_key]
                sigma_value = sig[pair_key]
                
                # Calculate pairwise score matrix
                score_matrix = self.calculate_pair_scores_matrix(
                    pos[type1], pos[type2], target_dist, sigma_value
                )
                
                # Mask diagonal (for same-type pairs)
                if type1 == type2:
                    np.fill_diagonal(score_matrix, np.inf)  # Use inf for clearer intent
                
                # Exclude certain pairs - more efficient implementation
                if excluded_pairs:
                    relevant_pairs = [(p[1], p[3]) for p in excluded_pairs 
                                    if p[0] == type1 and p[2] == type2]
                    for i, j in relevant_pairs:
                        if i < len(pos[type1]) and j < len(pos[type2]):
                            score_matrix[i, j] = np.inf
                
                # Identify minimal row/column pairs
                row_indices = np.argmin(score_matrix, axis=1)
                col_indices = np.argmin(score_matrix, axis=0)
                row_pairs = {(i, row_indices[i]) for i in range(len(row_indices))}
                col_pairs = {(col_indices[j], j) for j in range(len(col_indices))}
                unique_pairs = row_pairs.union(col_pairs)
                
                for i, j in unique_pairs:
                    pairwise_score += self.pair_weight * score_matrix[i, j]
                    # Log the chosen pair (e.g. (A1, B2))
                    chosen_pairs_log.append(f"({type1}{i}, {type2}{j})")
        
        # 3) Prior penalty
        if not use_sigma_distribution:
            prior_penalty = self.priors.neg_log_prior(sig, sig_range)
        else:
            prior_penalty = prior_penalty_from_distribution
        
        total_score = exclusion_score + pairwise_score + prior_penalty
        
        # Avoid opening/closing the file in every call by making this optional
        # or passing a file handle instead
        if hasattr(self, 'pairs_log_file') and self.pairs_log_file:
            with open(self.pairs_log_file, "a") as log_file:
                log_file.write("Chosen pairs:\n")
                for cpair in chosen_pairs_log:
                    log_file.write(f"{cpair}\n")
                log_file.write("\n")  # Blank line between scoring steps
        
        return total_score, exclusion_score, pairwise_score, prior_penalty

    def run_mc(self, n_steps: int = 50000, save_freq: int = 100,
                output_dir: str = "output_analysis/pairsampler_results/",
                position_move_prob: float = 0.9) -> Tuple[Dict[str, np.ndarray], List[Dict], str]:
            best_positions = None
            final_positions = None  # Also track final positions
            trajectory = []
            sigma_history = {key: [] for key in self.sigma}

            # Initialize debug file
            debug_file = os.path.join(output_dir, "debug_mcmc.txt")
            with open(debug_file, "w") as f:
                f.write("# MCMC Debugging Information\n")
                f.write("# This file tracks excluded volume violations and acceptance decisions\n")
                f.write("-" * 80 + "\n\n")

            current_score, curr_excl_initial, curr_pair_initial, curr_prior_initial = self.calculate_score(
                self.positions_ps, self.sigma, self.sigma_range, step=0
            )
            best_score = current_score  # Start with current as best

            # Log initial state
            with open(debug_file, "a") as f:
                f.write(f"INITIAL STATE\n")
                f.write(f"Initial Score: {current_score:.3f}, ExVol: {curr_excl_initial:.3f}, "
                        f"Pair: {curr_pair_initial:.3f}, Prior: {curr_prior_initial:.3f}\n\n")
                
                # Check for overlaps in initial state
                self._log_overlaps(f, self.positions_ps, "Initial")

            # Rest of setup is unchanged
            initial_temp = 5.0
            final_temp = 1.0
            cooling_factor = -np.log(final_temp / initial_temp) / n_steps
            temperatures = initial_temp * np.exp(-cooling_factor * np.arange(n_steps * 50))

            os.makedirs(output_dir, exist_ok=True)
            csv_log_file = os.path.join(output_dir, "all_info_mcmc.csv")
            all_log_file = os.path.join(output_dir, "all_steps_scores.csv")
            
            # Reset the pairs log file at the start of each run
            self.pairs_log_file = os.path.join(output_dir, "pairs_log.txt")
            with open(self.pairs_log_file, "w") as f:
                f.write("# MCMC Pair Selection Log\n")

            with open(csv_log_file, "w") as f:
                f.write("Step,Prior,Exvol_score,Pair_score,Score\n")
            
            with open(all_log_file, "w") as f:
                f.write("Step,Proposed_Score,Proposed_Exvol,Proposed_Pair,Proposed_Prior,Delta_E,accepted\n")

            accepted_moves = 0
            total_moves = 0
            temp_index = 0  # Only increment on accepted moves
            
            # Track statistics for debugging
            total_with_overlaps = 0
            accepted_with_overlaps = 0

            while accepted_moves < n_steps and total_moves < (n_steps * 50):
                total_moves += 1
                # Always use safe temperature access to prevent index errors
                temp = temperatures[min(temp_index, len(temperatures)-1)]
                move_type_is_position = (np.random.random() < position_move_prob)

                if move_type_is_position:
                    proposed_positions = self.propose_position_move(self.positions_ps)
                    proposed_sigma = self.sigma
                    move_description = "position move"
                else:
                    proposed_positions = self.positions_ps
                    proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)
                    move_description = f"sigma move for {pair_type}"

                # Calculate score with the current step number for proper logging
                proposed_score, curr_excl, curr_pair, curr_prior = self.calculate_score(
                    proposed_positions, proposed_sigma, self.sigma_range, step=total_moves
                )
                delta_e = proposed_score - current_score
                
                # Check for overlaps before Metropolis decision
                has_overlaps = self._check_overlaps(proposed_positions)
                if has_overlaps:
                    total_with_overlaps += 1

                # Metropolis acceptance criterion
                accepted = False
                if delta_e < 0 or np.random.random() < np.exp(-delta_e / temp):
                    accepted = True
                    # Move accepted - update state and increment temp_index
                    self.positions_ps = proposed_positions
                    self.sigma = proposed_sigma
                    current_score = proposed_score
                    accepted_moves += 1
                    temp_index += 1  # Only increment temperature on accepted moves
                    
                    if has_overlaps:
                        accepted_with_overlaps += 1

                    if current_score < best_score:
                        best_score = current_score
                        best_positions = {k: v.copy() for k, v in self.positions_ps.items()}

                    if accepted_moves % save_freq == 0:
                        for key in sigma_history:
                            sigma_history[key].append(self.sigma[key])

                        trajectory.append(
                            self.save_state(
                                accepted_moves, self.positions_ps, self.sigma, current_score,
                                prior_score=curr_prior,
                                pair_score=curr_pair,
                                exvol_score=curr_excl
                            )
                        )

                        with open(csv_log_file, "a") as f:
                            f.write(f"{accepted_moves},{curr_prior:.3f},{curr_excl:.3f},"
                                    f"{curr_pair:.3f},{current_score:.3f}\n")

                        accept_rate = accepted_moves / total_moves
                        print(f"Accepted Step {accepted_moves}, Score: {current_score:.2f}, T: {temp:.4f}, AcceptRate: {accept_rate:.2f}")
                
                # Write debug information for all steps with overlaps, and periodic summaries
                if has_overlaps or total_moves % 1000 == 0:
                    with open(debug_file, "a") as f:
                        f.write(f"Move {total_moves} ({move_description}): ")
                        f.write(f"Proposed Score={proposed_score:.1f}, ExVol={curr_excl:.3f}, ")
                        f.write(f"Pair={curr_pair:.1f}, Prior={curr_prior:.3f}, ")
                        f.write(f"Delta={delta_e:.1f}, Temp={temp:.3f}, ")
                        f.write(f"Accepted: {accepted}\n")
                        
                        if has_overlaps:
                            self._log_overlaps(f, proposed_positions, "Proposed")
                            f.write(f"OVERLAP STATE {'ACCEPTED' if accepted else 'REJECTED'}\n\n")
                
                # Log all steps to all_log_file
                with open(all_log_file, "a") as f:
                    f.write(f"{total_moves},{proposed_score:.1f},{curr_excl:.1f},{curr_pair:.1f},{curr_prior:.1f},{delta_e:.1f},{accepted_moves:.1f}\n")
                
                # Periodically write summary statistics
                if total_moves % 1000 == 0:
                    with open(debug_file, "a") as f:
                        overlap_rate = total_with_overlaps / total_moves * 100
                        if total_with_overlaps > 0:
                            overlap_accept_rate = accepted_with_overlaps / total_with_overlaps * 100
                        else:
                            overlap_accept_rate = 0
                            
                        f.write(f"\nSUMMARY at move {total_moves}:\n")
                        f.write(f"- Total moves with overlaps: {total_with_overlaps}/{total_moves} ({overlap_rate:.1f}%)\n")
                        f.write(f"- Moves with overlaps that were accepted: {accepted_with_overlaps}/{total_with_overlaps} ({overlap_accept_rate:.1f}%)\n")
                        f.write(f"- Overall acceptance rate: {accepted_moves/total_moves:.3f}\n\n")
                        f.write("-" * 50 + "\n\n")
                        
            # Store final positions
            final_positions = {k: v.copy() for k, v in self.positions_ps.items()}

            # Write final debug summary
            with open(debug_file, "a") as f:
                f.write("\nFINAL SUMMARY:\n")
                f.write(f"- Total moves: {total_moves}\n")
                f.write(f"- Accepted moves: {accepted_moves}\n")
                f.write(f"- Overall acceptance rate: {accepted_moves/total_moves:.3f}\n")
                f.write(f"- Total moves with overlaps: {total_with_overlaps}/{total_moves} ({total_with_overlaps/total_moves*100:.1f}%)\n")
                if total_with_overlaps > 0:
                    f.write(f"- Acceptance rate for moves with overlaps: {accepted_with_overlaps}/{total_with_overlaps} ({accepted_with_overlaps/total_with_overlaps*100:.1f}%)\n")
                f.write(f"- Final score: {current_score:.2f}\n")
                f.write(f"- Best score: {best_score:.2f}\n")
                
                # Check final state for overlaps
                self._log_overlaps(f, final_positions, "Final")
            
            # Clean up
            self.pairs_log_file = None  # Clear the file handle reference

            sigma_history_df = pd.DataFrame(sigma_history)
            sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)

            trajectory_file = os.path.join(output_dir, "trajectory.h5")
            final_file = self.save_trajectory(trajectory, trajectory_file)

            # Return both best and final positions for better analysis
            return final_positions, trajectory, final_file
            
    def _check_overlaps(self, positions):
        """Check if there are any overlaps in the given positions."""
        has_overlaps = False
        for type1, pos1 in positions.items():
            for type2, pos2 in positions.items():
                if type1 <= type2:  # Only check unique pairs
                    min_dist = self.params.radii[type1] + self.params.radii[type2]
                    distances = cdist(pos1, pos2)
                    
                    # Handle same-type particles
                    if type1 == type2:
                        mask = np.triu(np.ones_like(distances), k=1)
                        if np.any((distances < min_dist) & (mask > 0)):
                            has_overlaps = True
                            break
                    else:
                        if np.any(distances < min_dist):
                            has_overlaps = True
                            break
            if has_overlaps:
                break
        return has_overlaps
        
    def _log_overlaps(self, file_handle, positions, state_name):
        """Log detailed overlap information to the given file handle."""
        file_handle.write(f"{state_name} State Overlaps:\n")
        
        total_overlaps = 0
        for type1, pos1 in positions.items():
            for type2, pos2 in positions.items():
                if type1 <= type2:  # Only check unique pairs
                    min_dist = self.params.radii[type1] + self.params.radii[type2]
                    distances = cdist(pos1, pos2)
                    
                    # Handle same-type particles
                    if type1 == type2:
                        mask = np.triu(np.ones_like(distances), k=1)
                        viol_mask = (distances < min_dist) & (mask > 0)
                    else:
                        viol_mask = distances < min_dist
                    
                    if np.any(viol_mask):
                        viol_indices = np.where(viol_mask)
                        for i, j in zip(viol_indices[0], viol_indices[1]):
                            distance = distances[i, j]
                            overlap = min_dist - distance
                            if overlap > 0:
                                total_overlaps += 1
                                file_handle.write(f"  {type1}{i}-{type2}{j}: dist={distance:.3f}, "
                                                f"min_dist={min_dist:.3f}, overlap={overlap:.3f}\n")
                                
                                # Limit output to avoid huge files
                                if total_overlaps >= 20:
                                    file_handle.write(f"  ... additional overlaps not shown\n")
                                    return
        
        if total_overlaps == 0:
            file_handle.write("  No overlaps detected\n")
        file_handle.write("\n")
#----------------------------------------------------------------------
if __name__ == "__main__":
    pair_samp = PairSampler()
    pair_samp.run_mc(n_steps=10000, save_freq=100)
#----------------------------------------------------------------------
