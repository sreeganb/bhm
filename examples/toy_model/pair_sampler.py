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
import h5py
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
            #visualize_3d_configuration(self.positions_ps, self.params.radii)  # Visualize initial configuration
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
        debug_pairs=False,
        step: int = 0
    ) -> Tuple[float, float, float, float]:
        """Calculate the log posterior with detailed debugging of chosen pairs."""
        # 1) Excluded volume contribution
        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)
        
        # 2) Pairwise negative log-likelihood
        pairwise_score = 0.0
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
        
        # Convert excluded_pairs to a set of standardized tuples for faster lookup
        excluded_set = set()
        if excluded_pairs:
            for p in excluded_pairs:
                standardized = (p[0], p[1], p[2], p[3])
                excluded_set.add(standardized)
        
        # For tracking chosen pairs and detecting duplicates
        chosen_pairs_log = []
        chosen_pairs_set = set()  # For fast duplicate checking
        excluded_but_chosen = []  # Track excluded pairs that were still chosen
        duplicate_chosen = []     # Track pairs chosen multiple times
        
        # Log file for debugging
        pairs_log_file = "output_analysis/pair_debugging.txt"
        if debug_pairs:
            with open(pairs_log_file, "a") as log_file:
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"PAIR SCORING - STEP {step}\n")
                log_file.write(f"{'='*50}\n")
                log_file.write(f"Number of excluded pairs: {len(excluded_set)}\n\n")
        
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
                
                # Exclude certain pairs
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
                
                if debug_pairs:
                    with open(pairs_log_file, "a") as log_file:
                        log_file.write(f"\n--- {pair_key} pairs ---\n")
                
                for i, j in unique_pairs:
                    # Skip infinite scores (excluded or diagonal)
                    if np.isinf(score_matrix[i, j]):
                        continue
                        
                    # Construct pair identifier
                    pair_tuple = (type1, i, type2, j)
                    pair_str = f"{type1}({i})-{type2}({j})"
                    
                    # Check if this pair is excluded but still being scored (shouldn't happen)
                    if pair_tuple in excluded_set:
                        excluded_but_chosen.append(pair_str)
                    
                    # Check if this pair is already chosen (shouldn't happen either)
                    if pair_str in chosen_pairs_set:
                        duplicate_chosen.append(pair_str)
                    
                    # Add to chosen sets and calculate score
                    chosen_pairs_set.add(pair_str)
                    chosen_pairs_log.append(pair_str)
                    pairwise_score += self.pair_weight * score_matrix[i, j]
                    
                    # Log pair information
                    if debug_pairs:
                        with open(pairs_log_file, "a") as log_file:
                            dist = np.sqrt(np.sum((pos[type1][i] - pos[type2][j])**2))
                            score = score_matrix[i, j]
                            log_file.write(f"  {pair_str}: dist={dist:.2f}, target={target_dist:.2f}, score={score:.2f}\n")
        
        # 3) Prior penalty
        if not use_sigma_distribution:
            prior_penalty = self.priors.neg_log_prior(sig, sig_range)
        else:
            prior_penalty = prior_penalty_from_distribution
        
        total_score = exclusion_score + pairwise_score + prior_penalty
        
        # Write debugging summary with emphasis on potential issues
        if debug_pairs:
            with open(pairs_log_file, "a") as log_file:
                log_file.write(f"\n{'='*50}\n")
                log_file.write(f"PAIR SCORING SUMMARY - STEP {step}\n")
                log_file.write(f"{'='*50}\n")
                log_file.write(f"Total scored pairs: {len(chosen_pairs_log)}\n")
                
                # Report on excluded pairs that were wrongly scored
                if excluded_but_chosen:
                    log_file.write(f"\n>>> ERROR: Found {len(excluded_but_chosen)} excluded pairs that were scored! <<<\n")
                    for p in excluded_but_chosen:
                        log_file.write(f"  {p}\n")
                else:
                    log_file.write("\n✓ No excluded pairs were scored (correct)\n")
                
                # Report on duplicate pairs
                if duplicate_chosen:
                    log_file.write(f"\n>>> ERROR: Found {len(duplicate_chosen)} duplicate pairs! <<<\n")
                    for p in duplicate_chosen:
                        log_file.write(f"  {p}\n")
                else:
                    log_file.write("\n✓ No duplicate pairs found (correct)\n")
                
                # Write final score summary
                log_file.write(f"\nScore Components:\n")
                log_file.write(f"  Exclusion score: {exclusion_score:.2f}\n")
                log_file.write(f"  Pairwise score: {pairwise_score:.2f}\n")
                log_file.write(f"  Prior penalty: {prior_penalty:.2f}\n")
                log_file.write(f"  Total score: {total_score:.2f}\n")
                log_file.write(f"\n{'='*50}\n\n")
        
        return total_score, exclusion_score, pairwise_score, prior_penalty
    #----------------------------------------------------------------------
    def run_mc(self,
            n_steps: int = 50000,
            save_freq: int = 100,
            output_dir: str = "output_analysis/pairsampler_results/",
            position_move_prob: float = 0.9,
            debug: bool = True) -> Tuple[Dict[str, np.ndarray], str]:
        """
        Monte Carlo sampling with position and sigma moves.
        Adds optional debug prints and removes garbage collection.
        """

        # Setup output directory
        os.makedirs(output_dir, exist_ok=True)
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        with h5py.File(trajectory_file, 'w'):
            pass

        # Track sigma history in lists
        sigma_history = {key: [] for key in self.sigma}

        # Initial score
        current_score, curr_excl, curr_pair, curr_prior = self.calculate_score(
            self.positions_ps, self.sigma, self.sigma_range, debug_pairs = False, step=0
        )
        best_score = current_score

        # Temperature schedule
        initial_temp, final_temp = 5.0, 1.0
        cooling_factor = -np.log(final_temp / initial_temp) / n_steps
        temperatures = initial_temp * np.exp(-cooling_factor * np.arange(n_steps * 50))

        # Logging files
        csv_log_file = os.path.join(output_dir, "all_info_mcmc.csv")
        with open(csv_log_file, "w") as f:
            f.write("Step,Prior,Exvol_score,Pair_score,Score\n")

        self.pairs_log_file = os.path.join(output_dir, "pairs_log.txt")
        with open(self.pairs_log_file, "w") as f:
            f.write("# MCMC Pair Selection Log\n")

        if debug:
            print(f"Starting MCMC with up to {n_steps} accepted steps.")

        accepted_moves, total_moves, temp_index = 0, 0, 0
        while accepted_moves < n_steps and total_moves < (n_steps * 50):
            total_moves += 1
            temp = temperatures[min(temp_index, len(temperatures) - 1)]
            move_type_is_position = (np.random.random() < position_move_prob)

            # Propose move
            if move_type_is_position:
                proposed_positions = self.propose_position_move(self.positions_ps)
                proposed_sigma = self.sigma
            else:
                proposed_positions = self.positions_ps
                proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)

            # Score proposed move
            proposed_score, excl_val, pair_val, prior_val = self.calculate_score(
                proposed_positions, proposed_sigma, self.sigma_range, debug_pairs=False, step=total_moves
            )
            delta_e = proposed_score - current_score

            # Metropolis accept/reject
            if delta_e < 0 or np.random.random() < np.exp(-delta_e / temp):
                self.positions_ps = proposed_positions
                self.sigma = proposed_sigma
                current_score = proposed_score
                accepted_moves += 1
                temp_index += 1

                # Update best score
                if current_score < best_score:
                    best_score = current_score

                # Periodic saving
                if accepted_moves % save_freq == 0:
                    for key in sigma_history:
                        sigma_history[key].append(self.sigma[key])

                    self.save_state_to_disk(
                        accepted_moves, self.positions_ps, self.sigma, current_score,
                        prior_score=prior_val, pair_score=pair_val, exvol_score=excl_val,
                        traj_file=trajectory_file
                    )

                    with open(csv_log_file, "a") as f:
                        f.write(f"{accepted_moves},{prior_val:.3f},{excl_val:.3f},"
                                f"{pair_val:.3f},{current_score:.3f}\n")

                    if debug:
                        accept_rate = accepted_moves / total_moves
                        print(f"Accepted Step {accepted_moves}, "
                            f"Score: {current_score:.2f}, T: {temp:.4f}, "
                            f"AcceptRate: {accept_rate:.2f}")

        # Store final positions
        final_positions = {k: v.copy() for k, v in self.positions_ps.items()}
        self.pairs_log_file = None

        # Save sigma history
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)

        if debug:
            print("MCMC sampling complete.")
            print(f"Accepted moves: {accepted_moves} / {total_moves}")
            print(f"Best score: {best_score:.3f}")

        return final_positions, trajectory_file     
#----------------------------------------------------------------------
if __name__ == "__main__":
    pair_samp = PairSampler()
    pair_samp.run_mc(n_steps=10000, save_freq=100)
#----------------------------------------------------------------------
