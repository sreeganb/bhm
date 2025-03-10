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
        pos: Dict[str, np.ndarray],  # input positions
        sig: Dict[str, float],       # input sigma values
        sig_range: Dict[str, Tuple[float, float]] = None,  # input sigma ranges
        excluded_pairs=None,
        use_sigma_distribution=False,
        prior_penalty_from_distribution=0.0
    ) -> Tuple[float, float, float, float]:
        """Calculate the log posterior for the pair-level interactions."""
        # 1) Excluded volume contribution (from BaseMCSampler)
        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)

        # 2) Pairwise negative log-likelihood (pair-specific)
        pairwise_score = 0.0
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
        for type1, type2 in pair_types:
            pair_key = f"{type1}{type2}"
            if pair_key in self.params.pair_distances:
                target_dist = self.params.pair_distances[pair_key]
                sigma_value = sig[pair_key]  # Use the passed sigma

                # Calculate pairwise score matrix
                score_matrix = self.calculate_pair_scores_matrix(
                    pos[type1],
                    pos[type2],
                    target_dist,
                    sigma_value
                )

                # Mask diagonal (for same-type pairs)
                if type1 == type2:
                    mask = np.triu(np.ones_like(score_matrix), k=1)
                    score_matrix *= mask

                # excluded pairs should not be considered, so set their score to a large value to avoid selection
                if excluded_pairs:
                    for i in range(len(pos[type1])):
                        for j in range(len(pos[type2])):
                            if (type1, i, type2, j) in excluded_pairs:
                                #score_matrix[i, j] = 0
                                score_matrix[i, j] = 99999999.0  # Set to high value to avoid selection

                # Identify minimal row/column pairs and accumulate the score
                row_indices = np.argmin(score_matrix, axis=1)
                col_indices = np.argmin(score_matrix, axis=0)
                row_pairs = {(i, row_indices[i]) for i in range(len(row_indices))}
                col_pairs = {(col_indices[j], j) for j in range(len(col_indices))}
                unique_pairs = row_pairs.union(col_pairs)
                for i, j in unique_pairs:
                    pairwise_score += self.pair_weight * score_matrix[i, j]

        # 3) Prior penalty: if not using the sigma distribution, compute using Priors.
        prior_penalty = 0.0
        if not use_sigma_distribution:
            prior_penalty = self.priors.neg_log_prior(sig, sig_range)
        else:
            prior_penalty = prior_penalty_from_distribution

        total_score = exclusion_score + pairwise_score + prior_penalty
        return total_score, exclusion_score, pairwise_score, prior_penalty

    def run_mc(self, n_steps: int = 50000, save_freq: int = 100, output_dir : str = "output_analysis/pairsampler_results/") -> Tuple[Dict[str, np.ndarray], List[Dict], str]:
        best_positions = None  # Will set on first improvement
        trajectory = []
        sigma_history = {key: [] for key in self.sigma}  # Pre-allocate sigma history

        # Initial score calculation
        current_score, _, _, _ = self.calculate_score(self.positions_ps, self.sigma, self.sigma_range)
        best_score = float('inf')  # Initialize high to ensure first improvement updates

        # Pre-calculate temperature schedule
        initial_temp = 5.0
        final_temp = 1.0
        cooling_factor = -np.log(final_temp / initial_temp) / n_steps
        temperatures = initial_temp * np.exp(-cooling_factor * np.arange(n_steps))

        # Prepare output directories and files
        #output_dir = os.path.join("output_analysis", "pairsampler_results")
        os.makedirs(output_dir, exist_ok=True)
        csv_log_file = os.path.join(output_dir, "all_info_mcmc.csv")

        # Write CSV header
        with open(csv_log_file, "w") as f:
            f.write("Step,Prior,Exvol_score,Pair_score,Score,Accepted\n")

        accepted_moves = 0

        # Main MCMC loop
        for step in range(n_steps):
            temp = temperatures[step]
            move_type_is_position = (np.random.random() < 0.9)

            # Propose move: position move or sigma move
            if move_type_is_position:
                proposed_positions = self.propose_position_move(self.positions_ps)
                proposed_sigma = self.sigma  # use current sigma (by reference)
            else:
                proposed_positions = self.positions_ps  # positions remain unchanged
                proposed_sigma, pair_type = self.propose_sigma_move(self.sigma)

            # Calculate proposed score
            proposed_score, curr_excl, curr_pair, curr_prior = self.calculate_score(
                proposed_positions, proposed_sigma, self.sigma_range
            )
            delta_e = proposed_score - current_score

            # Metropolis criterion
            acceptance = 0
            if delta_e < 0 or np.random.random() < np.exp(-delta_e / temp):
                self.positions_ps = proposed_positions
                self.sigma = proposed_sigma
                current_score = proposed_score
                accepted_moves += 1
                acceptance = 1
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_ps.items()}

            # Batch logging and saving trajectory
            if step % save_freq == 0:
                accept_rate = accepted_moves / (step + 1)
                # Update sigma history
                for key in sigma_history:
                    sigma_history[key].append(self.sigma[key])

                # Append current state with detailed score breakdown.
                trajectory.append(
                    self.save_state(
                        step, self.positions_ps, self.sigma, current_score,
                        prior_score=curr_prior,
                        pair_score=curr_pair,
                        exvol_score=curr_excl
                    )
                )

                # Write progress to CSV log
                with open(csv_log_file, "a") as f:
                    f.write(f"{step},{curr_prior:.3f},{curr_excl:.3f},"
                            f"{curr_pair:.3f},{current_score:.3f},{acceptance:.1f}\n")
                
                print(f"Step {step}, Score: {current_score:.2f}, T: {temp:.4f}, AcceptRate: {accept_rate:.2f}")

        # Save sigma history as CSV in the output directory.
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)

        # Save the trajectory in HDF5 format in the output folder.
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        final_file = self.save_trajectory(trajectory, trajectory_file)
        
        return best_positions, trajectory, final_file
#----------------------------------------------------------------------
if __name__ == "__main__":
    pair_samp = PairSampler()
    pair_samp.run_mc(n_steps=10000, save_freq=100)
#----------------------------------------------------------------------
