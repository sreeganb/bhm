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
            #self.priors = Priors("jeffreys")  # Initialize priors
            self.priors = Priors()  # Initialize priors

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
        pos: Dict[str, np.ndarray],  # input positions
        sig: Dict[str, float],       # input sigma values
        sig_range: Dict[str, Tuple[float, float]],  # input sigma ranges
        excluded_pairs=None,
        use_sigma_distribution=False,
        prior_penalty_from_distribution=0.0,
        debug=False,  # Whether to log detailed pair scoring info
        debug_file="pair_score_debug.csv"  # File to write debug info to
    ) -> Tuple[float, float, float, float]:
        """Calculate the log posterior for the pair-level interactions."""
        # Setup debug file if requested
        if debug:
            debug_fh = open(debug_file, 'w')
            debug_fh.write("Pair Type, Particle1 Type, Particle1 Index, Particle2 Type, Particle2 Index, Distance, Target Distance, Sigma, Score\n")
        else:
            debug_fh = None
            
    def calculate_score(
        self,
        pos: Dict[str, np.ndarray],  # input positions
        sig: Dict[str, float],       # input sigma values
        sig_range: Dict[str, Tuple[float, float]],  # input sigma ranges
        excluded_pairs=None,
        use_sigma_distribution=False,
        prior_penalty_from_distribution=0.0,
        debug=False,  # Whether to log detailed pair scoring info
        debug_file="pair_score_debug.csv"  # File to write debug info to
    ) -> Tuple[float, float, float, float]:
        """Calculate the log posterior for the pair-level interactions."""
        # Setup debug file if requested
        if debug:
            debug_fh = open(debug_file, 'w')
            debug_fh.write("Pair Type, Particle1 Type, Particle1 Index, Particle2 Type, Particle2 Index, Distance, Target Distance, Sigma, Score\n")
        else:
            debug_fh = None
            
        # 1) Excluded volume contribution
        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)
                    
        # 2) Pairwise negative log-likelihood (pair-specific)
        pairwise_score = 0.0
        #pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C')]
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
                # ----------------------------------------------------------------------
                # Zero out excluded pairs if provided
                # ----------------------------------------------------------------------
                if excluded_pairs:
                    print(f"Exclusion pairs: {excluded_pairs}")
                    for i_idx in range(len(pos[type1])):
                        for j_idx in range(len(pos[type2])):
                            # Check if (type1, i_idx, type2, j_idx) or the reverse is in excluded_pairs
                            if ((type1, i_idx, type2, j_idx) in excluded_pairs or
                            (type2, j_idx, type1, i_idx) in excluded_pairs):
                                score_matrix[i_idx, j_idx] = np.inf
                #print(f"Pairwise score matrix for {pair_key}:\n{score_matrix}")
                # For debugging, calculate the actual distances between particles
                if debug_fh:
                    distance_matrix = cdist(pos[type1], pos[type2])

                # Mask the lower triangle and diagonal to avoid self-interactions and duplicates
                if type1 == type2:
                #----------------------------------------------------------------------
                # For same type particles, use loops to only consider upper triangle - diagonal elements 
                # for each row and column, get the minimum element
                # do not fill matrix with inf 
                #----------------------------------------------------------------------
                    np.fill_diagonal(score_matrix, np.inf)
                    m,n = score_matrix.shape
                    selected_scores = []
                    selected_indices = []
                    # for each row i, find the minimum element in columns > i
                    for i in range(m):
                        # only consider j > i in that row
                        if i < n-1:
                            row_segment = score_matrix[i, i+1:]
                            # If all inf, skip
                            if np.all(np.isinf(row_segment)):
                                continue
                            j_offset = np.argmin(row_segment)
                            j = i + 1 + j_offset
                            selected_scores.append(score_matrix[i, j])
                            selected_indices.append((i, j))
                    # detailed debug information for each selected pair 
                    if debug_fh:
                        for i, j in selected_indices:
                            distance = distance_matrix[i, j]
                            pair_score = self.pair_weight * score_matrix[i, j]
                            debug_fh.write(f"{pair_key},{type1},{i},{type2},{j},{distance:.6f},{target_dist:.6f},{sigma_value:.6f},{pair_score:.6f}\n")
                else:
                    # For different types, use argmin (indices), not min (values)
                    row_min_indices = np.argmin(score_matrix, axis=1)
                    col_min_indices = np.argmin(score_matrix, axis=0)

                    # Collect unique index pairs
                    U = set()
                    for i, min_idx in enumerate(row_min_indices):
                        U.add((i, min_idx))
                    for j, min_idx in enumerate(col_min_indices):
                        U.add((min_idx, j))

                    # Now gather the scores for those index pairs
                    selected_scores = [score_matrix[i, j] for (i, j) in U if i < score_matrix.shape[0] and j < score_matrix.shape[1]]

                # Sum up the selected scores
                pairwise_score_part = np.sum(selected_scores)
                pairwise_score += self.pair_weight * pairwise_score_part
                #----------------------------------------------------------------------
        # 3) Prior penalty: if not using the sigma distribution, compute using Priors
        if not use_sigma_distribution:
            prior_penalty = self.priors.neg_log_prior(sig, sig_range)
        else:
            prior_penalty = prior_penalty_from_distribution

        # Clean up debug file if opened
        if debug_fh:
            debug_fh.close()
            print(f"Particle pairing debug information written to {debug_file}")

        total_score = exclusion_score + pairwise_score + prior_penalty
        return total_score, exclusion_score, pairwise_score, prior_penalty                
    #----------------------------------------------------------------------
    def run_mc(
        self,
        n_steps: int = 50000,
        save_freq: int = 100,
        output_dir: str = "output_analysis/pairsampler_results/",
        debug: bool = True
    ) -> Tuple[Dict[str, np.ndarray], str]:
        """
        Monte Carlo sampling with position and sigma moves.
        Runs for a fixed number of total steps (accepted + rejected).
        """
        # Setup output directory
        os.makedirs(output_dir, exist_ok=True)
        trajectory_file = os.path.join(output_dir, "trajectory.h5")
        with h5py.File(trajectory_file, 'w') as f:
            pass  # Create empty file
        
        # Initialize tracking variables
        best_positions = None
        best_score = float('inf')
        sigma_history = {key: np.zeros(n_steps // save_freq + 1) for key in self.sigma}
        accepts = {'position': 0, 'sigma': 0}
        attempts = {'position': 0, 'sigma': 0}
        accepted_moves = 0  # Track total accepted moves
        
        # Initial score calculation
        current_score, curr_excl, curr_pair, curr_prior = self.calculate_score(
            self.positions_ps, self.sigma, self.sigma_range, 
            set(), False, 0.0
        )
        
        # Store initial sigma values
        for key in self.sigma:
            sigma_history[key][0] = self.sigma[key]
        
        # Setup logging files
        csv_log_file = os.path.join(output_dir, "all_info_mcmc.csv")
        with open(csv_log_file, "w") as f:
            f.write("Step,Prior,Exvol_score,Pair_score,Score,Accept_rate\n")
            f.write(f"0,{curr_prior:.3f},{curr_excl:.3f},{curr_pair:.3f},{current_score:.3f},0.0\n")

        if debug:
            self.pairs_log_file = os.path.join(output_dir, "pairs_log.txt")
            with open(self.pairs_log_file, "w") as f:
                f.write("# MCMC Pair Selection Log\n")
        
        # Main MCMC loop: run for exactly n_steps total moves
        move_types = ['position', 'sigma']
        move_probs = [0.9, 0.1]  # position, sigma
        
        # Simple cooling schedule based on total steps
        temp_start, temp_end = 5.0, 1.0
        temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
        
        print(f"Starting MCMC sampling for {n_steps} total steps...")
        
        for step in range(1, n_steps + 1):
            # Temperature schedule based on total steps
            temp = temp_start * (temp_decay ** step)
            
            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            attempts[move_type] += 1
            
            # Propose move
            pair_type = None  # Will track which sigma parameter was modified
            if move_type == 'position':
                proposed_positions = self.propose_position_move(
                    self.positions_ps, 
                    accepts['position'] / max(1, step)
                )
                proposed_sigma = self.sigma
            else:  # sigma move
                proposed_positions = self.positions_ps
                proposed_sigma, pair_type = self.propose_sigma_move(
                    self.sigma, 
                    accepts['sigma'] / max(1, step)
                )
            
            # Calculate new score
            proposed_score, prop_excl, prop_pair, prop_prior = self.calculate_score(
                proposed_positions, proposed_sigma, self.sigma_range,
                set(), False, 0.0
            )
            
            # Metropolis criterion with Jacobian correction for sigma moves
            delta = proposed_score - current_score
            
            # Add Jacobian correction for sigma moves (log(sigma'/sigma))
            if move_type == 'sigma' and pair_type is not None:
                # Add log(sigma'/sigma) to delta for proper detailed balance
                jacobian_term = np.log(proposed_sigma[pair_type] / self.sigma[pair_type])
                delta += jacobian_term
                if debug and step % save_freq == 0:
                    print(f"  Sigma move: {pair_type} {self.sigma[pair_type]:.4f}->{proposed_sigma[pair_type]:.4f}, "
                        f"Jacobian term: {jacobian_term:.4f}")
            
            # Acceptance probability
            accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
            
            if accept:
                # Update state
                self.positions_ps = proposed_positions
                self.sigma = proposed_sigma
                current_score = proposed_score
                curr_excl, curr_pair, curr_prior = prop_excl, prop_pair, prop_prior
                
                # Track best state
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in self.positions_ps.items()}
                
                accepts[move_type] += 1
                accepted_moves += 1
            
            # Save state every save_freq total steps
            if step % save_freq == 0:
                # Store sigma history
                save_idx = step // save_freq
                if save_idx < len(sigma_history[list(sigma_history.keys())[0]]):
                    for key in self.sigma:
                        sigma_history[key][save_idx] = self.sigma[key]
                
                # Save to disk using existing method
                self.save_state_to_disk(
                    step, self.positions_ps, self.sigma, current_score,
                    prior_score=curr_prior, pair_score=curr_pair,
                    exvol_score=curr_excl, tet_score=0.0,
                    traj_file=trajectory_file
                )
                
                # Log to CSV
                acceptance_rate = accepted_moves / step
                with open(csv_log_file, "a") as f:
                    f.write(f"{step},{curr_prior:.3f},{curr_excl:.3f},"
                            f"{curr_pair:.3f},{current_score:.3f},{acceptance_rate:.3f}\n")
                
                # Print progress
                print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={acceptance_rate:.2f}")
        
        # Save sigma history
        sigma_history_df = pd.DataFrame(sigma_history)
        sigma_history_df.to_csv(os.path.join(output_dir, "sigma_history.csv"), index=False)
        
        # Clean up
        if debug:
            self.pairs_log_file = None
        
        # Print final statistics
        print("\nSampling complete:")
        for move_type in move_types:
            rate = accepts[move_type] / max(1, attempts[move_type])
            print(f"- {move_type}: {rate:.2f} acceptance ({accepts[move_type]}/{attempts[move_type]})")
        print(f"- Total steps: {n_steps}")
        print(f"- Accepted moves: {accepted_moves}")
        print(f"- Overall acceptance rate: {accepted_moves / n_steps:.2f}")
        print(f"- Best score: {best_score:.2f}")
        
        # Return best positions if found, otherwise current positions
        return best_positions if best_positions is not None else self.positions_ps, trajectory_file
#----------------------------------------------------------------------
if __name__ == "__main__":
    pair_samp = PairSampler()
    pair_samp.run_mc(n_steps=10000, save_freq=100)
#----------------------------------------------------------------------
