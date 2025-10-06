#----------------------------------------------------------------------
# PairSampler class which will use the base sampler functions 
# and run the MCMC simulations for the pairwise sampling and 
# track sigma values
#----------------------------------------------------------------------
from typing import List, Dict, Tuple
import os
import numpy as np
import pandas as pd
from base_sampler import BaseMCSampler  # Import BaseMCSampler
from parameters import SystemParameters
from visualization import visualize_3d_configuration
from scipy.spatial.distance import cdist
from new_sigma_provider import GMMSigmaProvider
import h5py
import random
import pathlib
from scipy.special import logsumexp
#----------------------------------------------------------------------
class PairSampler(BaseMCSampler):
    """
    Sampler for pair-level interactions, inheriting from BaseMCSampler.
    Includes run_mc method.
    """
    def __init__(self, sampler_sequence: List[str], sequence_idx: int,
                 ex_weight: float = 1.0, pair_weight: float = 1.0, 
                 base_output_dir: str = None, specific_chain: int = None,
                 sigma_ranges: Dict[str, Tuple[float, float]] = None,
                 prior_type: str = "gamma",
                 pos_passed: Dict[str, np.ndarray] = None):
        """
        Initialize PairSampler with sampler sequence information.
        
        Args:
            sampler_sequence: List of sampler names in execution order
            sequence_idx: Current position in the sampler sequence
            ex_weight: Weight for exclusion volume term
            pair_weight: Weight for pairwise interactions
            base_output_dir: Base directory for analysis results
            specific_chain: Specific chain to load GMM from (if not first sampler)
            sigma_ranges: Ranges for sigma values
            prior_type: Prior type for first sampler ("uniform" or "jeffreys")
            pos_passed: Optional positions to use instead of initializing or loading
        """
        super().__init__()  # Call BaseMCSampler constructor
        
        # Store sampler sequence information
        self.sampler_sequence = sampler_sequence
        self.sequence_idx = sequence_idx
        self.specific_chain = specific_chain
        self.exclusion_weight = ex_weight
        self.pair_weight = pair_weight
        self.params = SystemParameters()  # Initialize system parameters
        
        # Initialize sigma provider
        self.sigma_provider = GMMSigmaProvider(
            sampler_sequence=sampler_sequence,
            sequence_idx=sequence_idx,
            base_output_dir=base_output_dir,
            specific_chain=specific_chain,
            sigma_ranges=sigma_ranges,
            prior_type=prior_type
        )
        
        # Initialize sigma values and ranges from provider
        self.sigma = self.sigma_provider.sample_sigma_values()
        self.sigma_range = self.sigma_provider.sigma_ranges
        
        # Initialize positions based on sequence and parameters
        if pos_passed is not None:
            # Use explicitly passed positions (highest priority)
            self.positions_ps = pos_passed
            print("PairSampler initialized with passed positions.")
        else:
            # Use sequence-based position loading
            self.positions_ps = self.get_positions()
        
        # Print initialization info
        sampler_name = sampler_sequence[sequence_idx]
        if self.sigma_provider.is_first_sampler:
            print(f"PairSampler ({sampler_name}) initialized as first sampler with {prior_type} prior.")
        else:
            print(f"PairSampler ({sampler_name}) initialized with GMM prior from previous sampler.")
        
        print(f"Initial sigma values: {self.sigma}")

    def get_positions(self) -> Dict[str, np.ndarray]:
        """
        Load positions from the previous sampler in the sequence or initialize if first.
        This method uses the class attributes and calls the standalone function.
        """
        return self.get_positions_from_previous_sampler(
            sampler_sequence=self.sampler_sequence,
            current_sequence_idx=self.sequence_idx,
            specific_chain=getattr(self, 'specific_chain', None),
            base_output_dir="output_analysis",
            params=self.params
        )

    def calculate_score(
        self,
        pos: Dict[str, np.ndarray],
        sig: Dict[str, float],
        sig_range: Dict[str, Tuple[float, float]] = None,
        excluded_pairs=None,
        debug=False,
        debug_file="pair_score_debug.csv"
    ) -> Tuple[float, float, float, float]:
        """Calculate the log posterior using union-of-argmin pairing strategy."""
        
        # 1) Excluded volume contribution
        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)
                    
        # 2) Pairwise negative log-likelihood
        pairwise_score = 0.0
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C')]
        
        for type1, type2 in pair_types:
            pair_key = f"{type1}{type2}"
            if pair_key not in self.params.pair_distances:
                continue
                
            target_dist = self.params.pair_distances[pair_key]
            sigma_value = sig[pair_key]

            # Calculate score matrix (WITH proper normalization term)
            score_matrix = self.calculate_pair_scores_matrix(
                pos[type1],
                pos[type2],
                target_dist,
                sigma_value
            )
            
            # Apply exclusions if provided
            if excluded_pairs:
                for i_idx in range(len(pos[type1])):
                    for j_idx in range(len(pos[type2])):
                        if ((type1, i_idx, type2, j_idx) in excluded_pairs or
                            (type2, j_idx, type1, i_idx) in excluded_pairs):
                            score_matrix[i_idx, j_idx] = np.inf

            # Use union-of-argmin strategy (boss's simple approach)
            if type1 == type2:
                # Same type: use upper triangle only
                np.fill_diagonal(score_matrix, np.inf)
                unique_pairs = set()
                
                # Row minimums (each particle finds best partner with index > itself)
                for i in range(len(score_matrix)):
                    row_segment = score_matrix[i, i+1:]
                    if len(row_segment) > 0 and not np.all(np.isinf(row_segment)):
                        j = i + 1 + np.argmin(row_segment)
                        unique_pairs.add((i, j))
                
                # Column minimums (check reverse direction)
                for j in range(len(score_matrix)):
                    col_segment = score_matrix[:j, j]
                    if len(col_segment) > 0 and not np.all(np.isinf(col_segment)):
                        i = np.argmin(col_segment)
                        unique_pairs.add((min(i,j), max(i,j)))  # Store in canonical order
            else:
                # Different types: union of row and column argmins
                unique_pairs = set()
                
                # Each type1 particle finds best type2 partner
                for i, j in enumerate(np.argmin(score_matrix, axis=1)):
                    if not np.isinf(score_matrix[i, j]):
                        unique_pairs.add((i, j))
                
                # Each type2 particle finds best type1 partner
                for j, i in enumerate(np.argmin(score_matrix, axis=0)):
                    if not np.isinf(score_matrix[i, j]):
                        unique_pairs.add((i, j))
            
            # Sum scores for unique pairs
            selected_scores = [score_matrix[i, j] for i, j in unique_pairs]
            pairwise_score += self.pair_weight * np.sum(selected_scores)
        
        # 3) Prior penalty using sigma provider (YOUR KEY IMPROVEMENT)
        prior_penalty = self.sigma_provider.calculate_negative_log_prior(sig)

        total_score = exclusion_score + pairwise_score + prior_penalty
        
        return total_score, exclusion_score, pairwise_score, prior_penalty

#    def calculate_score(
#        self,
#        pos: Dict[str, np.ndarray],
#        sig: Dict[str, float],
#        sig_range: Dict[str, Tuple[float, float]] = None,
#        excluded_pairs=None,
#        debug=False,
#        debug_file="pair_score_debug.csv"
#    ) -> Tuple[float, float, float, float]:
#        """Calculate the log posterior for the pair-level interactions."""
#        
#        # 1) Excluded volume contribution
#        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)
#                    
#        # 2) Pairwise negative log-likelihood
#        pairwise_score = 0.0
#        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C')]
#        
#        # Temperature parameter for Boltzmann selection (can be adjusted)
#        selection_temp = 3.0
#        
#        for type1, type2 in pair_types:
#            pair_key = f"{type1}{type2}"
#            if pair_key not in self.params.pair_distances:
#                continue
#                
#            target_dist = self.params.pair_distances[pair_key]
#            sigma_value = sig[pair_key]
#
#            # Calculate pairwise score matrix
#            score_matrix = self.calculate_pair_scores_matrix(
#                pos[type1],
#                pos[type2],
#                target_dist,
#                sigma_value
#            )
#            
#            # Apply exclusions if provided
#            if excluded_pairs:
#                for i_idx in range(len(pos[type1])):
#                    for j_idx in range(len(pos[type2])):
#                        if ((type1, i_idx, type2, j_idx) in excluded_pairs or
#                            (type2, j_idx, type1, i_idx) in excluded_pairs):
#                            score_matrix[i_idx, j_idx] = np.inf
#
#            # Handle same type vs different type particles
#            if type1 == type2:
#                # Same type: avoid self-interaction and use optimal pairing
#                np.fill_diagonal(score_matrix, np.inf)
#                selected_scores = self._solve_same_type_pairing_probabilistic(
#                    score_matrix, selection_temp
#                )
#            else:
#                # Different types: use Hungarian algorithm or probabilistic bipartite matching
#                selected_scores = self._solve_bipartite_pairing_probabilistic(
#                    score_matrix, selection_temp
#                )
#
#            # Sum selected scores
#            pairwise_score += self.pair_weight * np.sum(selected_scores)
#        
#        # 3) Prior penalty using sigma provider
#        prior_penalty = self.sigma_provider.calculate_negative_log_prior(sig)
#
#        total_score = exclusion_score + pairwise_score + prior_penalty
#        
#        return total_score, exclusion_score, pairwise_score, prior_penalty
#
#    def _solve_same_type_pairing_probabilistic(self, score_matrix, temp):
#        """
#        Solve same-type pairing with probabilistic selection but avoid double counting.
#        """
#        n = len(score_matrix)
#        used = set()
#        selected_scores = []
#        
#        # Create list of all valid pairs
#        valid_pairs = []
#        for i in range(n):
#            for j in range(i + 1, n):
#                if not np.isinf(score_matrix[i, j]):
#                    valid_pairs.append((i, j, score_matrix[i, j]))
#        
#        if not valid_pairs:
#            return []
#        
#        # Sort by score for better convergence
#        valid_pairs.sort(key=lambda x: x[2])
#        
#        # Probabilistic selection with exclusion
#        for _ in range(n // 2):  # Maximum possible pairs
#            # Filter available pairs (not already used)
#            available_pairs = [(i, j, score) for i, j, score in valid_pairs 
#                            if i not in used and j not in used]
#            
#            if not available_pairs:
#                break
#            
#            # Convert to Boltzmann probabilities
#            scores = np.array([score for _, _, score in available_pairs])
#            scores_shifted = scores - np.min(scores)
#            probs = np.exp(-scores_shifted / temp)
#            probs = probs / np.sum(probs)
#            
#            # Select pair probabilistically
#            idx = np.random.choice(len(available_pairs), p=probs)
#            i, j, score = available_pairs[idx]
#            
#            selected_scores.append(score)
#            used.add(i)
#            used.add(j)
#        
#        return selected_scores
#
#    def _solve_bipartite_pairing_probabilistic(self, score_matrix, temp):
#        """
#        Solve bipartite matching avoiding double counting using Hungarian-inspired approach.
#        """
#        n_type1, n_type2 = score_matrix.shape
#        
#        # Option 1: Use scipy's Hungarian algorithm for exact solution
#        try:
#            from scipy.optimize import linear_sum_assignment
#            row_ind, col_ind = linear_sum_assignment(score_matrix)
#            return [score_matrix[i, j] for i, j in zip(row_ind, col_ind) 
#                    if not np.isinf(score_matrix[i, j])]
#        except ImportError:
#            pass
#        
#        # Option 2: Greedy probabilistic assignment (fallback)
#        selected_scores = []
#        used_type2 = set()
#        
#        for i in range(n_type1):
#            # Find available partners for particle i
#            available_partners = []
#            for j in range(n_type2):
#                if j not in used_type2 and not np.isinf(score_matrix[i, j]):
#                    available_partners.append((j, score_matrix[i, j]))
#            
#            if not available_partners:
#                continue
#            
#            # Probabilistic selection among available partners
#            partners, scores = zip(*available_partners)
#            scores = np.array(scores)
#            scores_shifted = scores - np.min(scores)
#            probs = np.exp(-scores_shifted / temp)
#            probs = probs / np.sum(probs)
#            
#            # Select partner
#            if np.random.random() < 0.7:  # 70% greedy, 30% exploration
#                j_idx = np.argmin(scores)
#            else:
#                j_idx = np.random.choice(len(partners), p=probs)
#            
#            j = partners[j_idx]
#            selected_scores.append(score_matrix[i, j])
#            used_type2.add(j)
#        
#        return selected_scores
#    def calculate_score(
#        self,
#        pos: Dict[str, np.ndarray],  # input positions
#        sig: Dict[str, float],       # input sigma values
#        sig_range: Dict[str, Tuple[float, float]] = None,  # input sigma ranges
#        excluded_pairs=None,
#        debug=False,  # Whether to log detailed pair scoring info
#        debug_file="pair_score_debug.csv"  # File to write debug info to
#    ) -> Tuple[float, float, float, float]:
#        """Calculate the log posterior for the pair-level interactions."""
#        # Setup debug file if requested
#        if debug:
#            debug_fh = open(debug_file, 'w')
#            debug_fh.write("Pair Type, Particle1 Type, Particle1 Index, Particle2 Type, Particle2 Index, Distance, Target Distance, Sigma, Score\n")
#        else:
#            debug_fh = None
#            
#        # 1) Excluded volume contribution
#        exclusion_score = self.exclusion_weight * self.excluded_volume_nll(pos)
#                    
#        # 2) Pairwise negative log-likelihood (pair-specific)
#        pairwise_score = 0.0
#        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C')]
#        for type1, type2 in pair_types:
#            pair_key = f"{type1}{type2}"
#            if pair_key in self.params.pair_distances:
#                target_dist = self.params.pair_distances[pair_key]
#                sigma_value = sig[pair_key]  # Use the passed sigma
#
#                # Calculate pairwise score matrix
#                score_matrix = self.calculate_pair_scores_matrix(
#                    pos[type1],
#                    pos[type2],
#                    target_dist,
#                    sigma_value
#                )
#                
#                # Zero out excluded pairs if provided
#                if excluded_pairs:
#                    print(f"Exclusion pairs: {excluded_pairs}")
#                    for i_idx in range(len(pos[type1])):
#                        for j_idx in range(len(pos[type2])):
#                            # Check if (type1, i_idx, type2, j_idx) or the reverse is in excluded_pairs
#                            if ((type1, i_idx, type2, j_idx) in excluded_pairs or
#                            (type2, j_idx, type1, i_idx) in excluded_pairs):
#                                score_matrix[i_idx, j_idx] = np.inf
#                
#                # For debugging, calculate the actual distances between particles
#                if debug_fh:
#                    distance_matrix = cdist(pos[type1], pos[type2])
#
#                # Handle same type vs different type particles
#                if type1 == type2:
#                    # For same type particles, consider upper triangle only
#                    np.fill_diagonal(score_matrix, np.inf)
#                    m, n = score_matrix.shape
#                    selected_scores = []
#                    selected_indices = []
#                    
#                    # For each row i, find the minimum element in columns > i
#                    for i in range(m):
#                        if i < n-1:
#                            row_segment = score_matrix[i, i+1:]
#                            # If all inf, skip
#                            if np.all(np.isinf(row_segment)):
#                                continue
#                            j_offset = np.argmin(row_segment)
#                            j = i + 1 + j_offset
#                            selected_scores.append(score_matrix[i, j])
#                            selected_indices.append((i, j))
#                    
#                    # Debug information for each selected pair 
#                    if debug_fh:
#                        for i, j in selected_indices:
#                            distance = distance_matrix[i, j]
#                            pair_score = self.pair_weight * score_matrix[i, j]
#                            debug_fh.write(f"{pair_key},{type1},{i},{type2},{j},{distance:.6f},{target_dist:.6f},{sigma_value:.6f},{pair_score:.6f}\n")
#                else:
#                    # For different types, use argmin approach
#                    row_min_indices = np.argmin(score_matrix, axis=1)
#                    col_min_indices = np.argmin(score_matrix, axis=0)
#
#                    # Collect unique index pairs
#                    U = set()
#                    for i, min_idx in enumerate(row_min_indices):
#                        U.add((i, min_idx))
#                    for j, min_idx in enumerate(col_min_indices):
#                        U.add((min_idx, j))
#
#                    # Gather the scores for those index pairs
#                    selected_scores = [score_matrix[i, j] for (i, j) in U if i < score_matrix.shape[0] and j < score_matrix.shape[1]]
#
#                # Sum up the selected scores
#                pairwise_score_part = np.sum(selected_scores)
#                pairwise_score += self.pair_weight * pairwise_score_part
#        
#        # 3) Prior penalty using sigma provider
#        prior_penalty = self.sigma_provider.calculate_negative_log_prior(sig)
#
#        total_score = exclusion_score + pairwise_score + prior_penalty
#
#        debug = False  # Disable debug logging after first run
#        if debug:
#            # Create debug folder if it doesn't exist
#            debug_folder = "debug_info"
#            os.makedirs(debug_folder, exist_ok=True)
#            
#            # Use fixed filenames without timestamps
#            debug_file = f"{debug_folder}/debug_scores.csv"
#            sigma_file = f"{debug_folder}/sigma_values.csv"
#            
#            # Check if files exist to determine if headers are needed
#            scores_file_exists = os.path.exists(debug_file)
#            sigma_file_exists = os.path.exists(sigma_file)
#            
#            # Write scores
#            with open(debug_file, 'a') as f:
#                if not scores_file_exists:
#                    f.write("Total Score, Exclusion Score, Pairwise Score, Prior Penalty\n")
#                f.write(f"{total_score:.1f}, {exclusion_score:.1f}, {pairwise_score:.1f}, {prior_penalty:.1f}\n")
#            
#            # Write sigma values
#            with open(sigma_file, 'a') as f:
#                if not sigma_file_exists:
#                    f.write("Sigma Type, Sigma Value\n")
#                for key, value in sig.items():
#                    f.write(f"{key}, {value:.1f}\n")
#            
#        return total_score, exclusion_score, pairwise_score, prior_penalty

    def run_mc(
        self,
        n_steps: int = 50000,
        save_freq: int = 100,
        output_dir: str = None,
        debug: bool = True
    ) -> Tuple[Dict[str, np.ndarray], str]:
        """
        Monte Carlo sampling with position and sigma moves.
        Runs for a fixed number of total steps (accepted + rejected).
        """
        # Setup output directory with sampler sequence info
        if output_dir is None:
            sampler_name = self.sampler_sequence[self.sequence_idx]
            current_idx = self.sigma_provider.current_idx
            output_dir = f"output_analysis/{sampler_name.lower()}_results_{current_idx}/"
        
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
            set()
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
        temp_start, temp_end = 5.0, 0.1
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
            log_q_ratio = 0.0  # For asymmetric proposals
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
            
            if move_type == 'sigma' and pair_type is not None:
                # Calculate asymmetric proposal ratio for sigma moves
                #old_sigma = self.sigma[pair_type]
                #new_sigma = proposed_sigma[pair_type]
                #log_q_ratio = np.log(new_sigma) - np.log(old_sigma)
                log_q_ratio = 0.0
            else:
                log_q_ratio = 0.0

            # Calculate new score
            proposed_score, prop_excl, prop_pair, prop_prior = self.calculate_score(
                proposed_positions, proposed_sigma, self.sigma_range,
                set()
            )
            
            # Metropolis criterion with Jacobian correction for sigma moves
            delta = proposed_score - current_score
            log_accept = -delta / temp + log_q_ratio
            
            # Acceptance probability
            if log_accept >=0.0 or np.log(np.random.random()) < log_accept:
                accept = True
            else:
                accept = False
            #accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
            
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
    # Example usage with sampler sequence
    sampler_sequence = ['PairSampler', 'TetramerSampler', 'OctetSampler']
    sequence_idx = 0  # First sampler in sequence
    
    pair_samp = PairSampler(
        sampler_sequence=sampler_sequence,
        sequence_idx=sequence_idx,
        #prior_type="jeffreys",  # Use Jeffreys prior for first sampler
        prior_type="gamma",  # Use Jeffreys prior for first sampler
    )
    pair_samp.run_mc(n_steps=10000, save_freq=100)
#----------------------------------------------------------------------
