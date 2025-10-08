# samplers/base.py
import numpy as np
import os
import h5py
from typing import Dict, Tuple, Any, Optional, Callable
from core.state import SystemState
from core.sigma import GMMSigmaProvider
from core.io_utils import save_state_to_disk

def run_mcmc_sampling(
    state: SystemState,
    score_fn: Callable,
    propose_fn_dict: Dict[str, Callable],
    move_probs: Dict[str, float],
    n_steps: int = 1000,
    save_freq: int = 100,
    output_dir: str = "output",
    temp_start: float = 10.0,
    temp_end: float = 1.0,
    adapt_step_sizes: Optional[Callable] = None,
    debug: bool = False
) -> Tuple[SystemState, str]:
    """
    Generic MCMC sampling function that can be used by all samplers
    
    Args:
        state: Initial system state
        score_fn: Function to calculate score (neg_log_posterior)
        propose_fn_dict: Dict of proposal functions for each move type
        move_probs: Dict of probabilities for each move type
        n_steps: Number of MCMC steps to run
        save_freq: How often to save trajectory frames
        output_dir: Directory to save output
        temp_start: Initial temperature
        temp_end: Final temperature
        adapt_step_sizes: Function to adapt step sizes based on acceptance
        debug: Whether to print debug info
    
    Returns:
        Final state and path to trajectory file
    """
    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)
    trajectory_file = os.path.join(output_dir, "trajectory.h5")
    
    # Remove existing trajectory file if it exists
    if os.path.exists(trajectory_file):
        os.remove(trajectory_file)
    
    # Initialize tracking variables
    best_state = state.copy()
    best_score = float('inf')
    accepts = {move: 0 for move in propose_fn_dict}
    attempts = {move: 0 for move in propose_fn_dict}
    
    # Initialize sigma provider
    sig_provider = GMMSigmaProvider()
    
    # Calculate initial score
    prior_penalty = sig_provider.calculate_negative_log_prior(state)
    current_score, *score_components = score_fn(state, prior_penalty)
    
    # Unpack score components (adjust based on your score_fn return signature)
    # Assuming score_fn returns: (total_score, exvol_score, pair_score)
    # Add more components as needed: tet_score, oct_score, etc.
    exvol_score = score_components[0] if len(score_components) > 0 else 0.0
    pair_score = score_components[1] if len(score_components) > 1 else 0.0
    tet_score = score_components[2] if len(score_components) > 2 else 0.0
    oct_score = score_components[3] if len(score_components) > 3 else 0.0
    
    # Convert move_probs dict to arrays for random choice
    move_types = list(move_probs.keys())
    move_weights = np.array([move_probs[m] for m in move_types])
    move_weights /= move_weights.sum()  # Normalize
    
    # Temperature schedule
    temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
    temp = temp_start
    
    print(f"Starting MCMC sampling for {n_steps} steps...")
    print(f"Output will be saved to: {trajectory_file}")
    
    # Main MCMC loop
    for step in range(1, n_steps + 1):
        # Select move type
        move_type = np.random.choice(move_types, p=move_weights)
        attempts[move_type] += 1
        
        # Get proposal function for this move
        propose_fn = propose_fn_dict[move_type]
        
        # Make proposal
        proposed_state = state.copy()
        propose_fn(proposed_state, attempts[move_type] / max(1, accepts[move_type]))
        
        # Calculate new score with prior
        new_prior = sig_provider.calculate_negative_log_prior(proposed_state)
        proposed_score, *prop_components = score_fn(proposed_state, new_prior)
        
        # Unpack proposed score components
        prop_exvol = prop_components[0] if len(prop_components) > 0 else 0.0
        prop_pair = prop_components[1] if len(prop_components) > 1 else 0.0
        prop_tet = prop_components[2] if len(prop_components) > 2 else 0.0
        prop_oct = prop_components[3] if len(prop_components) > 3 else 0.0
        
        # Metropolis acceptance criterion
        delta = proposed_score - current_score
        accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
        
        if debug and step % 10 == 0:
            print(f"DEBUG: Move: {move_type}, Delta: {delta:.2f}, Temp: {temp:.2f}")
            print(f"DEBUG: Current score: {current_score:.2f} (exvol={exvol_score:.2f}, pair={pair_score:.2f})")
            print(f"DEBUG: Proposed score: {proposed_score:.2f} (exvol={prop_exvol:.2f}, pair={prop_pair:.2f})")
        
        if accept:
            # Accept the move
            state = proposed_state
            current_score = proposed_score
            exvol_score = prop_exvol
            pair_score = prop_pair
            tet_score = prop_tet
            oct_score = prop_oct
            prior_penalty = new_prior
            accepts[move_type] += 1
            
            # Update best state if improved
            if current_score < best_score:
                best_score = current_score
                best_state = state.copy()
        
        # Update temperature
        temp = temp_start * (temp_decay ** step)
        
        # Adapt step sizes if provided
        if adapt_step_sizes and step % 100 == 0:
            acceptance_rates = {k: accepts[k] / max(1, attempts[k]) for k in accepts}
            adapt_step_sizes(acceptance_rates)
        
        # Save trajectory using io_utils
        if step % save_freq == 0 or step == n_steps:
            save_state_to_disk(
                step=step,
                positions=state.positions,
                sigmas=state.sigma,
                score=current_score,
                prior_score=prior_penalty,
                pair_score=pair_score,
                exvol_score=exvol_score,
                tet_score=tet_score,
                oct_score=oct_score,
                types=getattr(state, 'types', None),  # If state has these attributes
                bead_numbers=getattr(state, 'bead_numbers', None),
                traj_file=trajectory_file
            )
            
            # Status update
            accept_rate = sum(accepts.values()) / max(1, sum(attempts.values()))
            print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={accept_rate:.2f}")
    
    # Print final acceptance rates
    print("\nSampling complete:")
    for move in move_types:
        rate = accepts[move] / max(1, attempts[move])
        print(f"- {move}: {rate:.2f} acceptance ({accepts[move]}/{attempts[move]})")
    print(f"\nTrajectory saved to: {trajectory_file}")
    
    return best_state, trajectory_file