# samplers/base.py
import numpy as np
import os
import h5py
from typing import Dict, Tuple, Any, Optional, Callable
from core.state import SystemState
from core.sigma import GMMSigmaProvider

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
    with h5py.File(trajectory_file, 'w') as f:
        pass  # Create empty file
    
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
    
    # Convert move_probs dict to arrays for random choice
    move_types = list(move_probs.keys())
    move_weights = np.array([move_probs[m] for m in move_types])
    move_weights /= move_weights.sum()  # Normalize
    
    # Temperature schedule
    temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
    temp = temp_start
    
    print(f"Starting MCMC sampling for {n_steps} steps...")
    
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
        new_prior = sig_provider.calculate_negative_log_prior(proposed_state)  # Fixed: added sig_provider
        proposed_score, *prop_components = score_fn(proposed_state, new_prior)
        
        # Metropolis acceptance criterion
        delta = proposed_score - current_score
        accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
        
        if debug and step % 10 == 0:
            print(f"DEBUG: Move: {move_type}, Delta: {delta:.2f}, Temp: {temp:.2f}")
            print(f"DEBUG: Score components: {score_components}")
        
        if accept:
            # Accept the move
            state = proposed_state
            current_score = proposed_score
            score_components = prop_components
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
        
        # Save trajectory
        if step % save_freq == 0 or step == n_steps:
            with h5py.File(trajectory_file, 'a') as f:
                frame_group = f.create_group(f"frame_{step}")
                for component, positions in state.positions.items():
                    frame_group.create_dataset(component, data=positions)
                # Save sigma values
                sigma_group = frame_group.create_group("sigma")
                for key, value in state.sigma.items():
                    sigma_group.attrs[key] = value
                # Save score
                frame_group.attrs["score"] = current_score
            
            # Status update
            accept_rate = sum(accepts.values()) / max(1, sum(attempts.values()))
            print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={accept_rate:.2f}")
    
    # Print final acceptance rates
    print("\nSampling complete:")
    for move in move_types:
        rate = accepts[move] / max(1, attempts[move])
        print(f"- {move}: {rate:.2f} acceptance ({accepts[move]}/{attempts[move]})")
    
    return best_state, trajectory_file