import numpy as np
import os
import h5py
from typing import Dict, Tuple, Any, Optional, Callable
from core.state import SystemState
from core.sigma import calculate_negative_log_prior

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
    trajectory_file: Optional[str] = None,
    chain_number: int = 1,  # Added chain_number parameter
    debug: bool = False,
    **kwargs  # Add **kwargs to capture any other arguments
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
        trajectory_file: Path to save trajectory file (will be generated if None)
        chain_number: Chain identifier for multi-chain runs
        debug: Whether to print debug info
        **kwargs: Additional arguments not used directly
    
    Returns:
        Final state and path to trajectory file
    """
    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create trajectory file path if not specified
    if trajectory_file is None:
        trajectory_file = os.path.join(output_dir, f"trajectory_chain{chain_number}.h5")
    
    # Create empty file
    with h5py.File(trajectory_file, 'w') as f:
        f.create_group('trajectory')
    
    # Initialize tracking variables
    best_state = state.copy()
    best_score = float('inf')
    accepts = {move: 0 for move in propose_fn_dict}
    attempts = {move: 0 for move in propose_fn_dict}
    
    # Calculate initial score
    prior_penalty = calculate_negative_log_prior(state)
    current_score, *score_components = score_fn(state, prior_penalty)
    
    # Convert move_probs dict to arrays for random choice
    move_types = list(move_probs.keys())
    move_weights = np.array([move_probs[m] for m in move_types])
    move_weights = move_weights / move_weights.sum()  # Normalize
    
    # Temperature schedule
    temp_decay = (temp_end / temp_start) ** (1.0 / n_steps)
    temp = temp_start
    
    print(f"Starting MCMC sampling for {n_steps} steps (chain {chain_number})...")
    
    # Save initial state
    with h5py.File(trajectory_file, 'a') as f:
        traj_grp = f['trajectory']
        state_grp = traj_grp.create_group(f"state_{0:05d}")
        
        # Save state attributes
        state_grp.attrs["step"] = 0
        state_grp.attrs["total_score"] = float(current_score)
        if score_components:
            state_grp.attrs["prior_score"] = float(score_components[0]) if len(score_components) > 0 else 0.0
            state_grp.attrs["pair_score"] = float(score_components[1]) if len(score_components) > 1 else 0.0
            state_grp.attrs["exvol_score"] = float(score_components[2]) if len(score_components) > 2 else 0.0
            state_grp.attrs["tet_score"] = float(score_components[3]) if len(score_components) > 3 else 0.0
            state_grp.attrs["oct_score"] = float(score_components[4]) if len(score_components) > 4 else 0.0
        
        # Save sigma values
        sigma_grp = state_grp.create_group("sigma")
        for key, value in state.sigma.items():
            sigma_grp.attrs[key] = float(value)
        
        # Save positions
        pos_grp = state_grp.create_group("positions")
        for comp, coords in state.positions.items():
            pos_grp.create_dataset(comp, data=coords.astype(np.float32), compression="gzip")
    
    # Main MCMC loop
    for step in range(1, n_steps + 1):
        # Select move type
        move_type = np.random.choice(move_types, p=move_weights)
        attempts[move_type] += 1
        
        # Get proposal function for this move
        propose_fn = propose_fn_dict[move_type]
        
        # Make proposal
        proposed_state = state.copy()
        acceptance_rate = accepts[move_type] / max(1, attempts[move_type])
        propose_fn(proposed_state, acceptance_rate)
        
        # Calculate new score
        new_prior = calculate_negative_log_prior(proposed_state)
        proposed_score, *prop_components = score_fn(proposed_state, new_prior)
        
        # Metropolis acceptance criterion
        delta = proposed_score - current_score

        # Add inside the MCMC loop, after calculating proposed_score:
        if step % 10 == 0 or abs(delta) > 1000:
            print(f"\nStep {step} details:")
            print(f"Move type: {move_type}")
            print(f"Prior penalty: {new_prior:.2f}")
            print(f"Current score: {current_score:.2f}, Proposed score: {proposed_score:.2f}")
            if len(prop_components) >= 3:
                print(f"Components: Exclusion={prop_components[0]:.2f}, Pair={prop_components[1]:.2f}, Prior penalty={prop_components[2]:.2f}")
            
            if move_type == 'sigma':
                print("Sigma values:")
                for k, v in proposed_state.sigma.items():
                    print(f"  {k}: {v:.4f}")
            print("")
        
        if abs(delta) < 1e-6 and step%10 == 0:
            print(f"DEBUG: Step {step}: No change in score, delta = {delta:.2f}")
            
        accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
        
        if debug and step % 10 == 0:
            print(f"DEBUG: Move: {move_type}, Delta: {delta:.2f}, Temp: {temp:.2f}")
        
        if accept:
            # Accept the move
            state = proposed_state
            current_score = proposed_score
            score_components = prop_components
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
                traj_grp = f['trajectory']
                state_grp = traj_grp.create_group(f"state_{step:05d}")
                
                # Save state attributes
                state_grp.attrs["step"] = step
                state_grp.attrs["total_score"] = float(current_score)
                if score_components:
                    state_grp.attrs["prior_score"] = float(score_components[0]) if len(score_components) > 0 else 0.0
                    state_grp.attrs["pair_score"] = float(score_components[1]) if len(score_components) > 1 else 0.0
                    state_grp.attrs["exvol_score"] = float(score_components[2]) if len(score_components) > 2 else 0.0
                    state_grp.attrs["tet_score"] = float(score_components[3]) if len(score_components) > 3 else 0.0
                    state_grp.attrs["oct_score"] = float(score_components[4]) if len(score_components) > 4 else 0.0
                
                # Save sigma values
                sigma_grp = state_grp.create_group("sigma")
                for key, value in state.sigma.items():
                    sigma_grp.attrs[key] = float(value)
                
                # Save positions
                pos_grp = state_grp.create_group("positions")
                for comp, coords in state.positions.items():
                    pos_grp.create_dataset(comp, data=coords.astype(np.float32), compression="gzip")
            
            # Status update
            accept_rate = sum(accepts.values()) / max(1, sum(attempts.values()))
            print(f"Chain {chain_number}, Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, Accept={accept_rate:.2f}")
    
    # Print final acceptance rates
    print(f"\nChain {chain_number} sampling complete:")
    for move in move_types:
        rate = accepts[move] / max(1, attempts[move])
        print(f"- {move}: {rate:.2f} acceptance ({accepts[move]}/{attempts[move]})")
    
    return best_state, trajectory_file