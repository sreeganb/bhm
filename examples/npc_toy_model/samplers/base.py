# samplers/base.py
import numpy as np
import os
import h5py
from typing import Dict, Tuple, Any, Optional, Callable
from core.state import SystemState
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
    equilibration_steps: int = 500,
    adapt_step_sizes: Optional[Callable] = None,
    debug: bool = False
) -> Tuple[SystemState, str]:
    """
    Generic MCMC sampling function that can be used by all samplers
    
    Args:
        state: Initial system state (with sigma_prior already attached)
        score_fn: Function to calculate score (neg_log_posterior)
                 Should compute prior internally from state.sigma_prior
        propose_fn_dict: Dict of proposal functions for each move type
                        Each function should take (state) and modify it in-place
        move_probs: Dict of probabilities for each move type
        n_steps: Number of MCMC steps to run
        save_freq: How often to save trajectory frames
        output_dir: Directory to save output
        temp_start: Initial temperature
        temp_end: Final temperature
        equilibration_steps: Number of steps at high temperature before annealing
        adapt_step_sizes: Function to adapt step sizes based on acceptance
        debug: Whether to print debug info
    
    Returns:
        Final state and path to trajectory file
    """
    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)
    trajectory_file = os.path.join(output_dir, "trajectory.h5")
    
    if os.path.exists(trajectory_file):
        os.remove(trajectory_file)
    
    # Initialize tracking
    best_state = state.copy()
    best_score = float('inf')
    accepts = {move: 0 for move in propose_fn_dict}
    attempts = {move: 0 for move in propose_fn_dict}
    
    # Verify sigma prior is attached (should be done by pipeline)
    sigma_prior = getattr(state, 'sigma_prior', None)
    if sigma_prior is None:
        raise ValueError("state.sigma_prior not found - pipeline must attach it before sampling")

    # Calculate initial score (score_fn computes prior internally!)
    # The score_fn signature is: score_fn(state, prior_penalty=0.0)
    # But we pass prior_penalty as a dummy - the function computes it internally
    current_score, *score_components = score_fn(state, prior_penalty=0.0)
    
    # Unpack score components (order depends on sampler)
    # pair sampler: (total, exvol, pair, prior)
    # tetramer sampler: (total, exvol, pair, tetramer, prior)
    # octet sampler: (total, exvol, pair, tetramer, octet, prior)
    exvol_score = score_components[0] if len(score_components) > 0 else 0.0
    pair_score = score_components[1] if len(score_components) > 1 else 0.0
    tet_score = score_components[2] if len(score_components) > 2 else 0.0
    oct_score = score_components[3] if len(score_components) > 3 else 0.0
    prior_penalty = score_components[-1] if len(score_components) > 0 else 0.0  # Last component is always prior
    
    # Setup move selection
    move_types = list(move_probs.keys())
    move_weights = np.array([move_probs[m] for m in move_types])
    move_weights /= move_weights.sum()
    
    # Temperature schedule
    annealing_steps = n_steps - equilibration_steps
    temp_decay = (temp_end / temp_start) ** (1.0 / max(1, annealing_steps))
    temp = temp_start
    
    print(f"Starting MCMC sampling for {n_steps} steps...")
    print(f"  - Equilibration: {equilibration_steps} steps at T={temp_start:.2f}")
    print(f"  - Annealing: {annealing_steps} steps from T={temp_start:.2f} to T={temp_end:.2f}")
    print(f"  - Using {'GMM' if sigma_prior.use_gmm else sigma_prior.prior_type} prior for sigma")
    print(f"  - Initial score: {current_score:.2f} (prior penalty: {prior_penalty:.2f})")
    print(f"Output will be saved to: {trajectory_file}")
    
    # Main MCMC loop
    for step in range(1, n_steps + 1):
        # Select and apply move
        move_type = np.random.choice(move_types, p=move_weights)
        attempts[move_type] += 1
        
        # Create proposed state and apply move
        proposed_state = state.copy()
        propose_fn_dict[move_type](proposed_state)
        
        # CRITICAL: Copy sigma_prior reference to proposed state
        # (same prior object applies to both current and proposed)
        proposed_state.sigma_prior = sigma_prior
        
        # Calculate new score (score_fn computes prior internally from state.sigma_prior)
        # The function will evaluate: -log_prior(proposed_state.sigma)
        proposed_score, *prop_components = score_fn(proposed_state, prior_penalty=0.0)
        
        # Unpack proposed scores
        prop_exvol = prop_components[0] if len(prop_components) > 0 else 0.0
        prop_pair = prop_components[1] if len(prop_components) > 1 else 0.0
        prop_tet = prop_components[2] if len(prop_components) > 2 else 0.0
        prop_oct = prop_components[3] if len(prop_components) > 3 else 0.0
        new_prior = prop_components[-1] if len(prop_components) > 0 else 0.0
        
        # Metropolis-Hastings acceptance criterion
        # delta = -log(posterior_proposed) - (-log(posterior_current))
        # Accept if posterior_proposed > posterior_current (i.e., delta < 0)
        delta = proposed_score - current_score
        accept = delta < 0 or np.random.random() < np.exp(-delta / temp)
        
        if debug and step % 10 == 0:
            sigma_str = ', '.join([f'{k}={v:.2f}' for k, v in state.sigma.items()])
            print(f"Step {step}: Move={move_type}, Delta={delta:.2f}, T={temp:.2f}, "
                  f"Prior={prior_penalty:.2f}, Sigma=[{sigma_str}]")
            if move_type == 'sigma':
                prop_sigma_str = ', '.join([f'{k}={v:.2f}' for k, v in proposed_state.sigma.items()])
                print(f"  Proposed sigma: [{prop_sigma_str}], New prior: {new_prior:.2f}")
        
        if accept:
            state = proposed_state
            current_score = proposed_score
            exvol_score, pair_score = prop_exvol, prop_pair
            tet_score, oct_score = prop_tet, prop_oct
            prior_penalty = new_prior
            accepts[move_type] += 1
            
            if current_score < best_score:
                best_score = current_score
                best_state = state.copy()
        
        # Update temperature (simulated annealing)
        if step > equilibration_steps:
            temp = temp_start * (temp_decay ** (step - equilibration_steps))
        else:
            temp = temp_start
        
        # Adaptive step size tuning
        if adapt_step_sizes and step % 100 == 0:
            acceptance_rates = {k: accepts[k] / max(1, attempts[k]) for k in accepts}
            adapt_step_sizes(acceptance_rates)
        
        # Save trajectory to disk
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
                types=getattr(state, 'types', None),
                bead_numbers=getattr(state, 'bead_numbers', None),
                traj_file=trajectory_file
            )
            
            accept_rate = sum(accepts.values()) / max(1, sum(attempts.values()))
            sigma_str = ', '.join([f'{k}={v:.2f}' for k, v in state.sigma.items()][:3])
            print(f"Step {step}/{n_steps}: Score={current_score:.2f}, T={temp:.2f}, "
                  f"Accept={accept_rate:.2%}, Sigma=[{sigma_str}]")
    
    # Print final statistics
    print("\n" + "="*60)
    print("Sampling complete!")
    print("="*60)
    print(f"Best score: {best_score:.2f}")
    print(f"Final score: {current_score:.2f}")
    print(f"\nAcceptance rates:")
    for move in move_types:
        rate = accepts[move] / max(1, attempts[move])
        print(f"  {move:>12s}: {rate:>6.2%} ({accepts[move]}/{attempts[move]})")
    
    print(f"\nFinal sigma values:")
    for k, v in state.sigma.items():
        print(f"  {k}: {v:.3f}")
    
    print(f"\nTrajectory saved to: {trajectory_file}")
    print("="*60 + "\n")
    
    return best_state, trajectory_file