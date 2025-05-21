# samplers/pair.py
import numpy as np
from typing import Dict, Any, Optional, Callable, Tuple
from core.state import SystemState
from core.scoring import calculate_excluded_volume, calculate_pair_scores
from samplers.base import run_mcmc_sampling
from core.sigma import calculate_negative_log_prior

def propose_position_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """
    Propose a move for a randomly selected component's position using a symmetric,
    zero-mean Gaussian proposal with reflection at boundaries.
    Updates state in-place.
    
    Args:
        state: Current system state to modify
        acceptance_rate: Current acceptance rate to adapt step size
    """
    def reflect_in_box(coord, box_size):
        """More efficient reflection logic"""
        # First fold into [0, 2*box_size] range
        coord = coord % (2 * box_size)
        # Then reflect if in the [box_size, 2*box_size] range
        if coord > box_size:
            coord = 2 * box_size - coord
        return coord
    
    # Choose particle type and index
    particle_type = np.random.choice(list(state.positions.keys()))
    if len(state.positions[particle_type]) == 0:
        # Skip if there are no particles of this type
        return
    
    # Get a random particle index
    idx = np.random.randint(len(state.positions[particle_type]))
    
    # Adaptive step size based on acceptance rate
    # Use a more sophisticated adaptive scheme similar to original
    base_step = 0.5  # Base step size
    target_acceptance = 0.4  # Target acceptance rate
    
    # Adjust step size based on current acceptance rate
    adjustment = np.clip(1.0 + 5.0 * (acceptance_rate - target_acceptance), 0.1, 2.75)
    step_size = base_step * adjustment
    
    # Apply displacement using zero-mean Gaussian
    displacement = np.random.normal(0, step_size, 3)
    
    # Apply the displacement
    new_position = state.positions[particle_type][idx] + displacement
    
    # Use reflection boundary conditions instead of periodic
    for i in range(3):
        new_position[i] = reflect_in_box(new_position[i], state.box_size)
    
    # Update the position in place
    state.positions[particle_type][idx] = new_position

#def propose_sigma_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """
    Propose a move for one sigma parameter using symmetrical sampling in log-space.
    Updates state in-place.
    
    Args:
        state: Current system state to modify
        acceptance_rate: Current acceptance rate to adapt step size
    """
#    # Choose a sigma key to modify
#    sigma_key = np.random.choice(list(state.sigma.keys()))
#    
#    # Current value in log-space
#    log_current = np.log(state.sigma[sigma_key])
#    
#    # Adjust step size based on acceptance rate
#    target_acceptance = 0.4  # Target acceptance rate
#    base_step_size = 0.0001  # Base step size in log space
#    
#    # More sophisticated step adaptation from first implementation
#    step_factor = 1.0 + 3.0 * np.clip(acceptance_rate - target_acceptance, -0.1, 0.3)
#    step_size = base_step_size * step_factor
#    
#    # Generate proposal in log space (symmetric sampling)
#    log_proposed = log_current + np.random.normal(0, step_size)
#    
#    # Apply the new value (transformed back from log space)
#    state.sigma[sigma_key] = np.exp(log_proposed)
#    
#    # Enforce bounds if sigma_range is defined
#    if hasattr(state, 'sigma_range') and sigma_key in state.sigma_range:
#        min_val, max_val = state.sigma_range[sigma_key]
#        state.sigma[sigma_key] = max(min_val, min(max_val, state.sigma[sigma_key]))
def propose_sigma_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """
    Propose a move for one sigma parameter using symmetrical sampling in log-space.
    Updates state in-place with better exploration properties.
    """
    # Choose a sigma key to modify
    sigma_key = np.random.choice(list(state.sigma.keys()))
    
    # Get bounds for this sigma
    min_val, max_val = 0.5, 10.0  # Default tighter bounds
    if hasattr(state, 'sigma_range') and sigma_key in state.sigma_range:
        min_val, max_val = state.sigma_range[sigma_key]
    
    # Current value in log-space
    current_sigma = state.sigma[sigma_key]
    log_current = np.log(current_sigma)
    log_min = np.log(min_val)
    log_max = np.log(max_val)
    
    # Use a much larger step size for better exploration
    # This is key - previous step size of 0.001 was much too small
    target_acceptance = 0.25  # Target slightly lower for better exploration
    base_step_size = 0.01  # 100x larger than previous value
    
    # Adjust step size based on acceptance rate
    step_factor = 1.0 + 2.0 * np.clip(acceptance_rate - target_acceptance, -0.2, 0.2)
    step_size = base_step_size * step_factor
    
    # Add a slight downward bias to counteract log-space transformation asymmetry
    # This helps prevent the natural drift toward larger values
    bias = -0.01
    
    # Try up to 5 times to get a valid proposal within bounds
    for attempt in range(5):
        # Generate proposal in log space with slight downward bias
        log_proposed = log_current + bias + np.random.normal(0, step_size)
        
        # Check if within bounds in log space
        if log_proposed >= log_min and log_proposed <= log_max:
            # Valid proposal, transform back to linear space
            state.sigma[sigma_key] = np.exp(log_proposed)
            return
    
    # If we couldn't generate a valid proposal after attempts, make a small move
    # This ensures we don't waste too many steps with rejected boundary proposals
    center = (log_min + log_max) / 2
    # Move slightly toward center of allowed range
    log_proposed = log_current + 0.1 * (center - log_current)
    state.sigma[sigma_key] = np.exp(log_proposed)
    
def neg_log_posterior(
    state: SystemState,
    prior_type: str = 'uniform',
    excluded_pairs: Optional[set] = None,
    prior_penalty: Optional[float] = None
) -> Tuple[float, float, float, float]:
    """
    Calculate negative log posterior for pair sampler
    
    Args:
        state: Current system state
        prior_type: Type of prior to use ('uniform', 'gaussian', etc.)
        excluded_pairs: Set of particle pairs to exclude from scoring
        prior_penalty: Pre-calculated prior penalty (if None, will be calculated)
    """
    # Use provided prior penalty if available, otherwise calculate it
    print(f"prior_penalty from neg log post: {prior_penalty}")
    if prior_penalty is None:
        prior_penalty = calculate_negative_log_prior(state, prior_type=prior_type)
    
    # Calculate excluded volume contribution
    exclusion_score = calculate_excluded_volume(state.positions)
    
    # Calculate pairwise score, excluding specified pairs
    pair_score = calculate_pair_scores(
        state.positions,
        state.sigma,
        excluded_pairs=excluded_pairs or set()
    )
    
    # Total score
    total_score = exclusion_score + pair_score + prior_penalty
    
    return total_score, exclusion_score, pair_score, prior_penalty

def run_pair_sampling(
    state: SystemState, 
    n_steps: int = 1000,
    output_dir: str = "output/pair_sampler",
    save_freq: int = 100,
    trajectory_file: Optional[str] = None,
    chain_number: int = 1,
    **kwargs
) -> Tuple[SystemState, str]:
    """Run pair-level MCMC sampling"""
    
    # Remove system_state from kwargs if it exists
    if 'system_state' in kwargs:
        kwargs.pop('system_state')
    
    # Extract prior_type from kwargs, defaulting to 'uniform'
    prior_type = kwargs.pop('prior_type', 'uniform')
    
    propose_fns = {
        'position': propose_position_move,
        'sigma': propose_sigma_move
    }
    
    # Define move probabilities
    move_probs = {
        'position': 0.7,
        'sigma': 0.3
    }
    
    # Create a wrapper for the scoring function to handle prior_type
    def score_wrapper(s, prior_penalty):
        return neg_log_posterior(s, prior_type=prior_type, prior_penalty=prior_penalty)
    
    # Run MCMC
    return run_mcmc_sampling(
        state=state,
        score_fn=score_wrapper,
        propose_fn_dict=propose_fns,
        move_probs=move_probs,
        n_steps=n_steps,
        save_freq=save_freq,
        output_dir=output_dir,
        trajectory_file=trajectory_file,
        chain_number=chain_number,
        **kwargs
    )