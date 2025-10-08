# samplers/pair.py
import numpy as np
from typing import Dict, Any, Optional, Callable, Tuple
from core.state import SystemState
from samplers.base import run_mcmc_sampling
from scoring.pair_score import PairNLL
from scoring.exvol_score import ExvolNLL

def propose_position_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """
    Single-particle Gaussian move with reflective walls.
    Updates state in-place.
    """
    box_half = state.box_size / 2.0
    
    # Get all particle types and counts
    type_names = list(state.positions.keys())
    counts = [len(state.positions[t]) for t in type_names]
    
    if sum(counts) == 0:
        return  # No particles to move
    
    # Flatten particle selection across all types
    cum_counts = np.cumsum(counts)
    flat_idx = np.random.randint(cum_counts[-1])
    
    # Find which type and local index
    type_name = None
    local_idx = 0
    for tn, count, cum in zip(type_names, counts, cum_counts):
        if flat_idx < cum:
            type_name = tn
            local_idx = flat_idx - (cum - count)
            break
    
    if type_name is None or len(state.positions[type_name]) == 0:
        return
    
    # Get radius for this particle type (from params if available)
    radius = getattr(state, 'radii', {}).get(type_name, 1.0)
    max_radius = max(getattr(state, 'radii', {type_name: 1.0}).values())
    
    # Step size scaled by particle size
    step_sigma = 2.0 * (max_radius / radius)
    
    # Generate proposal
    current_pos = state.positions[type_name][local_idx].copy()
    proposal = current_pos + np.random.normal(0.0, step_sigma, 3)
    
    # Apply reflective boundary conditions for each dimension
    for dim in range(3):
        while proposal[dim] > box_half or proposal[dim] < -box_half:
            if proposal[dim] > box_half:
                proposal[dim] = 2 * box_half - proposal[dim]
            else:
                proposal[dim] = -2 * box_half - proposal[dim]
    
    # Update state in-place
    state.positions[type_name][local_idx] = proposal

def propose_sigma_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """
    Non-adaptive Metropolis proposal that preserves detailed balance.
    Updates state in-place.
    
    - Selects a single pair_type uniformly at random.
    - Uses an additive Gaussian step in linear sigma with constant scale
      per parameter (independent of the current value/state).
    - Applies exact reflective boundary conditions on [low, high],
      which keeps proposals unbiased and symmetric within bounds.
    """
    if len(state.sigma) == 0:
        return
    
    # Choose which parameter to update (uniform)
    pair_type = np.random.choice(list(state.sigma.keys()))
    current_val = float(state.sigma[pair_type])
    
    # Bounds (strictly positive)
    if pair_type in state.sigma_range:
        low, high = state.sigma_range[pair_type]
    else:
        low, high = 1e-6, 20.0
    
    # Guard invalid current values by snapping into bounds
    if not np.isfinite(current_val) or current_val <= 0.0:
        current_val = np.clip((low + high) * 0.5 if np.isfinite(current_val) else (low + high) * 0.5, low, high)
    
    # Constant, state-independent proposal width (ensures symmetry q(x->y)=q(y->x))
    width = max(high - low, 1e-9)
    step_sd = 0.15 * width  # tune as needed; constant for this pair_type
    
    # Symmetric additive Gaussian proposal in sigma-space
    proposed = current_val + np.random.normal(0.0, step_sd)
    
    # Reflective boundary conditions preserve symmetry on [low, high]
    def reflect(x: float, a: float, b: float) -> float:
        w = b - a
        if w <= 0.0:
            return float(np.clip(x, a, b))
        # Repeated reflections until inside [a, b]
        while x < a or x > b:
            if x < a:
                x = a + (a - x)
            if x > b:
                x = b - (x - b)
        return float(x)
    
    proposed_val = reflect(proposed, low, high)
    
    # Update state in-place
    state.sigma[pair_type] = proposed_val

def neg_log_posterior(
    state: SystemState, 
    prior_penalty: float = 0.0,
    excluded_pairs: Optional[set] = None
) -> Tuple[float, float, float, float]:
    """Calculate negative log posterior for pair sampler"""
    # Calculate excluded volume contribution
    exs = ExvolNLL(state.positions, kappa=100.0)
    exclusion_score = exs.compute_score()
    
    # Calculate pairwise score, excluding specified pairs
    ps = PairNLL(state.positions, state.sigma)
    pair_score = ps.compute_score()
    
    # Total score
    total_score = exclusion_score + pair_score + prior_penalty
    
    return total_score, exclusion_score, pair_score, prior_penalty

def run_pair_sampling(
    state: SystemState,
    n_steps: int = 1000,
    output_dir: str = "output/pair_sampler",
    **kwargs
) -> Tuple[SystemState, str]:
    """Run pair-level MCMC sampling"""
    # Define proposal functions
    propose_fns = {
        'position': propose_position_move,
        'sigma': propose_sigma_move
    }
    
    # Define move probabilities
    move_probs = {
        'position': 0.7,
        'sigma': 0.3
    }
    
    # Run MCMC
    score_fn = lambda s, p: neg_log_posterior(s, p)
    
    return run_mcmc_sampling(
        state=state,
        score_fn=score_fn,
        propose_fn_dict=propose_fns,
        move_probs=move_probs,
        n_steps=n_steps,
        output_dir=output_dir,
        **kwargs
    )