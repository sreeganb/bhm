# samplers/pair.py
import numpy as np
from typing import Dict, Any, Optional, Callable, Tuple
from core.state import SystemState
from samplers.base import run_mcmc_sampling
from scoring.pair_score import PairNLL
from scoring.exvol_score import ExvolNLL

def propose_position_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """Propose position move and update state in-place"""
    # Adaptive step size based on acceptance rate
    step_size = 0.1 * min(2.0, max(0.1, acceptance_rate))
    
    # Choose particle type and index
    particle_type = np.random.choice(list(state.positions.keys()))
    if len(state.positions[particle_type]) > 0:
        idx = np.random.randint(len(state.positions[particle_type]))
        
        # Apply small displacement
        displacement = np.random.normal(0, step_size, 3)
        state.positions[particle_type][idx] += displacement
        
        # Apply periodic boundary conditions
        state.positions[particle_type][idx] %= state.box_size

def propose_sigma_move(state: SystemState, acceptance_rate: float = 0.5) -> None:
    """Propose sigma move and update state in-place"""
    # Adaptive step size based on acceptance rate
    step_factor = 1.0 + 0.1 * min(1.0, max(0.1, acceptance_rate))
    
    # Choose a sigma key to modify
    sigma_key = np.random.choice(list(state.sigma.keys()))
    
    # Apply log-normal perturbation
    if np.random.random() < 0.5:
        state.sigma[sigma_key] *= step_factor
    else:
        state.sigma[sigma_key] /= step_factor
    
    # Enforce bounds if sigma_range is defined
    if sigma_key in state.sigma_range:
        min_val, max_val = state.sigma_range[sigma_key]
        state.sigma[sigma_key] = max(min_val, min(max_val, state.sigma[sigma_key]))

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
