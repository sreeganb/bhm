# samplers/pair.py
import numpy as np
from typing import Dict, Any, Optional, Callable, Tuple
from core.state import SystemState
from samplers.base import run_mcmc_sampling
from scoring.pair_score import PairNLL
from scoring.exvol_score import ExvolNLL
from core.movers import propose_particle_move, propose_sigma_move

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
    # Define proposal functions - direct from core/movers.py
    propose_fns = {
        'position': propose_particle_move,
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