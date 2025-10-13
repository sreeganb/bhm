# samplers/tetramer.py
from typing import Optional, Tuple
from core.state import SystemState
from samplers.base import run_mcmc_sampling
from scoring.pair_score import PairNLL
from scoring.tetramer_score import TetramerNLL
from scoring.exvol_score import ExvolNLL
from core.movers import propose_particle_move, propose_sigma_move, propose_tetramer_move

def neg_log_posterior(
    state: SystemState, 
    prior_penalty: float = 0.0,
    excluded_pairs: Optional[set] = None
) -> Tuple[float, float, float, float, float]:
    """Calculate negative log posterior for tetramer sampler"""
    
    # Calculate excluded volume contribution
    exs = ExvolNLL(state.positions, kappa=100.0)
    exclusion_score = exs.compute_score()
    
    # Calculate pairwise score (background pairs not in tetramers)
    ps = PairNLL(state.positions, state.sigma)
    pair_score = ps.compute_score()
    
    # Calculate tetrameric score (structured ABCC units)
    ts = TetramerNLL(state)
    tetramer_score = ts.compute_score()

    # Total score = excluded volume + pairs + tetramers + prior
    total_score = exclusion_score + pair_score + tetramer_score + prior_penalty
    
    return total_score, exclusion_score, pair_score, tetramer_score, prior_penalty

def run_tetramer_sampling(
    state: SystemState,
    n_steps: int = 1000,
    output_dir: str = "output/tetramer_sampler",
    **kwargs
) -> Tuple[SystemState, str]:
    """
    Run tetramer-level MCMC sampling.
    
    Move hierarchy:
    - 60% tetramer moves (rigid body translation/rotation)
    - 25% individual particle moves (fine adjustments)
    - 15% sigma moves (parameter updates)
    
    This emphasizes structured tetramer moves while allowing
    fine-tuning through individual particle and parameter moves.
    """
    
    # Define proposal functions
    propose_fns = {
        'tetramer': propose_tetramer_move,   # Rigid body moves
        'position': propose_particle_move,   # Individual particle moves
        'sigma': propose_sigma_move          # Parameter moves
    }
    
    # Define move probabilities - tetramer moves dominate
    move_probs = {
        'tetramer': 0.60,   # Dominant: structured rigid body moves
        'position': 0.25,   # Secondary: fine particle adjustments
        'sigma': 0.15       # Tertiary: parameter optimization
    }
    
    # Create score function
    score_fn = lambda s, p: neg_log_posterior(s, p)
    
    # Run MCMC with tetramer-focused sampling
    return run_mcmc_sampling(
        state=state,
        score_fn=score_fn,
        propose_fn_dict=propose_fns,
        move_probs=move_probs,
        n_steps=n_steps,
        output_dir=output_dir,
        **kwargs
    )

def get_tetramers(state) -> list:
    """
    Convenience function to get tetramers from state.
    Delegates to the movers module implementation.
    """
    from core.movers import get_tetramers as _get_tetramers
    return _get_tetramers(state)