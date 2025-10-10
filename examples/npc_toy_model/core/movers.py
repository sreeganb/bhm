# core/movers.py
#==================================================================================
# Particle move proposals for MCMC sampling
#==================================================================================
import numpy as np
from core.parameters import SystemParameters

def propose_particle_move(state, use_pbc: bool = False):
    """
    Propose a random displacement for one randomly selected particle.
    
    Standard MCMC particle move (single-particle Metropolis):
    1. Pick ONE particle uniformly from ALL particles (not per-type sampling)
    2. Propose displacement uniformly in [-δ, +δ]³ where δ is step size
    3. Apply periodic boundary conditions (PBC) or reflective walls
    
    Args:
        state: SystemState object (modified in-place)
        use_pbc: If True, use periodic boundaries; else use reflective walls
        
    Returns:
        bool: True if move applied, False if rejected (only for hard violations)
    """
    params = SystemParameters()
    
    # Select ONE particle uniformly across ALL types
    types = list(state.positions.keys())
    counts = [state.positions[t].shape[0] for t in types]
    total_particles = sum(counts)
    
    # Weighted random choice by particle count
    type_weights = np.array(counts, dtype=float) / total_particles
    particle_type = np.random.choice(types, p=type_weights)
    particle_idx = np.random.randint(state.positions[particle_type].shape[0])
    
    # Get particle radius for boundary check
    radius = params.radii.get(particle_type, 0.0)
    
    # Choose step size randomly from a range for adaptive exploration
    # Larger particles get smaller base step (inverse scaling with radius)
    base_step_min = 1.0
    base_step_max = 5.0
    ref_radius = params.radii.get('A', 24.0)
    scale_factor = radius / ref_radius
    
    # Random step size in range, scaled by particle size
    step_size = np.random.uniform(base_step_min, base_step_max) / scale_factor
    
    # Propose displacement: uniform in cube [-step, +step]³
    displacement = np.random.uniform(-step_size, step_size, size=3)
    
    # Calculate new position
    old_pos = state.positions[particle_type][particle_idx]
    new_pos = old_pos + displacement
    
    # Handle boundary conditions
    box_size = state.box_size
    
    if use_pbc:
        # Periodic boundary conditions: wrap around
        new_pos = new_pos % box_size
        state.positions[particle_type][particle_idx] = new_pos
        return True
    else:
        # Reflective walls with soft constraint
        # Instead of hard rejection, reflect particles back into box
        for dim in range(3):
            # Lower boundary
            if new_pos[dim] < radius:
                new_pos[dim] = 2 * radius - new_pos[dim]  # Reflect
            # Upper boundary
            if new_pos[dim] > box_size - radius:
                new_pos[dim] = 2 * (box_size - radius) - new_pos[dim]  # Reflect
        
        # Final sanity check: still reject if completely outside after reflection
        if (np.any(new_pos < 0) or np.any(new_pos > box_size)):
            return False
        
        # Accept the move
        state.positions[particle_type][particle_idx] = new_pos
        return True

def propose_sigma_move(state): 
    """
    Propose a change to the interaction range (sigma) for one randomly selected pair type.
    
    Standard approach for scale parameters (common in Bayesian MCMC):
    - Log-normal random walk: propose in log-space, convert back
    - Symmetric in log-space → satisfies detailed balance
    - Step size chosen from range for adaptive exploration
    
    Args:
        state: SystemState object (modified in-place)
        
    Returns:
        None (modifies state.sigma in-place)
    """
    if len(state.sigma) == 0:
        return
    
    # Choose which sigma parameter to update
    pair_type = np.random.choice(list(state.sigma.keys()))
    
    # Current value in log-space
    current_sigma = state.sigma[pair_type]
    log_sigma_current = np.log(current_sigma)
    
    # Random step size for this move (adaptive exploration)
    log_step_size = np.random.uniform(0.05, 0.2)  # Tunable parameter
    
    # Symmetric proposal: log(σ') = log(σ) + N(0, δ²)
    # This is symmetric in log-space: q(σ'|σ) = q(σ|σ')
    log_sigma_proposed = log_sigma_current + np.random.normal(0.0, log_step_size)
    
    # Convert back to linear space
    proposed_sigma = np.exp(log_sigma_proposed)
    
    # Optional: enforce bounds if state has sigma_range defined
    if hasattr(state, 'sigma_range') and pair_type in state.sigma_range:
        lower, upper = state.sigma_range[pair_type]
        if proposed_sigma < lower or proposed_sigma > upper:
            # Out of bounds - reject by not modifying state
            return
    
    # Apply move in-place
    state.sigma[pair_type] = proposed_sigma