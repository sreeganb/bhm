# core/movers.py
#==================================================================================
# Particle move proposals for MCMC sampling
#==================================================================================
import numpy as np
from scipy.spatial.distance import cdist
from typing import Dict, List, Tuple
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
    base_step_min = 1.5
    base_step_max = 6.0
    # Reference radius for scaling (e.g., radius of type 'A' or default)
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
    log_step_size = np.random.uniform(0.01, 0.1)  # Tunable parameter
    
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

def get_tetramers(state) -> List[Tuple[int, ...]]:
    """Use Hungarian algorithm for optimal A-B matching, then greedy C selection."""
    try:
        from scipy.optimize import linear_sum_assignment
        
        positions = state.positions
        params = SystemParameters()
        
        # Validate input
        if not all(k in positions and len(positions[k]) > 0 for k in ['A', 'B', 'C']) or len(positions['C']) < 2:
            return []
            
        a_pos, b_pos, c_pos = positions['A'], positions['B'], positions['C']
        
        # Get target distances
        ab_target = params.pair_distances['AB']
        bc_target = params.pair_distances['BC']
        
        # Calculate AB cost matrix (deviation from target distance)
        dist_AB = cdist(a_pos, b_pos)
        cost_matrix = np.abs(dist_AB - ab_target)
        
        # Solve optimal assignment problem
        a_indices, b_indices = linear_sum_assignment(cost_matrix)
        
        # Pre-calculate BC distances
        dist_BC = cdist(b_pos, c_pos)
        
        # Now assign C particles greedily based on B assignments
        c_used = set()
        tetramers = []
        
        # Sort A-B pairs by their cost (best matches first)
        pair_costs = cost_matrix[a_indices, b_indices]
        sorted_pairs = np.argsort(pair_costs)
        
        for pair_idx in sorted_pairs:
            a_idx = a_indices[pair_idx]
            b_idx = b_indices[pair_idx]
            
            # Find available C particles
            available_c = [i for i in range(len(c_pos)) if i not in c_used]
            if len(available_c) < 2:
                break
            
            # Get best C pair for this B
            bc_dists = dist_BC[b_idx, available_c]
            c_scores = np.abs(bc_dists - bc_target)
            
            best_c_local = np.argsort(c_scores)[:2]
            best_c_indices = [available_c[i] for i in best_c_local]
            
            # Form tetramer
            tetramers.append((a_idx, b_idx, best_c_indices[0], best_c_indices[1]))
            c_used.update(best_c_indices)
            
            if len(tetramers) >= min(len(a_pos), len(b_pos), len(c_pos) // 2):
                break
        
        return tetramers
        
    except Exception as e:
        print(f"Error in Hungarian tetramer generation: {e}")
        return []

def propose_tetramer_move(state):
    """
    Tetramer move proposal with translation or rotation.
    60% probability for translation, 40% for rotation.
    Applies move in-place if particles stay within box boundaries.
    
    Args:
        state: SystemState object (modified in-place)
    """
    tetramers = get_tetramers(state)
    
    if not tetramers:
        return
    
    # Get box boundaries
    box_size = state.box_size
    
    # Select random tetramer
    tetramer = tetramers[np.random.randint(len(tetramers))]
    a_idx, b_idx, c_idx1, c_idx2 = tetramer
    
    # Extract tetramer particle information
    particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]
    
    # Get current coordinates
    coords = np.array([state.positions[part][idx] for part, idx in particles])
    centroid = np.mean(coords, axis=0)
    
    # Fixed step sizes
    trans_step = 2.0
    rot_step = 0.15
    
    # Choose move type: 60% translation, 40% rotation
    if np.random.random() < 0.6:
        # --- TRANSLATION MOVE ---
        displacement = np.random.normal(0.0, trans_step, 3)
        new_coords = coords + displacement
        
    else:
        # --- ROTATION MOVE ---
        # Generate random rotation axis using Marsaglia method
        while True:
            x1, x2 = np.random.uniform(-1, 1, 2)
            if x1*x1 + x2*x2 < 1:
                break
        
        sqrt_term = np.sqrt(1 - x1*x1 - x2*x2)
        rotation_axis = np.array([2*x1*sqrt_term, 2*x2*sqrt_term, 1 - 2*(x1*x1 + x2*x2)])
        
        # Generate rotation angle
        rotation_angle = np.random.normal(0.0, rot_step)
        
        # Build rotation matrix using quaternion
        half_angle = rotation_angle / 2.0
        qw = np.cos(half_angle)
        qx = rotation_axis[0] * np.sin(half_angle)
        qy = rotation_axis[1] * np.sin(half_angle)
        qz = rotation_axis[2] * np.sin(half_angle)
        
        rot_matrix = np.array([
            [1 - 2*(qy**2 + qz**2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
            [2*(qx*qy + qw*qz), 1 - 2*(qx**2 + qz**2), 2*(qy*qz - qw*qx)],
            [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx**2 + qy**2)]
        ])
        
        # Apply rotation around centroid
        new_coords = centroid + (rot_matrix @ (coords - centroid).T).T
    
    # Check if all particles stay within box bounds
    if np.all(new_coords >= 0) and np.all(new_coords <= box_size):
        # Apply move in-place
        for i, (part, idx) in enumerate(particles):
            state.positions[part][idx] = new_coords[i]