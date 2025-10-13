# core/movers.py
import random
import numpy as np
from scipy.spatial.distance import cdist
from typing import Dict, List, Optional, Tuple
from core.parameters import SystemParameters

# --------------------------
# Helpers
# --------------------------
def _reflect_scalar(x: float, a: float, b: float) -> float:
    """Exact reflective boundary on [a, b] (preserves symmetry)."""
    if not np.isfinite(x):
        return float(np.clip(x, a, b))
    w = b - a
    if w <= 0.0:
        return float(np.clip(x, a, b))
    # Reflect repeatedly until within bounds
    while x < a or x > b:
        if x < a:
            x = a + (a - x)
        if x > b:
            x = b - (x - b)
    return float(x)

def _flat_index_choice(positions: Dict[str, np.ndarray]) -> Tuple[str, int]:
    """Choose a single particle uniformly over all particles across types."""
    types = list(positions.keys())
    counts = [positions[t].shape[0] for t in types]
    total = sum(counts)
    if total == 0:
        raise ValueError("No particles to move: positions are empty")
    # Draw a flat index then map to (type, local_idx)
    k = np.random.randint(total)
    acc = 0
    for t, c in zip(types, counts):
        if k < acc + c:
            return t, k - acc
        acc += c
    # Fallback (should not happen)
    return types[-1], counts[-1] - 1

# --------------------------
# Position proposal (additive Gaussian, reflective walls)
# --------------------------
def propose_particle_move(state, accept_rate: float = 0.5):
    """
    Single-particle Gaussian move with reflective walls on [0, box_size]^3.

    - Select one particle uniformly across all types.
    - Use type-dependent step scale based on radii (larger particles move less).
    - Reflect at the walls to preserve proposal symmetry.
    """
    params = SystemParameters()
    box_size = getattr(state, "box_size", getattr(params, "box_size", 600.0))

    # Select particle uniformly across all types
    ptype, local_idx = _flat_index_choice(state.positions)

    # Step size heuristic using radii
    radii = getattr(params, "radii", {})
    radius = radii.get(ptype, 1.0)
    max_radius = max(radii.values()) if len(radii) > 0 else 1.0
    # Larger radius -> smaller step. Gaussian proposal.
    step_sigma = 2.0 * (max_radius / max(radius, 1e-6))

    current = state.positions[ptype][local_idx]
    proposal = current + np.random.normal(0.0, step_sigma, size=3)

    # Reflective walls on [0, box_size]
    for d in range(3):
        proposal[d] = _reflect_scalar(float(proposal[d]), 0.0, float(box_size))

    state.positions[ptype][local_idx] = proposal
    return True

# --------------------------
# Sigma proposal (additive Gaussian in linear sigma, reflective)
# --------------------------
def propose_sigma_move(
    state,
    accept_rate: Optional[float] = None
):
    """
    Non-adaptive Metropolis proposal that preserves detailed balance.

    - Selects a single pair_type uniformly at random.
    - Uses an additive Gaussian step in linear sigma with constant scale
      per parameter (independent of the current value/state).
    - Applies exact reflective boundary conditions on [low, high].

    Modifies state.sigma in-place.
    """
    if not hasattr(state, "sigma") or not isinstance(state.sigma, dict) or len(state.sigma) == 0:
        return False

    # Choose parameter uniformly
    pair_type = random.choice(list(state.sigma.keys()))
    current_val = float(state.sigma[pair_type])

    # Bounds: from prior if present, else from state.sigma_range, else defaults
    low, high = 1e-6, 20.0
    if hasattr(state, "sigma_prior") and getattr(state.sigma_prior, "sigma_ranges", None):
        low, high = state.sigma_prior.sigma_ranges.get(pair_type, (low, high))
    elif hasattr(state, "sigma_range") and isinstance(state.sigma_range, dict):
        low, high = state.sigma_range.get(pair_type, (low, high))

    if not np.isfinite(current_val) or current_val <= 0.0:
        # Snap invalid current value inside the bounds
        current_val = float(np.clip(0.5 * (low + high), low, high))

    # Constant, state-independent proposal width for symmetry
    width = max(high - low, 1e-9)
    step_sd = 0.15 * width  # tune globally if needed

    # Symmetric additive Gaussian proposal in sigma-space
    proposed = current_val + np.random.normal(0.0, step_sd)

    # Reflect to [low, high]
    proposed_val = _reflect_scalar(proposed, float(low), float(high))

    # Apply
    state.sigma[pair_type] = float(proposed_val)
    return True

# --------------------------
# Tetramer helpers and moves (tetramer move matches your previous method)
# --------------------------
def get_tetramers(state) -> List[Tuple[int, ...]]:
    """Use Hungarian algorithm for optimal A-B matching, then greedy C selection."""
    try:
        from scipy.optimize import linear_sum_assignment

        positions = state.positions
        params = SystemParameters()

        if not all(k in positions and len(positions[k]) > 0 for k in ['A', 'B', 'C']) or len(positions['C']) < 2:
            return []

        a_pos, b_pos, c_pos = positions['A'], positions['B'], positions['C']

        ab_target = params.pair_distances['AB']
        bc_target = params.pair_distances['BC']

        dist_AB = cdist(a_pos, b_pos)
        cost_matrix = np.abs(dist_AB - ab_target)

        a_indices, b_indices = linear_sum_assignment(cost_matrix)

        dist_BC = cdist(b_pos, c_pos)

        c_used = set()
        tetramers = []

        pair_costs = cost_matrix[a_indices, b_indices]
        sorted_pairs = np.argsort(pair_costs)

        for pair_idx in sorted_pairs:
            a_idx = a_indices[pair_idx]
            b_idx = b_indices[pair_idx]

            available_c = [i for i in range(len(c_pos)) if i not in c_used]
            if len(available_c) < 2:
                break

            bc_dists = dist_BC[b_idx, available_c]
            c_scores = np.abs(bc_dists - bc_target)

            best_c_local = np.argsort(c_scores)[:2]
            best_c_indices = [available_c[i] for i in best_c_local]

            tetramers.append((a_idx, b_idx, best_c_indices[0], best_c_indices[1]))
            c_used.update(best_c_indices)

            if len(tetramers) >= min(len(a_pos), len(b_pos), len(c_pos) // 2):
                break

        return tetramers

    except Exception as e:
        print(f"Error in Hungarian tetramer generation: {e}")
        return []

def propose_tetramer_move(state, acceptance_rate: float = 0.5):
    """
    Optimized tetramer move proposal with decoupled translation/rotation.
    60% probability for translation, 30% for rotation, 10% for mixed moves.
    Clips coordinates to stay within box boundaries minus particle radius.

    Modifies state.positions in-place.
    """
    tetramers = get_tetramers(state)
    if not tetramers:
        return False

    params = SystemParameters()
    box_size = getattr(params, 'box_size', 800.0)
    half_box = box_size / 2.0

    # Choose random tetramer
    tetramer = tetramers[np.random.randint(len(tetramers))]
    a_idx, b_idx, c_idx1, c_idx2 = tetramer

    # Pre-extract
    particles = [('A', a_idx), ('B', b_idx), ('C', c_idx1), ('C', c_idx2)]
    coords = np.array([state.positions[part][idx] for part, idx in particles])
    centroid = np.mean(coords, axis=0)

    # Compute buffer based on particle radii (fallback to 5.0)
    radii = getattr(params, 'radii', {})
    tet_radii = [radii.get(p, getattr(params, 'particle_radius', 5.0)) for p, _ in particles]
    buffer_radius = max(tet_radii) if len(tet_radii) else getattr(params, 'particle_radius', 5.0)

    max_coord = half_box - 2.0 * buffer_radius
    min_coord = -max_coord

    # Adaptive size factor (bounded)
    distances_from_center = np.linalg.norm(coords - centroid, axis=1)
    tetramer_radius = float(np.max(distances_from_center)) if distances_from_center.size else 1.0
    size_factor = max(0.5, min(2.0, tetramer_radius))

    # Fixed step sizes
    trans_step = 0.25
    rot_step = 0.2

    # Choose move type: 60% translation, 30% rotation, 10% mixed
    r = np.random.random()

    if r < 0.6:
        # Translation
        displacement = np.random.normal(0.0, trans_step * size_factor, 3)
        for i, (part, idx) in enumerate(particles):
            final_pos = coords[i] + displacement
            final_pos = np.clip(final_pos, min_coord, max_coord)
            state.positions[part][idx] = final_pos.astype(float)

    elif r < 0.9:
        # Rotation (Marsaglia)
        while True:
            x1, x2 = np.random.uniform(-1, 1, 2)
            if x1 * x1 + x2 * x2 < 1:
                break
        sqrt_term = np.sqrt(1 - x1 * x1 - x2 * x2)
        rotation_axis = np.array([2 * x1 * sqrt_term, 2 * x2 * sqrt_term, 1 - 2 * (x1 * x1 + x2 * x2)])
        rotation_axis /= max(np.linalg.norm(rotation_axis), 1e-12)

        rotation_angle = np.random.normal(0.0, rot_step)
        half_angle = rotation_angle / 2.0
        qw = np.cos(half_angle)
        sin_half = np.sin(half_angle)
        qx, qy, qz = rotation_axis * sin_half

        rot_matrix = np.array([
            [1 - 2 * (qy ** 2 + qz ** 2), 2 * (qx * qy - qw * qz), 2 * (qx * qz + qw * qy)],
            [2 * (qx * qy + qw * qz), 1 - 2 * (qx ** 2 + qz ** 2), 2 * (qy * qz - qw * qx)],
            [2 * (qx * qz - qw * qy), 2 * (qy * qz + qw * qx), 1 - 2 * (qx ** 2 + qy ** 2)]
        ])

        for i, (part, idx) in enumerate(particles):
            vec = coords[i] - centroid
            rotated = rot_matrix @ vec
            final_pos = centroid + rotated
            final_pos = np.clip(final_pos, min_coord, max_coord)
            state.positions[part][idx] = final_pos.astype(float)

    else:
        # Mixed (translation + rotation) with reduced steps
        trans_scale = 0.7
        rot_scale = 0.7

        displacement = np.random.normal(0.0, trans_step * size_factor * trans_scale, 3)

        while True:
            x1, x2 = np.random.uniform(-1, 1, 2)
            if x1 * x1 + x2 * x2 < 1:
                break
        sqrt_term = np.sqrt(1 - x1 * x1 - x2 * x2)
        rotation_axis = np.array([2 * x1 * sqrt_term, 2 * x2 * sqrt_term, 1 - 2 * (x1 * x1 + x2 * x2)])
        rotation_axis /= max(np.linalg.norm(rotation_axis), 1e-12)

        rotation_angle = np.random.normal(0.0, rot_step * rot_scale)
        half_angle = rotation_angle / 2.0
        qw = np.cos(half_angle)
        sin_half = np.sin(half_angle)
        qx, qy, qz = rotation_axis * sin_half

        rot_matrix = np.array([
            [1 - 2 * (qy ** 2 + qz ** 2), 2 * (qx * qy - qw * qz), 2 * (qx * qz + qw * qy)],
            [2 * (qx * qy + qw * qz), 1 - 2 * (qx ** 2 + qz ** 2), 2 * (qy * qz - qw * qx)],
            [2 * (qx * qz - qw * qy), 2 * (qy * qz + qw * qx), 1 - 2 * (qx ** 2 + qy ** 2)]
        ])

        for i, (part, idx) in enumerate(particles):
            vec = coords[i] - centroid
            rotated = rot_matrix @ vec
            final_pos = centroid + rotated + displacement
            final_pos = np.clip(final_pos, min_coord, max_coord)
            state.positions[part][idx] = final_pos.astype(float)

    return True