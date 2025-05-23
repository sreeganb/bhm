#--------------------------------------------------------------------------------
# Scoring functions for pair, tetramer, and octamer interaction models
# These are optimized using vectorized operations and designed to be modular
#--------------------------------------------------------------------------------
import numpy as np
from typing import Dict, List, Tuple, Set, Optional, Any, Union
import numba as nb
from functools import lru_cache
from scipy.spatial.distance import cdist
import time

#************************************
# Imports from files here
#************************************
from core.parameters import SystemParameters

#--------------------------------------------------------------------------------
# Base scoring components
#--------------------------------------------------------------------------------
def calculate_excluded_volume(
    positions: Dict[str, np.ndarray],
    params: Optional[Dict[str, Any]] = None,
    debug: bool = True
) -> float:
    """
    Calculate excluded volume score with a quadratic penalty for overlapping particles.

    Args:
        positions: Dictionary mapping particle types to position arrays (shape: Nx3)
        params: Optional parameters dictionary containing 'radii'
        debug: Whether to print debug information and write to file

    Returns:
        Excluded volume score (higher for overlapping particles)
    """
    if params is None:
        params = SystemParameters()  # Assumes this function exists
    
    ex_score = 0.0
    radii = params.radii
    sigma = 0.1  # Penalty strength, matching excluded_volume_nll
    
    # Initialize debug file
    debug_fh = None
    if debug:
        import tempfile
        import os
        debug_dir = os.path.join(os.getcwd(), "debug")
        os.makedirs(debug_dir, exist_ok=True)
        debug_file = os.path.join(debug_dir, f"excluded_volume_debug_{int(time.time())}.csv")
        debug_fh = open(debug_file, "w")
        debug_fh.write("type1,type2,particle1_idx,particle2_idx,distance,min_distance,overlap,score_contribution\n")
        print(f"Excluded volume debug: writing to {debug_file}")
    
    # Track statistics for debugging
    total_pairs = 0
    overlapping_pairs = 0
    type_pair_stats = {}
    
    for type1, pos1 in positions.items():
        for type2, pos2 in positions.items():
            if type1 <= type2:
                min_dist = radii[type1] + radii[type2]
                distances = cdist(pos1, pos2)
                pair_key = f"{type1}-{type2}"
                
                if pair_key not in type_pair_stats:
                    type_pair_stats[pair_key] = {
                        "total": 0,
                        "overlapping": 0,
                        "score": 0.0
                    }
                
                # handle self-interactions
                if type1 == type2:
                    # use upper triangle to avoid double-counting
                    mask = np.tril(np.ones_like(distances), k=-1)
                    viol_mask = (distances < min_dist) & (mask > 0)
                else:
                    viol_mask = distances < min_dist
                
                # Count all pairs considered
                if type1 == type2:
                    # For same type, count only upper triangle pairs
                    num_pairs = (len(pos1) * (len(pos1) - 1)) // 2
                else:
                    num_pairs = len(pos1) * len(pos2)
                
                total_pairs += num_pairs
                type_pair_stats[pair_key]["total"] = num_pairs
                
                if np.any(viol_mask):
                    # Get indices of violating pairs
                    viol_indices = np.where(viol_mask)
                    overlaps = min_dist - distances[viol_mask]
                    scores = (overlaps**2) / (sigma**2)
                    pair_score = np.sum(scores)
                    ex_score += pair_score
                    
                    # Update debugging counters
                    overlapping_pairs += len(viol_indices[0])
                    type_pair_stats[pair_key]["overlapping"] = len(viol_indices[0])
                    type_pair_stats[pair_key]["score"] = pair_score
                    
                    # Write detailed information to debug file
                    if debug_fh:
                        for idx in range(len(viol_indices[0])):
                            i, j = viol_indices[0][idx], viol_indices[1][idx]
                            overlap = overlaps[idx]
                            score = scores[idx]
                            dist = distances[i, j]
                            debug_fh.write(f"{type1},{type2},{i},{j},{dist:.6f},{min_dist:.6f},{overlap:.6f},{score:.6f}\n")
    
    # Print summary if debugging
    if debug:
        print(f"\n===== EXCLUDED VOLUME SUMMARY =====")
        print(f"Total pairs evaluated: {total_pairs}")
        print(f"Overlapping pairs: {overlapping_pairs} ({100.0 * overlapping_pairs / max(1, total_pairs):.2f}%)")
        print(f"Total excluded volume score: {ex_score:.6f}")
        print("\nPer-type statistics:")
        for pair_key, stats in type_pair_stats.items():
            if stats["total"] > 0:
                overlap_pct = 100.0 * stats["overlapping"] / stats["total"]
                print(f"  {pair_key}: {stats['overlapping']}/{stats['total']} overlaps ({overlap_pct:.2f}%), score: {stats['score']:.6f}")
        
        # Close the debug file
        if debug_fh:
            debug_fh.close()
    
    return ex_score

def calculate_pair_scores_matrix(
    pos1: np.ndarray,
    pos2: np.ndarray,
    target_dist: float,
    sigma: float
) -> np.ndarray:
    """
    Calculate pairwise score matrix between two sets of positions.
    
    Args:
        pos1: Positions of first particle type (shape: Nx3)
        pos2: Positions of second particle type (shape: Mx3)
        target_dist: Target distance between particles
        sigma: Standard deviation for the distance
        
    Returns:
        Score matrix of shape (N, M)
    """
    # Calculate distance matrix
    delta = pos1[:, np.newaxis, :] - pos2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(delta**2, axis=2))
    
    # Calculate Gaussian log-likelihood
    score_matrix = ((distances - target_dist)**2) / (2 * sigma**2) + np.log(2 * np.pi * sigma**2)
    
    return score_matrix, distances

def calculate_pair_scores(
    positions: Dict[str, np.ndarray],
    sigma: Dict[str, float],
    excluded_pairs: Optional[Set[Tuple]] = None,
    params: Optional[Dict[str, Any]] = None,
    debug: bool = True
) -> float:
    """
    Calculate pairwise interaction scores with optional excluded pairs.
    
    Args:
        positions: Dictionary mapping particle types to position arrays
        sigma: Dictionary mapping pair types to sigma values
        excluded_pairs: Set of (type1, idx1, type2, idx2) tuples to exclude
        params: Optional parameters dictionary
        debug: Whether to print debug information
        
    Returns:
        Total pairwise score
    """
    if params is None:
        params = SystemParameters()  # Get default parameters
    
    # Get pair distances from parameters
    pair_distances = params.pair_distances
    
    # Default pair types to consider
    pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C')]
    
    # Initialize score and excluded pairs set
    total_score = 0.0
    if excluded_pairs is None:
        excluded_pairs = set()
    
    # Debug file handling
    debug_fh = None
    if debug:
        import tempfile
        debug_file = tempfile.NamedTemporaryFile(delete=False, prefix="pair_scores_", suffix=".csv").name
        debug_fh = open(debug_file, "w")
        debug_fh.write("pair_key,type1,idx1,type2,idx2,distance,target,sigma,score\n")
    
    # Process each pair type
    for type1, type2 in pair_types:
        pair_key = f"{type1}{type2}"
        
        # Skip if particle types don't exist or we don't have a target distance
        if (type1 not in positions or type2 not in positions or 
            len(positions[type1]) == 0 or len(positions[type2]) == 0 or
            pair_key not in pair_distances):
            continue
        
        # Get target distance and sigma
        target_dist = pair_distances[pair_key]
        sigma_value = sigma[pair_key]
        
        # Calculate score matrix and distances
        score_matrix, distances = calculate_pair_scores_matrix(
            positions[type1],
            positions[type2],
            target_dist,
            sigma_value
        )
        
        # Process differently based on whether types are the same
        if type1 == type2:
            # For same type, handle special case
            np.fill_diagonal(score_matrix, np.inf)  # Exclude self-interactions
            m, n = score_matrix.shape
            selected_scores = []
            selected_indices = []
            
            # For each row i, find minimum element in columns j > i
            for i in range(m):
                if i < n-1:  # Check if there are columns to the right
                    row_segment = score_matrix[i, i+1:]
                    # Skip if all values are infinity
                    if np.all(np.isinf(row_segment)):
                        continue
                    
                    # Find minimum element
                    j_offset = np.argmin(row_segment)
                    j = i + 1 + j_offset
                    
                    # Check if this pair is excluded
                    if ((type1, i, type2, j) not in excluded_pairs and
                        (type2, j, type1, i) not in excluded_pairs):
                        selected_scores.append(score_matrix[i, j])
                        selected_indices.append((i, j))
            
            # Debug output
            if debug_fh:
                for i, j in selected_indices:
                    pair_score = score_matrix[i, j]
                    dist = distances[i, j]
                    debug_fh.write(f"{pair_key},{type1},{i},{type2},{j},{dist:.1f},{target_dist:.1f},{sigma_value:.1f},{pair_score:.1f}\n")

        else:
            # For different types, find mutual best matches
            row_min_indices = np.argmin(score_matrix, axis=1)
            col_min_indices = np.argmin(score_matrix, axis=0)
            
            # Collect unique index pairs
            unique_pairs = set()
            for i, j in enumerate(row_min_indices):
                if ((type1, i, type2, j) not in excluded_pairs and
                    (type2, j, type1, i) not in excluded_pairs):
                    unique_pairs.add((i, j))
                    
            for j, i in enumerate(col_min_indices):
                if ((type1, i, type2, j) not in excluded_pairs and
                    (type2, j, type1, i) not in excluded_pairs):
                    unique_pairs.add((i, j))
            
            # Get scores for unique index pairs
            selected_scores = [score_matrix[i, j] for i, j in unique_pairs 
                               if i < score_matrix.shape[0] and j < score_matrix.shape[1]]
            
            # Debug output
            if debug_fh:
                for i, j in unique_pairs:
                    if i < score_matrix.shape[0] and j < score_matrix.shape[1]:
                        pair_score = score_matrix[i, j]
                        dist = distances[i, j]
                        debug_fh.write(f"{pair_key},{type1},{i},{type2},{j},{dist:.1f},{target_dist:.1f},{sigma_value:.1f},{pair_score:.1f}\n")

        # Sum selected scores and add to total
        if selected_scores:
            pairwise_score_part = np.sum(selected_scores)
            total_score += pairwise_score_part
    
    # Clean up debug file if opened
    if debug_fh:
        debug_fh.close()
        print(f"Particle pairing debug information written to {debug_file}")
    
    return total_score

#--------------------------------------------------------------------------------
# Tetramer scoring functions
#--------------------------------------------------------------------------------
def get_tetramers(
    positions: Dict[str, np.ndarray],
    params: Optional[Dict[str, Any]] = None,
    cutoff_factor: float = 1.5
) -> List[Tuple[int, int, int, int]]:
    """
    Identify tetramers in the system based on spatial proximity.
    Tetramers consist of A-B-C-C units where B connects to two C particles.
    
    Args:
        positions: Dictionary mapping particle types to position arrays
        params: Optional parameters dictionary
        cutoff_factor: Distance cutoff factor relative to target distances
        
    Returns:
        List of tetramers as (a_idx, b_idx, c1_idx, c2_idx) tuples
    """
    if params is None:
        params = SystemParameters()  # Assumes this function exists
    
    pair_distances = params.pair_distances
    
    # Verify required particle types exist
    for particle_type in ['A', 'B', 'C']:
        if particle_type not in positions or len(positions[particle_type]) == 0:
            return []  # Cannot form tetramers if missing any component
    
    # Calculate cutoff distances
    ab_cutoff = cutoff_factor * pair_distances['AB']
    bc_cutoff = cutoff_factor * pair_distances['BC']
    cc_cutoff = cutoff_factor * pair_distances['CC']
    
    # Calculate distance matrices
    ab_distances = calculate_pairwise_distances(
        positions, [('A', 'B')]
    )[('A', 'B')]
    
    bc_distances = calculate_pairwise_distances(
        positions, [('B', 'C')]
    )[('B', 'C')]
    
    cc_distances = calculate_pairwise_distances(
        positions, [('C', 'C')]
    )[('C', 'C')]
    
    # Identify potential AB connections
    ab_connections = []
    for a_idx in range(len(positions['A'])):
        # Get indices of B particles within cutoff
        b_indices = np.where(ab_distances[a_idx] < ab_cutoff)[0]
        for b_idx in b_indices:
            ab_connections.append((a_idx, b_idx))
    
    # Find tetramers
    tetramers = []
    
    # For each AB connection, find C pairs connected to B
    for a_idx, b_idx in ab_connections:
        # Get C particles connected to this B
        c_indices = np.where(bc_distances[b_idx] < bc_cutoff)[0]
        
        # Check all pairs of C particles
        for i, c1_idx in enumerate(c_indices):
            for c2_idx in c_indices[i+1:]:
                # Verify C-C distance
                c_dist = cc_distances[c1_idx, c2_idx]
                if c_dist < cc_cutoff:
                    tetramers.append((a_idx, b_idx, c1_idx, c2_idx))
    
    return tetramers

def calculate_tetramer_scores(
    positions: Dict[str, np.ndarray],
    tetramers: List[Tuple[int, int, int, int]],
    sigma: Dict[str, float],
    params: Optional[Dict[str, Any]] = None,
    debug: bool = False
) -> np.ndarray:
    """
    Calculate scores for all tetramers with vectorized operations.
    
    Args:
        positions: Dictionary mapping particle types to position arrays
        tetramers: List of (a_idx, b_idx, c1_idx, c2_idx) tuples
        sigma: Dictionary mapping pair types to sigma values
        params: Optional parameters dictionary
        debug: Whether to print debug information
        
    Returns:
        Array of scores, one per tetramer
    """
    if not tetramers:
        return np.array([], dtype=np.float32)
    
    if params is None:
        params = SystemParameters()  # Assumes this function exists
        
    pair_distances = params.pair_distances
    
    n_tetramers = len(tetramers)
    
    # Extract all tetramer indices efficiently
    a_indices = np.array([t[0] for t in tetramers], dtype=np.int32)
    b_indices = np.array([t[1] for t in tetramers], dtype=np.int32)
    c1_indices = np.array([t[2] for t in tetramers], dtype=np.int32)
    c2_indices = np.array([t[3] for t in tetramers], dtype=np.int32)
    
    # Get all positions in single vectorized operations
    pos_a = positions['A'][a_indices]
    pos_b = positions['B'][b_indices]
    pos_c1 = positions['C'][c1_indices]
    pos_c2 = positions['C'][c2_indices]
    
    # Calculate all distances, handling periodic boundary conditions
    ab_delta = pos_a - pos_b
    bc1_delta = pos_b - pos_c1
    bc2_delta = pos_b - pos_c2
    cc_delta = pos_c1 - pos_c2
        
    # Calculate all distances at once
    ab_dists = np.sqrt(np.sum(ab_delta**2, axis=1))
    bc1_dists = np.sqrt(np.sum(bc1_delta**2, axis=1))
    bc2_dists = np.sqrt(np.sum(bc2_delta**2, axis=1))
    cc_dists = np.sqrt(np.sum(cc_delta**2, axis=1))
    
    # Cache target distances for performance
    ab_target = pair_distances['AB']
    bc_target = pair_distances['BC']
    cc_target = pair_distances['CC']
    
    # Calculate scores individually for detailed logging
    ab_scores = ((ab_dists - ab_target)**2) / (2 * sigma['AB']**2) + np.log(2 * np.pi * sigma['AB'])
    bc1_scores = ((bc1_dists - bc_target)**2) / (2 * sigma['BC']**2) + np.log(2 * np.pi * sigma['BC'])
    bc2_scores = ((bc2_dists - bc_target)**2) / (2 * sigma['BC']**2) + np.log(2 * np.pi * sigma['BC'])
    cc_scores = ((cc_dists - cc_target)**2) / (2 * sigma['CC']**2) + np.log(2 * np.pi * sigma['CC'])
    
    # Calculate total scores
    scores = ab_scores + bc1_scores + bc2_scores + cc_scores
    
    # Debug output if requested
    if debug:
        for i in range(min(5, n_tetramers)):  # Show first 5 for debugging
            print(f"Tetramer {i}: Score={scores[i]:.2f}")
            print(f"  A-B: {ab_dists[i]:.2f} vs {ab_target:.2f}, Score: {ab_scores[i]:.2f}")
            print(f"  B-C1: {bc1_dists[i]:.2f} vs {bc_target:.2f}, Score: {bc1_scores[i]:.2f}")
            print(f"  B-C2: {bc2_dists[i]:.2f} vs {bc_target:.2f}, Score: {bc2_scores[i]:.2f}")
            print(f"  C-C: {cc_dists[i]:.2f} vs {cc_target:.2f}, Score: {cc_scores[i]:.2f}")
    
    return scores

def get_tetramer_pairs(
    tetramers: List[Tuple[int, int, int, int]]
) -> Set[Tuple[str, int, str, int]]:
    """
    Get the set of pairs that are part of tetramers, for exclusion from pair scoring.
    
    Args:
        tetramers: List of (a_idx, b_idx, c1_idx, c2_idx) tuples
        
    Returns:
        Set of (type1, idx1, type2, idx2) tuples representing pairs in tetramers
    """
    if not tetramers:
        return set()
    
    tetramer_pairs = set()
    
    # Process all tetramers
    for a_idx, b_idx, c1_idx, c2_idx in tetramers:
        # Add pairs to excluded set (these are scored separately as tetramer pairs)
        tetramer_pairs.add(('A', a_idx, 'B', b_idx))
        tetramer_pairs.add(('B', b_idx, 'A', a_idx))  # Add reverse pair too
        
        tetramer_pairs.add(('B', b_idx, 'C', c1_idx))
        tetramer_pairs.add(('C', c1_idx, 'B', b_idx))
        
        tetramer_pairs.add(('B', b_idx, 'C', c2_idx))
        tetramer_pairs.add(('C', c2_idx, 'B', b_idx))
        
        tetramer_pairs.add(('C', c1_idx, 'C', c2_idx))
        tetramer_pairs.add(('C', c2_idx, 'C', c1_idx))
    
    return tetramer_pairs

#--------------------------------------------------------------------------------
# Octamer scoring functions 
#--------------------------------------------------------------------------------
def get_octets(
    positions: Dict[str, np.ndarray],
    tetramers: Optional[List[Tuple[int, int, int, int]]] = None,
    params: Optional[Dict[str, Any]] = None,
    temperature: float = 0.9,
    distance_cutoff_factor: float = 2.0
) -> List[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]]]:
    """
    Group tetramers into octets (pairs of tetramers) with temperature-based selection.
    
    Args:
        positions: Dictionary mapping particle types to position arrays
        tetramers: List of tetramers, if already calculated
        params: Optional parameters dictionary
        temperature: Temperature parameter for probabilistic selection (0-1)
        distance_cutoff_factor: How far tetramers can be to form an octet
        
    Returns:
        List of (tetramer1, tetramer2) pairs forming octets
    """
    if params is None:
        params = SystemParameters()
    
    # Calculate tetramers if not provided
    if tetramers is None:
        tetramers = get_tetramers(positions, params)
    
    if len(tetramers) < 2:
        return []
    
    # Calculate tetramer centers using vectorized operations
    centers = np.zeros((len(tetramers), 3))
    for i, (a_idx, b_idx, c_idx1, c_idx2) in enumerate(tetramers):
        coords = np.vstack([
            positions['A'][a_idx],
            positions['B'][b_idx],
            positions['C'][c_idx1],
            positions['C'][c_idx2]
        ])
        centers[i] = np.mean(coords, axis=0)
    
    # Calculate typical distance scale from parameters
    distance_scale = params['pair_distances']['AB']
    max_distance = distance_cutoff_factor * distance_scale
    
    # Form octets by pairing tetramers
    octets = []
    available = list(range(len(tetramers)))
    
    while len(available) >= 2:
        # Pick first tetramer randomly
        idx1 = np.random.choice(available)
        available.remove(idx1)
        
        # Calculate distances to all other tetramers
        diffs = centers[available] - centers[idx1]
        
        # Calculate distances and sort
        distances = np.sqrt(np.sum(diffs**2, axis=1))
        sorted_indices = np.argsort(distances)
        
        # Filter by maximum distance
        valid_indices = [i for i, d in enumerate(distances) if d <= max_distance]
        
        if not valid_indices:
            # No valid partners for this tetramer
            continue
            
        # Apply temperature-based selection
        if temperature > 0:
            # Calculate selection probabilities (closer tetramers more likely)
            valid_distances = distances[valid_indices]
            probs = np.exp(-valid_distances / (temperature * distance_scale))
            probs /= probs.sum()
            
            # Select partner based on probabilities
            selected_idx = np.random.choice(valid_indices, p=probs)
        else:
            # At zero temperature, always pick closest
            selected_idx = sorted_indices[0] if distances[sorted_indices[0]] <= max_distance else None
            
        if selected_idx is not None:
            # Get the actual index in available list
            idx2 = available[selected_idx]
            available.remove(idx2)
            
            # Form the octet
            octets.append((tetramers[idx1], tetramers[idx2]))
    
    return octets

def calculate_octet_scores(
    positions: Dict[str, np.ndarray],
    octets: List[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]]],
    sigma: Dict[str, float],
    params: Optional[Dict[str, Any]] = None,
    debug: bool = False
) -> np.ndarray:
    """
    Calculate scores for octets, treating each component tetramer individually.
    
    Args:
        positions: Dictionary mapping particle types to position arrays
        octets: List of (tetramer1, tetramer2) pairs
        sigma: Dictionary mapping pair types to sigma values
        params: Optional parameters dictionary
        debug: Whether to print debug information
        
    Returns:
        Array of scores, one per octet
    """
    if not octets:
        return np.array([], dtype=np.float32)
    
    # Extract tetramers from octets
    all_tetramers = []
    for tet1, tet2 in octets:
        all_tetramers.extend([tet1, tet2])
    
    # Calculate scores for all tetramers
    tetramer_scores = calculate_tetramer_scores(positions, all_tetramers, sigma, params, debug=False)
    
    # Group scores by octet
    octet_scores = np.zeros(len(octets), dtype=np.float32)
    for i in range(len(octets)):
        # Sum scores of both tetramers in the octet
        octet_scores[i] = tetramer_scores[2*i] + tetramer_scores[2*i + 1]
        
        # Special octet interactions - placeholder for inter-tetramer interactions
        # This could be expanded with specific octet interaction terms
    
    if debug and len(octets) > 0:
        print(f"Calculated scores for {len(octets)} octets:")
        for i in range(min(5, len(octets))):
            print(f"  Octet {i}: Score={octet_scores[i]:.2f}")
    
    return octet_scores

def get_octet_pairs(
    octets: List[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]]]
) -> Set[Tuple[str, int, str, int]]:
    """
    Get all pairs that are part of octets, for exclusion from pair scoring.
    This combines pairs from all tetramers that are part of octets.
    
    Args:
        octets: List of (tetramer1, tetramer2) pairs
        
    Returns:
        Set of (type1, idx1, type2, idx2) tuples representing pairs in octets
    """
    # Extract individual tetramers
    tetramers = []
    for tet1, tet2 in octets:
        tetramers.extend([tet1, tet2])
    
    # Get pairs using the tetramer_pairs function
    return get_tetramer_pairs(tetramers)

#--------------------------------------------------------------------------------
# Combined scoring function
#--------------------------------------------------------------------------------
def calculate_total_score(
    positions: Dict[str, np.ndarray],
    sigma: Dict[str, float],
    tetramers: Optional[List[Tuple]] = None,
    octets: Optional[List[Tuple]] = None,
    sigma_prior_penalty: float = 0.0,
    exclusion_weight: float = 1.0,
    pair_weight: float = 1.0,
    tetramer_weight: float = 1.0,
    octet_weight: float = 1.0,
    params: Optional[Dict[str, Any]] = None,
    debug: bool = False
) -> Tuple[float, float, float, float, float]:
    """
    Calculate combined score for the system with all components.
    
    Args:
        positions: Dictionary mapping particle types to position arrays
        sigma: Dictionary mapping pair types to sigma values
        tetramers: List of tetramers, if already calculated
        octets: List of octets, if already calculated
        sigma_prior_penalty: Prior penalty for sigma values
        exclusion_weight: Weight for excluded volume score
        pair_weight: Weight for pairwise score
        tetramer_weight: Weight for tetramer score
        octet_weight: Weight for octet score
        params: Optional parameters dictionary
        debug: Whether to print debug information
        
    Returns:
        Tuple of (total_score, exclusion_score, pair_score, tetramer_score, octet_score)
    """
    if params is None:
        params = SystemParameters()  # Assumes this function exists
    
    # Calculate excluded volume score
    exclusion_score = exclusion_weight * calculate_excluded_volume(positions, params)
    
    # Identify tetramers if not provided
    if tetramers is None:
        tetramers = get_tetramers(positions, params)
        
    # Identify octets if not provided
    if octets is None and tetramers:
        octets = get_octets(positions, tetramers, params)
    elif octets is None:
        octets = []
    
    # Get pairs to exclude from pair scoring (pairs in tetramers and octets)
    excluded_pairs = set()
    if tetramers:
        excluded_pairs.update(get_tetramer_pairs(tetramers))
    
    # Calculate pairwise score excluding tetramer/octet pairs
    pair_score = pair_weight * calculate_pair_scores(
        positions, sigma, excluded_pairs, params, debug=debug
    )
    
    # Calculate tetramer score
    tetramer_score = 0.0
    if tetramers:
        scores = calculate_tetramer_scores(positions, tetramers, sigma, params, debug=debug)
        tetramer_score = tetramer_weight * np.sum(scores)
    
    # Calculate octet score
    octet_score = 0.0
    if octets:
        scores = calculate_octet_scores(positions, octets, sigma, params, debug=debug)
        octet_score = octet_weight * np.sum(scores)
    
    # Calculate total score
    total_score = exclusion_score + pair_score + tetramer_score + octet_score + sigma_prior_penalty
    
    # Debug output if requested
    if debug:
        print(f"\n===== SCORE SUMMARY =====")
        print(f"Exclusion Score: {exclusion_score:.2f}")
        print(f"Pair Score: {pair_score:.2f}")
        print(f"Tetramer Score: {tetramer_score:.2f}")
        print(f"Octet Score: {octet_score:.2f}")
        print(f"Prior Penalty: {sigma_prior_penalty:.2f}")
        print(f"Total Score: {total_score:.2f}")
    
    return total_score, exclusion_score, pair_score, tetramer_score, octet_score
