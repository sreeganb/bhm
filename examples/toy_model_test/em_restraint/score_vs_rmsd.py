#!/usr/bin/env python3
"""
Simplified perturbation analysis with combined scoring and rejection sampling.
Total Score = CCC Score + Excluded Volume Score
Only accepts perturbations with Total Score < threshold
Includes advanced perturbation methods with overlap resolution.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import h5py
from scipy.spatial.distance import cdist

# Import the scoring function
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from working_example import calculate_ccc_score, parse_density, create_dummy_map_from_model

#--------------------------------------------------------------------------------
# System Definition
array_A = np.array([
    [ 63.  ,   0.  ,   0.  ],
    [ 44.55,  44.55,   0.  ],
    [  0.  ,  63.  ,   0.  ],
    [-44.55,  44.55,   0.  ],
    [-63.  ,   0.  ,   0.  ],
    [-44.55, -44.55,   0.  ],
    [ -0.  , -63.  ,   0.  ],
    [ 44.55, -44.55,   0.  ]
])

array_B = np.array([
    [ 63.  ,   0.  , -38.5 ],
    [ 44.55,  44.55, -38.5 ],
    [  0.  ,  63.  , -38.5 ],
    [-44.55,  44.55, -38.5 ],
    [-63.  ,   0.  , -38.5 ],
    [-44.55, -44.55, -38.5 ],
    [ -0.  , -63.  , -38.5 ],
    [ 44.55, -44.55, -38.5 ]
])

array_C = np.array([
    [ 47.00,   0.00, -68.50], [ 79.00,   0.00, -68.50],
    [ 55.86,  55.86, -68.50], [ 33.23,  33.23, -68.50],
    [  0.00,  47.00, -68.50], [  0.00,  79.00, -68.50],
    [-55.86,  55.86, -68.50], [-33.23,  33.23, -68.50],
    [-47.00,   0.00, -68.50], [-79.00,   0.00, -68.50],
    [-55.86, -55.86, -68.50], [-33.23, -33.23, -68.50],
    [  0.00, -47.00, -68.50], [  0.00, -79.00, -68.50],
    [ 55.86, -55.86, -68.50], [ 33.23, -33.23, -68.50],
])

# Particle properties
radius = {'A': 24.0, 'B': 14.0, 'C': 16.0}

#--------------------------------------------------------------------------------
# Advanced Overlap Detection and Resolution Functions
def check_overlaps(coords, radii, overlap_tolerance=0.9):
    """
    Check if any particles overlap beyond tolerance.
    
    Args:
        coords: (N, 3) array of particle coordinates
        radii: (N,) array of particle radii
        overlap_tolerance: fraction of sum of radii that's allowed (0.8 = 20% overlap allowed)
    
    Returns:
        bool: True if configuration has overlaps beyond tolerance
        list: pairs of overlapping particles
    """
    distances = cdist(coords, coords)
    overlapping_pairs = []
    
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            min_allowed_dist = (radii[i] + radii[j]) * overlap_tolerance
            if distances[i, j] < min_allowed_dist:
                overlapping_pairs.append((i, j, distances[i, j], min_allowed_dist))
    
    return len(overlapping_pairs) > 0, overlapping_pairs

def resolve_overlaps(coords, radii, overlap_tolerance=0.8, max_iterations=100):
    """
    Resolve overlaps by pushing overlapping particles apart.
    
    Args:
        coords: (N, 3) array of particle coordinates
        radii: (N,) array of particle radii
        overlap_tolerance: fraction of sum of radii that's allowed
        max_iterations: maximum iterations to resolve overlaps
    
    Returns:
        coords: Modified coordinates with overlaps resolved
        success: Whether all overlaps were resolved
    """
    coords_modified = coords.copy()
    
    for iteration in range(max_iterations):
        has_overlaps, overlapping_pairs = check_overlaps(coords_modified, radii, overlap_tolerance)
        
        if not has_overlaps:
            return coords_modified, True
        
        # Apply repulsive forces to overlapping pairs
        forces = np.zeros_like(coords_modified)
        
        for i, j, current_dist, min_dist in overlapping_pairs:
            if current_dist > 0:  # Avoid division by zero
                # Vector from i to j
                direction = coords_modified[j] - coords_modified[i]
                direction = direction / np.linalg.norm(direction)
                
                # Push apart proportional to overlap
                overlap_amount = min_dist - current_dist
                push_magnitude = overlap_amount * 0.5  # Each particle moves half
                
                forces[i] -= direction * push_magnitude
                forces[j] += direction * push_magnitude
        
        # Apply forces with damping
        coords_modified += forces * 0.8
    
    return coords_modified, False

def smart_perturbation_with_expansion(coords, radii, perturbation_scale=5.0, 
                                     expansion_factor=1.2, overlap_tolerance=0.8):
    """
    Create perturbation that naturally pushes particles outward to avoid overlaps.
    
    Args:
        coords: (N, 3) array of particle coordinates
        radii: (N,) array of particle radii
        perturbation_scale: Scale of random perturbations in Angstroms
        expansion_factor: How much to expand the structure (1.2 = 20% expansion)
        overlap_tolerance: Minimum allowed distance as fraction of sum of radii
    
    Returns:
        perturbed_coords: New coordinates with perturbation applied
        expansion_applied: Amount of expansion that was applied
    """
    N = len(coords)
    perturbed = coords.copy()
    
    # Step 1: Apply random perturbations
    random_shifts = np.random.randn(N, 3) * perturbation_scale
    perturbed += random_shifts
    
    # Step 2: Check for overlaps
    has_overlaps, overlapping_pairs = check_overlaps(perturbed, radii, overlap_tolerance)
    
    if has_overlaps:
        # Step 3: Apply radial expansion from center of mass
        center_of_mass = np.mean(perturbed, axis=0)
        
        # Progressively expand until no overlaps
        current_expansion = 1.0
        max_expansion = 2.0  # Don't expand more than 2x
        
        while has_overlaps and current_expansion < max_expansion:
            # Expand radially from center
            vectors_from_center = perturbed - center_of_mass
            perturbed = center_of_mass + vectors_from_center * expansion_factor
            current_expansion *= expansion_factor
            
            # Check overlaps again
            has_overlaps, overlapping_pairs = check_overlaps(perturbed, radii, overlap_tolerance)
        
        # Step 4: If still overlapping, use force-based resolution
        if has_overlaps:
            perturbed, _ = resolve_overlaps(perturbed, radii, overlap_tolerance)
        
        expansion_applied = current_expansion
    else:
        expansion_applied = 1.0
    
    return perturbed, expansion_applied

def rotate_z(array, angle):
    """Helper: Rotate around Z-axis."""
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    rotation_matrix = np.array([[cos_a, -sin_a, 0], [sin_a, cos_a, 0], [0, 0, 1]])
    return array @ rotation_matrix.T

def layered_perturbation(array_A, array_B, array_C, radius_A, radius_B, radius_C,
                         perturbation_type='smart', intensity=10.0, 
                         enforce_spacing=True, overlap_tolerance=0.8):
    """
    Apply layer-specific perturbations with inter-layer clash prevention.
    
    Args:
        array_A, array_B, array_C: Coordinate arrays for each layer
        radius_A, radius_B, radius_C: Radii for particles in each layer
        perturbation_type: 'smart', 'radial_expansion', 'twist', or 'breathing'
        intensity: Perturbation intensity (interpretation depends on type)
        enforce_spacing: Whether to enforce non-overlapping constraints
        overlap_tolerance: Fraction of radii sum allowed for overlap
    
    Returns:
        Perturbed arrays and perturbation metrics
    """
    A_pert = array_A.copy()
    B_pert = array_B.copy()
    C_pert = array_C.copy()
    
    if perturbation_type == 'smart':
        # Combine all arrays for overlap checking
        combined = np.vstack([A_pert, B_pert, C_pert])
        radii = np.concatenate([
            np.full(len(A_pert), radius_A),
            np.full(len(B_pert), radius_B),
            np.full(len(C_pert), radius_C)
        ])
        
        # Apply smart perturbation
        perturbed, expansion = smart_perturbation_with_expansion(
            combined, radii, intensity, 1.2, overlap_tolerance
        )
        
        # Split back
        n_A, n_B = len(array_A), len(array_B)
        A_pert = perturbed[:n_A]
        B_pert = perturbed[n_A:n_A+n_B]
        C_pert = perturbed[n_A+n_B:]
        
    elif perturbation_type == 'radial_expansion':
        # Expand each ring radially with different factors
        for arr, radius in [(A_pert, radius_A), (B_pert, radius_B), (C_pert, radius_C)]:
            center = np.mean(arr[:, :2], axis=0)  # Center in XY plane
            for i in range(len(arr)):
                xy_vec = arr[i, :2] - center
                if np.linalg.norm(xy_vec) > 0:
                    xy_vec = xy_vec / np.linalg.norm(xy_vec)
                    # Push outward proportional to intensity and inverse to radius
                    push = xy_vec * intensity * (1.0 + 1.0/radius)
                    arr[i, :2] += push
    
    elif perturbation_type == 'twist':
        # Twist the rings relative to each other
        angle_A = np.random.uniform(-intensity, intensity) * np.pi/180
        angle_B = np.random.uniform(-intensity, intensity) * np.pi/180
        angle_C = np.random.uniform(-intensity, intensity) * np.pi/180
        
        A_pert = rotate_z(A_pert, angle_A)
        B_pert = rotate_z(B_pert, angle_B)
        C_pert = rotate_z(C_pert, angle_C)
        
    elif perturbation_type == 'breathing':
        # Rings move apart/together along Z-axis
        z_shift_B = np.random.uniform(-intensity, intensity)
        z_shift_C = np.random.uniform(-intensity, intensity)
        
        B_pert[:, 2] += z_shift_B
        C_pert[:, 2] += z_shift_C
    
    # Enforce spacing if requested
    if enforce_spacing and perturbation_type != 'smart':  # smart already handles this
        combined = np.vstack([A_pert, B_pert, C_pert])
        radii = np.concatenate([
            np.full(len(A_pert), radius_A),
            np.full(len(B_pert), radius_B),
            np.full(len(C_pert), radius_C)
        ])
        
        combined, success = resolve_overlaps(combined, radii, overlap_tolerance)
        
        n_A, n_B = len(array_A), len(array_B)
        A_pert = combined[:n_A]
        B_pert = combined[n_A:n_A+n_B]
        C_pert = combined[n_A+n_B:]
    
    return A_pert, B_pert, C_pert

#--------------------------------------------------------------------------------
# Core Functions (keeping your existing ones)
def combine_arrays(array_A, array_B, array_C):
    """Combine coordinate arrays."""
    return np.vstack([array_A, array_B, array_C])

def create_radius_array(n_A, n_B, n_C):
    """Create radius array matching coordinates."""
    return np.concatenate([
        np.full(n_A, radius['A']),
        np.full(n_B, radius['B']),
        np.full(n_C, radius['C'])
    ])

def split_to_dict(coords, n_A, n_B, n_C):
    """Split combined coordinates into dictionary."""
    return {
        'A': coords[:n_A],
        'B': coords[n_A:n_A+n_B],
        'C': coords[n_A+n_B:]
    }

def excluded_volume_score(positions, radii, sigma=0.1):
    """Calculate excluded volume penalty."""
    score = 0.0
    for type1 in positions.keys():
        for type2 in positions.keys():
            if type1 <= type2:
                pos1, pos2 = positions[type1], positions[type2]
                if len(pos1) == 0 or len(pos2) == 0:
                    continue
                
                min_dist = radii[type1] + radii[type2]
                distances = cdist(pos1, pos2)
                
                if type1 == type2:
                    mask = np.triu(np.ones_like(distances), k=1)
                    viol_mask = (distances < min_dist) & (mask > 0)
                else:
                    viol_mask = distances < min_dist
                
                if np.any(viol_mask):
                    overlaps = min_dist - distances[viol_mask]
                    score += np.sum((overlaps**2) / (sigma**2))
    return score

def calculate_total_score(coords, radii_array, target_density, resolution):
    """Calculate total score = CCC + ExVol."""
    # CCC Score (lower is better, so we transform it)
    ccc_raw = calculate_ccc_score(coords, radii_array, target_density, resolution, backend='cpu')
    ccc_score = 100 * (1 - ccc_raw)  # Transform so higher is worse
    
    # Excluded Volume Score (higher is worse)
    positions = split_to_dict(coords, len(array_A), len(array_B), len(array_C))
    exvol_score = excluded_volume_score(positions, radius)
    
    total_score = ccc_score + exvol_score
    
    return total_score, ccc_score, exvol_score

def kabsch_align(P, Q):
    """Kabsch alignment algorithm."""
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    
    H = P_centered.T @ Q_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    P_aligned = (R @ P_centered.T).T + centroid_Q
    return P_aligned

def calculate_rmsd(ref_coords, pert_coords):
    """Calculate RMSD after optimal alignment."""
    aligned_coords = kabsch_align(pert_coords, ref_coords)
    diff = ref_coords - aligned_coords
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

#--------------------------------------------------------------------------------
# UPDATED: Advanced Perturbation Functions (replacing old simple ones)
def generate_advanced_perturbation(array_A, array_B, array_C, pert_type, intensity):
    """Generate perturbation using advanced methods with overlap resolution."""
    
    # Available perturbation types with new advanced methods
    if pert_type in ['smart', 'radial_expansion', 'twist', 'breathing']:
        # Use advanced layered perturbation
        A_pert, B_pert, C_pert = layered_perturbation(
            array_A, array_B, array_C,
            radius['A'], radius['B'], radius['C'],
            perturbation_type=pert_type,
            intensity=intensity,
            enforce_spacing=True,
            overlap_tolerance=0.8
        )
    else:
        # Fallback to simple perturbations for backward compatibility
        A_pert, B_pert, C_pert = array_A.copy(), array_B.copy(), array_C.copy()
        
        if pert_type == 'rotate':
            angle = np.random.uniform(-intensity * np.pi/180, intensity * np.pi/180)
            A_pert = rotate_z(A_pert, angle)
            B_pert = rotate_z(B_pert, angle)
            C_pert = rotate_z(C_pert, angle)
            
        elif pert_type == 'shift':
            shift = np.random.uniform(-intensity, intensity, 3)
            A_pert += shift
            B_pert += shift
            C_pert += shift
            
        elif pert_type == 'random':
            # Apply random displacements to individual particles
            for arr in [A_pert, B_pert, C_pert]:
                n_particles = max(1, len(arr) // 2)
                indices = np.random.choice(len(arr), n_particles, replace=False)
                for idx in indices:
                    random_shift = np.random.uniform(-intensity, intensity, 3)
                    arr[idx] += random_shift
    
    return A_pert, B_pert, C_pert

def generate_valid_perturbation(array_A, array_B, array_C, pert_type, intensity, 
                               radii_array, target_density, resolution, 
                               score_threshold, max_attempts=50):
    """Generate perturbation with score below threshold using rejection sampling."""
    
    for attempt in range(max_attempts):
        # Generate advanced perturbation
        A_pert, B_pert, C_pert = generate_advanced_perturbation(
            array_A, array_B, array_C, pert_type, intensity
        )
        coords_pert = combine_arrays(A_pert, B_pert, C_pert)
        
        # Check for overlaps
        has_overlaps, overlap_pairs = check_overlaps(coords_pert, radii_array, overlap_tolerance=0.8)
        
        # Calculate total score
        total_score, ccc_score, exvol_score = calculate_total_score(
            coords_pert, radii_array, target_density, resolution
        )
        
        # Accept if below threshold
        if total_score < score_threshold:
            return (A_pert, B_pert, C_pert), total_score, ccc_score, exvol_score, attempt + 1, has_overlaps, len(overlap_pairs)
    
    # If no valid perturbation found, return None
    return None, None, None, None, max_attempts, True, 0

#--------------------------------------------------------------------------------
# UPDATED: Data Storage (enhanced with overlap information)
def save_to_h5(data_list, output_file):
    """Save all data to H5 file."""
    
    with h5py.File(output_file, 'w') as f:
        # Create coordinate groups
        grp_A = f.create_group('array_A')
        grp_B = f.create_group('array_B')
        grp_C = f.create_group('array_C')
        metadata = f.create_group('metadata')
        
        # Store each frame
        for i, data in enumerate(data_list):
            frame_name = f'frame_{i:04d}'
            
            # Store coordinates
            grp_A.create_dataset(frame_name, data=data['coords_A'], compression='gzip')
            grp_B.create_dataset(frame_name, data=data['coords_B'], compression='gzip')
            grp_C.create_dataset(frame_name, data=data['coords_C'], compression='gzip')
            
            # Store metadata
            frame_meta = metadata.create_group(frame_name)
            frame_meta.attrs['perturbation_type'] = data['pert_type'].encode('utf-8')
            frame_meta.attrs['intensity'] = float(data['intensity'])
            frame_meta.attrs['total_score'] = float(data['total_score'])
            frame_meta.attrs['ccc_score'] = float(data['ccc_score'])
            frame_meta.attrs['exvol_score'] = float(data['exvol_score'])
            frame_meta.attrs['rmsd'] = float(data['rmsd'])
            frame_meta.attrs['attempts'] = int(data['attempts'])
            frame_meta.attrs['has_overlaps'] = bool(data['has_overlaps'])
            frame_meta.attrs['n_overlaps'] = int(data['n_overlaps'])
        
        # Store global arrays for easy access
        n_frames = len(data_list)
        metadata.create_dataset('all_pert_types', data=[d['pert_type'].encode('utf-8') for d in data_list])
        metadata.create_dataset('all_intensities', data=[d['intensity'] for d in data_list])
        metadata.create_dataset('all_total_scores', data=[d['total_score'] for d in data_list])
        metadata.create_dataset('all_ccc_scores', data=[d['ccc_score'] for d in data_list])
        metadata.create_dataset('all_exvol_scores', data=[d['exvol_score'] for d in data_list])
        metadata.create_dataset('all_rmsds', data=[d['rmsd'] for d in data_list])
        metadata.create_dataset('all_attempts', data=[d['attempts'] for d in data_list])
        metadata.create_dataset('all_has_overlaps', data=[d['has_overlaps'] for d in data_list])
        metadata.create_dataset('all_n_overlaps', data=[d['n_overlaps'] for d in data_list])
        
        # Global attributes
        f.attrs['num_frames'] = n_frames
        f.attrs['description'] = 'Advanced perturbation analysis with overlap resolution'.encode('utf-8')
        f.attrs['score_threshold'] = 1000.0
        f.attrs['overlap_tolerance'] = 0.8
        f.attrs['radii_A'] = radius['A']
        f.attrs['radii_B'] = radius['B']
        f.attrs['radii_C'] = radius['C']

#--------------------------------------------------------------------------------
# UPDATED: Main Analysis
def main():
    print("Advanced Perturbation Analysis with Overlap Resolution")
    print("=" * 60)
    
    # Parameters
    n_perturbations = 100
    score_threshold = 1000.0
    output_dir = "output_data"
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup target density
    target_map_file = "target_map.mrc"
    resolution = 15.0
    voxel_size = 4.0
    box_size = 800.0
    
    # Prepare reference structure
    ref_coords = combine_arrays(array_A, array_B, array_C)
    radii_array = create_radius_array(len(array_A), len(array_B), len(array_C))
    
    # Create/load target map
    if not os.path.exists(target_map_file):
        print("📊 Creating target density map...")
        create_dummy_map_from_model(ref_coords, radii_array, resolution, voxel_size, box_size, target_map_file)
    
    target_density = parse_density(target_map_file)
    
    # Calculate reference scores
    ref_total, ref_ccc, ref_exvol = calculate_total_score(ref_coords, radii_array, target_density, resolution)
    ref_overlaps, ref_overlap_pairs = check_overlaps(ref_coords, radii_array, overlap_tolerance=0.8)
    
    print(f"  Reference Scores:")
    print(f"   Total: {ref_total:.2f}, CCC: {ref_ccc:.4f}, ExVol: {ref_exvol:.2f}")
    print(f"   Has overlaps: {ref_overlaps}, N_overlaps: {len(ref_overlap_pairs)}")
    print(f"   Score threshold: {score_threshold}")
    
    # Enhanced perturbation types including new advanced methods
#    perturbation_types = ['smart', 'radial_expansion', 'twist', 'breathing', 'rotate', 'shift', 'random']
    perturbation_types = ['smart', 'radial_expansion', 'twist', 'breathing', 'shift', 'random']

    intensity_ranges = {
        'smart': (1, 15),          # Smart perturbation with expansion
        'radial_expansion': (1, 10), # Ring expansion
        'twist': (1, 30),          # Ring twisting (degrees)
        'breathing': (1, 20),      # Z-axis motion
#        'rotate': (1, 45),         # Simple rotation
        'shift': (1, 20),          # Simple translation
        'random': (1, 10)          # Random particle motion
    }
    
    valid_data = []
    rejected_count = 0
    total_attempts = 0
    perturbation_stats = {pt: {'count': 0, 'rejections': 0} for pt in perturbation_types}
    
    print(f"\n Generating {n_perturbations} valid perturbations...")
    print(f"   Available perturbation types: {perturbation_types}")
    
    for i in range(n_perturbations):
        # Random perturbation parameters
        pert_type = np.random.choice(perturbation_types)
        min_int, max_int = intensity_ranges[pert_type]
        intensity = np.random.uniform(min_int, max_int)
        
        perturbation_stats[pert_type]['count'] += 1
        
        # Generate valid perturbation
        result = generate_valid_perturbation(
            array_A, array_B, array_C, pert_type, intensity,
            radii_array, target_density, resolution, score_threshold
        )
        
        coords_tuple, total_score, ccc_score, exvol_score, attempts, has_overlaps, n_overlaps = result
        total_attempts += attempts
        
        if coords_tuple is None:
            rejected_count += 1
            perturbation_stats[pert_type]['rejections'] += 1
            print(f"Frame {i+1}: No valid {pert_type} perturbation found after {attempts} attempts")
            continue
        
        # Calculate RMSD
        A_pert, B_pert, C_pert = coords_tuple
        coords_pert = combine_arrays(A_pert, B_pert, C_pert)
        rmsd = calculate_rmsd(ref_coords, coords_pert)
        
        # Store valid data
        data_entry = {
            'coords_A': A_pert,
            'coords_B': B_pert,
            'coords_C': C_pert,
            'pert_type': pert_type,
            'intensity': intensity,
            'total_score': total_score,
            'ccc_score': ccc_score,
            'exvol_score': exvol_score,
            'rmsd': rmsd,
            'attempts': attempts,
            'has_overlaps': has_overlaps,
            'n_overlaps': n_overlaps
        }
        valid_data.append(data_entry)
        
        if (i + 1) % 20 == 0:
            print(f"Processed {i+1}/{n_perturbations} perturbations")
    
    print(f"\n Results Summary:")
    print(f"   Valid perturbations: {len(valid_data)}")
    print(f"   Rejected perturbations: {rejected_count}")
    print(f"   Total attempts: {total_attempts}")
    print(f"   Average attempts per perturbation: {total_attempts/n_perturbations:.1f}")
    
    # Perturbation type statistics
    print(f"\n Perturbation Type Statistics:")
    for pt in perturbation_types:
        stats = perturbation_stats[pt]
        success_rate = (stats['count'] - stats['rejections']) / max(stats['count'], 1) * 100
        print(f"   {pt}: {stats['count']} attempts, {stats['rejections']} rejections, {success_rate:.1f}% success")
    
    if len(valid_data) == 0:
        print("No valid perturbations generated. Consider increasing score threshold.")
        return
    
    # Save to H5 file
    output_file = os.path.join(output_dir, "advanced_perturbation_analysis.h5")
    save_to_h5(valid_data, output_file)
    print(f"Data saved to: {output_file}")
    
    # Quick analysis
    total_scores = [d['total_score'] for d in valid_data]
    ccc_scores = [d['ccc_score'] for d in valid_data]
    exvol_scores = [d['exvol_score'] for d in valid_data]
    rmsds = [d['rmsd'] for d in valid_data]
    overlaps = [d['n_overlaps'] for d in valid_data]
    
    print(f"\n Score Ranges:")
    print(f"   Total: {min(total_scores):.2f} - {max(total_scores):.2f}")
    print(f"   CCC: {min(ccc_scores):.4f} - {max(ccc_scores):.4f}")
    print(f"   ExVol: {min(exvol_scores):.2f} - {max(exvol_scores):.2f}")
    print(f"   RMSD: {min(rmsds):.2f} - {max(rmsds):.2f} Å")
    print(f"   Overlaps: {min(overlaps)} - {max(overlaps)}")
    
    # Enhanced correlation analysis
    import pandas as pd
    df = pd.DataFrame({
        'Total_Score': total_scores,
        'CCC_Score': ccc_scores,
        'ExVol_Score': exvol_scores,
        'RMSD': rmsds,
        'N_Overlaps': overlaps,
        'Pert_Type': [d['pert_type'] for d in valid_data]
    })
    
    print(f"\n Correlations with RMSD:")
    print(f"   RMSD vs Total Score: {df['RMSD'].corr(df['Total_Score']):.3f}")
    print(f"   RMSD vs CCC Score: {df['RMSD'].corr(df['CCC_Score']):.3f}")
    print(f"   RMSD vs ExVol Score: {df['RMSD'].corr(df['ExVol_Score']):.3f}")
    print(f"   RMSD vs N_Overlaps: {df['RMSD'].corr(df['N_Overlaps']):.3f}")
    
    # Count structures with overlaps
    structures_with_overlaps = sum(1 for d in valid_data if d['has_overlaps'])
    print(f"\n Structures with overlaps: {structures_with_overlaps}/{len(valid_data)}")
    
    print(f"\n Advanced analysis complete!")

if __name__ == "__main__":
    main()