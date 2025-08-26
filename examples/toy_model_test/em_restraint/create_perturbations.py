#!/usr/bin/env python3
"""
Perturbation generation + scoring (CCC + Excluded Volume) with diagnostics.

Goals:
- Generate a broad INTERNAL (aligned) RMSD distribution (not just rigid motions)
- Keep excluded volume scores bounded without dwarfing CCC
- Provide raw vs aligned RMSD (raw includes rigid-body motion; aligned removes it)
- Deterministic clash resolution consistent with scoring threshold
- Transparent logging of each accepted perturbation

Adjust key constants below to tune behavior.
"""

import numpy as np
import os
import sys
import h5py
from scipy.spatial.distance import cdist

# -----------------------------------------------------------------------------
# Import external scoring helpers (expected in same directory)
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from working_example import calculate_ccc_score, parse_density, create_dummy_map_from_model

# -----------------------------------------------------------------------------
# System Definition
# -----------------------------------------------------------------------------
array_A = np.array([
    [ 63.00,   0.00,   0.00], [ 44.55,  44.55,   0.00], [  0.00,  63.00,   0.00],
    [-44.55,  44.55,   0.00], [-63.00,   0.00,   0.00], [-44.55, -44.55,   0.00],
    [ -0.00, -63.00,   0.00], [ 44.55, -44.55,   0.00]
])
array_B = np.array([
    [ 63.00,   0.00, -38.50], [ 44.55,  44.55, -38.50], [  0.00,  63.00, -38.50],
    [-44.55,  44.55, -38.50], [-63.00,   0.00, -38.50], [-44.55, -44.55, -38.50],
    [ -0.00, -63.00, -38.50], [ 44.55, -44.55, -38.50]
])
array_C = np.array([
    [ 47.00,   0.00, -68.50], [ 79.00,   0.00, -68.50], [ 55.86,  55.86, -68.50],
    [ 33.23,  33.23, -68.50], [  0.00,  47.00, -68.50], [  0.00,  79.00, -68.50],
    [-55.86,  55.86, -68.50], [-33.23,  33.23, -68.50], [-47.00,   0.00, -68.50],
    [-79.00,   0.00, -68.50], [-55.86, -55.86, -68.50], [-33.23, -33.23, -68.50],
    [  0.00, -47.00, -68.50], [  0.00, -79.00, -68.50], [ 55.86, -55.86, -68.50],
    [ 33.23, -33.23, -68.50]
])

radius = {'A': 24.0, 'B': 14.0, 'C': 16.0}

# -----------------------------------------------------------------------------
# Tunable Parameters
# -----------------------------------------------------------------------------
N_TARGET                 = 200          # Number of accepted perturbations
JITTER_RANGE             = (0.1, 25.0)  # Std dev range for Gaussian per-particle noise (Å)
INCREMENTAL_STEPS        = 5            # Split jitter into steps (keeps overlaps saner)
OVERLAP_TOLERANCE        = 1.0          # Distance threshold factor for allowed contacts (1.0 = radii sum)
SIGMA_EXVOL              = 2.0          # Softness for excluded volume penalty
MAX_EXVOL_SCORE          = 5000.0       # Reject if above this (still allows CCC to matter)
MAX_CLASH_RELAX_ITERS    = 120          # Max clash resolution iterations
PAIR_BUFFER              = 0.25         # Extra separation beyond threshold when resolving (Å)
MAX_ATTEMPTS_FACTOR      = 6            # Max attempts = N_TARGET * this
LOG_EVERY_ACCEPT         = 1            # Print each accepted (set higher to reduce log spam)
STORE_RAW_RMSD           = True

# -----------------------------------------------------------------------------
# Basic Helpers
# -----------------------------------------------------------------------------
def combine_arrays(aA, aB, aC):
    return np.vstack([aA, aB, aC])

def create_radius_array(n_A, n_B, n_C):
    return np.concatenate([
        np.full(n_A, radius['A']),
        np.full(n_B, radius['B']),
        np.full(n_C, radius['C'])
    ])

# -----------------------------------------------------------------------------
# Overlap / Separation
# -----------------------------------------------------------------------------
def find_clashes(coords, radii, tolerance=OVERLAP_TOLERANCE):
    d = cdist(coords, coords)
    np.fill_diagonal(d, np.inf)
    min_allowed = (radii[:, None] + radii) * tolerance
    mask = d < min_allowed
    return mask, d, min_allowed

def pair_separation_relax(coords, radii, tolerance=OVERLAP_TOLERANCE,
                          max_iter=MAX_CLASH_RELAX_ITERS, buffer=PAIR_BUFFER):
    coords = coords.copy()
    for _ in range(max_iter):
        mask, d, min_allowed = find_clashes(coords, radii, tolerance)
        if not np.any(mask):
            break
        i_idx, j_idx = np.where(mask)
        processed = set()
        for i, j in zip(i_idx, j_idx):
            if i >= j: continue
            if (i, j) in processed: continue
            dist = d[i, j]
            if dist < 1e-8:
                direction = np.random.normal(size=3)
                direction /= np.linalg.norm(direction)
            else:
                direction = (coords[j] - coords[i]) / dist
            needed = (min_allowed[i, j] + buffer) - dist
            if needed <= 0: 
                continue
            shift = 0.5 * needed
            coords[i] -= direction * shift
            coords[j] += direction * shift
            processed.add((i, j))
    return coords

def count_overlap_pairs(coords, radii, tolerance=OVERLAP_TOLERANCE):
    mask, _, _ = find_clashes(coords, radii, tolerance)
    # Count unique pairs i<j
    i_idx, j_idx = np.where(mask)
    return sum(1 for i, j in zip(i_idx, j_idx) if i < j)

# -----------------------------------------------------------------------------
# Perturbation Generation
# -----------------------------------------------------------------------------
def incremental_gaussian_jitter(coords, std_dev, steps=INCREMENTAL_STEPS):
    """
    Apply Gaussian noise in several increments so final per-particle displacement
    std ≈ std_dev (variance adds linearly).
    """
    step_sigma = std_dev / np.sqrt(steps)
    out = coords.copy()
    for _ in range(steps):
        out += np.random.normal(scale=step_sigma, size=out.shape)
    return out

# -----------------------------------------------------------------------------
# Excluded Volume Scoring (consistent tolerance) 
# -----------------------------------------------------------------------------
def excluded_volume_score_by_type(positions, radii_dict, sigma=SIGMA_EXVOL, tolerance=OVERLAP_TOLERANCE):
    score = 0.0
    types = list(positions.keys())
    for i, t1 in enumerate(types):
        pos1 = positions[t1]
        r1 = radii_dict[t1]
        for j, t2 in enumerate(types[i:]):
            t2 = types[i + j]
            pos2 = positions[t2]
            r2 = radii_dict[t2]
            if len(pos1) == 0 or len(pos2) == 0:
                continue
            min_dist = (r1 + r2) * tolerance
            d = cdist(pos1, pos2)
            if t1 == t2:
                if len(pos1) < 2: 
                    continue
                triu = np.triu_indices(len(pos1), k=1)
                dd = d[triu]
                mask = dd < min_dist
                if np.any(mask):
                    overlap = (min_dist - dd[mask])
                    score += np.sum((overlap / sigma) ** 2)
            else:
                mask = d < min_dist
                if np.any(mask):
                    overlap = (min_dist - d[mask])
                    score += np.sum((overlap / sigma) ** 2)
    return float(score)

# -----------------------------------------------------------------------------
# CCC + Total Score
# -----------------------------------------------------------------------------
def calculate_total_score(coords, radii_array, target_density, resolution):
    ccc_raw = calculate_ccc_score(coords, radii_array, target_density, resolution, backend='cpu')
    ccc_score = 500.0 * (1.0 - ccc_raw)  # Lower is better if ccc_raw close to 1
    n_A, n_B, n_C = 8, 8, 16
    positions = {
        'A': coords[:n_A],
        'B': coords[n_A:n_A+n_B],
        'C': coords[n_A+n_B:]
    }
    exvol_score = excluded_volume_score_by_type(positions, radius)
    total = ccc_score + exvol_score
    return total, ccc_score, exvol_score

# -----------------------------------------------------------------------------
# RMSD (raw vs aligned)
# -----------------------------------------------------------------------------
def kabsch_align(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    P0, Q0 = P - cP, Q - cQ
    H = P0.T @ Q0
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    return (R @ P0.T).T + cQ

def aligned_rmsd(ref, pert):
    aligned = kabsch_align(pert, ref)
    return np.sqrt(np.mean(np.sum((ref - aligned)**2, axis=1)))

def raw_rmsd(ref, pert):
    return np.sqrt(np.mean(np.sum((ref - pert)**2, axis=1)))

# -----------------------------------------------------------------------------
# H5 Saving
# -----------------------------------------------------------------------------
def save_to_h5(results, output_file):
    with h5py.File(output_file, 'w') as f:
        gA = f.create_group('array_A')
        gB = f.create_group('array_B')
        gC = f.create_group('array_C')
        gM = f.create_group('metadata')
        for i, r in enumerate(results):
            name = f"frame_{i:04d}"
            gA.create_dataset(name, data=r['coords_A'])
            gB.create_dataset(name, data=r['coords_B'])
            gC.create_dataset(name, data=r['coords_C'])
            meta = gM.create_group(name)
            for k, v in r.items():
                if k.startswith("coords_"):
                    continue
                if isinstance(v, str):
                    meta.attrs[k] = v.encode('utf-8')
                else:
                    meta.attrs[k] = v
        f.attrs['num_frames'] = len(results)
        f.attrs['description'] = "Perturbations with diagnostics".encode('utf-8')
        f.attrs['overlap_tolerance'] = OVERLAP_TOLERANCE
        f.attrs['sigma_excluded_volume'] = SIGMA_EXVOL
        f.attrs['max_exvol_score'] = MAX_EXVOL_SCORE

# -----------------------------------------------------------------------------
# Main Perturbation Loop
# -----------------------------------------------------------------------------
def main():
    print("Perturbation Generation + Scoring")
    print("=" * 60)

    # Output setup
    out_dir = "output_data_revised"
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "perturbation_analysis_revised.h5")

    # Reference system
    ref_coords = combine_arrays(array_A, array_B, array_C)
    n_A, n_B, n_C = len(array_A), len(array_B), len(array_C)
    radii_array = create_radius_array(n_A, n_B, n_C)

    # Density map
    target_map_file = "target_map.mrc"
    resolution, voxel_size, box_size = 15.0, 4.0, 800.0
    if not os.path.exists(target_map_file):
        print("Creating target density map...")
        create_dummy_map_from_model(ref_coords, radii_array, resolution,
                                    voxel_size, box_size, target_map_file)
    target_density = parse_density(target_map_file)

    ref_total, ref_ccc, ref_exvol = calculate_total_score(ref_coords, radii_array,
                                                          target_density, resolution)
    print(f"Reference Scores: Total={ref_total:.2f}  CCC={ref_ccc:.2f}  ExVol={ref_exvol:.2f}")
    print(f"Settings: overlap_tol={OVERLAP_TOLERANCE} sigma_exvol={SIGMA_EXVOL} max_exvol={MAX_EXVOL_SCORE}")
    print(f"Jitter std range: {JITTER_RANGE} (Å) incremental_steps={INCREMENTAL_STEPS}\n")

    results = []
    attempts = 0
    max_attempts = N_TARGET * MAX_ATTEMPTS_FACTOR

    while len(results) < N_TARGET and attempts < max_attempts:
        attempts += 1
        intensity = np.random.uniform(*JITTER_RANGE)
        # Generate perturbation
        pert_coords = incremental_gaussian_jitter(ref_coords, intensity)
        # Resolve overlaps
        pert_coords = pair_separation_relax(pert_coords, radii_array, tolerance=OVERLAP_TOLERANCE)
        # Score
        total_score, ccc_score, exvol_score = calculate_total_score(
            pert_coords, radii_array, target_density, resolution
        )
        # Decide acceptance
        if exvol_score > MAX_EXVOL_SCORE:
            continue

        # RMSDs
        rmsd_raw_val = raw_rmsd(ref_coords, pert_coords)
        rmsd_aligned_val = aligned_rmsd(ref_coords, pert_coords)

        # Overlap stats (under scoring tolerance)
        overlap_pairs = count_overlap_pairs(pert_coords, radii_array, tolerance=OVERLAP_TOLERANCE)

        rec = {
            'coords_A': pert_coords[:n_A],
            'coords_B': pert_coords[n_A:n_A+n_B],
            'coords_C': pert_coords[n_A+n_B:],
            'total_score': total_score,
            'ccc_score': ccc_score,
            'exvol_score': exvol_score,
            'rmsd': rmsd_aligned_val,
            'rmsd_aligned': rmsd_aligned_val,
            'rmsd_raw': rmsd_raw_val if STORE_RAW_RMSD else -1.0,
            'intensity_requested': intensity,
            'overlap_pairs': overlap_pairs,
            'attempt_index': attempts
        }
        results.append(rec)

        if len(results) % LOG_EVERY_ACCEPT == 0:
            print(f"Accepted {len(results):03d}/{N_TARGET} | attempt {attempts:04d} | "
                  f"req_int={intensity:5.1f} | rawRMSD={rmsd_raw_val:6.2f} | "
                  f"alignRMSD={rmsd_aligned_val:6.2f} | CCC={ccc_score:6.2f} | "
                  f"ExVol={exvol_score:7.1f} | overlap_pairs={overlap_pairs}")

    # Summary
    print("\nSummary")
    print("-" * 60)
    print(f"Accepted: {len(results)}  Attempts: {attempts}  Acceptance Rate: {len(results)/max(1,attempts):.2%}")

    if results:
        aligned_arr = np.array([r['rmsd_aligned'] for r in results])
        raw_arr     = np.array([r['rmsd_raw'] for r in results])
        totals      = np.array([r['total_score'] for r in results])
        exvols      = np.array([r['exvol_score'] for r in results])
        cccs        = np.array([r['ccc_score'] for r in results])

        print(f"Aligned RMSD: mean {aligned_arr.mean():.2f} ± {aligned_arr.std():.2f} "
              f"(min {aligned_arr.min():.2f}, max {aligned_arr.max():.2f})")
        print(f"Raw RMSD:     mean {raw_arr.mean():.2f} ± {raw_arr.std():.2f} "
              f"(min {raw_arr.min():.2f}, max {raw_arr.max():.2f})")
        print(f"Total Score:  mean {totals.mean():.1f} ± {totals.std():.1f}")
        print(f"CCC Score:    mean {cccs.mean():.1f} ± {cccs.std():.1f}")
        print(f"ExVol Score:  mean {exvols.mean():.1f} ± {exvols.std():.1f}")
        if len(results) > 1:
            corr = np.corrcoef(aligned_arr, totals)[0,1]
            print(f"Correlation (Aligned RMSD vs Total): {corr:.3f}")

    # Save
    save_to_h5(results, out_file)
    print(f"\nSaved H5: {out_file}")
    print("Done.")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
