#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging
from pathlib import Path
import argparse
from types import SimpleNamespace # Added for the mock model

# Math and file libraries
import math
import numpy as np
import scipy
import mrcfile

# This script does not require ProDy for the coarse-grained example,
# but the original functions keep it as a dependency.
# import prody 

# --- Backend Selection ---
BACKEND = 'cpu'
try:
    import cupy as cp
    import cupyx.scipy.ndimage
    BACKEND = 'gpu'
    logging.warning('CuPy is available, using GPU backend.')
except ImportError:
    logging.info('CuPy is unavailable, using NumPy backend for CPU.')


# =====================================================================
# MODIFIED SCORING FUNCTION FOR MCMC SIMULATION
# =====================================================================
def calculate_ccc_score(sphere_coords, sphere_radii, target_density_map, resolution, backend='cpu'):
    """
    Calculates the CCC score for a coarse-grained model against a target map.

    Args:
        sphere_coords (np.ndarray): An array of shape (N, 3) for sphere coordinates.
        sphere_radii (np.ndarray): An array of shape (N,) for sphere radii.
        target_density_map (mrcfile object): The pre-loaded target map.
        resolution (float): The resolution to use for blurring the model.
        backend (str): 'cpu' or 'gpu'.

    Returns:
        float: The calculated cross-correlation coefficient (CCC).
    """
    # 1. Create a "mock" model that provides coordinates and masses (derived from radii).
    mock_model = SimpleNamespace()
    mock_model.getCoords = lambda: sphere_coords
    # Use the volume (radius^3) as the weight for each sphere.
    mock_model.getMasses = lambda: sphere_radii**3

    # 2. Run the comparison using the script's existing functions
    ccc = -1.0 # Default bad score
    use_gpu = (BACKEND == 'gpu' and backend == 'gpu')

    if use_gpu:
        ccc = compare_data_gpu(target_density_map, mock_model, resolution)
    else:
        ccc = compare_data_cpu(target_density_map, mock_model, resolution)

    return ccc


# =====================================================================
# ORIGINAL SCRIPT FUNCTIONS (Unchanged)
# =====================================================================

def resolution_to_sigma(resolution: float, pixel_size: float) -> float:
    rsig_ = resolution / (4 * math.sqrt(2. * math.log(2.))) / pixel_size
    return rsig_

def parse_density(fname: str) -> mrcfile.mrcfile.MrcFile:
    density = mrcfile.open(fname, permissive=True)
    try:
        assert np.allclose(density.voxel_size.x, density.voxel_size.y) and np.allclose(density.voxel_size.x, density.voxel_size.z)
    except AssertionError:
        logging.error('Non-uniform grids are not supported')
        density.close()
        sys.exit()
    return density

def bins_from_density(density) -> list:
    binsx = (np.linspace(0, density.header.nx, density.header.nx + 1) - density.header.nx/2) * density.voxel_size.x
    binsy = (np.linspace(0, density.header.ny, density.header.ny + 1) - density.header.ny/2) * density.voxel_size.y
    binsz = (np.linspace(0, density.header.nz, density.header.nz + 1) - density.header.nz/2) * density.voxel_size.z
    return (binsx, binsy, binsz)

def compare_data_cpu(density, model, resolution: float) -> float:
    bins = bins_from_density(density)
    coords = model.getCoords()
    weights = model.getMasses()
    projection = calc_projection_cpu(coords, weights, bins, resolution)
    ccc = pairwise_correlation_cpu(projection.flatten(), density.data.flatten())  # REMOVED [0]
    return ccc

def calc_projection_cpu(coords, weights, bins, resolution: float):
    img_, _ = np.histogramdd(coords, weights=weights, bins=bins)
    img_ = np.swapaxes(img_, 0, 2)
    voxel_size = bins[0][1] - bins[0][0]
    sigma  = resolution_to_sigma(resolution, voxel_size)
    img_sm = scipy.ndimage.gaussian_filter(img_, sigma, truncate=4).astype(np.float32)
    return img_sm

def compare_data_gpu(density, model, resolution: float) -> float:
    bins = bins_from_density(density)
    coords = model.getCoords()
    weights = model.getMasses()
    bins_gpu = cp.asarray(bins)
    coords_gpu = cp.asarray(coords)
    weights_gpu = cp.asarray(weights)
    projection = calc_projection_gpu(coords_gpu, weights_gpu, bins_gpu, resolution)
    density_data_gpu = cp.asarray(density.data)
    ccc = pairwise_correlation_gpu(projection.flatten(), density_data_gpu.flatten())
    return float(cp.asnumpy(ccc))  # REMOVED [0] and simplified conversion

def calc_projection_gpu(coords, weights, bins, resolution: float):
    img_, _ = cp.histogramdd(coords, weights=weights, bins=bins)
    img_ = cp.swapaxes(img_, 0, 2)
    voxel_size = bins[0][1] - bins[0][0]
    sigma  = resolution_to_sigma(resolution, float(voxel_size))
    img_sm = cupyx.scipy.ndimage.gaussian_filter(img_, sigma, truncate=4).astype(np.float32)
    return img_sm

def pairwise_correlation_cpu(A, B):
    am = A - np.mean(A)
    bm = B - np.mean(B)
    return np.sum(am * bm) / (np.sqrt(np.sum(am**2)) * np.sqrt(np.sum(bm**2)))

def pairwise_correlation_gpu(A, B):
    am = A - cp.mean(A)
    bm = B - cp.mean(B)
    return cp.sum(am * bm) / (cp.sqrt(cp.sum(am**2)) * cp.sqrt(cp.sum(bm**2)))

# =====================================================================
# EXAMPLE RUNNER
# =====================================================================
def create_dummy_map_from_model(coords, radii, resolution, voxel_size, 
                                box_size, filename="target_map.mrc"):
    """Generates and saves a simulated .mrc map from a coarse-grained model.
       Coordinate frame: bins span [-box_size/2, +box_size/2] uniformly; 
       origin (0,0,0) at box center.
       Header origin is set to -box_size/2 so voxel (0,0,0) corresponds 
       to (-L/2 + Δ/2) + shift semantics.
    """
    print(f"Creating a dummy target map '{filename}'...")
    coords = np.asarray(coords, dtype=np.float32)
    radii  = np.asarray(radii, dtype=np.float32)
    weights = radii**3

    # Diagnostics
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    print(f"Coordinate bounds: X=[{min_coords[0]:.1f},{max_coords[0]:.1f}] "
          f"Y=[{min_coords[1]:.1f},{max_coords[1]:.1f}] "
          f"Z=[{min_coords[2]:.1f},{max_coords[2]:.1f}]")
    half = box_size / 2.0
    print(f"Intended physical cube: [-{half:.1f}, +{half:.1f}] Å per axis")

    # Grid
    grid_dim = int(round(box_size / voxel_size))
    assert math.isclose(grid_dim * voxel_size, box_size, rel_tol=1e-6), "box_size not divisible by voxel_size"
    print(f"Grid: {grid_dim}^3 voxels; voxel_size={voxel_size} Å")

    # Independent bin arrays (edges)
    bins_x = np.linspace(-half, half, grid_dim + 1, dtype=np.float64)
    bins_y = np.linspace(-half, half, grid_dim + 1, dtype=np.float64)
    bins_z = np.linspace(-half, half, grid_dim + 1, dtype=np.float64)
    bins = (bins_x, bins_y, bins_z)

    # Histogram + blur
    raw, _ = np.histogramdd(coords, bins=bins, weights=weights)
    # Optional axis swap: keep consistent with rest of pipeline
    raw = np.swapaxes(raw, 0, 2)  # matches your existing calc_projection_cpu
    sigma_vox = resolution_to_sigma(resolution, voxel_size)
    density = scipy.ndimage.gaussian_filter(raw, sigma_vox, truncate=4).astype(np.float32)

    print(f"Density stats before write: min={density.min():.4g} max={density.max():.4g} "
          f"mean={density.mean():.4g} std={density.std():.4g} nonzero={np.count_nonzero(density)}")

    if density.std() < 1e-10:
        print("ERROR: zero variance density")
        return None

    # Write MRC with correct origin and labels
    with mrcfile.new(filename, overwrite=True) as mrc:
        mrc.set_data(density)              # shape (nz, ny, nx) given swap; confirm orientation downstream
        mrc.voxel_size = voxel_size        # sets cella.{x,y,z} = nx*voxel
        # Set origin so that voxel index 0 corresponds to physical -half
        mrc.header.origin.x = -half
        mrc.header.origin.y = -half
        mrc.header.origin.z = -half
        # (Optional) also ensure nxstart,y,z start = 0 (already default)
        mrc.update_header_stats()          # recompute dmin/dmax/dmean
        # Add informative labels
        labels = [
            f"Simulated map; box={box_size:.1f}A; vox={voxel_size:.3f}A",
            f"Resolution(param)={resolution:.2f}A; sigma(vox)={sigma_vox:.3f}",
            f"Origin set to (-{half:.1f}, -{half:.1f}, -{half:.1f}) A",
            f"grid_dim={grid_dim}"
        ]
        for i, lab in enumerate(labels):
            mrc.header.label[i] = lab.encode('ascii', 'replace')[:80].ljust(80, b'\x00')
        mrc.header.nlabl = len(labels)

    # Report header interpretation
    with mrcfile.open(filename, permissive=True) as chk:
        nx, ny, nz = chk.header.nx, chk.header.ny, chk.header.nz
        print(f"Written header: nx,ny,nz = {nx},{ny},{nz}")
        print(f"Cell dims (Å): {chk.header.cella.x:.1f},{chk.header.cella.y:.1f},{chk.header.cella.z:.1f}")
        print(f"Voxel size (Å): {chk.voxel_size.x:.3f}")
        print(f"Origin (Å): ({chk.header.origin.x:.1f},{chk.header.origin.y:.1f},{chk.header.origin.z:.1f})")
        print(f"Center (Å): ({chk.header.origin.x + 0.5*box_size:.1f}, "
              f"{chk.header.origin.y + 0.5*box_size:.1f}, "
              f"{chk.header.origin.z + 0.5*box_size:.1f}) expected ≈ (0,0,0)")
    return density

if __name__ == "__main__":
    # 1. DEFINE YOUR TOY MODEL DATA
    # In a real MCMC, this data would come from your RMF3 parser and be updated in a loop
    print("--- Setting up the toy model ---")
    array_A  = np.array([[ 63.,0.,0.],[ 44.55,44.55,0.],[0.,63.,0.],[-44.55,44.55,0.],[-63.,0.,0.],[-44.55,-44.55,0.],[ -0.,-63.,0.],[ 44.55,-44.55,0.]])
    array_B = np.array([[ 63.,0.,-38.5],[ 44.55,44.55,-38.5],[0.,63.,-38.5],[-44.55,44.55,-38.5],[-63.,0.,-38.5],[-44.55,-44.55,-38.5],[-0.,-63.,-38.5],[ 44.55,-44.55,-38.5]]) 
    array_C = np.array([[ 47.,0.,-68.5],[ 79.,0.,-68.5],[ 55.86,55.86,-68.5],[ 33.23,33.23,-68.5],[0.,47.,-68.5],[0.,79.,-68.5],[-55.86,55.86,-68.5],[-33.23,33.23,-68.5],[-47.,0.,-68.5],[-79.,0.,-68.5],[-55.86,-55.86,-68.5],[-33.23,-33.23,-68.5],[0.,-47.,-68.5],[0.,-79.,-68.5],[ 55.86,-55.86,-68.5],[ 33.23,-33.23,-68.5]])

    radius_A, radius_B, radius_C = 24.0, 14.0, 16.0
    copy_numbers = {'A': 8, 'B': 8, 'C': 16}

    # Consolidate all coordinates and radii into single arrays
    ideal_coords = np.vstack([array_A, array_B, array_C])
    ideal_radii = np.repeat([radius_A, radius_B, radius_C], repeats=[copy_numbers['A'], copy_numbers['B'], copy_numbers['C']])
    print(f"Model has {len(ideal_coords)} total particles.")

    # 2. CREATE A DUMMY "EXPERIMENTAL" MAP FROM THE IDEAL MODEL
    RESOLUTION = 50.0  # Resolution of the map in Angstroms
    VOXEL_SIZE = 5.0   # Voxel size of the map in Angstroms/pixel
    BOX_SIZE = 800.0   # Box size of the map in Angstroms
    TARGET_MAP_FILE = "target_map.mrc"
    
    create_dummy_map_from_model(ideal_coords, ideal_radii, RESOLUTION, VOXEL_SIZE, BOX_SIZE, filename=TARGET_MAP_FILE)

    mrcfile.validate(TARGET_MAP_FILE)
    print(f"Dummy target map '{TARGET_MAP_FILE}' created and validated.")
    
    #-------------------------------------------------------------------------------
    # give a perturbed system and then get the density map for this 
    array_A_perturbed = array_A + np.random.normal(0, 12.0, array_A.shape)
    array_B_perturbed = array_B + np.random.normal(0, 10.0, array_B.shape)
    array_C_perturbed = array_C + np.random.normal(0, 15.0, array_C.shape)
    pert_coords = np.vstack([array_A_perturbed, array_B_perturbed, array_C_perturbed])
    create_dummy_map_from_model(pert_coords, ideal_radii, RESOLUTION, VOXEL_SIZE, BOX_SIZE, filename="perturbed_map.mrc")

    #-------------------------------------------------------------------------------
    
    # 3. LOAD THE TARGET MAP (this would be done once before an MCMC loop)
    target_map = parse_density(TARGET_MAP_FILE)

    # 4. SCORE THE IDEAL MODEL (this should be ~1.0)
    print("\n--- Scoring the ideal (ground truth) model ---")
    score_ideal = calculate_ccc_score(ideal_coords, ideal_radii, target_map, RESOLUTION)
    print(f"CCC score for ideal model: {score_ideal:.6f}")

    # 5. SCORE A PERTURBED MODEL (to simulate an MCMC step)
    print("\n--- Scoring a perturbed model ---")
    # Add some random noise to the coordinates
    perturbed_coords = ideal_coords + np.random.normal(0, 2.0, ideal_coords.shape) # 2.0 Angstrom standard deviation
    score_perturbed = calculate_ccc_score(perturbed_coords, ideal_radii, target_map, RESOLUTION)
    print(f"CCC score for perturbed model: {score_perturbed:.6f}")
    
    target_map.close()
