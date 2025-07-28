#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The script to calculate aproximate model -> denisty projection and Pearson's r
"""

__author__ = "Arthur Zalevsky"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Arthur Zalevsky"
__email__ = "aozalevsky@gmailc.com"
__status__ = "alpha"

import sys
import logging
from pathlib import Path
import argparse

import prody
import mrcfile

import math
import scipy

BACKEND = 'cpu'
import numpy as np
try:
    import cupy as cp
    import cupyx.scipy.ndimage
    BACKEND = 'gpu'
except ImportError:
    logging.warning('CuPy is unavailable, switching to NumPy backend')


def resolution_to_sigma(resolution: float, pixel_size: float) -> float:
    '''Estimate sigma for gaussian smoothing from resolution'''
    # Adapted from IMP.em2d
    rsig_ = resolution / (4 * math.sqrt(2. * math.log(2.))) / pixel_size  # sigma
    return rsig_

def parse_density(fname: str) -> mrcfile.mrcfile.MrcFile:
    '''Read EMDB density in map format'''
    density = mrcfile.open(fname)
    try:
        assert density.voxel_size.x == density.voxel_size.y == density.voxel_size.z
    except AssertionError as e:
        logging.error('Non-uniform grids are not supported')
        density.close()
        sys.exit()
    return density

def parse_model(fname: str) -> prody.atomic.atomgroup.AtomGroup:
    '''Parse model. Currently support only atomic PDB and PDBx/mmCIF formats'''
    path = Path(fname)

    if path.suffixes[-1] == '.pdb' or path.suffixes[-2:] == ['.pdb', '.gz']:
        parser = prody.parsePDB
    elif path.suffixes[-1] == '.cif' or path.suffixes[-2:] == ['.cif', '.gz']:
        parser = prody.parseMMCIF
    else:
        logging.error(f'Wrong file format for {fname}')

    try:
        pdb = parser(fname)
    except Exception as e:
        logging.error(f'Model parsing failed for {fname}')
        logging.error(e)

    return pdb

def bins_from_density(density) -> list:
    '''Setup grid bins from density map'''
    # Setup grid
    binsx = (np.linspace(0, density.header.nx, density.header.nx + 1) - density.nstart.x) * density.voxel_size.x
    binsy = (np.linspace(0, density.header.ny, density.header.ny + 1) - density.nstart.y) * density.voxel_size.y
    binsz = (np.linspace(0, density.header.nz, density.header.nz + 1) - density.nstart.z) * density.voxel_size.z

    return (binsx, binsy, binsz)

def compare_data_cpu(density, model, resolution: float, min_level: float|None=None) -> float:
    # Setup grid

    bins = bins_from_density(density)

    # Calculate projections
    coords = model.getCoords()
    weights = model.getMasses()
    projection = calc_projection_cpu(coords, weights, bins, resolution)
    # Compare data
    density_data = density.data
    ccc = pairwise_correlation_cpu(projection.flat, density.data.flat)[0]

    return ccc

def calc_projection_cpu(coords: list, weights: list, bins: list, resolution: float):
    # digitize particels
    img_, bins_ = np.histogramdd(coords, weights=weights, bins=bins)
    # swap axes to match map layout
    img_ = np.swapaxes(img_, 0, 2)
    # get voxel size from grid
    voxel_size = bins[0][1] - bins[0][0]
    # use desired resolution and voxel size to estimate sigma
    sigma  = resolution_to_sigma(resolution, voxel_size)
    # apply gaussian smoothing
    img_sm = scipy.ndimage.gaussian_filter(img_, sigma, truncate=4).astype(np.float32)
    return img_sm

def compare_data_gpu(density, model, resolution: float, min_level: float|None=None) -> float:
    # Setup grid
    bins = bins_from_density(density)
    # Calculate projections
    coords = model.getCoords()
    weights = model.getMasses()
    # Copy data to gpu
    bins = cp.asarray(bins)
    coords = cp.asarray(coords)
    weights = cp.asarray(weights)

    projection = calc_projection_gpu(coords, weights, bins, resolution)

    # Compare data
    density_data = cp.asarray(density.data)
    ccc = pairwise_correlation_gpu(projection.flatten(), density_data.flatten())
    # Copy data back from gpu
    ccc = cp.asnumpy(ccc)[0]

    return ccc

def calc_projection_gpu(coords: list, weights: list, bins: list, resolution: float):
    # digitize particels
    img_, bins_ = cp.histogramdd(coords, weights=weights, bins=bins)
    # swap axes to match map layout
    img_ = cp.swapaxes(img_, 0, 2)
    # get voxel size from grid
    voxel_size = bins[0][1] - bins[0][0]
    # use desired resolution and voxel size to estimate sigma
    sigma  = resolution_to_sigma(resolution, voxel_size)
    # apply gaussian smoothing
    img_sm = cupyx.scipy.ndimage.gaussian_filter(img_, sigma, truncate=4).astype(np.float32)
    return img_sm

def pairwise_correlation_cpu(A, B):
    # https://stackoverflow.com/a/71847068
    am = A - np.mean(A, axis=0, keepdims=True)
    bm = B - np.mean(B, axis=0, keepdims=True)
    return am.T @ bm /  (np.sqrt(
        np.sum(am**2, axis=0,
               keepdims=True)).T * np.sqrt(
        np.sum(bm**2, axis=0, keepdims=True)))

def pairwise_correlation_gpu(A, B):
    # Adapted from https://stackoverflow.com/a/71847068
    am = A - cp.mean(A, axis=0, keepdims=True)
    bm = B - cp.mean(B, axis=0, keepdims=True)
    return am.T @ bm /  (np.sqrt(
        cp.sum(am**2, axis=0,
               keepdims=True)).T * cp.sqrt(
        cp.sum(bm**2, axis=0, keepdims=True)))

def main():
    """Main function."""

    parser = argparse.ArgumentParser(
                    prog='Approximate Cross-Correlation Coefficient',
                    description='Calculates map-model cross-correlation coefficient',
                    )

    parser.add_argument('--density', type=str, required=True,
                        help='map or map.gz model from EMDDB')
    parser.add_argument('--resolution', type=float, required=True,
                        help='Estimated resolution of the map')
    parser.add_argument('--model', type=str, required=True,
                        help='Atomic model in any of the pdb, pdb.gz, cif, cif.gz formats')
    parser.add_argument('--backend', type=str, required=False, default='cpu',
                        choices=['cpu', 'gpu'],
                        help='Computational backend')


    args = parser.parse_args()

    density = parse_density(args.density)
    resolution = args.resolution
    model = parse_model(args.model)
    backend = args.backend

    ccc = None

    if BACKEND ==  backend == 'cpu':
        ccc = compare_data_cpu(density, model, resolution)
    elif BACKEND == 'gpu' and backend == 'cpu':
        ccc = compare_data_cpu(density, model, resolution)
    elif BACKEND == backend == 'gpu':
        ccc = compare_data_gpu(density, model, resolution)
    else:
        logging.error(f'Backend {backend} is not available')

    if ccc is not None:
        print(f'Map - Model CCC {backend}: {ccc:.3f}')

    density.close()

if __name__ == "__main__":
    sys.exit(main())
