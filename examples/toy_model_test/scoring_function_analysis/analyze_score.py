#!/usr/bin/env python3
"""
Simple scoring script - no fancy imports, just add the path and import directly.
"""

import sys
import os
import numpy as np

# Just add the parent directory to Python's path - that's it!
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

# Now import normally - Python will find the files
import parameters
import pair_sampler
import tetramer_sampler
import octet_sampler

from parameters import SystemParameters
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from octet_sampler import OctetSampler

#-----------------------------------------------------
# Instantiate the parameters object
#-----------------------------------------------------
params = SystemParameters()
# get the ideal coordinates
ideal_coords = params.latest_ideal()
print(f"Ideal coordinates: {ideal_coords}")
#-----------------------------------------------------
# create instance of the pair sampler
#-----------------------------------------------------
def pair_sampler_score(coords, sigmas, sigma_ranges):
    sampler_sequence = ["pair_sampler"]
    sampler_idx = 0
    sampler = PairSampler(sampler_sequence, sampler_idx)
    tot_score, ex_score, pair_score, prior_score = sampler.calculate_score(pos = coords, sig = sigmas, sig_range=sigma_ranges)
    return tot_score, ex_score, pair_score, prior_score

def tetramer_sampler_score(coords, sigmas, sigma_ranges):
    sampler_sequence = ["tetramer_sampler"]
    sampler_idx = 0
    sampler = TetramerSampler(sampler_sequence, sampler_idx, positions_ts=coords)
    tetramers = sampler.get_tetramers(positions =coords)
    tot_score, ex_score, pair_score, tet_score = sampler.neg_log_posterior(positions = coords,
                                                                             tetramers = tetramers, 
                                                                           sig = sigmas)
    return tot_score, ex_score, pair_score, tet_score

# Create a dictionary of sigma values AA, AB, BC
sigmas = {
    'AA': 0.9,
    'AB': 0.7,
    'BC': 0.1,
}
# this is how I want the ranges to be defined sigma_ranges: [str, Tuple[float, float]]]
sigma_ranges = {
    'AA': (0.1, 2.2),
    'AB': (0.1, 2.3),
    'BC': (0.2, 2.4),
}
tot_score, ex_score, pair_score, prior_score = pair_sampler_score(ideal_coords, sigmas, sigma_ranges)

print(f"Total Score: {tot_score}, Excluded Volume Score: {ex_score}, Pair Score: {pair_score}, Prior Score: {prior_score}")

tot_score, ex_score, pair_score, tet_score = tetramer_sampler_score(ideal_coords, sigmas, sigma_ranges)
print(f"Total Score: {tot_score}, Excluded Volume Score: {ex_score}, Pair Score: {pair_score}, Tetramer Score: {tet_score}")
