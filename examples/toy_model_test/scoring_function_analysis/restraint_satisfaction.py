#------------------------------------------------------------------
# Script to compute the number of restraints that are satisfied
# for a given set of coordinates.
#-------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd 
import os
import matplotlib.pyplot as plt
#-------------------------------------------------------------------
# Just add the parent directory to Python's path - that's it!
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

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
def find_and_plot_satisfied_restraints(coords):
    """
    Coords is a numpy array of shape (n, 3) where n is the number of atoms.
    The restraints in this system are defined as AA distances being smaller than 
    or equal to 48.22, AB distances being smaller than or equal to 38.5
    and BC distances being smaller than or equal to 34.0. 
    """
    # Calculate all vs all distances between A particles 
    a1 = coords['A'][0]
    a2 = coords['A'][1]
    # distance between A1 and A2
    d_a1_a2 = np.linalg.norm(a1 - a2)
    print(f"Distance between A1 and A2: {d_a1_a2}")

# call this function
find_and_plot_satisfied_restraints(ideal_coords)