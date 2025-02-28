#----------------------------------------------------------------------
# Running the hierarchical sampling and analyzing the results
#----------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Type

from parameters import SystemParameters
from base_sampler import BaseMCSampler
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from octamer_sampler import OctamerSampler
from visualization import visualize_3d_configuration, visualize_trajectory_plotly
from graph_results import SigmaAnalyzer

#----------------------------------------------------------------------
# First define a list of samplers
# Define the number of chains to run for each sampler
#----------------------------------------------------------------------
sampler_classes = [PairSampler, TetramerSampler, OctamerSampler]
n_chains = 1  # Number of independent chains to run for each sampler
# Define the number of steps for each chain
n_steps = 200000  # Number of Monte Carlo steps per chain
# Define the frequency at which to save trajectory data
save_freq = 500  # Save every 500 steps

sampl = sampler_classes[0].__name__  # Access the name of the first sampler class
for i in range(n_chains):
    sampl.run_mc(n_steps=n_steps, save_freq=save_freq)  # Run the MCMC simulation
    # label the trajectory names based on the chain number 
     
