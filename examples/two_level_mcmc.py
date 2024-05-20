"""MCMC protocol for the Bayesian hierarchical modeling
   This involves two states of parameter estimation
   In the first stage it is a simple MCMC step, whereas the 
   second stage involves the joint posterior estimation of the 
   local and global parameters"""

from __future__ import print_function, division
import IMP
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.analysis
import IMP.pmi.io
import IMP.pmi.alphabets
import IMP.rmf
import IMP.isd
import IMP.pmi.dof
import os
try:
    from pathlib import Path
except ImportError:  # Use bundled pathlib on Python 2 without pathlib
    from IMP._compat_pathlib import Path
import glob
from operator import itemgetter
from collections import defaultdict
import numpy as np
import itertools
import warnings
import math

import pickle

