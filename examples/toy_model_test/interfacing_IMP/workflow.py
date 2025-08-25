import numpy as np
from prob_sampler import ProbabilisticSampler
from imp_to_prob import IMPToProbabilistic
from sampling_interface import SamplingInterface
import IMP 
import IMP.pmi
from IMP.pmi.restraints.basic import DistanceRestraint

# Your existing IMP setup
root_hier = ...  # Your IMP hierarchy
restraints = [
    DistanceRestraint(root_hier, ...),
    # More restraints
]

# Create bridge
bridge = IMPToProbabilistic(root_hier.get_model(), restraints)

# Create sampler
sampler = ProbabilisticSampler(bridge)

# Create interface
interface = SamplingInterface(root_hier.get_model(), bridge, sampler)

# Run advanced sampling
trace = interface.run_pymc_sampling(draws=5000)

# Get best sample
best_idx = np.argmin(trace.posterior['coords'].mean(axis=0))
best_coords = trace.posterior['coords'].values.reshape(-1, bridge.n_particles, 3)[best_idx]

# Update IMP with best configuration
interface.update_imp_coordinates(best_coords)

# Continue with IMP operations