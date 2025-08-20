import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.atom
import IMP.algebra
import IMP.core
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.dof
import numpy as np
import pymc as pm
import os
import sys
from typing import List, Callable, Dict, Tuple

# Optional JAX/NumPyro imports with proper error handling
HAS_JAX = False

# Force JAX to use CPU before importing
os.environ["JAX_PLATFORM_NAME"] = "cpu"
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # Hide GPUs from JAX
os.environ["JAX_PLATFORMS"] = "cpu"

try:
    import jax
    # Additional config to ensure CPU usage
    if hasattr(jax, 'config'):
        jax.config.update("jax_platform_name", "cpu")
        jax.config.update("jax_enable_x64", True)  # Use 64-bit for better precision
    elif hasattr(jax, 'default_backend'):
        # For older JAX versions
        jax.default_backend = 'cpu'
    
    import jax.numpy as jnp
    import numpyro
    import numpyro.distributions as dist
    from numpyro.infer import MCMC, NUTS, SA
    
    # Verify we're using CPU
    try:
        test_array = jnp.array([1.0])
        backend = jax.default_backend()
        print(f"JAX backend: {backend}")
        if 'gpu' in backend.lower() or 'cuda' in backend.lower():
            print("Warning: JAX is using GPU, forcing CPU-only mode")
            raise RuntimeError("GPU detected, need CPU-only mode")
    except Exception as e:
        print(f"JAX backend check: {e}")
    
    HAS_JAX = True
    print("JAX/NumPyro available (CPU mode enforced)")
except Exception as e:
    print(f"JAX/NumPyro not available: {e}")
    print("Will proceed with simple Monte Carlo")

#--------------------------------------------------------------------------
# System Builder Class - Corrected for IMP PMI API
#--------------------------------------------------------------------------
class System_builder:
    """Build IMP system with multiple protein types"""
    
    def __init__(self, ntype, copy_numbers, chain_ids, sequences_type, colors):
        self.ntype = ntype
        self.copy_numbers = copy_numbers
        self.chain_ids = chain_ids
        self.sequences_type = sequences_type
        self.colors = colors
        
    def build_system(self):
        """Build the IMP system using correct PMI API"""
        mdl = IMP.Model()
        system = IMP.pmi.topology.System(mdl)
        
        # Create a single state (most common case)
        state = system.create_state()
        
        parts = {}
        chain_idx = 0
        
        for i in range(self.ntype):
            for j in range(self.copy_numbers[i]):
                name = f"prot_{i+1}"
                mol_name = f"{name}_{j+1}"
                
                # Create molecule using the state
                molecule = state.create_molecule(
                    name=mol_name,
                    sequence=self.sequences_type[i],
                    chain_id=self.chain_ids[chain_idx]
                )
                
                # Add representation to the molecule
                molecule.add_representation(
                    molecule,
                    resolutions=[1],
                    bead_ca_centers=True,
                    color=self.colors[i]
                )
                
                parts[mol_name] = molecule
                chain_idx += 1
        
        return mdl, system, parts

#--------------------------------------------------------------------------
# Scoring Function Class - Fixed
#--------------------------------------------------------------------------
class ScoringFunction:
    """Manage restraints and scoring for the IMP system"""
    
    def __init__(self, model, hierarchy, parts):
        self.model = model
        self.hierarchy = hierarchy
        self.parts = parts
        self.restraints = []
        self.restraint_objects = []  # Keep actual restraint objects for scoring function
        
    def add_excluded_volume_restraint(self):
        """Add excluded volume restraint"""
        ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects=[self.hierarchy],
            resolution=1
        )
        ev.add_to_model()
        self.restraints.append(ev)
        
        # Get the actual IMP restraint object
        if hasattr(ev, 'rs'):
            self.restraint_objects.extend(ev.rs.get_restraints())
        elif hasattr(ev, 'get_restraint'):
            self.restraint_objects.append(ev.get_restraint())
        
        print("Added excluded volume restraint")
        return ev
        
    def add_distance_restraint(self, mol1_spec, mol2_spec, distance_min, distance_max, kappa):
        """Add distance restraint between two molecules"""
        mol1_name, mol1_res, mol1_copy = mol1_spec
        mol2_name, mol2_res, mol2_copy = mol2_spec
        
        # Get molecules from parts
        mol1_key = f"{mol1_name}_{mol1_copy+1}"  # Convert 0-based to 1-based
        mol2_key = f"{mol2_name}_{mol2_copy+1}"
        
        if mol1_key in self.parts and mol2_key in self.parts:
            dr = IMP.pmi.restraints.basic.DistanceRestraint(
                root_hier=self.hierarchy,
                tuple_selection1=(mol1_res, mol1_res, mol1_key),
                tuple_selection2=(mol2_res, mol2_res, mol2_key),
                distancemin=distance_min,
                distancemax=distance_max,
                kappa=kappa
            )
            dr.add_to_model()
            self.restraints.append(dr)
            
            # Get the actual IMP restraint object
            if hasattr(dr, 'rs'):
                self.restraint_objects.extend(dr.rs.get_restraints())
            elif hasattr(dr, 'get_restraint'):
                self.restraint_objects.append(dr.get_restraint())
            
            print(f"Added distance restraint: {mol1_key} - {mol2_key}")
            return dr
        else:
            print(f"Warning: Could not find molecules {mol1_key} or {mol2_key}")
            return None
    
    def get_total_score(self):
        """Get total score from all restraints"""
        return sum(r.evaluate() for r in self.restraints)
    
    def get_scoring_function(self):
        """Get IMP scoring function from restraints"""
        # Create a RestraintSet and add all restraint objects
        rs = IMP.RestraintSet(self.model, "all_restraints")
        for r in self.restraint_objects:
            rs.add_restraint(r)
        
        # Create scoring function from the RestraintSet
        return IMP.core.RestraintsScoringFunction([rs])

#--------------------------------------------------------------------------
# Bridge Class between IMP and Probabilistic Samplers 
#--------------------------------------------------------------------------
class IMPToProbabilistic:
    """Bridge IMP system with probabilistic sampling frameworks"""
    
    def __init__(self, model, hierarchy, restraints):
        self.model = model
        self.hierarchy = hierarchy
        self.restraints = restraints
        self.particles = []
        
        # Extract all particles
        for leaf in IMP.atom.get_leaves(hierarchy):
            if IMP.core.XYZR.get_is_setup(leaf):
                self.particles.append(IMP.core.XYZR(leaf))
        
        print(f"Found {len(self.particles)} particles to sample")
    
    def get_coordinates_array(self):
        """Extract coordinates as numpy array"""
        coords = []
        for p in self.particles:
            xyz = p.get_coordinates()
            coords.append([xyz[0], xyz[1], xyz[2]])
        return np.array(coords)
    
    def set_coordinates_from_array(self, coords_array):
        """Set coordinates from numpy array"""
        # Handle both numpy and JAX arrays
        if hasattr(coords_array, '__array__'):
            coords_array = np.array(coords_array)
        
        coords_array = coords_array.reshape(-1, 3)
        
        for i, p in enumerate(self.particles):
            if i < len(coords_array):
                p.set_coordinates(IMP.algebra.Vector3D(
                    float(coords_array[i][0]), 
                    float(coords_array[i][1]), 
                    float(coords_array[i][2])
                ))
    
    def evaluate_score(self, coords):
        """Evaluate score for given coordinates"""
        self.set_coordinates_from_array(coords)
        score = sum(r.evaluate() for r in self.restraints)
        return score
    
    def build_likelihood_function(self):
        """Build likelihood function for PyMC"""
        def likelihood(coords):
            score = self.evaluate_score(coords)
            return -score  # Negative because we want to minimize score
        return likelihood
    
    def build_jax_likelihood(self):
        """Build JAX-compatible likelihood for NumPyro - avoiding array conversion issues"""
        if not HAS_JAX:
            raise RuntimeError("JAX not available")
        
        # Pre-compute some values to avoid tracing issues
        n_particles = len(self.particles)
        
        def likelihood(coords_flat):
            # NumPyro will pass JAX arrays, we need to handle them carefully
            # We'll use a custom operation that doesn't trigger array conversion
            
            # For now, return a simple distance-based score that JAX can handle
            # This is a placeholder - in production you'd implement restraints in JAX
            coords = coords_flat.reshape(-1, 3)
            
            # Simple pairwise distance penalty (JAX-compatible)
            score = 0.0
            for i in range(min(5, len(coords))):  # Limit for speed
                for j in range(i+1, min(5, len(coords))):
                    dist = jnp.linalg.norm(coords[i] - coords[j])
                    # Soft penalty for being too close
                    score = score + jnp.where(dist < 5.0, (5.0 - dist)**2, 0.0)
                    # Soft penalty for being too far
                    score = score + jnp.where(dist > 50.0, (dist - 50.0)**2, 0.0)
            
            return -score  # Negative log likelihood
        
        return likelihood

#--------------------------------------------------------------------------
# Probabilistic Sampler 
#--------------------------------------------------------------------------
class ProbabilisticSampler:
    """Interface with PyMC and NumPyro"""
    
    def __init__(self, bridge: IMPToProbabilistic, prior_std=10.0):
        self.bridge = bridge
        self.initial_coords = bridge.get_coordinates_array()
        self.n_particles = len(self.initial_coords)
        self.prior_std = prior_std
        
    def build_pymc_model(self):
        """Build PyMC model - simplified to avoid tracing issues"""
        initial_coords = self.initial_coords
        likelihood_fn = self.bridge.build_likelihood_function()
        
        with pm.Model() as model:
            # Prior on coordinates - centered around initial positions
            coords_flat = pm.Normal('coords',
                                   mu=initial_coords.flatten(),
                                   sigma=self.prior_std,
                                   shape=initial_coords.size)
            
            # We'll evaluate likelihood outside of PyMC graph
            # This is a workaround for black-box likelihood
            
        return model, likelihood_fn
    
    def build_numpyro_model(self):
        """Build NumPyro model with JAX-compatible likelihood"""
        if not HAS_JAX:
            raise RuntimeError("NumPyro/JAX not available")
            
        initial_coords = self.initial_coords
        likelihood_fn = self.bridge.build_jax_likelihood()
        
        def model():
            # Prior on coordinates
            coords_flat = numpyro.sample(
                'coords',
                dist.Normal(initial_coords.flatten(), self.prior_std)
            )
            
            # Likelihood (using JAX-compatible version)
            log_lik = likelihood_fn(coords_flat)
            numpyro.factor('likelihood', log_lik)
            
            return coords_flat
        
        return model

#--------------------------------------------------------------------------
# Sampling Interface - Fixed with better Monte Carlo
#--------------------------------------------------------------------------
class SamplingInterface:
    """Run sampling and integrate results with IMP"""
    
    def __init__(self, bridge: IMPToProbabilistic, sampler: ProbabilisticSampler):
        self.bridge = bridge
        self.sampler = sampler
        
    def run_pymc_sampling(self, draws=1000, chains=2, tune=500):
        """Run MCMC with PyMC - using black-box likelihood"""
        print(f"Running PyMC sampling: {draws} draws, {chains} chains")
        
        # For black-box likelihood, we'll use a custom sampler
        # PyMC's built-in samplers expect differentiable likelihoods
        print("PyMC with black-box likelihood is complex - using fallback")
        raise NotImplementedError("Use simple_mc or numpyro instead")
    
    def run_numpyro_sampling(self, num_samples=10000, num_warmup=500):
        """Run sampling with NumPyro using JAX-compatible likelihood"""
        if not HAS_JAX:
            raise RuntimeError("NumPyro/JAX not available")
            
        print(f"Running NumPyro sampling: {num_samples} samples")
        print("Note: Using simplified JAX-compatible likelihood")
        
        try:
            # Double-check we're using CPU
            import jax
            print(f"JAX devices: {jax.devices()}")
            
            model = self.sampler.build_numpyro_model()
            
            # Use Random Walk Metropolis for non-differentiable likelihood
            from numpyro.infer import HMC
            kernel = HMC(model)  # Try HMC with the JAX-compatible likelihood
            
            mcmc = MCMC(kernel, num_warmup=num_warmup, num_samples=num_samples)
            
            # Create RNG key - this is where the previous error occurred
            rng_key = jax.random.PRNGKey(0)
            mcmc.run(rng_key)
            
            samples = mcmc.get_samples()
            
            # Now evaluate the actual IMP score for each sample
            print("Evaluating IMP scores for NumPyro samples...")
            coords_samples = samples['coords'].reshape(-1, self.sampler.n_particles, 3)
            actual_scores = []
            
            for i, coords in enumerate(coords_samples):
                if i % 100 == 0:
                    print(f"  Evaluating sample {i}/{len(coords_samples)}")
                score = self.bridge.evaluate_score(np.array(coords))
                actual_scores.append(score)
            
            return coords_samples, np.array(actual_scores)
            
        except Exception as e:
            print(f"NumPyro sampling error: {e}")
            # Re-raise to trigger fallback
            raise
    
    def run_simple_mc(self, n_samples=1000, step_size=2.0, temperature=10.0):
        """Improved Monte Carlo sampler with adaptive step size"""
        print(f"Running simple Monte Carlo: {n_samples} samples")
        print(f"Initial step size: {step_size}, Temperature: {temperature}")
        
        coords_samples = []
        all_scores = []
        accept_count = 0
        
        current_coords = self.sampler.initial_coords.copy()
        current_score = self.bridge.evaluate_score(current_coords)
        
        best_coords = current_coords.copy()
        best_score = current_score
        
        print(f"Initial score: {current_score:.2f}")
        
        # Adaptive step size
        target_accept_rate = 0.5
        adapt_interval = 100
        
        for i in range(n_samples):
            # Propose move
            proposal = current_coords + np.random.normal(0, step_size, current_coords.shape)
            
            # Evaluate
            proposal_score = self.bridge.evaluate_score(proposal)
            
            # Metropolis criterion
            if proposal_score < current_score:
                accept = True
            else:
                delta = proposal_score - current_score
                accept_prob = np.exp(-delta / temperature)
                accept = np.random.random() < accept_prob
            
            if accept:
                current_coords = proposal.copy()
                current_score = proposal_score
                accept_count += 1
                
                if current_score < best_score:
                    best_score = current_score
                    best_coords = current_coords.copy()
            
            coords_samples.append(current_coords.copy())
            all_scores.append(current_score)
            
            # Adapt step size
            if (i + 1) % adapt_interval == 0:
                current_accept_rate = accept_count / (i + 1)
                if current_accept_rate < target_accept_rate - 0.1:
                    step_size *= 0.9  # Decrease step size
                elif current_accept_rate > target_accept_rate + 0.1:
                    step_size *= 1.1  # Increase step size
                
                if i % 100 == 0:
                    print(f"MC step {i}: score={current_score:.2f}, best={best_score:.2f}, "
                          f"accept_rate={current_accept_rate:.3f}, step_size={step_size:.3f}")
        
        final_accept_rate = accept_count / n_samples
        print(f"Final acceptance rate: {final_accept_rate:.3f}")
        print(f"Final step size: {step_size:.3f}")
        
        # Set to best configuration
        self.bridge.set_coordinates_from_array(best_coords)
        
        return np.array(coords_samples), np.array(all_scores)
    
    def analyze_samples(self, samples_data, output_best=True):
        """Analyze samples and optionally output best configuration"""
        if isinstance(samples_data, tuple) and len(samples_data) == 2:
            # Direct output from sampling (coords, scores)
            coords_samples, scores = samples_data
            if len(coords_samples.shape) == 2:
                coords_samples = coords_samples.reshape(-1, self.sampler.n_particles, 3)
        else:
            raise ValueError("Unexpected sample format")
        
        best_idx = np.argmin(scores)
        
        print(f"\nScore statistics:")
        print(f"  Best score: {scores[best_idx]:.2f}")
        print(f"  Mean score: {np.mean(scores):.2f}")
        print(f"  Median score: {np.median(scores):.2f}")
        print(f"  Std score: {np.std(scores):.2f}")
        print(f"  Min score: {np.min(scores):.2f}")
        print(f"  Max score: {np.max(scores):.2f}")
        
        if output_best:
            # Set to best configuration
            best_coords = coords_samples[best_idx]
            self.bridge.set_coordinates_from_array(best_coords)
            print("IMP model updated with best sampled configuration")
        
        return coords_samples, scores

#--------------------------------------------------------------------------
# Main Workflow - Fixed
#--------------------------------------------------------------------------
def run_probabilistic_sampling_workflow():
    """Complete workflow integrating IMP system with probabilistic sampling"""
    
    # 1. Build your IMP system
    print("="*60)
    print("Building IMP System")
    print("="*60)
    
    ntype = 3
    copy_numbers = [8, 8, 16]
    chain_ids = []
    for i in range(np.array(copy_numbers).sum()):
        chain_ids.append(chr(ord('A') + i))
    sequences_type = ['AA', 'L', 'A']
    colors = ['red', 'blue', 'green']
    
    builder = System_builder(ntype, copy_numbers, chain_ids, sequences_type, colors)
    mdl, system, parts = builder.build_system()
    
    hierarchy = system.build()
    IMP.atom.show_with_representations(hierarchy)
    
    # 2. Set up initial configuration and degrees of freedom
    print("\n" + "="*60)
    print("Setting up initial configuration and degrees of freedom")
    print("="*60)
    
    # Set up degrees of freedom for optimization
    dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
    
    # Add rigid body or flexible movers for each molecule
    for mol_name, mol in parts.items():
        # Option 1: Rigid body movement for entire molecule
        dof.create_rigid_body(mol, max_trans=5.0, max_rot=0.3)
        
        # Option 2 (alternative): Flexible beads
        #dof.create_flexible_beads(mol, max_trans=3.0)
    
    IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=20.0)
    
    # Save initial configuration
    output_initial = IMP.pmi.output.Output()
    output_initial.init_rmf("initial_config.rmf3", [hierarchy])
    output_initial.write_rmf("initial_config.rmf3")
    
    # 3. Set up scoring function with restraints
    print("\n" + "="*60)
    print("Setting up restraints")
    print("="*60)
    
    sf = ScoringFunction(mdl, hierarchy, parts)
    
    # Add excluded volume
    sf.add_excluded_volume_restraint()
    
    # Add some distance restraints
    sf.add_distance_restraint(
        mol1_spec=("prot_1", 1, 0),
        mol2_spec=("prot_2", 1, 0),
        distance_min=10.0,
        distance_max=30.0,
        kappa=1.0
    )
    
    sf.add_distance_restraint(
        mol1_spec=("prot_2", 1, 0),
        mol2_spec=("prot_3", 1, 0),
        distance_min=15.0,
        distance_max=35.0,
        kappa=0.5
    )
    
    sf.add_distance_restraint(
        mol1_spec=("prot_3", 1, 0),
        mol2_spec=("prot_1", 1, 0),
        distance_min=12.0,
        distance_max=32.0,
        kappa=0.8
    )
    sf.add_distance_restraint(
        mol1_spec=("prot_5", 1, 0),
        mol2_spec=("prot_6", 1, 0),
        distance_min=10.0,
        distance_max=30.0,
        kappa=1.0
    )
    initial_score = sf.get_total_score()
    print(f"Total initial score: {initial_score:.2f}")
    
    # 4. Create probabilistic sampling bridge
    print("\n" + "="*60)
    print("Setting up probabilistic sampling")
    print("="*60)
    
    bridge = IMPToProbabilistic(mdl, hierarchy, sf.restraints)
    sampler = ProbabilisticSampler(bridge, prior_std=5.0)
    interface = SamplingInterface(bridge, sampler)
    
    # 5. Run sampling
    print("\n" + "="*60)
    print("Running advanced sampling")
    print("="*60)
    
    sampling_method = ""
    samples_data = None
    
    # Try NumPyro if available
    if HAS_JAX:  # Enable NumPyro by default if available
        try:
            print("Trying NumPyro with JAX-compatible likelihood...")
            samples_data = interface.run_numpyro_sampling(num_samples=5000, num_warmup=200)
            sampling_method = "NumPyro"
        except Exception as e:
            print(f"NumPyro sampling failed: {e}")
            import traceback
            traceback.print_exc()
            samples_data = None
    
    # Fall back to simple Monte Carlo
    if samples_data is None:
        print("Using optimized Monte Carlo sampler")
        samples_data = interface.run_simple_mc(
            n_samples=10000, 
            step_size=3.0,  # Larger initial step
            temperature=10.0  # Adjust temperature for your system
        )
        sampling_method = "Monte Carlo"
    
    # Analyze results
    coords_samples, scores = interface.analyze_samples(samples_data, output_best=True)
    print(f"\nSuccessfully sampled using: {sampling_method}")
    
    # 6. Save best configuration
    print("\n" + "="*60)
    print("Saving results")
    print("="*60)
    
    output_best = IMP.pmi.output.Output()
    output_best.init_rmf("best_sampled_config.rmf3", [hierarchy])
    output_best.write_rmf("best_sampled_config.rmf3")
    
    score_after_sampling = sf.get_total_score()
    print(f"Score after sampling: {initial_score:.2f} -> {score_after_sampling:.2f}")
    
    # 7. Local refinement with IMP
    print("\n" + "="*60)
    print("Local refinement with IMP")
    print("="*60)
    
    try:
        # Get proper scoring function
        scoring_function = sf.get_scoring_function()
        
        # Try Monte Carlo optimization if Conjugate Gradients fails
        try:
            optimizer = IMP.core.ConjugateGradients(mdl)
            optimizer.set_scoring_function(scoring_function)
            
            print(f"Starting conjugate gradient optimization from score: {score_after_sampling:.2f}")
            optimizer.optimize(100)
            final_score = sf.get_total_score()
            print(f"Score after CG optimization: {score_after_sampling:.2f} -> {final_score:.2f}")
        except Exception as e:
            print(f"Conjugate gradients failed ({e}), trying Monte Carlo optimization...")
            
            # Use Monte Carlo optimizer instead
            mc = IMP.core.MonteCarlo(mdl)
            mc.set_scoring_function(scoring_function)
            mc.set_kt(1.0)  # Temperature
            
            # Add movers if available
            if 'dof' in locals():
                movers = dof.get_movers()
                for mover in movers:
                    mc.add_mover(mover)
                
                print(f"Starting Monte Carlo optimization with {len(movers)} movers")
                mc.optimize(1000)
                final_score = sf.get_total_score()
                print(f"Score after MC optimization: {score_after_sampling:.2f} -> {final_score:.2f}")
            else:
                print("No movers available for optimization")
                final_score = score_after_sampling
                
    except Exception as e:
        print(f"Local optimization error (non-critical): {e}")
        final_score = score_after_sampling
    
    # Save final optimized configuration
    output_final = IMP.pmi.output.Output()
    output_final.init_rmf("final_optimized_config.rmf3", [hierarchy])
    output_final.write_rmf("final_optimized_config.rmf3")
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Initial score:        {initial_score:.2f}")
    print(f"Best sampled score:   {np.min(scores):.2f}")
    print(f"After MC sampling:    {score_after_sampling:.2f}")
    print(f"After optimization:   {final_score:.2f}")
    print(f"Total improvement:    {initial_score - final_score:.2f}")
    
    return hierarchy, coords_samples, scores

# Run the complete workflow
if __name__ == "__main__":
    try:
        hierarchy, samples, scores = run_probabilistic_sampling_workflow()
        print("\n" + "="*60)
        print("Workflow completed successfully!")
        print("="*60)
        print("\nOutput files:")
        print("  - initial_config.rmf3: Initial shuffled configuration")
        print("  - best_sampled_config.rmf3: Best configuration from sampling")
        print("  - final_optimized_config.rmf3: Locally optimized configuration")
    except Exception as e:
        print(f"\nError in workflow: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)