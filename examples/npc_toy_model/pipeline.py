# pipeline.py - Pipeline to run multiple samplers in sequence with analysis
from typing import List, Dict, Any, Optional
from core.state import SystemState
from core.sigma import SigmaPrior, create_sigma_prior
from core.system import SystemBuilder
import os
import json
import time
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from analysis.mcmc_diagnostics import run_mcmc_diagnostics
import random
from core.parameters import SystemParameters

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline with analysis between stages"""
    
    def __init__(self, initial_state: SystemState, 
                 base_seed: int = 1234, 
                 init_mode: str = "random",
                 prior_type: str = "uniform"):
        self.initial_state = initial_state
        self.stages = []
        self.base_seed = int(base_seed)
        self.init_mode = init_mode
        self.prior_type = prior_type
                    
    def add_stage(self, sampler_function, n_steps: int = 1000, 
                  name: str = None, **sampler_kwargs):
        """Add a sampler stage to the pipeline"""
        if name is None:
            name = sampler_function.__name__.replace('run_', '')
        
        self.stages.append({
            'function': sampler_function,
            'n_steps': n_steps,
            'name': name,
            'kwargs': sampler_kwargs
        })
        return self
    
    def _run_single_chain(self, args):
        """Run a single MCMC chain - used for multiprocessing"""
        state, stage, chain_id, stage_output, n_steps, chain_seed = args
        self._seed_everything(chain_seed)

        chain_output = os.path.join(stage_output, f"chain_{chain_id}")
        os.makedirs(chain_output, exist_ok=True)

        final_state, trajectory_file = stage['function'](
            state=state,  # keep the prepared state (sigma & prior stay attached)
            n_steps=n_steps,
            output_dir=chain_output,
            **stage['kwargs']
        )

        return {
            'chain_id': chain_id,
            'final_state': final_state,
            'trajectory_file': trajectory_file,
            'final_sigma': {k: float(v) for k, v in final_state.sigma.items()}
        }
    def _run_parallel_chains(self, states, stage, stage_output, n_chains):
        """Run multiple MCMC chains in parallel for a stage"""
        print(f"  Running {n_chains} parallel chains...")
        
        args_list = [
            (states[i], stage, i+1, stage_output, stage['n_steps'], 
            self.base_seed + 1000 * len(self.stages) + i + 1)
            for i in range(n_chains)
        ]
        
        max_workers = min(n_chains, mp.cpu_count())
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            return list(executor.map(self._run_single_chain, args_list))

    def _seed_everything(self, seed: int) -> None:
        """Seed all RNGs for reproducibility"""
        np.random.seed(seed)
        random.seed(seed)
        try:
            import torch
            torch.manual_seed(seed)
            if torch.cuda.is_available():
                torch.cuda.manual_seed_all(seed)
        except:
            pass
    
    def _create_initial_states(self, stage_idx: int, n_chains: int, 
                            prev_stage_output: str = None) -> List[SystemState]:
        """Create initial states for a stage using SystemBuilder"""
        states = []
        
        for chain_id in range(n_chains):
            if stage_idx == 0:
                # First stage: use ideal or random based on init_mode
                source = self.init_mode if self.init_mode in ["ideal", "random"] else "random"
                builder = SystemBuilder(
                    params=self._extract_params_from_state(self.initial_state),
                    sampler_sequence=self.initial_state.sampler_sequence,
                    current_sampler=self.stages[stage_idx]['name'],
                    source=source
                )
            else:
                # Subsequent stages: load from random trajectory frames
                chain_dirs = [d for d in os.listdir(prev_stage_output) 
                            if d.startswith('chain_') and os.path.isdir(os.path.join(prev_stage_output, d))]
                selected_chain = np.random.choice(chain_dirs)
                traj_file = os.path.join(prev_stage_output, selected_chain, "trajectory.h5")
                
                builder = SystemBuilder(
                    params=self._extract_params_from_state(self.initial_state),
                    sampler_sequence=self.initial_state.sampler_sequence,
                    current_sampler=self.stages[stage_idx]['name'],
                    source="trajectory",
                    trajectory_file=traj_file,
                    frame=-1  # Last frame
                )
            
            state = builder.build()
            
            # Ensure sigma dict exists (may be empty from builder)
            if not hasattr(state, 'sigma') or state.sigma is None:
                state.sigma = {}
            
            states.append(state)
        
        return states
    
    def _extract_params_from_state(self, state: SystemState) -> 'SystemParameters':
        """Extract SystemParameters from SystemState (simplified)"""
        from core.parameters import SystemParameters
        params = SystemParameters()
        params.box_size = state.box_size
        # Add other parameter extraction as needed
        return params
    
    def _find_gmm_file(self, stage_idx: int, output_base: str) -> Optional[str]:
        """
        Find GMM file from previous stage for initializing sigma prior.
        
        Args:
            stage_idx: Current stage index
            output_base: Base output directory
            
        Returns:
            Path to GMM file or None if not found (first stage or missing file)
        """
        if stage_idx == 0:
            return None
        
        # Look for GMM file from previous stage
        prev_stage_name = self.stages[stage_idx - 1]['name']
        prev_stage_dir = os.path.join(output_base, f"stage_{stage_idx}_{prev_stage_name}")
        
        # Check for consolidated GMM file
        gmm_file = os.path.join(prev_stage_dir, "gmm_all_chains.json")
        if os.path.exists(gmm_file):
            return gmm_file
        
        # If not found, look in analysis subdirectory
        gmm_file_alt = os.path.join(prev_stage_dir, "analysis", "gmm_all_chains.json")
        if os.path.exists(gmm_file_alt):
            return gmm_file_alt
        
        print(f"  Warning: GMM file not found for stage {stage_idx}. Will use simple prior.")
        return None

    def _initialize_sigma_for_stage(
        self,
        states: List[SystemState],
        stage_idx: int,
        prev_stage_output: Optional[str] = None
    ) -> None:
        """Initialize sigma values and attach sigma_prior to each state"""
        print(f"  Initializing sigma for stage {stage_idx+1}...")

        default_ranges = {
            'AA': (1.0, 5.0),
            'AB': (1.5, 8.0),
            'BC': (2.0, 12.0),
        }
        pair_types = list(default_ranges.keys())

        for chain_id, state in enumerate(states):
            # make sure sigma dict and ranges exist
            if not state.sigma:
                state.sigma = {pt: 1.0 for pt in pair_types}
            if not state.sigma_range:
                state.sigma_range = default_ranges.copy()

            # choose prior source
            gmm_file = None
            if stage_idx > 0 and prev_stage_output:
                chain_dirs = [
                    d for d in os.listdir(prev_stage_output)
                    if d.startswith("chain_") and os.path.isdir(os.path.join(prev_stage_output, d))
                ]
                if chain_dirs:
                    selected = np.random.choice(chain_dirs)
                    candidate = os.path.join(prev_stage_output, selected, "gmm_posterior.json")
                    if os.path.exists(candidate):
                        gmm_file = candidate

            # attach prior and draw sigma
            state.sigma_prior = create_sigma_prior(state=state, gmm_file=gmm_file, prior_type=self.prior_type)

            rng = np.random.default_rng(self.base_seed + stage_idx * 1000 + chain_id)
            state.sigma.update(state.sigma_prior.initialize_sigma(rng=rng))

            if chain_id < 3:
                sigma_str = ", ".join(f"{k}={v:.3f}" for k, v in state.sigma.items())
                prior_str = "GMM" if gmm_file else self.prior_type
                print(f"    Chain {chain_id+1}: {sigma_str} (prior: {prior_str})")

        if len(states) > 3:
            print(f"    ... and {len(states)-3} more chains")
    
    def run(self, output_base: str = "output/pipeline", 
            n_chains: int = 4) -> Dict[str, Any]:
        """Run the pipeline with parallel chains and analysis between stages"""
        os.makedirs(output_base, exist_ok=True)
        
        # Save config
        config = {
            'stages': [{'name': s['name'], 'n_steps': s['n_steps'], 
                       'kwargs': {k: str(v) for k, v in s['kwargs'].items()}} for s in self.stages],
            'n_chains': n_chains,
            'sampler_sequence': self.initial_state.sampler_sequence,
            'prior_type': self.prior_type,
            'base_seed': self.base_seed
        }
        with open(os.path.join(output_base, "pipeline_config.json"), 'w') as f:
            json.dump(config, f, indent=2)
        
        all_results = {}
        prev_stage_output = None
        
        for i, stage in enumerate(self.stages):
            stage_name = stage['name']
            stage_output = os.path.join(output_base, f"stage_{i+1}_{stage_name}")
            os.makedirs(stage_output, exist_ok=True)
            
            print(f"\n=== Stage {i+1}/{len(self.stages)}: {stage_name} ===")
            print(f"Running {stage['n_steps']} steps with {n_chains} parallel chains")
            
            start_time = time.time()
            
            # Create initial states using unified SystemBuilder
            starting_states = self._create_initial_states(i, n_chains, prev_stage_output)
            
            # Initialize sigma and get GMM file paths
            self._initialize_sigma_for_stage(starting_states, i, prev_stage_output)
            stage_results = self._run_parallel_chains(starting_states, stage, stage_output, n_chains)
                      
            elapsed = time.time() - start_time
            print(f"  Parallel sampling completed in {elapsed:.1f} seconds")
            
            # Run analysis
            analysis_results = self._run_analysis(stage_output, stage_name)
            
            # Store results
            all_results[f"stage_{i+1}_{stage_name}"] = {
                'chain_results': stage_results,
                'analysis': analysis_results,
                'elapsed_time': elapsed
            }
            
            prev_stage_output = stage_output
            
            # Print summary
            print(f"  Stage {i+1} complete. Final sigma ranges:")
            for chain_result in stage_results[:3]:
                chain_id = chain_result['chain_id']
                sigma_vals = chain_result['final_sigma']
                print(f"    Chain {chain_id}: {', '.join([f'{k}={v:.3f}' for k, v in list(sigma_vals.items())[:3]])}")
            if len(stage_results) > 3:
                print(f"    ... and {len(stage_results)-3} more chains")
        
        # Save summary (without final_state objects for JSON serialization)
        with open(os.path.join(output_base, "results_summary.json"), 'w') as f:
            json_results = {
                stage_key: {
                    'chain_results': [{k: v for k, v in result.items() if k != 'final_state'} 
                                    for result in stage_data['chain_results']],
                    'analysis': stage_data['analysis'],
                    'elapsed_time': stage_data['elapsed_time']
                }
                for stage_key, stage_data in all_results.items()
            }
            json.dump(json_results, f, indent=2)
        
        print(f"\n=== Pipeline Complete ===")
        print(f"Results saved to: {output_base}")
        
        return all_results
    
    def _run_analysis(self, stage_output, stage_name):
        """Run MCMC analysis"""
        print(f"  Running convergence analysis for {stage_name}...")
        
        try:
            chain_dirs = sorted([d for d in os.listdir(stage_output) 
                               if d.startswith('chain_') and os.path.isdir(os.path.join(stage_output, d))])
            
            print(f"  Found {len(chain_dirs)} chains to analyze")
            
            # Convert H5 to RMF3
            rmf_files = []
            for chain_dir in chain_dirs:
                h5_file = os.path.join(stage_output, chain_dir, "trajectory.h5")
                rmf_file = os.path.join(stage_output, chain_dir, "trajectory.rmf3")
                
                if os.path.exists(h5_file):
                    print(f"    Converting {chain_dir}/trajectory.h5 to RMF3...")
                    try:
                        from analysis import h5_to_rmf3
                        h5_to_rmf3.convert_hdf5_to_rmf3(h5_file, rmf_file)
                        rmf_files.append(rmf_file)
                        print(f"      Converted to {os.path.basename(rmf_file)}")
                    except Exception as e:
                        print(f"      Warning: RMF3 conversion failed: {e}")
                        
            # Save chain-specific GMMs for next stage
            self._save_chain_gmms(stage_output, chain_dirs)
                        
            # Run diagnostics
            analysis_results = run_mcmc_diagnostics(stage_output, stage_name, chain_dirs, rmf_files)
            
            # Save and print results
            analysis_file = os.path.join(stage_output, "analysis_results.json")
            with open(analysis_file, 'w') as f:
                json.dump(analysis_results, f, indent=2)
            
            if 'rhat' in analysis_results:
                print(f"\n  R-hat values for {stage_name}:")
                for param, rhat_value in analysis_results['rhat'].items():
                    status = "✓" if rhat_value < 1.1 else "✗"
                    print(f"    {status} {param}: {rhat_value:.3f}")
            
            return analysis_results
            
        except Exception as e:
            print(f"  Warning: Analysis failed with error: {e}")
            return {}

    def _save_chain_gmms(self, stage_output: str, chain_dirs: List[str]):
        """Save GMM parameters for each chain for use in next stage"""
        print(f"  Saving chain-specific GMMs...")
        
        for chain_dir in chain_dirs:
            h5_file = os.path.join(stage_output, chain_dir, "trajectory.h5")
            gmm_file = os.path.join(stage_output, chain_dir, "gmm_posterior.json")
            
            if not os.path.exists(h5_file):
                continue
            
            try:
                # Extract sigma values from this chain's trajectory
                sigma_samples = self._extract_sigma_from_h5(h5_file, burn_in=0.5)
                
                # Fit GMM to each pair type
                from sklearn.mixture import GaussianMixture
                gmm_params = {}
                
                for pair_type, samples in sigma_samples.items():
                    if len(samples) < 10:
                        continue
                    
                    # Fit GMM (1-3 components)
                    best_gmm = None
                    best_bic = np.inf
                    
                    for n_comp in range(1, min(4, len(samples)//3)):
                        gmm = GaussianMixture(n_components=n_comp, random_state=42)
                        gmm.fit(samples.reshape(-1, 1))
                        bic = gmm.bic(samples.reshape(-1, 1))
                        
                        if bic < best_bic:
                            best_bic = bic
                            best_gmm = gmm
                    
                    if best_gmm:
                        gmm_params[pair_type] = {
                            'n_components': best_gmm.n_components,
                            'weights': best_gmm.weights_.tolist(),
                            'means': best_gmm.means_.tolist(),
                            'covariances': best_gmm.covariances_.tolist()
                        }
                
                # Save to JSON
                with open(gmm_file, 'w') as f:
                    json.dump(gmm_params, f, indent=2)
                
                print(f"    Saved GMM for {chain_dir}")
                
            except Exception as e:
                print(f"    Warning: Failed to save GMM for {chain_dir}: {e}")

    def _extract_sigma_from_h5(self, h5_file: str, burn_in: float = 0.5) -> Dict[str, np.ndarray]:
        """Extract sigma samples from H5 trajectory after burn-in"""
        sigma_samples = {}
        import h5py
        
        with h5py.File(h5_file, 'r') as f:
            if 'trajectory' not in f:
                return sigma_samples
            
            traj = f['trajectory']
            state_names = sorted([k for k in traj.keys() if k.startswith('state_')])
            
            # Apply burn-in
            start_idx = int(len(state_names) * burn_in)
            state_names = state_names[start_idx:]
            
            for state_name in state_names:
                state = traj[state_name]
                if 'sigma' in state:
                    sigma_grp = state['sigma']
                    for key in sigma_grp.attrs.keys():
                        if key not in sigma_samples:
                            sigma_samples[key] = []
                        sigma_samples[key].append(float(sigma_grp.attrs[key]))
        
        # Convert to numpy arrays
        return {k: np.array(v) for k, v in sigma_samples.items()}