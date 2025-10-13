# pipeline.py - Pipeline to run multiple samplers in sequence with analysis
from typing import List, Dict, Any, Optional, Tuple
import os
import json
import time
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import random
import h5py
from sklearn.mixture import GaussianMixture
import warnings

from analysis.mcmc_diagnostics import run_mcmc_diagnostics
from core.state import SystemState
from core.sigma import create_sigma_prior, get_default_sigma_ranges, initialize_sigma_dict, DEFAULT_SIGMA_RANGES
from core.system import SystemBuilder
from core.parameters import SystemParameters

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline with analysis between stages"""
    
    def __init__(self, initial_state: SystemState, base_seed: int = 1234, 
                 init_mode: str = "random", prior_type: str = "uniform"):
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

    def _load_from_h5_last_frame(self, h5_file: str) -> Tuple[Dict[str, float], bool]:
        """
        Load sigma dict from the last frame of trajectory.h5.
        Returns (sigma_dict, success).
        """
        result = {}
        if not h5_file or not os.path.exists(h5_file):
            return result, False
        
        try:
            with h5py.File(h5_file, 'r') as f:
                if 'trajectory' not in f:
                    return result, False
                traj = f['trajectory']
                state_names = sorted([k for k in traj.keys() if k.startswith('state_')])
                if not state_names:
                    return result, False
                
                last_state = traj[state_names[-1]]
                if 'sigma' not in last_state:
                    return result, False
                
                sigma_grp = last_state['sigma']
                for key in sigma_grp.attrs.keys():
                    result[str(key)] = float(sigma_grp.attrs[key])
                
                return result, True
        except Exception:
            return {}, False

    def _extract_sigma_samples(self, h5_file: str, burn_in: float = 0.5) -> Dict[str, np.ndarray]:
        """Extract all sigma samples from H5 trajectory after burn-in"""
        sigma_samples = {}
        
        try:
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
            
            return {k: np.array(v) for k, v in sigma_samples.items()}
        except Exception:
            return {}

    def _create_initial_states(self, stage_idx: int, n_chains: int, 
                               prev_stage_output: str = None) -> List[SystemState]:
        """Create initial states for a stage using SystemBuilder"""
        states = []
        params = SystemParameters()
        params.box_size = self.initial_state.box_size
        
        for chain_id in range(n_chains):
            if stage_idx == 0:
                # First stage: random or ideal initialization
                source = self.init_mode if self.init_mode in ["ideal", "random"] else "random"
                builder = SystemBuilder(
                    params=params,
                    sampler_sequence=self.initial_state.sampler_sequence,
                    current_sampler=self.stages[stage_idx]['name'],
                    source=source
                )
                traj_file = None
            else:
                # Later stages: pick random chain from previous stage
                chain_dirs = [d for d in os.listdir(prev_stage_output) 
                              if d.startswith('chain_') and 
                              os.path.isdir(os.path.join(prev_stage_output, d))]
                selected_chain = np.random.choice(chain_dirs)
                traj_file = os.path.join(prev_stage_output, selected_chain, "trajectory.h5")
                
                builder = SystemBuilder(
                    params=params,
                    sampler_sequence=self.initial_state.sampler_sequence,
                    current_sampler=self.stages[stage_idx]['name'],
                    source="trajectory",
                    trajectory_file=traj_file,
                    frame=-1
                )
            
            state = builder.build()
            state._init_traj_file = traj_file  # Remember source for sigma loading
            states.append(state)
        
        return states

    def _initialize_sigma_for_stage(
        self,
        states: List[SystemState],
        stage_idx: int,
        prev_stage_output: Optional[str] = None
    ) -> None:
        """Initialize sigma values and attach sigma_prior to each state."""
        print(f"  Initializing sigma for stage {stage_idx+1}...")

        pair_types = list(get_default_sigma_ranges().keys())
        sigma_ranges = get_default_sigma_ranges()

        for chain_id, state in enumerate(states):
            state.sigma_range = sigma_ranges
            gmm_file = None
            loaded_from_last_frame = False

            # Stage > 0: try to load sigma from the same trajectory that provided positions
            if stage_idx > 0 and getattr(state, "_init_traj_file", None):
                chain_dir = os.path.dirname(state._init_traj_file)
                candidate_gmm = os.path.join(chain_dir, "gmm_posterior.json")
                if os.path.exists(candidate_gmm):
                    gmm_file = candidate_gmm

                # Load sigma from last frame
                last_sigma, success = self._load_from_h5_last_frame(state._init_traj_file)
                if success:
                    # Validate and clip
                    cleaned = {}
                    for pt in pair_types:
                        if pt in last_sigma and np.isfinite(last_sigma[pt]):
                            low, high = sigma_ranges[pt]
                            cleaned[pt] = float(np.clip(last_sigma[pt], low, high))
                    if len(cleaned) == len(pair_types):
                        state.sigma = cleaned
                        loaded_from_last_frame = True

            # Create prior
            state.sigma_prior = create_sigma_prior(
                pair_types=pair_types,
                sigma_ranges=sigma_ranges,
                gmm_file=gmm_file,
                prior_type=self.prior_type
            )

            # Fallback: sample from prior if not loaded
            if not loaded_from_last_frame:
                rng = np.random.default_rng(self.base_seed + stage_idx * 1000 + chain_id)
                try:
                    state.sigma = state.sigma_prior.initialize_sigma(rng)
                except Exception:
                    state.sigma = initialize_sigma_dict(
                        pair_types=pair_types,
                        sigma_ranges=sigma_ranges,
                        rng=rng,
                        prior_type=self.prior_type
                    )

            # Print summary
            if chain_id < 3:
                sigma_str = ", ".join(f"{k}={state.sigma[k]:.3f}" for k in pair_types)
                source_str = "last_frame" if loaded_from_last_frame else ("GMM" if gmm_file else self.prior_type)
                print(f"    Chain {chain_id+1}: {sigma_str} (source: {source_str})")

        if len(states) > 3:
            print(f"    ... and {len(states)-3} more chains")

    def _run_single_chain(self, args):
        """Run a single MCMC chain"""
        state, stage, chain_id, stage_output, n_steps, chain_seed = args
        self._seed_everything(chain_seed)

        chain_output = os.path.join(stage_output, f"chain_{chain_id}")
        os.makedirs(chain_output, exist_ok=True)

        final_state, trajectory_file = stage['function'](
            state=state,
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
        """Run multiple MCMC chains in parallel"""
        print(f"  Running {n_chains} parallel chains...")
        
        args_list = [
            (states[i], stage, i+1, stage_output, stage['n_steps'], 
             self.base_seed + 1000 * len(self.stages) + i + 1)
            for i in range(n_chains)
        ]
        
        max_workers = min(n_chains, mp.cpu_count())
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            return list(executor.map(self._run_single_chain, args_list))

    def _fit_robust_gmm(self, samples: np.ndarray) -> Dict:
        """Fit GMM with 1-3 components using BIC, handling low-variance cases"""
        samples = samples.reshape(-1, 1)
        n_samples = len(samples)
        
        # Check for sufficient variation
        unique_samples = np.unique(samples)
        sample_std = np.std(samples)
        
        if len(unique_samples) < 3 or sample_std < 1e-6:
            # Insufficient variation: single Gaussian
            mean_val = float(np.mean(samples))
            std_val = max(0.1, float(sample_std))
            return {
                'n_components': 1,
                'weights': [1.0],
                'means': [[mean_val]],
                'covariances': [[std_val**2]]
            }
        
        # Fit with adaptive component count
        max_components = min(3, max(1, n_samples // 10))
        best_gmm = None
        best_bic = np.inf
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            
            for n_comp in range(1, max_components + 1):
                try:
                    gmm = GaussianMixture(n_components=n_comp, random_state=42, 
                                         max_iter=100, tol=1e-4)
                    gmm.fit(samples)
                    
                    if len(np.unique(gmm.predict(samples))) >= n_comp:
                        bic = gmm.bic(samples)
                        if bic < best_bic:
                            best_bic = bic
                            best_gmm = gmm
                except Exception:
                    continue
        
        if best_gmm is None:
            mean_val = float(np.mean(samples))
            std_val = max(0.1, float(np.std(samples)))
            return {
                'n_components': 1,
                'weights': [1.0],
                'means': [[mean_val]],
                'covariances': [[std_val**2]]
            }
        
        return {
            'n_components': best_gmm.n_components,
            'weights': best_gmm.weights_.tolist(),
            'means': best_gmm.means_.tolist(),
            'covariances': best_gmm.covariances_.tolist()
        }

    def _save_chain_gmms(self, stage_output: str, chain_dirs: List[str]):
        """Save robust GMM parameters for each chain"""
        print(f"  Saving chain-specific GMMs...")
        
        for chain_dir in chain_dirs:
            h5_file = os.path.join(stage_output, chain_dir, "trajectory.h5")
            gmm_file = os.path.join(stage_output, chain_dir, "gmm_posterior.json")
            
            if not os.path.exists(h5_file):
                continue
            
            try:
                sigma_samples = self._extract_sigma_samples(h5_file, burn_in=0.5)
                gmm_params = {}
                
                for pair_type, samples in sigma_samples.items():
                    if len(samples) >= 5:
                        gmm_params[pair_type] = self._fit_robust_gmm(samples)
                
                # Fallback: ensure all pair types have entries
                for pair_type in ['AA', 'AB', 'BC']:
                    if pair_type not in gmm_params:
                        low, high = DEFAULT_SIGMA_RANGES.get(pair_type, (0.5, 5.0))
                        mid = (low + high) / 2
                        spread = (high - low) / 6
                        gmm_params[pair_type] = {
                            'n_components': 1,
                            'weights': [1.0],
                            'means': [[mid]],
                            'covariances': [[spread**2]]
                        }
                
                with open(gmm_file, 'w') as f:
                    json.dump(gmm_params, f, indent=2)
                
                print(f"    Saved GMM for {chain_dir}")
                
            except Exception as e:
                print(f"    Warning: Failed to save GMM for {chain_dir}: {e}")

    def _run_analysis(self, stage_output, stage_name):
        """Run MCMC convergence analysis and save GMMs"""
        print(f"  Running convergence analysis for {stage_name}...")
        
        try:
            chain_dirs = sorted([d for d in os.listdir(stage_output) 
                               if d.startswith('chain_') and 
                               os.path.isdir(os.path.join(stage_output, d))])
            
            print(f"  Found {len(chain_dirs)} chains to analyze")
            
            # Convert H5 to RMF3
            rmf_files = []
            for chain_dir in chain_dirs:
                h5_file = os.path.join(stage_output, chain_dir, "trajectory.h5")
                rmf_file = os.path.join(stage_output, chain_dir, "trajectory.rmf3")
                
                if os.path.exists(h5_file):
                    try:
                        from analysis import h5_to_rmf3
                        h5_to_rmf3.convert_hdf5_to_rmf3(h5_file, rmf_file)
                        rmf_files.append(rmf_file)
                    except Exception as e:
                        print(f"      Warning: RMF3 conversion failed for {chain_dir}: {e}")
            
            # Save GMMs for next stage
            self._save_chain_gmms(stage_output, chain_dirs)
            
            # Run diagnostics
            analysis_results = run_mcmc_diagnostics(stage_output, stage_name, 
                                                   chain_dirs, rmf_files)
            
            # Save results
            with open(os.path.join(stage_output, "analysis_results.json"), 'w') as f:
                json.dump(analysis_results, f, indent=2)
            
            # Print R-hat
            if 'rhat' in analysis_results:
                print(f"\n  R-hat values for {stage_name}:")
                for param, rhat_value in analysis_results['rhat'].items():
                    status = "✓" if rhat_value < 1.1 else "✗"
                    print(f"    {status} {param}: {rhat_value:.3f}")
            
            return analysis_results
            
        except Exception as e:
            print(f"  Warning: Analysis failed: {e}")
            return {}

    def run(self, output_base: str = "output", n_chains: int = 4) -> Dict[str, Any]:
        """Run the full pipeline"""
        os.makedirs(output_base, exist_ok=True)
        
        # Save config
        config = {
            'stages': [{'name': s['name'], 'n_steps': s['n_steps']} for s in self.stages],
            'n_chains': n_chains,
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
            
            # Create states, initialize sigma, run chains
            starting_states = self._create_initial_states(i, n_chains, prev_stage_output)
            self._initialize_sigma_for_stage(starting_states, i, prev_stage_output)
            stage_results = self._run_parallel_chains(starting_states, stage, 
                                                     stage_output, n_chains)
            
            elapsed = time.time() - start_time
            print(f"  Parallel sampling completed in {elapsed:.1f} seconds")
            
            # Analyze
            analysis_results = self._run_analysis(stage_output, stage_name)
            
            all_results[f"stage_{i+1}_{stage_name}"] = {
                'chain_results': stage_results,
                'analysis': analysis_results,
                'elapsed_time': elapsed
            }
            
            prev_stage_output = stage_output
            
            # Print summary
            print(f"  Stage {i+1} complete. Final sigma values:")
            for res in stage_results[:3]:
                sigma_str = ", ".join(f"{k}={v:.3f}" for k, v in res['final_sigma'].items())
                print(f"    Chain {res['chain_id']}: {sigma_str}")
            if len(stage_results) > 3:
                print(f"    ... and {len(stage_results)-3} more chains")
        
        print(f"\n=== Pipeline Complete ===")
        print(f"Results saved to: {output_base}")
        
        return all_results