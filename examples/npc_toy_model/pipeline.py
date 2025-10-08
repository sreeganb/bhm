from typing import List, Dict, Any, Tuple
from core.state import SystemState
from core.sigma import GMMSigmaProvider
import os
import json
import time
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from analysis.mcmc_diagnostics import run_mcmc_diagnostics
import random
import h5py

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline with analysis between stages"""
    
    def __init__(self, initial_state: SystemState, base_seed: int = 1234, init_mode: str = "random"):
        self.initial_state = initial_state
        self.stages = []
        self.base_seed = int(base_seed)
        self.init_mode = init_mode
                    
    def add_stage(self, sampler_function, n_steps: int = 1000, name: str = None, **sampler_kwargs):
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
            state=state.copy(),
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
 
    def _randomize_initial_positions(self, state: SystemState, rng: np.random.Generator, mode: str = "random") -> None:
        """Set or perturb positions per chain"""
        box = getattr(state, "box_size", 0.0)
        if not box or box <= 0 or not state.positions:
            return

        if mode == "random":
            for comp, coords in state.positions.items():
                state.positions[comp] = rng.uniform(0.0, box, size=(len(coords), 3)).astype(coords.dtype)
        elif mode == "jitter":
            sigma = 0.02 * box
            for comp, coords in state.positions.items():
                prop = coords + rng.normal(0.0, sigma, size=coords.shape)
                state.positions[comp] = np.clip(prop, 0.0, box).astype(coords.dtype)

    def _load_random_trajectory_frame(self, stage_output: str, n_chains: int) -> List[SystemState]:
        """Load last frame from randomly selected trajectories from previous stage"""
        chain_dirs = [d for d in os.listdir(stage_output) 
                     if d.startswith('chain_') and os.path.isdir(os.path.join(stage_output, d))]
        
        starting_states = []
        for i in range(n_chains):
            # Randomly select a chain directory
            selected_chain = np.random.choice(chain_dirs)
            h5_file = os.path.join(stage_output, selected_chain, "trajectory.h5")
            
            if os.path.exists(h5_file):
                state = self._load_last_frame(h5_file)
                starting_states.append(state)
            else:
                # Fallback to copying initial state
                starting_states.append(self.initial_state.copy())
        
        return starting_states

    def _load_last_frame(self, h5_file: str) -> SystemState:
        """Load the last frame from an HDF5 trajectory file"""
        state = self.initial_state.copy()
        
        try:
            with h5py.File(h5_file, 'r') as f:
                if 'trajectory' not in f:
                    return state
                
                traj = f['trajectory']
                state_names = sorted([k for k in traj.keys() if k.startswith('state_')])
                
                if not state_names:
                    return state
                
                # Get last state
                last_state = traj[state_names[-1]]
                
                # Load positions
                if 'positions' in last_state:
                    positions = last_state['positions']
                    for comp_name in positions.keys():
                        state.positions[comp_name] = np.array(positions[comp_name][:])
                
                # Load sigma values
                if 'sigma' in last_state:
                    sigma_grp = last_state['sigma']
                    for key in sigma_grp.attrs.keys():
                        state.sigma[key] = float(sigma_grp.attrs[key])
                        
        except Exception as e:
            print(f"    Warning: Failed to load {h5_file}: {e}")
        
        return state
    
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
    
    def _initialize_sigma_for_stage(self, states, stage_name):
        """Initialize sigma values for each state based on the current stage"""
        sampler_name_map = {
            'pair_sampling': 'PairSampler',
            'tetramer_sampling': 'TetramerSampler', 
            'octet_sampling': 'OctetSampler'
        }
        
        sampler_name = sampler_name_map.get(stage_name, 'PairSampler')
        sig_provider = GMMSigmaProvider(sampler_name=sampler_name, 
                                        prior_type='uniform',
                                        sequence_position=0)
        
        for chain_id, state in enumerate(states):
            sig_provider.initialize_sigma(
                state=state,
                sigma_source="uniform",
                sampler_name=sampler_name,
                spread_chains=True,
                chain_id=chain_id,
                n_chains=len(states)
            )
    
    def run(self, output_base: str = "output/pipeline", n_chains: int = 4) -> Dict[str, Any]:
        """Run the pipeline with parallel chains and analysis between stages"""
        os.makedirs(output_base, exist_ok=True)
        
        # Save config
        config = {
            'stages': [{'name': s['name'], 'n_steps': s['n_steps'], 
                       'kwargs': {k: str(v) for k, v in s['kwargs'].items()}} for s in self.stages],
            'n_chains': n_chains
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
            
            # Determine starting states
            if i == 0:
                # First stage: randomized initial states
                starting_states = []
                for chain_id in range(n_chains):
                    seed = self.base_seed + 10_000 + chain_id
                    rng = np.random.default_rng(seed)
                    st = self.initial_state.copy()
                    self._randomize_initial_positions(st, rng, mode=self.init_mode)
                    starting_states.append(st)
            else:
                # Subsequent stages: random trajectory frames from previous stage
                starting_states = self._load_random_trajectory_frame(prev_stage_output, n_chains)
            
            # Initialize sigma
            self._initialize_sigma_for_stage(starting_states, stage_name)
            
            # Run parallel chains
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