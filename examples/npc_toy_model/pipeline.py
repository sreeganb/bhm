# pipeline.py - Pipeline to run multiple samplers in sequence with analysis
from typing import List, Dict, Any, Tuple
from core.state import SystemState
from core.sigma import GMMSigmaProvider
from core.system import SystemBuilder
import os
import json
import time
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from analysis.mcmc_diagnostics import run_mcmc_diagnostics
import random
import h5py
from core.io_utils import load_trajectory_from_disk

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline with analysis between stages"""
    
    def __init__(self, initial_state: SystemState, 
                 base_seed: int = 1234, 
                 init_mode: str = "random"):
        self.initial_state = initial_state
        self.stages = []
        self.base_seed = int(base_seed)
        self.init_mode = init_mode
                    
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
            states.append(state)
        
        return states
    
    def _extract_params_from_state(self, state: SystemState) -> 'SystemParameters':
        """Extract SystemParameters from SystemState (simplified)"""
        from core.parameters import SystemParameters
        params = SystemParameters()
        params.box_size = state.box_size
        # Add other parameter extraction as needed
        return params
    
    def _initialize_sigma_for_stage(self, states: List[SystemState], stage_idx: int):
        """Initialize sigma values for each state based on the current stage"""
        sampler_name_map = {
            'pair_sampling': 'PairSampler',
            'tetramer_sampling': 'TetramerSampler', 
            'octet_sampling': 'OctetSampler',
            'full_sampling': 'FullSampler'
        }
        
        stage_name = self.stages[stage_idx]['name']
        sampler_name = sampler_name_map.get(stage_name, 'PairSampler')
        
        # Use stage_idx as sequence_position for proper prior selection
        sig_provider = GMMSigmaProvider(
            sampler_name=sampler_name, 
            prior_type='uniform',
            sequence_position=stage_idx  # Fixed: use actual stage index
        )
        
        for chain_id, state in enumerate(states):
            sig_provider.initialize_sigma(
                state=state,
                sigma_source="auto",  # Let provider decide based on sequence_position
                sampler_name=sampler_name
            )
    
    def run(self, output_base: str = "output/pipeline", 
            n_chains: int = 4) -> Dict[str, Any]:
        """Run the pipeline with parallel chains and analysis between stages"""
        os.makedirs(output_base, exist_ok=True)
        
        # Save config
        config = {
            'stages': [{'name': s['name'], 'n_steps': s['n_steps'], 
                       'kwargs': {k: str(v) for k, v in s['kwargs'].items()}} for s in self.stages],
            'n_chains': n_chains,
            'sampler_sequence': self.initial_state.sampler_sequence
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
            
            # Initialize sigma with proper sequence position
            self._initialize_sigma_for_stage(starting_states, i)
            
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