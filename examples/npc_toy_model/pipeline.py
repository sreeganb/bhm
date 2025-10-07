# pipeline.py
from typing import List, Dict, Any, Tuple
from core.state import SystemState
from core.sigma import GMMSigmaProvider
import os
import json
import time
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import importlib

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline with analysis between stages"""
    
    def __init__(self, initial_state: SystemState):
        self.initial_state = initial_state
        self.stages = []
        
        # Create base analysis directory structure
        base_analysis_dir = "output_analysis"
        os.makedirs(base_analysis_dir, exist_ok=True)
        
        # Create placeholder directories for each sampler type
        for sampler_type in ["pairsampler_results", "tetramersampler_results", "octetsampler_results"]:
            sampler_dir = os.path.join(base_analysis_dir, sampler_type)
            os.makedirs(sampler_dir, exist_ok=True)
    
    def add_stage(
        self, 
        sampler_function, 
        n_steps: int = 1000, 
        name: str = None,
        **sampler_kwargs
    ):
        """Add a sampler stage to the pipeline"""
        if name is None:
            name = sampler_function.__name__.replace('run_', '')
        
        self.stages.append({
            'function': sampler_function,
            'n_steps': n_steps,
            'name': name,
            'kwargs': sampler_kwargs
        })
        
        return self  # Enable method chaining
    
    def _run_single_chain(self, args):
        """Run a single MCMC chain - used for multiprocessing"""
        state, stage, chain_id, stage_output, n_steps = args
        
        # Create chain-specific output directory
        chain_output = os.path.join(stage_output, f"chain_{chain_id}")
        os.makedirs(chain_output, exist_ok=True)
        
        # Run the sampler
        sampler_fn = stage['function']
        final_state, trajectory_file = sampler_fn(
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
        
        # Prepare arguments for each chain
        args_list = []
        for chain_id in range(1, n_chains + 1):
            # Use the appropriate state for each chain
            state = states[chain_id - 1] if len(states) > 1 else states[0]
            args_list.append((state, stage, chain_id, stage_output, stage['n_steps']))
        
        # Run chains in parallel
        max_workers = min(n_chains, mp.cpu_count())
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(self._run_single_chain, args_list))
        
        return results
    
    def _run_analysis(self, stage_output, stage_name):
        """Run MCMC analysis using analysis.py"""
        print(f"  Running convergence analysis for {stage_name}...")
        
        try:
            # Import analysis module
            analysis = importlib.import_module('analysis')
            
            # Run R-hat, cross-correlation, and other MCMC diagnostics
            analysis_results = analysis.run_mcmc_diagnostics(
                stage_output, 
                stage_name=stage_name
            )
            
            # Save analysis results
            analysis_file = os.path.join(stage_output, "analysis_results.json")
            with open(analysis_file, 'w') as f:
                json.dump(analysis_results, f, indent=2)
            
            print(f"  Analysis complete. R-hat values saved to {analysis_file}")
            return analysis_results
            
        except ImportError:
            print(f"  Warning: analysis.py not found. Skipping convergence analysis.")
            return {}
        except Exception as e:
            print(f"  Warning: Analysis failed with error: {e}")
            return {}
    
    def _select_starting_states(self, previous_results, n_chains):
        """Randomly select starting states from previous stage results"""
        if not previous_results:
            return [self.initial_state] * n_chains
        
        # Extract final states from previous results
        available_states = [result['final_state'] for result in previous_results]
        
        # Randomly select states for the new chains
        selected_states = []
        for _ in range(n_chains):
            selected_state = np.random.choice(available_states).copy()
            selected_states.append(selected_state)
        
        return selected_states
    
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
        
        for state in states:
            # Initialize sigma values for this sampler type
            sig_provider.initialize_sigma(
                state=state,
                sigma_source="gmm",
                sampler_name=sampler_name
            )
    
    def run(
        self, 
        output_base: str = "output/pipeline", 
        n_chains: int = 4
    ) -> Dict[str, Any]:
        """Run the pipeline with parallel chains and analysis between stages"""
        os.makedirs(output_base, exist_ok=True)
        
        # Save pipeline configuration
        config = {
            'stages': [
                {
                    'name': stage['name'],
                    'n_steps': stage['n_steps'],
                    'kwargs': {k: str(v) for k, v in stage['kwargs'].items()}
                }
                for stage in self.stages
            ],
            'n_chains': n_chains
        }
        
        with open(os.path.join(output_base, "pipeline_config.json"), 'w') as f:
            json.dump(config, f, indent=2)
        
        all_results = {}
        previous_results = []
        
        # Run each stage in sequence
        for i, stage in enumerate(self.stages):
            stage_name = stage['name']
            stage_output = os.path.join(output_base, f"stage_{i+1}_{stage_name}")
            os.makedirs(stage_output, exist_ok=True)
            
            print(f"\n=== Stage {i+1}/{len(self.stages)}: {stage_name} ===")
            print(f"Running {stage['n_steps']} steps with {n_chains} parallel chains")
            
            start_time = time.time()
            
            # Determine starting states for this stage
            if i == 0:
                # First stage: use initial state for all chains
                starting_states = [self.initial_state] * n_chains
                # Initialize sigma for first stage
                self._initialize_sigma_for_stage(starting_states, stage_name)
            else:
                # Subsequent stages: randomly select from previous results
                starting_states = self._select_starting_states(previous_results, n_chains)
                # Initialize sigma for current stage
                self._initialize_sigma_for_stage(starting_states, stage_name)
            
            # Run parallel chains for this stage
            stage_results = self._run_parallel_chains(
                starting_states, stage, stage_output, n_chains
            )
            
            elapsed = time.time() - start_time
            print(f"  Parallel sampling completed in {elapsed:.1f} seconds")
            
            # Run MCMC analysis
            analysis_results = self._run_analysis(stage_output, stage_name)
            
            # Store results for this stage
            all_results[f"stage_{i+1}_{stage_name}"] = {
                'chain_results': stage_results,
                'analysis': analysis_results,
                'elapsed_time': elapsed
            }
            
            # Update previous_results for next stage
            previous_results = stage_results
            
            # Print summary
            print(f"  Stage {i+1} complete. Final sigma ranges:")
            for chain_result in stage_results[:3]:  # Show first 3 chains
                chain_id = chain_result['chain_id']
                sigma_vals = chain_result['final_sigma']
                print(f"    Chain {chain_id}: {', '.join([f'{k}={v:.3f}' for k, v in list(sigma_vals.items())[:3]])}")
            if len(stage_results) > 3:
                print(f"    ... and {len(stage_results)-3} more chains")
        
        # Save overall results summary
        with open(os.path.join(output_base, "results_summary.json"), 'w') as f:
            # Remove final_state objects for JSON serialization
            json_results = {}
            for stage_key, stage_data in all_results.items():
                json_results[stage_key] = {
                    'chain_results': [
                        {k: v for k, v in result.items() if k != 'final_state'}
                        for result in stage_data['chain_results']
                    ],
                    'analysis': stage_data['analysis'],
                    'elapsed_time': stage_data['elapsed_time']
                }
            json.dump(json_results, f, indent=2)
        
        print(f"\n=== Pipeline Complete ===")
        print(f"Results saved to: {output_base}")
        
        return all_results