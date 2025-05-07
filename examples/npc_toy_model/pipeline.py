# pipeline.py
from typing import List, Dict, Any, Tuple
from core.state import SystemState
import os
import json
import time

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline"""
    
    def __init__(self, initial_state: SystemState):
        self.initial_state = initial_state
        self.stages = []
    
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
    
    def run(
        self, 
        output_base: str = "output/pipeline", 
        n_chains: int = 1
    ) -> Dict[str, Any]:
        """Run the pipeline for multiple chains"""
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
        
        results = {}
        
        # Run each chain
        for chain in range(1, n_chains + 1):
            chain_output = os.path.join(output_base, f"chain_{chain}")
            os.makedirs(chain_output, exist_ok=True)
            
            # Start with initial state for each chain
            current_state = self.initial_state.copy()
            chain_results = []
            
            print(f"\nRunning chain {chain}/{n_chains}...")
            
            # Run each stage in sequence
            for i, stage in enumerate(self.stages):
                stage_name = stage['name']
                stage_output = os.path.join(chain_output, f"{i+1}_{stage_name}")
                
                print(f"  Stage {i+1}/{len(self.stages)}: {stage_name} ({stage['n_steps']} steps)")
                start_time = time.time()
                
                # Run the sampler
                sampler_fn = stage['function']
                final_state, trajectory_file = sampler_fn(
                    state=current_state,
                    n_steps=stage['n_steps'],
                    output_dir=stage_output,
                    **stage['kwargs']
                )
                
                elapsed = time.time() - start_time
                print(f"  Completed in {elapsed:.1f} seconds")
                
                # Use the final state as input to the next stage
                current_state = final_state
                
                # Record results
                chain_results.append({
                    'stage': stage_name,
                    'trajectory': trajectory_file,
                    'final_sigma': {k: float(v) for k, v in current_state.sigma.items()}
                })
            
            results[f"chain_{chain}"] = chain_results
            
            print(f"Chain {chain} complete. Final sigma values:")
            for k, v in current_state.sigma.items():
                print(f"  {k}: {v:.4f}")
        
        # Save overall results summary
        with open(os.path.join(output_base, "results_summary.json"), 'w') as f:
            json.dump(results, f, indent=2)
        
        return results