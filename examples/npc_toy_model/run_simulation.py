# run_simulation.py
from core.parameters import SystemParameters
from core.system import SystemBuilder
# Updated imports to use new utility functions
from core.sigma import initialize_sigma, calculate_negative_log_prior
from samplers.pair import run_pair_sampling
from samplers.tetramer import run_tetramer_sampling
from samplers.octet import run_octet_sampling
from pipeline import SamplerPipeline

def setup_system(params, sampler_sequence, current_sampler_idx, source=None, trajectory_file=None):
    """
    Setup the system state based on sampler position in sequence
    
    Args:
        params: SystemParameters instance
        sampler_sequence: List of sampler names in execution order
        current_sampler_idx: Index of current sampler in the sequence
        source: Optional override for position source ("ideal", "random", "trajectory")
        trajectory_file: Optional explicit trajectory file path
    """
    current_sampler = sampler_sequence[current_sampler_idx]
    
    # Automatically determine source if not specified
    if source is None:
        if current_sampler_idx == 0:
            # First sampler uses random positions by default
            source = "random"
        else:
            # Non-first samplers use trajectory from previous
            source = "trajectory"
    
    # Create and return system state
    builder = SystemBuilder(
        params=params,
        sampler_sequence=sampler_sequence,
        current_sampler=current_sampler,
        source=source,
        trajectory_file=trajectory_file
    )
    
    return builder.build()

def create_pipeline(sampler_configs, params=None, prior_type='jeffreys'):
    """
    Create a complete sampling pipeline from a sequence of sampler configurations
    """
    if params is None:
        params = SystemParameters()
    
    # Map sampler names to their functions
    sampler_functions = {
        "PairSampler": run_pair_sampling,
        "TetramerSampler": run_tetramer_sampling,
        "OctetSampler": run_octet_sampling
    }
    
    # Extract just the sequence of sampler names
    sampler_sequence = [config["name"] for config in sampler_configs]
    
    # Setup system state for the first sampler
    first_sampler_idx = 0
    first_config = sampler_configs[first_sampler_idx]
    initial_state = setup_system(
        params=params,
        sampler_sequence=sampler_sequence,
        current_sampler_idx=first_sampler_idx,
        source=first_config.get("source", "random")
    )
    
    # Initialize sigma values for the first sampler
    initialize_sigma(
        state=initial_state,
        sampler_sequence=sampler_sequence,
        current_sampler_idx=first_sampler_idx,
        sigma_source="default",
        prior_type=first_config.get("prior_type", prior_type)
    )
    
    # Create pipeline with the initial state
    pipeline = SamplerPipeline(initial_state)
    
    # Add each sampler to the pipeline
    for idx, config in enumerate(sampler_configs):
        sampler_name = config["name"]
        sampler_func = sampler_functions[sampler_name]
        
        # Only pass a system_state to the first stage
        # For subsequent stages, pipeline will automatically use previous stage's output
        system_state = None
        if idx == 0:
            system_state = initial_state
        
        # Add this stage to the pipeline
        pipeline.add_stage(
            sampler_func,
            system_state=system_state,  # Only non-None for the first sampler
            n_steps=config.get("n_steps", 1000),
            save_freq=config.get("save_freq", 100),
            temp_start=config.get("temp_start", 10.0),
            temp_end=config.get("temp_end", 1.0),
            prior_type=config.get("prior_type", prior_type)
        )
    
    return pipeline

def main():
    # Define the sequence of samplers with their configurations
    sampler_configs = [
        {
            "name": "PairSampler",
            "n_steps": 50000,
            "save_freq": 100,
            "temp_start": 5.0,
            "temp_end": 0.1,
            "prior_type": "jeffreys"  # Can specify prior type per sampler
        },
#        {
#            "name": "TetramerSampler",
#            "n_steps": 2000,
#            "save_freq": 100,
#            "temp_start": 5.0,
#            "temp_end": 1.0
#        },
#        {
#            "name": "OctetSampler",
#            "n_steps": 3000,
#            "save_freq": 100,
#            "temp_start": 3.0,
#            "temp_end": 1.0
#        }
    ]
    
    # Create the pipeline with these samplers
    pipeline = create_pipeline(sampler_configs)
    
    # Run the pipeline
    results = pipeline.run(
        output_base="output/final_simulation",
        n_chains=2
    )
    
    print("Simulation complete!")

# Example of a different sampler sequence
def alternative_sequence_example():
    # You can easily create any sequence you want
    sampler_configs = [
        {"name": "TetramerSampler", "n_steps": 1500, "temp_start": 8.0},
        {"name": "OctetSampler", "n_steps": 2000, "temp_start": 4.0},
        {"name": "PairSampler", "n_steps": 1000, "temp_start": 2.0, "prior_type": "jeffreys"},
        {"name": "TetramerSampler", "n_steps": 2000, "temp_start": 1.5},
        {"name": "PairSampler", "n_steps": 1000, "temp_start": 1.0, "prior_type": "halfcauchy"}
    ]
    
    pipeline = create_pipeline(sampler_configs)
    # Run as needed...
    pipeline.run(output_base="output/alternative_simulation", n_chains=4)

if __name__ == "__main__":
    main()
    # Or run alternative_sequence_example() for a different sequence
