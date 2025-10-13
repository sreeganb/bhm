# run_simulation.py
from core.parameters import SystemParameters
from core.system import setup_system
from samplers.pair import run_pair_sampling
from samplers.tetramer import run_tetramer_sampling
from samplers.octet import run_octet_sampling
from pipeline import SamplerPipeline

def main():
    # Initialize system parameters
    params = SystemParameters()
    
    # Define sampler sequence
    sampler_sequence = ["pair_sampling", "tetramer_sampling", "octet_sampling"]
    
    # Setup initial system state with sequence
    system_state = setup_system(
        params=params,
        source="random",  # Will be overridden by pipeline for first stage
        sampler_sequence=sampler_sequence,
        current_sampler=sampler_sequence[0]  # First sampler
    )
    
    # Create and run pipeline
    pipeline = SamplerPipeline(system_state, prior_type="inv_gamma")
    
    # Add sampling stages (automatically uses sequence for initialization)
    pipeline.add_stage(
        run_pair_sampling,
        n_steps=50000,
        save_freq=50,
        temp_start=10.0,
        temp_end=1.0
    )
    
    # Uncomment for additional stages
    pipeline.add_stage(
        run_tetramer_sampling,
        n_steps=40000,
        save_freq=40,
        temp_start=10.0,
        temp_end=1.0
    )
    
    # Uncomment for additional stages
    pipeline.add_stage(
        run_octet_sampling,
        n_steps=20000,
        save_freq=40,
        temp_start=10.0,
        temp_end=1.0
    )    
    
    # pipeline.add_stage(
    #     run_octet_sampling,
    #     n_steps=3000,
    #     save_freq=100,
    #     temp_start=3.0,
    #     temp_end=1.0
    # )
    
    # Run the pipeline with multiple chains
    results = pipeline.run(
        output_base="output",
        n_chains=8  # Run 8 parallel chains for each stage
    )
    
    print("Simulation complete!")

if __name__ == "__main__":
    main()
