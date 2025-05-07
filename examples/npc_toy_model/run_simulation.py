# run_simulation.py
from parameters import SystemParametersTorch
from core.system import setup_system
from core.sigma import initialize_sigma
from samplers.pair import run_pair_sampling
from samplers.tetramer import run_tetramer_sampling
from samplers.octet import run_octet_sampling
from pipeline import SamplerPipeline

def main():
    # Initialize system parameters
    params = SystemParametersTorch()
    
    # Setup initial system state
    system_state = setup_system(
        params=params,
        source="ideal"  # or "random" or "trajectory"
        # trajectory_file="path/to/trajectory.h5"
    )
    
    # Initialize sigma values
    initialize_sigma(
        state=system_state,
        sigma_source="gmm",
        sampler_name="PairSampler"
    )
    
    # Create and run pipeline
    pipeline = SamplerPipeline(system_state)
    
    # Add sampling stages
    pipeline.add_stage(
        run_pair_sampling,
        n_steps=1000,
        save_freq=100,
        temp_start=10.0,
        temp_end=1.0
    )
    
    pipeline.add_stage(
        run_tetramer_sampling,
        n_steps=2000,
        save_freq=100,
        temp_start=5.0,
        temp_end=1.0
    )
    
    pipeline.add_stage(
        run_octet_sampling,
        n_steps=3000,
        save_freq=100,
        temp_start=3.0,
        temp_end=1.0
    )
    
    # Run the pipeline with multiple chains
    results = pipeline.run(
        output_base="output/final_simulation",
        n_chains=8
    )
    
    print("Simulation complete!")

if __name__ == "__main__":
    main()