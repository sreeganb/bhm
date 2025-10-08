# run_simulation.py
from core.parameters import SystemParameters
from core.system import setup_system
from samplers.pair import run_pair_sampling
#from samplers.tetramer import run_tetramer_sampling
#from samplers.octet import run_octet_sampling
from pipeline import SamplerPipeline

def main():
    # Initialize system parameters
    params = SystemParameters()  # Fixed: removed 'Torch' suffix
    
    # Setup initial system state
    system_state = setup_system(
        params=params,
        source="ideal"  # or "random" or "trajectory"
        # trajectory_file="path/to/trajectory.h5"
    )
    
    # Create and run pipeline
    pipeline = SamplerPipeline(system_state)
    
    # Add sampling stages
    pipeline.add_stage(
        run_pair_sampling,
        n_steps=100,
        save_freq=1,
        temp_start=5.0,
        temp_end=1.0
    )
    
#    pipeline.add_stage(
#        run_tetramer_sampling,
#        n_steps=2000,
#        save_freq=100,
#        temp_start=5.0,
#        temp_end=1.0
#    )
    
#    pipeline.add_stage(
#        run_octet_sampling,
#        n_steps=3000,
#        save_freq=100,
#        temp_start=3.0,
#        temp_end=1.0
#    )
    
    # Run the pipeline with multiple chains
    results = pipeline.run(
        output_base="output",
        n_chains=8  # Run 8 parallel chains for each stage
    )
    
    print("Simulation complete!")

if __name__ == "__main__":
    main()
