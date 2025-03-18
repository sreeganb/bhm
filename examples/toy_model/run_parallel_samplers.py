import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from datetime import datetime
from typing import Dict, Type
from functools import partial
from tqdm import tqdm

from parameters import SystemParameters
from base_sampler import BaseMCSampler
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from visualization import visualize_3d_configuration

SAMPLER_MAP = {
    "PairSampler": PairSampler,
    "TetramerSampler": TetramerSampler,
}

def run_single_chain(
    chain_idx: int,
    sampler_class: Type[BaseMCSampler],
    config: Dict,
    base_output_dir: str,
    timestamp: str
) -> Dict:
    """Run a single MCMC chain in its own process."""
    # Create chain-specific output directory
    chain_dir = os.path.join(
        base_output_dir,
        f"{sampler_class.__name__.lower()}_{timestamp}",
        f"chain_{chain_idx}"
    )
    os.makedirs(chain_dir, exist_ok=True)
    
    # Initialize sampler
    sampler = sampler_class(config.get("use_sigma_dist", True))
    
    # Run the chain with chain-specific output directory
    #best_positions, trajectory, traj_file = sampler.run_mc(
    #    n_steps=config["n_steps"],
    #    save_freq=config["save_freq"],
    #    output_dir=chain_dir  # Pass unique directory for each chain
    #)
    best_positions, traj_file = sampler.run_mc(
        n_steps=config["n_steps"],
        save_freq=config["save_freq"],
        output_dir=chain_dir  # Pass unique directory for each chain
    )
    
    # Save configuration visualization (remove save_path parameter)
    params = SystemParameters()
    plt.figure(figsize=(10, 8))
    visualize_3d_configuration(
        best_positions, 
        params.radii, 
        f"Chain {chain_idx} Final Configuration"
    )
    fig_path = os.path.join(chain_dir, "final_configuration.png")
    plt.savefig(fig_path)
    plt.close()
    
    # Don't try to save trajectory visualization if function doesn't return a figure
    try:
        # Only call visualization if needed - note it may display directly
        # For debugging and visualization purposes
        print(f"Generating trajectory visualization for chain {chain_idx}")
        
# This part of the code snippet is handling file operations related to saving the trajectory data
# generated during the sampling process. Here's a breakdown of what it does:
        # Use plotly directly to create a visualization
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        # Create a simple test trajectory plot
        fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])
        for type_name in best_positions:
            # Add points for current positions
            fig.add_trace(
                go.Scatter3d(
                    x=best_positions[type_name][:, 0],
                    y=best_positions[type_name][:, 1],
                    z=best_positions[type_name][:, 2],
                    mode='markers',
                    marker=dict(
                        size=5,
                        color={
                            'A': 'red',
                            'B': 'blue',
                            'C': 'green'
                        }.get(type_name, 'gray'),
                    ),
                    name=type_name
                )
            )
        
        fig.update_layout(
            title=f"Final configuration for chain {chain_idx}",
            scene=dict(
                xaxis=dict(range=[0, params.box_size]),
                yaxis=dict(range=[0, params.box_size]),
                zaxis=dict(range=[0, params.box_size])
            )
        )
        
        # Save this figure instead
        traj_fig_path = os.path.join(chain_dir, "trajectory.html")
        fig.write_html(traj_fig_path)
        
    except Exception as e:
        print(f"Warning: Could not save trajectory visualization: {e}")
    
    return {
        'chain_idx': chain_idx,
        'output_dir': chain_dir,
        'best_positions': best_positions,
        'trajectory_file': traj_file,
        'final_sigma': sampler.sigma
    }

def run_parallel_sampling(
    sampler_class: Type[BaseMCSampler],
    config: Dict,
    output_folder: str = "output_analysis",
    n_processes: int = None
) -> Dict[int, Dict]:
    """Run multiple MCMC chains in parallel."""
    n_chains = config["n_chains"]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Set environment variable to disable HDF5 file locking (backup solution)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    # Prepare parallel execution
    n_processes = n_processes or mp.cpu_count()
    pool = mp.Pool(processes=min(n_processes, n_chains))
    
    # Create partial function with fixed arguments
    run_chain = partial(
        run_single_chain,
        sampler_class=sampler_class,
        config=config,
        base_output_dir=output_folder,
        timestamp=timestamp
    )
    
    # Run chains in parallel with progress bar
    print(f"Starting {n_chains} chains using {min(n_processes, n_chains)} processes...")
    results = {}
    
    try:
        with tqdm(total=n_chains, desc="Running chains") as pbar:
            for result in pool.imap_unordered(run_chain, range(n_chains)):
                chain_idx = result['chain_idx']
                results[chain_idx] = result
                pbar.update(1)
                pbar.set_postfix({'Chain': chain_idx})
    finally:
        pool.close()
        pool.join()
    
    return results

def analyze_chains(results: Dict[int, Dict], output_folder: str, sampler_name: str):
    """Analyze results across all chains and organize output files."""
    
    # check if sampler specific folder exists
    sampler_folder = os.path.join(output_folder, f"{sampler_name.lower()}_results")
    old_sampler_folder = os.path.join(output_folder, f"{sampler_name.lower()}_results_old")
    if os.path.exists(sampler_folder):
        if os.path.exists(old_sampler_folder):
            import shutil
            shutil.rmtree(old_sampler_folder)  # Remove old_sampler_folder if it exists        # move this folder and rename it as old_sampler_folder
        os.rename(sampler_folder, old_sampler_folder)
    os.makedirs(sampler_folder, exist_ok=True)
    
    # Copy/move trajectory files with consistent naming
    for chain_idx, result in results.items():
        # Define consistent filename format
        chain_name = f"chain_{chain_idx+1}"
        new_traj_file = os.path.join(sampler_folder, f"trajectory_{chain_name}.h5")
        
        # Copy trajectory file to centralized location
        src_file = result['trajectory_file']
        if os.path.exists(src_file):
            import shutil
            shutil.copy2(src_file, new_traj_file)
            print(f"Copied trajectory for {chain_name} to {new_traj_file}")
            
            # Update result with new centralized file path
            results[chain_idx]['centralized_trajectory'] = new_traj_file
    
    # Collect sigma values across chains
    sigma_values = pd.DataFrame([
        {**{'chain': f"chain_{idx+1}"}, **result['final_sigma']}
        for idx, result in results.items()
    ])
    
    # Save combined results
    sigma_values.to_csv(os.path.join(sampler_folder, "final_sigma_values.csv"), index=False)
    
    # Print summary statistics
    print("\nSummary of sigma values across chains:")
    print(sigma_values.describe())
    
    return sigma_values

def main():
    run_config = {
        "pair_sampler": {
            "run": False,
            #"run": True,
            "n_chains": 50,
            "n_steps": 1000000,
            "save_freq": 500,
            "use_sigma_dist": False, # Example of sampler-specific parameter
        },
        "tetramer_sampler": {
            "run": True,
            #"run": False,
            "n_chains": 8,
            "n_steps": 100000,
            "save_freq": 500,
            "use_sigma_dist": True,
        }
    }
    
    output_folder = "output_analysis"
    sampler_results = {}
    
    # Number of processes to use (None = use CPU count)
    #n_processes = None
    n_processes = 8
    
    #sampler_sequence = ["pair_sampler", "tetramer_sampler"] # Define the sequence of samplers to run
    #sampler_sequence = ["pair_sampler"] # Example: Run only pair sampler
    sampler_sequence = ["tetramer_sampler"] # Example: Run only tetramer sampler
    #sampler_sequence = ["pair_sampler", "tetramer_sampler", "octamer_sampler"] # Example: Run all, assuming octamer_sampler config is present
    
    for sampler_key in sampler_sequence:
        config = run_config[sampler_key]
        if config["run"]:
            sampler_name = sampler_key.replace("_sampler", "").capitalize() + "Sampler"
            sampler_class = SAMPLER_MAP[sampler_name]
            
            print(f"\nStarting parallel sampling with {sampler_name}...")
            results = run_parallel_sampling(
                sampler_class=sampler_class,
                config=config,
                output_folder=output_folder,
                n_processes=n_processes
            )
            sampler_results[sampler_key] = results
            
            # After all chains complete, organize files in the format from run_samplers.py
            sigma_values = analyze_chains(results, output_folder, sampler_name)
                        
            print(f"\n{sampler_name} parallel sampling complete.")
            print(f"Results saved in: {output_folder}")

if __name__ == "__main__":
    main()
