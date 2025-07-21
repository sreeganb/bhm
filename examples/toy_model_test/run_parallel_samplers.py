import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from datetime import datetime
from typing import Dict, Type, Any, List, Tuple
from functools import partial
from tqdm import tqdm
import shutil
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import shutil
import subprocess
import time

from parameters import SystemParameters
from base_sampler import BaseMCSampler
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from visualization import visualize_3d_configuration
from octet_sampler import OctetSampler

# Map of sampler names to their classes
SAMPLER_MAP = {
    "pair": PairSampler,
    "tetramer": TetramerSampler,
    "octet": OctetSampler
}

# Default configuration
DEFAULT_CONFIG = {
    "n_chains": 4,
    "n_steps": 10000,
    "save_freq": 100,
    "use_sigma_dist": True,
}

class SamplerSequenceManager:
    """Manages the sequence of samplers and their indexing."""
    
    def __init__(self, sampler_sequence: List[str]):
        self.sampler_sequence = sampler_sequence
        self.sampler_counts = self._count_samplers()
        self.current_indices = {name: 0 for name in set(sampler_sequence)}
    
    def _count_samplers(self) -> Dict[str, int]:
        """Count total occurrences of each sampler in the sequence."""
        counts = {}
        for sampler in self.sampler_sequence:
            counts[sampler] = counts.get(sampler, 0) + 1
        return counts
    
    def get_sampler_info(self, sequence_idx: int) -> Tuple[str, int, int]:
        """
        Get sampler information for a given position in the sequence.
        
        Returns:
            Tuple of (sampler_name, current_index, total_count)
        """
        sampler_name = self.sampler_sequence[sequence_idx]
        self.current_indices[sampler_name] += 1
        current_idx = self.current_indices[sampler_name]
        total_count = self.sampler_counts[sampler_name]
        
        return sampler_name, current_idx, total_count
    
    def get_directory_name(self, sampler_name: str, current_idx: int) -> str:
        """Generate directory name for the sampler."""
        return f"{sampler_name}sampler_results_{current_idx}"

def run_analysis(sampler_key: str, sampler_index: int, sampler_sequence: List[str], sequence_position: int):
    """Run the fit_gmm.py script for the specified sampler with sequence information."""
    
    print(f"\nRunning analysis for {sampler_key} sampler (occurrence #{sampler_index}) at sequence position {sequence_position}...")
    try:
        # Convert sampler sequence to comma-separated string
        sequence_str = ",".join(sampler_sequence)
        
        # Build command with the new arguments
        cmd = [
            "python", 
            "fit_gmm.py",
            "--sequence", sequence_str,
            "--sampler", sampler_key,
            "--position", str(sequence_position),
            "--burnin", "0.4"
        ]
        
        print(f"Executing command: {' '.join(cmd)}")
        
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Analysis failed with exit code {process.returncode}")
            print(f"Error output: {stderr}")
            if stdout:
                print(f"Standard output: {stdout}")
            return False
        
        if stdout:
            print(f"Analysis output: {stdout}")
        
        print(f"Analysis completed for {sampler_key} sampler (occurrence #{sampler_index})")
        return True
        
    except Exception as e:
        print(f"Error running analysis: {e}")
        return False

def visualize_final_configuration(positions: Dict[str, np.ndarray], 
                               params: SystemParameters, 
                               chain_idx: int, 
                               output_dir: str) -> None:
    """Create and save a static and interactive visualization of final positions"""
    # Static matplotlib plot
    plt.figure(figsize=(10, 8))
    visualize_3d_configuration(positions, params.radii, f"Chain {chain_idx} Final Configuration")
    plt.savefig(os.path.join(output_dir, "final_configuration.png"))
    plt.close()

def run_single_chain(chain_idx: int,
                    sampler_class: Type[BaseMCSampler],
                    config: Dict[str, Any],
                    base_output_dir: str,
                    timestamp: str,
                    sampler_sequence: List[str],
                    sequence_position: int) -> Dict[str, Any]:
    """Run a single MCMC chain in its own process."""
    
    # Setup chain-specific directory
    chain_dir = os.path.join(
        base_output_dir,
        f"{sampler_class.__name__.lower()}_{timestamp}",
        f"chain_{chain_idx}"
    )
    os.makedirs(chain_dir, exist_ok=True)
    
    # Initialize sampler with sequence information
    sampler = sampler_class(
        sampler_sequence=sampler_sequence,
        sequence_idx=sequence_position,
        # ex_weight=1.0,
        # pair_weight=1.0,
        # base_output_dir=None,
        # specific_chain=None,
        # sigma_ranges=None,
        prior_type="inverse_gamma"
        # pos_passed=None 
    )
    
    try:
        best_positions, traj_file = sampler.run_mc(
            n_steps=config["n_steps"],
            save_freq=config["save_freq"],
            output_dir=chain_dir
        )
        
        # Generate visualizations
        params = SystemParameters()
            
        return {
            'chain_idx': chain_idx,
            'output_dir': chain_dir,
            'best_positions': best_positions,
            'trajectory_file': traj_file,
            'final_sigma': sampler.sigma
        }
    
    except Exception as e:
        print(f"Error in chain {chain_idx}: {e}")
        return {
            'chain_idx': chain_idx,
            'output_dir': chain_dir,
            'error': str(e)
        }

def run_parallel_sampling(sampler_class: Type[BaseMCSampler],
                         config: Dict[str, Any],
                         output_folder: str = "output_analysis",
                         n_processes: int = None,
                         sampler_sequence: List[str] = None,
                         sequence_position: int = 0) -> Dict[int, Dict[str, Any]]:
    """Run multiple MCMC chains in parallel."""
    
    # Setup
    n_chains = config["n_chains"]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # Prevent HDF5 file locking issues
    os.makedirs(output_folder, exist_ok=True)
    
    # Prepare parallel execution
    n_processes = min(n_processes or mp.cpu_count(), n_chains)
    
    # Create partial function with fixed arguments
    run_chain = partial(
        run_single_chain,
        sampler_class=sampler_class,
        config=config,
        base_output_dir=output_folder,
        timestamp=timestamp,
        sampler_sequence=sampler_sequence or [],
        sequence_position=sequence_position
    )
    
    # Run chains in parallel with progress bar
    print(f"Starting {n_chains} chains using {n_processes} processes...")
    results = {}
    
    with mp.Pool(processes=n_processes) as pool:
        with tqdm(total=n_chains, desc="Running chains") as pbar:
            for result in pool.imap_unordered(run_chain, range(n_chains)):
                chain_idx = result['chain_idx']
                results[chain_idx] = result
                pbar.update(1)
                pbar.set_postfix({'Chain': chain_idx})
    
    return results

def analyze_chains(results: Dict[int, Dict], 
                  output_folder: str, 
                  directory_name: str) -> pd.DataFrame:
    """Organize result files and generate summary statistics."""
    
    # Setup output directory
    sampler_folder = os.path.join(output_folder, directory_name)
    old_sampler_folder = os.path.join(output_folder, f"{directory_name}_old")
    
    # Archive previous results if they exist
    if os.path.exists(sampler_folder):
        if os.path.exists(old_sampler_folder):
            shutil.rmtree(old_sampler_folder)
        os.rename(sampler_folder, old_sampler_folder)
    os.makedirs(sampler_folder, exist_ok=True)
    
    # Process chains with successful runs
    valid_results = {idx: res for idx, res in results.items() if 'error' not in res}
    if not valid_results:
        print("No valid chains completed. Check logs for errors.")
        return pd.DataFrame()
    
    # Copy trajectory files to centralized location
    for chain_idx, result in valid_results.items():
        chain_name = f"chain_{chain_idx+1}"
        new_traj_file = os.path.join(sampler_folder, f"trajectory_{chain_name}.h5")
        
        src_file = result.get('trajectory_file')
        if src_file and os.path.exists(src_file):
            shutil.copy2(src_file, new_traj_file)
            result['centralized_trajectory'] = new_traj_file
    
    # Collect sigma values across chains
    sigma_values = pd.DataFrame([
        {**{'chain': f"chain_{idx+1}"}, **result['final_sigma']}
        for idx, result in valid_results.items()
    ])
    
    # Save and summarize results
    sigma_values.to_csv(os.path.join(sampler_folder, "final_sigma_values.csv"), index=False)
    
    print("\nSummary of sigma values across chains:")
    summary = sigma_values.describe()
    print(summary)
    summary.to_csv(os.path.join(sampler_folder, "sigma_statistics.csv"))
    
    return sigma_values

def main():
    """Main execution function with configuration."""
    # Define sampler configurations
    base_config = DEFAULT_CONFIG.copy()
    
    # Define the sequence of samplers to run and their steps
    #sampler_sequence = ["pair", "tetramer", "octet", "pair", "octet", "tetramer"]
    #mcmc_steps = [5000000, 200000, 100000, 3000000, 150000, 100000]
    
    # Alternative sequences for testing
    sampler_sequence = ["pair", "tetramer", "pair", "octet"]
    mcmc_steps = [100000, 50000, 50000, 50000]
    #sampler_sequence = ["pair", "tetramer", "octet", "pair", "tetramer", "octet", "pair"]
    #mcmc_steps = [500000, 500000, 300000, 300000, 200000, 200000, 200000]
    # sampler_sequence = ["tetramer", "octet", "pair"]
    # mcmc_steps = [500000, 500000, 1000000]
    
    # Setup
    output_folder = "output_analysis"
    n_processes = 8  # Set to None to use all available CPUs
    
    # Initialize sequence manager
    sequence_manager = SamplerSequenceManager(sampler_sequence)
    
    print(f"\nSampler Sequence Plan:")
    print(f"Total sequence: {sampler_sequence}")
    print(f"Sampler counts: {sequence_manager.sampler_counts}")
    print(f"{'='*80}")
    
    for seq_idx, sampler_key in enumerate(sampler_sequence):
        # Get sampler information
        sampler_name, current_idx, total_count = sequence_manager.get_sampler_info(seq_idx)
        directory_name = sequence_manager.get_directory_name(sampler_name, current_idx)
        
        # Determine if this is the first sampler in the entire sequence
        is_first = seq_idx == 0
        
        # Create specific configuration for this sampler
        config = {
            **base_config,
            "run": True,
            "n_chains": 8,
            "n_steps": mcmc_steps[seq_idx],
            "save_freq": 200,
            "use_sigma_dist": False if is_first else True
        }
        
        # Get sampler class from name
        sampler_class = SAMPLER_MAP.get(sampler_key)
        
        if not sampler_class:
            print(f"Error: Unknown sampler '{sampler_key}'. Skipping.")
            continue
        
        print(f"\n{'='*80}")
        print(f"STEP {seq_idx+1}/{len(sampler_sequence)}: Starting {sampler_class.__name__}")
        print(f"Sampler: {sampler_name} (instance {current_idx}/{total_count})")
        print(f"Directory: {directory_name}")
        print(f"Configuration: {config['n_steps']} steps, {config['n_chains']} chains, use_sigma_dist={config['use_sigma_dist']}")
        print(f"{'='*80}")
        
        # Run sampling and analyze results
        results = run_parallel_sampling(
            sampler_class=sampler_class,
            config=config,
            output_folder=output_folder,
            n_processes=n_processes,
            sampler_sequence=sampler_sequence,
            sequence_position=seq_idx
        )
        
        sigma_values = analyze_chains(results, output_folder, directory_name)
        
        print(f"\n{sampler_class.__name__} sampling complete.")
        print(f"Results saved in: {output_folder}/{directory_name}/")
        
        # Run the fit_gmm analysis for this sampler - UPDATED CALL
        analysis_success = run_analysis(
            sampler_key=sampler_key, 
            sampler_index=current_idx,
            sampler_sequence=sampler_sequence,
            sequence_position=seq_idx
        )
        
        if analysis_success:
            print(f"Analysis for {sampler_key} (occurrence #{current_idx}) at position {seq_idx} completed successfully.")
        else:
            print(f"WARNING: Analysis for {sampler_key} (occurrence #{current_idx}) at position {seq_idx} may have failed.")
        
        # Add a small delay before next sampler to ensure files are properly written
        time.sleep(2)
    
    print(f"\n{'='*80}")
    print("ALL SAMPLERS COMPLETED!")
    print(f"Final directory structure in {output_folder}:")
    
    # List final directories
    if os.path.exists(output_folder):
        for item in sorted(os.listdir(output_folder)):
            if os.path.isdir(os.path.join(output_folder, item)) and "sampler_results" in item:
                print(f"  - {item}/")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()
