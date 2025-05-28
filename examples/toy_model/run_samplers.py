import os

from typing import Dict, List, Type
from collections import defaultdict

from parameters import SystemParameters
#from base_sampler import BaseMCSampler
from pair_sampler import PairSampler
from tetramer_sampler import TetramerSampler
from octet_sampler import OctetSampler
from visualization import visualize_3d_configuration, visualize_trajectory_plotly
import cProfile
import pstats
import io

#from graph_results import SigmaAnalyzer

SAMPLER_MAP = {
    "PairSampler": PairSampler,
    "TetramerSampler": TetramerSampler,
    "OctetSampler": OctetSampler, # Uncomment when OctamerSampler is implemented
}

def run_hierarchical_sampling(
    sampler_class: Type[PairSampler],
    n_chains: int = 1,
    n_steps: int = 300000,
    save_freq: int = 300,
    use_sigma_distribution: bool = False,
    output_folder: str = "output_analysis"
):
    """
    Runs hierarchical sampling for a given sampler class.
    Args:
        sampler_class: The class of the sampler to use (e.g., PairSampler, TetramerSampler).
        n_chains: Number of independent chains to run.
        n_steps: Number of Monte Carlo steps per chain.
        save_freq: Frequency at which to save trajectory data.
    Returns:
        A dict with results for each chain.
    """
    if sampler_class == PairSampler:
        output_folder = os.path.join(output_folder, "pairsampler_results")
    elif sampler_class == TetramerSampler:
        output_folder = os.path.join(output_folder, "tetramersampler_results")
    elif sampler_class == OctetSampler:
        output_folder = os.path.join(output_folder, "octetsampler_results")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    params = SystemParameters() # Consider making SystemParameters configurable if needed.
    results = {}
    levels = []
    sampler_name = sampler_class.__name__

    for chain_idx in range(n_chains):
        level = f"{sampler_name}_chain_{chain_idx+1}"
        levels.append(level)
        print(f"\nRunning {level} sampling ({sampler_name})...")

        sampler = sampler_class(use_sigma_distribution) # Initialize the sampler with default parameters
        
        #best_positions, trajectory, traj_file = sampler.run_mc(
        #    n_steps=n_steps,
        #    save_freq=save_freq,
        #)
        best_positions, traj_file = sampler.run_mc(
            n_steps=n_steps,
            save_freq=save_freq,
        )

        # --- File Handling ---
        base_filename = f"{sampler_name.lower()}_chain_{chain_idx+1}" # Consistent filename base

        # Trajectory file
        new_traj_file = os.path.join(output_folder, f"trajectory_{level}.h5")
        os.rename(traj_file, new_traj_file)
        traj_file = new_traj_file

        # Sigma history file - Assuming sampler.sigma_history and .sigma_history_tetramer are consistent or adaptable
#        sigma_history_attr = "sigma_history" if sampler_name == "PairSampler" else "sigma_history_tetramer" # Adapt attribute name
#        sigma_filename_base = "sigma_history" if sampler_name == "PairSampler" else "sigma_history_tetramer" # Adapt filename base
#        sigma_fil = os.path.join(output_folder, f"{sigma_filename_base}.csv") # Read intermediate file
#        new_sigmafil = os.path.join(output_folder, f"sigma_history_{level}.csv") # Rename
#        os.rename(sigma_fil, new_sigmafil)

        print(f"{level} sampling complete ({sampler_name}):")
        print(f"Final sigma values: {sampler.sigma}")

        #visualize_3d_configuration(best_positions, params.radii, f"Final {level} configuration")
        #visualize_trajectory_plotly(trajectory, params.box_size, params.radii)

        results[level] = {
            'best_positions': best_positions,
            'trajectory_file': traj_file,
            'final_sigma': sampler.sigma,
            #'sigma_history_path': csv_path # No longer needed as path is predictable
        }
    return results

def main():
    run_config = {
        "pair_sampler": {
            #"run": False,
            "run": True,
            "n_chains": 1,
            "n_steps": 10000,
            "save_freq": 5,
            "use_sigma_dist": False, # Example of sampler-specific parameter
        },
        "tetramer_sampler": {
            "run": True,
            #"run": False,
            "n_chains": 1,
            "n_steps": 100000,
            "save_freq": 50,
            "use_sigma_dist": True, # Example of sampler-specific parameter
        },
         "octet_sampler": { # Example configuration for OctamerSampler
             "run": True,
             #"run": False,
             "n_chains": 1,
             "n_steps": 50000,
             "save_freq": 500,
             "use_sigma_dist": True #, # Example of sampler-specific parameter
             #"read_pos": False,
         }
    }

    output_folder = "output_analysis"
    sampler_results = {}
    last_sigma = None # To hold sigma from previous level if hierarchical

    #sampler_sequence = ["pair_sampler", "tetramer_sampler"] # Define the sequence of samplers to run
    #sampler_sequence = ["pair_sampler"] # Example: Run only pair sampler
    sampler_sequence = ["tetramer_sampler"] # Example: Run only tetramer sampler
    #sampler_sequence = ["octet_sampler"] # Example: Run only tetramer sampler
    #sampler_sequence = ["pair_sampler", "tetramer_sampler", "octamer_sampler"] # Example: Run all, assuming octamer_sampler config is present

    # Profile the output 
    profil = True
    if profil == True:
        profile_output = "profile_output.prof"
    for sampler_key in sampler_sequence:
        config = run_config[sampler_key]
        if config["run"]:
            sampler_name = sampler_key.replace("_sampler", "").capitalize() # e.g., "PairSampler"
            sampler_name = sampler_name + "Sampler"
            sampler_class = SAMPLER_MAP[sampler_name]

            print(f"\nStarting sampling with {sampler_name}...")
            if profil == True:
                # Use runctx instead of run to properly access local variables
                statement = 'results = run_hierarchical_sampling(sampler_class=sampler_class, n_chains=config["n_chains"], n_steps=config["n_steps"], save_freq=config["save_freq"], use_sigma_distribution=config["use_sigma_dist"], output_folder=output_folder)'
                cProfile.runctx(statement, globals(), locals(), profile_output)
            else:
                results = run_hierarchical_sampling(
                    sampler_class=sampler_class,
                    n_chains=config["n_chains"],
                    n_steps=config["n_steps"],
                    save_freq=config["save_freq"],
                    use_sigma_distribution=config["use_sigma_dist"], 
                    output_folder=output_folder
                )
                sampler_results[sampler_key] = results

                print(f"{sampler_name} hierarchical sampling complete.")
        else:
            print(f"{sampler_key} run is skipped as per configuration.")

    if profil == True:
        # Print sorted stats to a readable file
        with open("profiled_output.txt", "w") as f:
            p = pstats.Stats(profile_output, stream=f)
            p.sort_stats('cumulative').print_stats(50)  # Top 50 functions by cumulative time
            f.write("\n\n")
            p.sort_stats('time').print_stats(50)  # Top 50 functions by internal time
    
    print("\nNext steps:")
    print(f"• Check '{output_folder}/sigma_history_*' CSV files for sigma evolution.")
    print("• Now call the SigmaAnalyzer for post-processing...")

    # Analyze results - using TetramerSampler analyzer for all for now, adjust if needed.
#    analyzer = SigmaAnalyzer(output_folder=output_folder, sampler="TetramerSampler") # Or adapt analyzer if needed per sampler type
#    analyzer.generate_report()
#    analyzer.plot_convergence()
#    analyzer.plot_trace_combined()
#    analyzer.constraint_violations()

if __name__ == "__main__":
    main()
