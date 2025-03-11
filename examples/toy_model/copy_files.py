import os
import shutil

def organize_results(base_output_folder: str, sampler_folder_name: str, timestamp_folder: str):
    """
    Renames the old results folder, creates a new one, and copies trajectory files
    with a consistent naming scheme.

    Args:
        base_output_folder: The main output directory (e.g., "output_analysis").
        sampler_folder_name: The base name of the sampler folder (e.g., "pairsampler").
        timestamp_folder: The timestamped folder containing individual chain results.
    """

    sampler_results_folder = os.path.join(base_output_folder, f"{sampler_folder_name}_results")
    old_sampler_results_folder = os.path.join(base_output_folder, f"{sampler_folder_name}_results_old")
    timestamped_chains_folder = os.path.join(base_output_folder, timestamp_folder)

    # 1. Rename existing results folder (if it exists)
    if os.path.exists(sampler_results_folder):
        if os.path.exists(old_sampler_results_folder):
            shutil.rmtree(old_sampler_results_folder)  # Remove any existing _old folder
        os.rename(sampler_results_folder, old_sampler_results_folder)
        print(f"Renamed '{sampler_results_folder}' to '{old_sampler_results_folder}'")

    # 2. Create a new results folder
    os.makedirs(sampler_results_folder, exist_ok=True)
    print(f"Created new directory: '{sampler_results_folder}'")

    # 3. Iterate through chain folders and copy/rename trajectory files
    if os.path.exists(timestamped_chains_folder):
        for chain_folder_name in os.listdir(timestamped_chains_folder):
            chain_folder_path = os.path.join(timestamped_chains_folder, chain_folder_name)
            # check that it is actually a folder
            if os.path.isdir(chain_folder_path):

                # Extract chain number.  More robust to variations in naming.
                try:
                    chain_number = int(chain_folder_name.split("_")[-1])  # Extract number
                except ValueError:
                    print(f"Skipping '{chain_folder_name}' - could not extract chain number.")
                    continue  # Skip if we can't get a valid number

                trajectory_file_path = os.path.join(chain_folder_path, "trajectory.h5")

                if os.path.exists(trajectory_file_path):
                    new_trajectory_file_name = f"trajectory_{sampler_folder_name.capitalize()}Sampler_chain_{chain_number}.h5"
                    new_trajectory_file_path = os.path.join(sampler_results_folder, new_trajectory_file_name)

                    shutil.copy2(trajectory_file_path, new_trajectory_file_path)  # copy2 preserves metadata
                    print(f"Copied '{trajectory_file_path}' to '{new_trajectory_file_path}'")
                else:
                    print(f"Warning: 'trajectory.h5' not found in '{chain_folder_path}'")
    else:
        print("Warning: timestamped chains folder:", timestamped_chains_folder, " not found")

if __name__ == '__main__':
    base_output_folder = "output_analysis"  #  Replace with your base output folder
    #sampler_folder_name = "pairsampler"    #  Replace with your sampler name (without "Sampler")
    sampler_folder_name = "tetramersampler"    #  Replace with your sampler name (without "Sampler")
    timestamp_folder = "tetramersampler_20250310_112337"   # Replace with the actual timestamp folder

    organize_results(base_output_folder, sampler_folder_name, timestamp_folder)
