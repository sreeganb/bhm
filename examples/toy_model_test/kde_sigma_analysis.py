#============================================================
# Analysis script across the different sampler sequences and chains within 
# each sequence. Plot the distribution of sigma values for 
# each pair type and how they evolved
# across the different samplers.
#============================================================
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import os
import pandas as pd
import h5py
#--------------------------------------------------------------
# Function to read the sampler sequence from a JSON file
#--------------------------------------------------------------
def read_sampler_sequence(file_path):
    import json
    with open(file_path, 'r') as f:
        data = json.load(f)
    print(f"Sampler sequence: {data['sequence']}")
    print(f"Number of chains: {data['n_chains']}")
    # Return the sequence and number of chains for further processing
    return data['sequence'], data['n_chains']
#--------------------------------------------------------------
# Extract sigma time series from HDF5 file
#--------------------------------------------------------------
def load_sigma_timeseries(filename: str) -> pd.DataFrame:
    """
    Loads only the sigma parameters (AA, AB, BC) from an HDF5 trajectory file.
    Returns a pandas DataFrame with columns for step and sigma values.
    """
    # Initialize lists to collect data
    steps = []
    sigma_AA = []
    sigma_AB = []
    sigma_BC = []
    
    if not os.path.isfile(filename):
        print(f"File not found: {filename}")
        return pd.DataFrame()

    with h5py.File(filename, 'r') as f:
        # Check if trajectory group exists
        if 'trajectory' not in f:
            print(f"No 'trajectory' group in file: {filename}")
            return pd.DataFrame()
        
        traj_grp = f['trajectory']
        for state_name in traj_grp:
            state_grp = traj_grp[state_name]
            step = state_grp.attrs.get("step", 0)
            
            # Extract only sigma values
            if 'sigma' in state_grp:
                sigma_grp = state_grp['sigma']
                
                # Add step number and sigma values to our lists
                steps.append(step)
                sigma_AA.append(sigma_grp.attrs.get('AA', float('nan')))
                sigma_AB.append(sigma_grp.attrs.get('AB', float('nan')))
                sigma_BC.append(sigma_grp.attrs.get('BC', float('nan')))
    
    # Create DataFrame and sort by step number
    df = pd.DataFrame({
        'step': steps,
        'sigma_AA': sigma_AA,
        'sigma_AB': sigma_AB,
        'sigma_BC': sigma_BC
    })
    
    # Sort by step number
    df = df.sort_values('step').reset_index(drop=True)
    print(f"Loaded sigma time series with {len(df)} steps from {filename}")
    
    return df
#--------------------------------------------------------------
# KDE fit got a DataFrame of sigma values   
#--------------------------------------------------------------
def fit_kde_to_sigma(df: pd.DataFrame,  bandwidth: float = 0.5):
    """
    Fit a Kernel Density Estimate (KDE) for all the sigma pair types in the DataFrame.
    Returns a dictionary of scipy.stats.gaussian_kde objects for each pair type.
    """
    kdes = {}
    for pair_type in ['AA', 'AB', 'BC']:
        sigma_col = f'sigma_{pair_type}'
        if sigma_col not in df.columns:
            print(f"Column {sigma_col} not found in DataFrame.")
            continue
        kde = stats.gaussian_kde(df[sigma_col].dropna(), bw_method=bandwidth)
        kdes[pair_type] = kde
    return kdes
#--------------------------------------------------------------
# For each sampler sequence
#---------------------------------------------------------------
def read_data_from_folders(sequence, n_chains):
    sampler_ids = []
    # check all occurences of the sampler: "pair" and mark down the indices
    ind_pair = 0
    ind_tet = 0
    ind_oct = 0
    for sampler in sequence:
        kind = 0
        if sampler == "pair":
            ind_pair += 1
            kind = ind_pair
        elif sampler == "tetramer":
            ind_tet += 1
            kind = ind_tet
        elif sampler == "octet":
            ind_oct += 1
            kind = ind_oct
        else:
            raise ValueError(f"Unknown sampler type: {sampler}")
        sampler_ids.append(kind)
    print(f"Sampler IDs: {sampler_ids}")
    
    # Dictionary to store KDEs for each sampler and sigma type
    all_kdes = {
        'AA': {},
        'AB': {},
        'BC': {}
    }
    
    # Collect all KDEs first
    for seq, id in zip(sequence, sampler_ids):
        folder_name = f"./output_analysis/{seq}sampler_results_{id}/"
        print(f"Reading data from folder: {folder_name}")
        # Check if the folder exists
        if not os.path.exists(folder_name):
            print(f"Folder {folder_name} does not exist. Skipping.")
            continue
        
        parent_dir = os.getcwd()
        os.chdir(folder_name)
        
        # Pool in data from all chains
        df_combined = pd.DataFrame()
        for chain in range(n_chains):
            filename = f"trajectory_chain_{chain+1}.h5"
            df_temp = load_sigma_timeseries(filename)
            if df_temp.empty:
                print(f"No data found in {filename}. Skipping.")
                continue
            # Keep pooling in all the data across chains
            df_combined = pd.concat([df_combined, df_temp], ignore_index=True)

        # After reading all chains go back to the main directory
        os.chdir(parent_dir)
        print(f"size of the combined dataframe: {df_combined.shape}")
        
        # Fit KDEs for this sampler
        kdes = fit_kde_to_sigma(df_combined)
        print(f"KDEs fitted for {seq} sampler with ID {id}.")
        
        # Store the KDEs with a label
        sampler_label = f"{seq}_{id}"
        for pair_type, kde in kdes.items():
            all_kdes[pair_type][sampler_label] = kde
    
    # Create a single figure with three subplots (one for each sigma type)
    plt.figure(figsize=(15, 12))
    
    # Color map for visually distinguishing samplers
    colors = plt.cm.tab10(np.linspace(0, 1, len(sequence)))
    
    # Plot each sigma type in a separate subplot
    for i, pair_type in enumerate(['AA', 'AB', 'BC']):
        plt.subplot(3, 1, i+1)
        
        # Plot all samplers for this sigma type
        x = np.linspace(0, 20, 1000)
        for j, (sampler_label, kde) in enumerate(all_kdes[pair_type].items()):
            color = colors[j % len(colors)]
            plt.plot(x, kde(x), label=sampler_label, color=color)
            plt.fill_between(x, kde(x), alpha=0.1, color=color)
        
        plt.title(f'Sigma {pair_type} KDE for All Samplers')
        plt.xlabel('Sigma Value')
        plt.ylabel('Density')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig('combined_kde_plots.png')
    plt.show()
    plt.close()
    
    return sampler_ids


if __name__ == "__main__":
    sequence, n_chains = read_sampler_sequence('sampler_sequence.json')
    sampler_ids = read_data_from_folders(sequence, n_chains)
