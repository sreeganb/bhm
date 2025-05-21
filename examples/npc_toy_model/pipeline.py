# pipeline.py
from typing import List, Dict, Any, Tuple
from core.state import SystemState
import os
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from sklearn.mixture import GaussianMixture
from collections import defaultdict
import h5py

class SamplerPipeline:
    """Run a sequence of samplers in a pipeline"""
    
    def __init__(self, initial_state: SystemState):
        self.initial_state = initial_state
        self.stages = []
    
    def add_stage(
        self, 
        sampler_function, 
        n_steps: int = 1000, 
        save_freq: int = 100,
        name: str = None,
        **sampler_kwargs
    ):
        """
        Add a sampler stage to the pipeline
        
        Args:
            sampler_function: Function to run the sampler
            n_steps: Number of MCMC steps
            save_freq: How often to save state to disk
            name: Name of this stage (defaults to function name)
            **sampler_kwargs: Additional kwargs for the sampler
        """
        if name is None:
            name = sampler_function.__name__.replace('run_', '')
        
        self.stages.append({
            'function': sampler_function,
            'n_steps': n_steps,
            'save_freq': save_freq,
            'name': name,
            'kwargs': sampler_kwargs
        })
        
        return self  # Enable method chaining
    
    def fit_gmm(self, data: np.ndarray, max_components: int = 5):
        """Fits a Gaussian Mixture Model to the data"""
        data = data.reshape(-1, 1)
        best_gmm = None
        best_bic = np.inf

        for n_components in range(1, max_components + 1):
            gmm = GaussianMixture(n_components=n_components, covariance_type='full', random_state=42)
            gmm.fit(data)
            bic = gmm.bic(data)
            if bic < best_bic:
                best_bic = bic
                best_gmm = gmm

        return best_gmm
    
    def analyze_trajectory_and_fit_gmm(self, trajectory_file: str, output_dir: str, stage_name: str, chain_id: int, burnin: float = 0.4):
        """
        Analyze a trajectory file and fit GMMs to sigma parameters
        
        Args:
            trajectory_file: Path to HDF5 trajectory file
            output_dir: Directory to save GMM files and plots
            stage_name: Name of the sampling stage
            chain_id: Chain number
            burnin: Fraction of initial samples to discard
            
        Returns:
            Dict mapping sigma types to GMM models
        """
        print(f"Analyzing trajectory and fitting GMMs for {stage_name} (chain {chain_id})...")
        
        # Load trajectory data
        all_states = self.load_trajectory_from_hdf5(trajectory_file)
        if not all_states:
            print(f"Warning: No valid states found in {trajectory_file}")
            return {}
        
        # Apply burnin
        num_samples = len(all_states)
        burnin_samples = int(num_samples * burnin)
        all_states = all_states[burnin_samples:]
        
        if not all_states:
            print(f"Warning: No states remain after burnin in {trajectory_file}")
            return {}
        
        # Extract sigma values
        sigma_data = defaultdict(list)
        for state in all_states:
            for sigma_type, sigma_value in state["sigma"].items():
                sigma_data[sigma_type].append(sigma_value)
        
        # Fit GMMs
        gmm_models = {}
        for sigma_type, values in sigma_data.items():
            if len(values) > 10:  # Need sufficient data
                data = np.array(values)
                gmm = self.fit_gmm(data)
                gmm_models[sigma_type] = gmm
                
                # Save GMM parameters to JSON
                gmm_params = {
                    "n_components": gmm.n_components,
                    "means": gmm.means_.flatten().tolist(),
                    "covariances": gmm.covariances_.flatten().tolist(),
                    "weights": gmm.weights_.tolist()
                }
                
                json_filename = os.path.join(output_dir, f"gmm_fit_{sigma_type}_chain{chain_id}.json")
                with open(json_filename, "w") as f:
                    json.dump(gmm_params, f, indent=4)
                
                # Create plots for this sigma parameter
                self.plot_gmm_fit(data, gmm, sigma_type, stage_name, chain_id, output_dir)
        
        return gmm_models
    
    def plot_gmm_fit(self, data, gmm, sigma_type, stage_name, chain_id, output_dir):
        """Create plots of the GMM fit for a sigma parameter"""
        plt.figure(figsize=(10, 6))
        
        # Histogram of data
        plt.hist(data, bins=30, density=True, alpha=0.6, label='Data')
        
        # GMM fit line
        x = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)
        log_prob = gmm.score_samples(x)
        plt.plot(x, np.exp(log_prob), 'r-', label=f'GMM ({gmm.n_components} components)')
        
        # Individual Gaussian components
        for i in range(gmm.n_components):
            pdf = np.exp(self._gaussian_pdf(x, gmm.means_[i, 0], gmm.covariances_[i, 0, 0]))
            plt.plot(x, gmm.weights_[i] * pdf, '--', label=f'Component {i+1}')
        
        plt.title(f'GMM Fit for {sigma_type} ({stage_name}, Chain {chain_id})')
        plt.xlabel('Sigma Value')
        plt.ylabel('Density')
        plt.legend()
        
        plt.savefig(os.path.join(output_dir, f"gmm_fit_{sigma_type}_chain{chain_id}.png"))
        plt.close()
    
    def _gaussian_pdf(self, x, mean, variance):
        """Calculate log probability for a 1D Gaussian"""
        return -0.5 * (np.log(2 * np.pi * variance) + (x - mean)**2 / variance)
    
    def load_trajectory_from_hdf5(self, filename: str) -> list:
        """Load trajectory from HDF5 file"""
        all_states = []
        if not os.path.isfile(filename):
            print(f"File not found: {filename}")
            return all_states

        try:
            with h5py.File(filename, 'r') as f:
                # Check if trajectory group exists
                if 'trajectory' not in f:
                    print(f"No 'trajectory' group in file: {filename}")
                    return all_states
                
                traj_grp = f['trajectory']
                for state_name in traj_grp:
                    state_grp = traj_grp[state_name]
                    state = {
                        "step": state_grp.attrs.get("step", 0),
                        "total_score": state_grp.attrs.get("total_score", 0.0),
                        "prior_score": state_grp.attrs.get("prior_score", 0.0),
                        "pair_score": state_grp.attrs.get("pair_score", 0.0),
                        "exvol_score": state_grp.attrs.get("exvol_score", 0.0),
                        "sigma": {},
                        "positions": {},
                    }
                    
                    # Read sigma values
                    if 'sigma' in state_grp:
                        sigma_grp = state_grp['sigma']
                        for key in sigma_grp.attrs:
                            state["sigma"][key] = float(sigma_grp.attrs[key])

                    # Read positions
                    if 'positions' in state_grp:
                        pos_grp = state_grp['positions']
                        for type_name in pos_grp:
                            state["positions"][type_name] = pos_grp[type_name][()]
                    
                    all_states.append(state)

            # Sort by step number
            all_states.sort(key=lambda s: s["step"])
            return all_states
            
        except Exception as e:
            print(f"Error loading trajectory file {filename}: {e}")
            return all_states
    
    def run(
        self, 
        output_base: str = "output/pipeline", 
        n_chains: int = 1,
        fit_gmm: bool = True,
        burnin: float = 0.4
    ) -> Dict[str, Any]:
        """Run the pipeline for multiple chains"""
        os.makedirs(output_base, exist_ok=True)
        
        # Save pipeline configuration
        config = {
            'stages': [
                {
                    'name': stage['name'],
                    'n_steps': stage['n_steps'],
                    'save_freq': stage.get('save_freq', 100),
                    'kwargs': {k: str(v) for k, v in stage['kwargs'].items() 
                               if not callable(v)}
                }
                for stage in self.stages
            ],
            'n_chains': n_chains,
            'fit_gmm': fit_gmm,
            'burnin': burnin
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
                os.makedirs(stage_output, exist_ok=True)
                
                print(f"  Stage {i+1}/{len(self.stages)}: {stage_name} ({stage['n_steps']} steps)")
                start_time = time.time()
                
                # Create trajectory filename with chain number
                traj_file = os.path.join(stage_output, f"trajectory_chain{chain}.h5")
                
                # Run the sampler, passing chain and trajectory information
                sampler_fn = stage['function']
                final_state, trajectory_file = sampler_fn(
                    state=current_state,
                    n_steps=stage['n_steps'],
                    output_dir=stage_output,
                    chain_number=chain,
                    save_freq=stage.get('save_freq', 100),
                    trajectory_file=traj_file,
                    **stage['kwargs']
                )
                
                elapsed = time.time() - start_time
                print(f"  Completed in {elapsed:.1f} seconds")
                
                # Fit GMMs if enabled and not the last stage
                gmm_models = {}
                if fit_gmm and i < len(self.stages) - 1:
                    gmm_models = self.analyze_trajectory_and_fit_gmm(
                        trajectory_file, 
                        stage_output,
                        stage_name,
                        chain,
                        burnin
                    )
                
                # Use the final state as input to the next stage
                current_state = final_state
                
                # Record results
                chain_results.append({
                    'stage': stage_name,
                    'trajectory': trajectory_file,
                    'final_sigma': {k: float(v) for k, v in current_state.sigma.items()},
                    'gmm_fitted': bool(gmm_models)
                })
            
            results[f"chain_{chain}"] = chain_results
            
            print(f"Chain {chain} complete. Final sigma values:")
            for k, v in current_state.sigma.items():
                print(f"  {k}: {v:.4f}")
        
        # Save overall results summary
        with open(os.path.join(output_base, "results_summary.json"), 'w') as f:
            json.dump(results, f, indent=2)
        
        return results