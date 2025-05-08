import os
import pickle
import json
import numpy as np
from typing import Dict, Tuple, Optional, Sequence, Any, Union
from .state import SystemState

class SigmaBuilder:
    """
    Class-based initializer for sigma values in SystemState, supporting
    default, explicit values, or GMM-derived sigmas.  When using GMM,
    it also stores the loaded GMM parameters in state.metadata['gmm_params'].
    """
    def __init__(
        self,
        state: SystemState,
        sigma_source: str = "default",
        sigma_values: Optional[Dict[str, float]] = None,
        sigma_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
        gmm_file: Optional[str] = None,
        sigma_keys: Optional[Sequence[str]] = None,
        gmm_folder: str = "output_analysis/gmm_parameters",
        sampler_name: Optional[str] = None
    ):
        self.state = state
        self.source = sigma_source
        self.values = sigma_values or {}
        self.ranges = sigma_ranges
        self.gmm_file = gmm_file
        self.gmm_folder = gmm_folder
        self.sampler_name = sampler_name
        # determine which keys to initialize
        self.keys = list(sigma_keys) if sigma_keys is not None else list(self.state.sigma_range.keys())

    def build(self) -> None:
        """
        Initialize state.sigma, state.sigma_range, state.use_sigma_distribution,
        and store any GMM parameters in state.metadata['gmm_params'].
        """
        # set ranges
        if self.ranges is None:
            # default ranges
            self.ranges = {k: (1.0, 10.0) for k in self.keys}
        self.state.sigma_range = {k: self.ranges[k] for k in self.keys}

        if self.source == "default":
            # midpoint of range
            mid = {k: 0.5 * (self.ranges[k][0] + self.ranges[k][1]) for k in self.keys}
            self.state.update_sigma(mid)
            self.state.use_sigma_distribution = False

        elif self.source == "values":
            if not self.values:
                raise ValueError("sigma_values must be provided when sigma_source='values'")
            self.state.update_sigma({k: self.values[k] for k in self.keys})
            self.state.use_sigma_distribution = False

        elif self.source == "gmm":
            # mark distribution usage
            self.state.use_sigma_distribution = True
            # locate gmm file if not explicit
            gmm_path = self.gmm_file or self._auto_select_gmm()
            # load parameters
            gmm_params = self._load_gmm_params(gmm_path)
            # store full gmm_params in metadata
            self.state.metadata['gmm_params'] = gmm_params
            # extract sigma
            sig_dict = {}
            for k in self.keys:
                if k in gmm_params:
                    params = gmm_params[k]
                    sig_dict[k] = self._extract_sigma_from_params(params)
                else:
                    # fallback mid-range
                    sig_dict[k] = 0.5 * (self.ranges[k][0] + self.ranges[k][1])
            self.state.update_sigma(sig_dict)

        else:
            raise ValueError(f"Invalid sigma_source: {self.source}")

    def _load_gmm_params(self, path: str) -> Dict[str, Any]:
        """
        Load GMM parameters from a file, supporting both pickle and JSON formats.
        """
        if path.endswith('.json'):
            # Extract parameter key from filename (e.g., 'AA' from gmm_fit_AA_chain_1.json)
            param_key = os.path.basename(path).split('_')[2]
            with open(path, 'r') as f:
                json_data = json.load(f)
                
            # Convert JSON format to a structure compatible with our code
            return {param_key: {
                'means_': np.array([[m] for m in json_data['means']]),
                'covariances_': np.array([[c] for c in json_data['covariances']]),
                'weights_': np.array(json_data['weights']),
                'n_components': json_data['n_components']
            }}
        else:
            # Original pickle format
            with open(path, 'rb') as f:
                return pickle.load(f)

    def _extract_sigma_from_params(self, params: Dict[str, Any]) -> float:
        """
        Extract sigma value from GMM parameters, handling both pickle-loaded
        sklearn models and our custom JSON format.
        """
        if hasattr(params, 'weights_'):  # sklearn GMM model from pickle
            comp = np.argmax(params.weights_)
            return float(params.means_[comp][0])
        else:  # Our JSON-derived format
            comp = np.argmax(params['weights_'])
            return float(params['means_'][comp][0])

    def _auto_select_gmm(self) -> str:
        """
        Pick the latest GMM file for sampler_name in gmm_folder.
        """
        if not self.sampler_name:
            raise ValueError("sampler_name must be provided to auto-select GMM file")
        
        # Find all files/directories that match the pattern
        gmm_files = []
        for root, _, files in os.walk(self.gmm_folder):
            for file in files:
                if file.startswith(f"gmm_fit_{self.sampler_name}_chain_") and (file.endswith('.json') or file.endswith('.pkl')):
                    gmm_files.append(os.path.join(root, file))
        
        if not gmm_files:
            # Fallback to previous directory-based logic
            chains = [d for d in os.listdir(self.gmm_folder)
                    if d.startswith(self.sampler_name)]
            if not chains:
                raise FileNotFoundError(f"No GMM chains for {self.sampler_name} in {self.gmm_folder}")
            # pick highest numerical suffix
            best = sorted(chains, key=lambda d: int(d.split('_')[-1]))[-1]
            path = os.path.join(self.gmm_folder, best)
        else:
            # Sort by chain number to get the latest one
            path = sorted(
                gmm_files,
                key=lambda f: int(os.path.basename(f).split('_chain_')[1].split('.')[0])
            )[-1]
            
        if not os.path.exists(path):
            raise FileNotFoundError(f"GMM file not found: {path}")
        return path