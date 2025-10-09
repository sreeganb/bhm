# io_utils.py - Utilities for saving state information
import numpy as np
import h5py
from typing import Dict, Any, Optional

def save_state_to_disk(
    step: int, 
    positions: Dict[str, np.ndarray], 
    sigmas: Dict[str, float], 
    score: float,
    prior_score: float = 0.0, 
    pair_score: float = 0.0, 
    exvol_score: float = 0.0, 
    tet_score: float = 0.0, 
    oct_score: float = 0.0,
    types: Optional[Dict] = None, 
    bead_numbers: Optional[Dict] = None, 
    traj_file: Optional[str] = None
) -> None:
    """
    Save state directly to an HDF5 file in a memory-efficient manner.
    
    Parameters:
        step: The current step number.
        positions: Dictionary mapping component names to position arrays.
        sigmas: Dictionary of sigma values.
        score: The total score (will be stored as 'total_score').
        prior_score, pair_score, exvol_score, tet_score, oct_score: Additional scores.
        types: Dictionary mapping bead indices to type names.
        bead_numbers: Dictionary mapping bead indices to bead numbers.
        traj_file: Path to the HDF5 file to write to.
    """
    # If types or bead_numbers are not provided, default to empty dictionaries.
    if types is None:
        types = {}
    if bead_numbers is None:
        bead_numbers = {}
    if traj_file is None:
        return

    try:
        with h5py.File(traj_file, 'a') as f:
            # Create or get the trajectory group.
            if 'trajectory' not in f:
                traj_grp = f.create_group('trajectory')
            else:
                traj_grp = f['trajectory']
            
            # Create a new state group named "state_XXXXX" where XXXXX is the step number zero-padded.
            state_name = f"state_{step:05d}"
            state_grp = traj_grp.create_group(state_name)
            
            # Save state attributes (same keys as before).
            state_grp.attrs["step"] = step
            state_grp.attrs["total_score"] = float(score)
            state_grp.attrs["prior_score"] = float(prior_score)
            state_grp.attrs["pair_score"] = float(pair_score)
            state_grp.attrs["exvol_score"] = float(exvol_score)
            state_grp.attrs["tet_score"] = float(tet_score)
            state_grp.attrs["oct_score"] = float(oct_score)
            
            # Save sigma as a subgroup.
            sigma_grp = state_grp.create_group("sigma")
            for key, value in sigmas.items():
                sigma_grp.attrs[key] = float(value)
            
            # Save positions as datasets (with gzip compression).
            pos_grp = state_grp.create_group("positions")
            for comp, coords in positions.items():
                pos_grp.create_dataset(comp, data=coords.astype(np.float32), compression="gzip")
            
            # Save types and bead_numbers as datasets.
            types_keys = list(types.keys())
            types_vals = [types[k] for k in types_keys]
            state_grp.create_dataset("types_keys", data=np.array(types_keys, dtype="S"))
            state_grp.create_dataset("types_vals", data=np.array(types_vals, dtype="S"))
            
            bead_keys = list(bead_numbers.keys())
            bead_vals = [bead_numbers[k] for k in bead_keys]
            state_grp.create_dataset("bead_keys", data=np.array(bead_keys))
            state_grp.create_dataset("bead_vals", data=np.array(bead_vals))
    
    except Exception as e:
        print(f"Warning: Failed to save state to HDF5: {e}")