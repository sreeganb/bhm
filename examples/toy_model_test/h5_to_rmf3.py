import h5py
import numpy as np
import IMP
import IMP.core
import IMP.atom
import IMP.algebra
import IMP.rmf  
import RMF
import sys

def convert_hdf5_to_rmf3(hdf5_file, rmf3_file):
    """
    Convert a custom HDF5 trajectory file to an RMF3 file compatible with ChimeraX.
    
    Assumptions:
    - Positions are stored per component, concatenated in order of sorted component keys.
    - Particle indices are 0 to n_particles-1, matching sorted types_keys (as strings '0', '1', ...).
    - Sigmas are radii per particle type.
    - Bead numbers are used to name particles as '{type}_{bead_number}'.
    - Radii and types are consistent across all states (only read from first state).
    - Coordinates update per state, but structure (particles) remains the same.
    
    Parameters:
    hdf5_file: Path to the input HDF5 file.
    rmf3_file: Path to the output RMF3 file.
    """
    with h5py.File(hdf5_file, 'r') as f:
        traj_grp = f['trajectory']
        states = sorted([k for k in traj_grp.keys() if k.startswith('state_')])
        if not states:
            print("No states found in trajectory.")
            return

        # Read first state to setup model and particles
        first_state_grp = traj_grp[states[0]]
        pos_grp = first_state_grp['positions']
        sigma_grp = first_state_grp['sigma']
        radii = {k: sigma_grp.attrs[k] for k in sigma_grp.attrs}  # radii per type

        # Read types
        types_keys = [k.decode('utf-8') for k in first_state_grp['types_keys'][:]]
        types_vals = [v.decode('utf-8') for v in first_state_grp['types_vals'][:]]
        types = dict(zip(types_keys, types_vals))

        # Read bead numbers (optional for naming)
        bead_keys = [str(k) for k in first_state_grp['bead_keys'][:]]  # Ensure str
        bead_vals = first_state_grp['bead_vals'][:]
        bead_numbers = dict(zip(bead_keys, bead_vals))

        # Sort keys as integers for consistent order
        sorted_indices = sorted(types_keys, key=int)
        n_particles = len(sorted_indices)

        # Get components and collect initial coordinates
        components = sorted(pos_grp.keys())
        init_coords = []
        for comp in components:
            init_coords.append(pos_grp[comp][:])
        all_init_coords = np.vstack(init_coords)
        if len(all_init_coords) != n_particles:
            raise ValueError("Number of positions does not match number of types.")

        # Create IMP Model
        model = IMP.Model()

        # Create root hierarchy
        p_root = IMP.Particle(model)
        root_h = IMP.atom.Hierarchy.setup_particle(p_root)
        p_root.set_name("root")

        # Create particles
        particles = []
        for idx_str in sorted_indices:
            i = int(idx_str)
            type_ = types[idx_str]
            radius_value = radii.get(type_, 1.0)  # Default radius if not found
            bead_num = bead_numbers.get(idx_str, i)  # Use bead number or index

            p = IMP.Particle(model)
            p.set_name(f"{type_}_{bead_num}")

            # Setup XYZR
            xyzr = IMP.core.XYZR.setup_particle(p)
            xyzr.set_coordinates(all_init_coords[i])
            xyzr.set_radius(radius_value)
            xyzr.set_coordinates_are_optimized(False)  # Not simulating

            # Optional mass
            IMP.atom.Mass.setup_particle(p, 1.0)

            # Setup hierarchy and add to root
            h = IMP.atom.Hierarchy.setup_particle(p)
            root_h.add_child(h)

            particles.append(p)

        # Create RMF file
        rmf = RMF.create_rmf_file(rmf3_file)
        rmf.set_description("Converted trajectory from HDF5 to RMF3.")

        # Add hierarchy to RMF
        IMP.rmf.add_hierarchy(rmf, root_h)

        # Save first frame
        IMP.rmf.save_frame(rmf, states[0])

        # Process remaining states
        for state_name in states[1:]:
            state_grp = traj_grp[state_name]
            pos_grp = state_grp['positions']

            # Collect coordinates for this state
            coords = []
            for comp in components:
                coords.append(pos_grp[comp][:])
            all_coords = np.vstack(coords)

            # Update coordinates
            for j, idx_str in enumerate(sorted_indices):
                i = int(idx_str)
                IMP.core.XYZ(particles[j]).set_coordinates(all_coords[i])

            # Update model (if needed)
            model.update()

            # Save frame
            IMP.rmf.save_frame(rmf, state_name)

        # Close RMF file
        rmf.close()

        print(f"Conversion complete. RMF3 file saved: {rmf3_file}")
# debugging function
def inspect_hdf5_structure(hdf5_file):
    """Inspect the structure of an HDF5 file to understand its layout."""
    def print_structure(name, obj):
        print(f"{name}: {type(obj)}")
        if hasattr(obj, 'attrs') and len(obj.attrs) > 0:
            print(f"  Attributes: {dict(obj.attrs)}")
        if hasattr(obj, 'shape'):
            print(f"  Shape: {obj.shape}, Dtype: {obj.dtype}")
    
    print(f"=== HDF5 Structure: {hdf5_file} ===")
    with h5py.File(hdf5_file, 'r') as f:
        f.visititems(print_structure)
def convert_fullsampler_hdf5_to_rmf3(hdf5_file, rmf3_file):
    """
    Convert FullSampler HDF5 trajectory to RMF3 format.
    Works with the structure: trajectory/state_XXXXX/positions/{A,B,C}
    """
    with h5py.File(hdf5_file, 'r') as f:
        traj_grp = f['trajectory']
        states = sorted([k for k in traj_grp.keys() if k.startswith('state_')])
        
        if not states:
            print("No states found in trajectory.")
            return
        
        print(f"Found {len(states)} states in trajectory")
        
        # Read first state to get structure
        first_state = traj_grp[states[0]]
        pos_grp = first_state['positions']
        
        # Define particle properties based on your system
        particle_types = []
        initial_coords = []
        radii = {'A': 24.0, 'B': 14.0, 'C': 16.0}  # Your actual radii
        
        # Process particles in order: A, B, C
        for particle_type in ['A', 'B', 'C']:
            if particle_type in pos_grp:
                coords = pos_grp[particle_type][:]
                n_particles = len(coords)
                particle_types.extend([particle_type] * n_particles)
                initial_coords.extend(coords)
                print(f"Found {n_particles} particles of type {particle_type}")
        
        initial_coords = np.array(initial_coords)
        n_total = len(initial_coords)
        
        print(f"Total particles: {n_total}")
        
        # Create RMF file FIRST
        rmf = RMF.create_rmf_file(rmf3_file)
        rmf.set_description("FullSampler trajectory converted from HDF5 to RMF3")
        
        # Create RMF node structure
        root_node = rmf.get_root_node()
        model_node = root_node.add_child("model", RMF.REPRESENTATION)
        
        # Add particle factory
        pf = RMF.ParticleFactory(rmf)
        
        # Create RMF particles
        rmf_particles = []
        particle_counter = {'A': 0, 'B': 0, 'C': 0}
        
        for i, ptype in enumerate(particle_types):
            # Create particle node
            particle_name = f"{ptype}_{particle_counter[ptype]}"
            particle_node = model_node.add_child(particle_name, RMF.REPRESENTATION)
            particle_counter[ptype] += 1
            
            # Set particle data
            particle_data = pf.get(particle_node)
            
            # Convert numpy array to RMF Vector3 - FIX HERE
            coord = initial_coords[i]
            rmf_coord = RMF.Vector3(float(coord[0]), float(coord[1]), float(coord[2]))
            particle_data.set_coordinates(rmf_coord)
            particle_data.set_radius(float(radii[ptype]))
            particle_data.set_mass(1.0)
            
            # Set particle name
            particle_node.set_name(particle_name)
            
            rmf_particles.append((particle_node, particle_data))
        
        # Add first frame
        rmf.add_frame("frame_0", RMF.FRAME)
        
        # Process remaining frames
        print(f"Converting {len(states)} frames...")
        for frame_idx, state_name in enumerate(states):
            if frame_idx % 10 == 0:
                print(f"Processing frame {frame_idx+1}/{len(states)}")
            
            # Add frame if not the first one
            if frame_idx > 0:
                rmf.add_frame(f"frame_{frame_idx}", RMF.FRAME)
            
            state_grp = traj_grp[state_name]
            pos_grp = state_grp['positions']
            
            # Collect coordinates for this state in the same order
            all_coords = []
            for particle_type in ['A', 'B', 'C']:
                if particle_type in pos_grp:
                    coords = pos_grp[particle_type][:]
                    all_coords.extend(coords)
            
            all_coords = np.array(all_coords)
            
            # Update coordinates for this frame - FIX HERE TOO
            for j, (particle_node, particle_data) in enumerate(rmf_particles):
                coord = all_coords[j]
                rmf_coord = RMF.Vector3(float(coord[0]), float(coord[1]), float(coord[2]))
                particle_data.set_coordinates(rmf_coord)
        
        # Flush and close
        rmf.flush()
        rmf.close()
        
        print(f"Conversion complete! RMF3 file saved: {rmf3_file}")
        print(f"Open with: chimerax {rmf3_file}")
                
# Example usage:
if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--inspect":
        # Usage: python h5_to_rmf3.py --inspect
        hdf5_file = "output_analysis/fullsampler_results_1/trajectory_chain_5.h5"
        inspect_hdf5_structure(hdf5_file)
    elif len(sys.argv) != 3:
        print("Usage: python h5_to_rmf3.py <input_hdf5_file> <output_rmf3_file>")
        print("   or: python h5_to_rmf3.py --inspect  (to inspect the HDF5 structure)")
        sys.exit(1)
    else:
        input_hdf5 = sys.argv[1]
        output_rmf3 = sys.argv[2]
        # Use the FullSampler-specific converter
        convert_fullsampler_hdf5_to_rmf3(input_hdf5, output_rmf3)