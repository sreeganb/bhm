import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.atom
import IMP.algebra
import IMP.core
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.test
import numpy as np
import IMP.core
import IMP.pmi.restraints.basic
                                  
#--------------------------------------------------------------------------
class System_builder:
    def __init__(self, ntype, copy_numbers, chain_ids, sequences_type, colors):
        self.ntype = ntype # number of particle types
        self.copy_numbers = copy_numbers # copy number for each particle type
        self.chain_ids = chain_ids # chain id for each copy of each particle type
        self.colors = colors # color for each particle type
        self.seq_type = sequences_type # sequence for each particle type
        
    def build_system(self):
        mdl = IMP.Model()
        system = IMP.pmi.topology.System(mdl)
        state = system.create_state()
        # create a list of the particles
        parts = []
        
        for i in range(self.ntype):
            name = f"prot_{i+1}"
            particle = state.create_molecule(name=name, sequence=self.seq_type[i], chain_id=self.chain_ids[i])
            color = self.colors[i]
            particle.add_representation(resolutions=[len(self.seq_type[i])], color=color)
            parts.append(particle)
            for j in range(self.copy_numbers[i]-1):
                cloned_particle = particle.create_clone(chain_id=self.chain_ids[j+1])
                parts.append(cloned_particle)

        return mdl, system, parts
#--------------------------------------------------------------------------
# Class for the scoring function, at this point just some excluded volume
#--------------------------------------------------------------------------
class scoring_function():
    def __init__(self, system, model, hierarchy, output_objects, parts, ntype, copy_numbers):
        self.output_objects = output_objects
        self.system = system
        self.parts = parts
        self.ntype = ntype
        self.copy_numbers = copy_numbers
        self.m = model
        self.hier = hierarchy

    def add_excluded_volume_restraint(self):
        evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects = self.parts,
            resolution = 1)
        evr.add_to_model()
        self.output_objects.append(evr)
        print("added excluded volume restraint", evr)
        
        return evr, self.output_objects
    
    def add_pair_distance_restraint(self, aa_distance=40.0, ab_distance=30.0, bc_distance=25.0):
        """Add distance restraints between different particle types with specified thresholds
        
        Args:
            aa_distance: Maximum distance for A-A interactions (default 40.0 Å)
            ab_distance: Maximum distance for A-B interactions (default 30.0 Å) 
            bc_distance: Maximum distance for B-C interactions (default 25.0 Å)
        """
        pairs = []
        
        # A-A type interactions (prot_1 with prot_1)
        print("Adding A-A interactions...")
        for i in range(self.copy_numbers[0]):
            for j in range(i + 1, self.copy_numbers[0]):
                tuple1 = (1, 1, "prot_1", i)
                tuple2 = (1, 1, "prot_1", j)
                
                dr = IMP.pmi.restraints.basic.DistanceRestraint(
                    self.hier, 
                    tuple1, 
                    tuple2, 
                    distancemin=0.0,
                    distancemax=aa_distance,
                    kappa=10.0,
                    resolution=1.0,
                    label=f"dist_AA_prot_1_{i}_{j}"
                )
                dr.add_to_model()
                self.output_objects.append(dr)
                pairs.append(("A-A", i, j))
                print(f"Added A-A distance restraint: prot_1 copy {i} - copy {j}, max distance: {aa_distance} Å")
        
        # A-B type interactions (prot_1 with prot_2)
        print("Adding A-B interactions...")
        for i in range(self.copy_numbers[0]):  # A particles (prot_1)
            for j in range(self.copy_numbers[1]):  # B particles (prot_2)
                tuple1 = (1, 1, "prot_1", i)
                tuple2 = (1, 1, "prot_2", j)
                
                dr = IMP.pmi.restraints.basic.DistanceRestraint(
                    self.hier,
                    tuple1,
                    tuple2,
                    distancemin=0.0,
                    distancemax=ab_distance,
                    kappa=10.0,
                    resolution=1.0,
                    label=f"dist_AB_prot_1_{i}_prot_2_{j}"
                )
                dr.add_to_model()
                self.output_objects.append(dr)
                pairs.append(("A-B", i, j))
                print(f"Added A-B distance restraint: prot_1 copy {i} - prot_2 copy {j}, max distance: {ab_distance} Å")
        
        # B-C type interactions (prot_2 with prot_3)
        print("Adding B-C interactions...")
        for i in range(self.copy_numbers[1]):  # B particles (prot_2)
            for j in range(self.copy_numbers[2]):  # C particles (prot_3)
                tuple1 = (1, 1, "prot_2", i)
                tuple2 = (1, 1, "prot_3", j)
                
                dr = IMP.pmi.restraints.basic.DistanceRestraint(
                    self.hier,
                    tuple1,
                    tuple2,
                    distancemin=0.0,
                    distancemax=bc_distance,
                    kappa=10.0,
                    resolution=1.0,
                    label=f"dist_BC_prot_2_{i}_prot_3_{j}"
                )
                dr.add_to_model()
                self.output_objects.append(dr)
                pairs.append(("B-C", i, j))
                print(f"Added B-C distance restraint: prot_2 copy {i} - prot_3 copy {j}, max distance: {bc_distance} Å")
        
        # Summary
        aa_count = len([p for p in pairs if p[0] == "A-A"])
        ab_count = len([p for p in pairs if p[0] == "A-B"])
        bc_count = len([p for p in pairs if p[0] == "B-C"])
        
        print(f"\nSummary of distance restraints added:")
        print(f"A-A restraints: {aa_count} (max distance: {aa_distance} Å)")
        print(f"A-B restraints: {ab_count} (max distance: {ab_distance} Å)")
        print(f"B-C restraints: {bc_count} (max distance: {bc_distance} Å)")
        print(f"Total restraints: {len(pairs)}")
        
        return pairs

if __name__ == "__main__":
    ntype = 3
    copy_numbers = [8, 8, 16]
    chain_ids = []
    for i in range(np.array(copy_numbers).sum()):
        chain_ids.append(chr(ord('A') + i))
    sequences_type = ['AA', 'L', 'A']
    colors = ['red', 'blue', 'green']
    
    builder = System_builder(ntype, copy_numbers, chain_ids, sequences_type, colors)
    mdl, system, parts = builder.build_system()
    
    hierarchy = system.build()
    IMP.atom.show_with_representations(hierarchy)
    
    output= IMP.pmi.output.Output()
    output.init_rmf("initial_particles.rmf3", [hierarchy])  
    output.write_rmf("initial_particles.rmf3")  # Write the RMF file
    
    IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=20.0)
    #IMP.atom.show_with_representations(hierarchy)

    # Set up RMF output
    output = IMP.pmi.output.Output()
    output.init_rmf("spherical_particles.rmf3", [hierarchy])  # Initialize the RMF file
    output.write_rmf("spherical_particles.rmf3")  # Write the RMF file
    
    output_objects = []

    sf = scoring_function(system, mdl, hierarchy, output_objects, parts, ntype, copy_numbers)
    sf.add_excluded_volume_restraint()
    sf.add_pair_distance_restraint()
    
    dof_s1 = IMP.pmi.dof.DegreesOfFreedom(mdl)
    # create flexible beads DOF
    print("what is inside parts:", parts)
    for particle in parts:
        dof_s1.create_rigid_body(particle)
    print("degrees of freedom", dof_s1.get_movers())
    
    rex = IMP.pmi.macros.ReplicaExchange(mdl, root_hier=hierarchy,
                                            monte_carlo_sample_objects=dof_s1.get_movers(),
                                            output_objects=output_objects,
                                            number_of_frames=1000)
    rex.execute_macro()