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
    def __init__(self, system, output_objects, parts, ntype, copy_numbers):
        self.output_objects = output_objects
        self.system = system
        self.parts = parts
        self.ntype = ntype
        self.copy_numbers = copy_numbers
        
    def add_excluded_volume_restraint(self):
        evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects = self.parts,
            resolution = 1)
        evr.add_to_model()
        self.output_objects.append(evr)
        print("added excluded volume restraint", evr)
        
        return evr, self.output_objects
    
    def add_pair_distance_restraint(self):
        """all vs all, so we have A-A, A-B and B-C restraints that we shall add.
        create particle pairs accordingly"""
        pairs = []
        # A-A type interactions
        for i in range(copy_numbers[0]):
            for j in range(i + 1, copy_numbers[0]):
                part1 = parts[i]
                part2 = parts[j]
                pairs.append((part1, part2))
                print(f"pair added is : ", part1, part2)
        # A-B type interactions
        for i in range(copy_numbers[0]):
            for j in range(copy_numbers[1]):
                part1 = parts[i]
                part2 = parts[j + copy_numbers[0]]
                pairs.append((part1, part2))
                print(f"pair added is : ", part1, part2)
        # B-C type interactions
        for i in range(copy_numbers[1]):
            for j in range(copy_numbers[2]):
                part1 = parts[i + copy_numbers[0]]
                part2 = parts[j + copy_numbers[0] + copy_numbers[1]]
                pairs.append((part1, part2))
                print(f"pair added is : ", part1, part2)

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

    sf = scoring_function(system, output_objects, parts, ntype, copy_numbers)
    sf.add_excluded_volume_restraint()
    sf.add_pair_distance_restraint()
    
    dof_s1 = IMP.pmi.dof.DegreesOfFreedom(mdl)
    # create flexible beads DOF
    print("what is inside parts:", parts)
    for particle in parts:
        dof_s1.create_rigid_body(particle)
    print("degrees of freedom", dof_s1.get_movers())
    
#    rex = IMP.pmi.macros.ReplicaExchange(mdl, root_hier=hierarchy,
#                                            monte_carlo_sample_objects=dof_s1.get_movers(),
#                                            output_objects=output_objects,
#                                            number_of_frames=100)
#    rex.execute_macro()