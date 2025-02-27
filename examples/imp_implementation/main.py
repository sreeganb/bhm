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
                particle.create_clone(chain_id=self.chain_ids[j+1])
                parts.append(particle)
                
        return mdl, system, parts
#--------------------------------------------------------------------------
# Class for the scoring function, at this point just some excluded volume
#--------------------------------------------------------------------------
class scoring_function():
    def __init__(self, system, output_objects, parts):
        self.output_objects = output_objects
        self.system = system
        self.parts = parts
        
    def add_excluded_volume_restraint(self):
        evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects = self.parts,
            resolution = 1)
        self.output_objects.append(evr)
        print("added excluded volume restraint", evr)
        
        return evr, self.output_objects
    
    def add_pair_distance_restraint(self, part1, part2):
        
        
        #rex = IMP.pmi.macros.ReplicaExchange(mdl, root_hier = hierarchy, 
        #                                     monte_carlo_sample_objects=dof.get_movers(), 
        #                                     output_objects=self.output_objects, 
        #                                     number_of_frames=100)
        #rex.execute_macro()
    
    # For the particles in the system hierarchy, calculate the score
#    def calculate_score_for_hierarchy(self, hierarchy):
#        score = 0.0
#        for restraint in self.restraints:
#            score += restraint.get_score()
#        return score
    
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
        
    dof_s1 = IMP.pmi.dof.DegreesOfFreedom(mdl)
    print("degrees of freedom", dof_s1)
    IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=20.0)
    #IMP.atom.show_with_representations(hierarchy)

    # Set up RMF output
    output = IMP.pmi.output.Output()
    output.init_rmf("spherical_particles.rmf3", [hierarchy])  # Initialize the RMF file
    output.write_rmf("spherical_particles.rmf3")  # Write the RMF file
    
    output_objects = []

    sf = scoring_function(system, output_objects, parts)
    sf.add_excluded_volume_restraint()
    
    


#for j in range(len(counts)):
#    for i in range(counts[j]):
#        m1A.create_clone(chain_id=chain_ids[i+1])

# Add spherical particles for each type
#for p_type, radius, count in zip(ids, radii, counts):
#    for i in range(count):
#        # create copies after the first particle is created
#        # clones and copies : clones are identical to the original particle, copies are not
#        name = f"{p_type}_{i+1}"
#        particle = state.create_molecule(name=name, sequence="X", chain_id=p_type)
#        color = color_map.get(p_type, "gray")  # Default color if not found
#        particle.add_representation(resolutions=[1], color=color)

# Assign random positions/radii to each molecule
#bbox = IMP.algebra.BoundingBox3D(
#    IMP.algebra.Vector3D(-100, -100, -100),
#    IMP.algebra.Vector3D(100, 100, 100)
#)

#for molecule in IMP.atom.get_by_type(hierarchy, IMP.atom.MOLECULE_TYPE):
#    p_type = molecule.get_name().split('_')[0]
#    radius = radii[ids.index(p_type)] 

    # Get the child particle representation
#    particle = molecule.get_child(0).get_particle()
#    xyzr = IMP.core.XYZR.setup_particle(particle)
#    xyzr.set_radius(radius)
#    xyzr.set_coordinates(IMP.algebra.get_random_vector_in(bbox))