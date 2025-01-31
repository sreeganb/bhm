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

mdl = IMP.Model()
system = IMP.pmi.topology.System(mdl)
state = system.create_state()
b1 = IMP.atom.Hierarchy.create_molecule(bead_radius = 20.0, Model=mdl)
print(b1)

m1A = state.create_molecule(name='ProtA', sequence='X', chain_id='A')
m1A.add_representation(resolutions=[1], color='red')

# Clone molecule B with a sequence too.
m2A = m1A.create_clone(chain_id='B')

# Build the hierarchy
hierarchy = system.build()
print(hierarchy.get_children())

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

# Shuffle configuration to reduce overlaps
IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=20.0)

# Set up RMF output
output = IMP.pmi.output.Output()
output.init_rmf("spherical_particles.rmf3", [hierarchy])  # Initialize the RMF file
output.write_rmf("spherical_particles.rmf3")  # Write the RMF file

# Optional: Print the particle hierarchy
#IMP.atom.show_with_representations(hierarchy)

#--------------------------------------------------------------------------
# Class for the scoring function, at this point just some excluded volume
#--------------------------------------------------------------------------
class scoring_function():
    def __init__(self, system):
        self.system = system
        self.restraints = []
        self.output_objects = []
        self.rmf_output_objects = []

    def add_restraint(self, restraint, mdl, dof, hierarchy):
        self.restraints.append(restraint)
        restraint.add_to_model()
        self.output_objects.append(restraint)
        self.rmf_output_objects.append(restraint)
        rex = IMP.pmi.macros.ReplicaExchange(mdl, root_hier = hierarchy, monte_carlo_sample_objects=dof.get_movers(), 
                                             output_objects=self.output_objects, number_of_frames=100)
        rex.execute_macro()
        

    def calculate_score(self):
        return sum(r.get_score() for r in self.restraints)
    
    # For the particles in the system hierarchy, calculate the score
    def calculate_score_for_hierarchy(self, hierarchy):
        score = 0.0
        for restraint in self.restraints:
            score += restraint.get_score()
        return score

dof_s1 = IMP.pmi.dof.DegreesOfFreedom(mdl)
IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=10.0)
IMP.atom.show_with_representations(hierarchy)

sf = scoring_function(system)
sf.add_restraint(IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(hierarchy), mdl, dof_s1, hierarchy)