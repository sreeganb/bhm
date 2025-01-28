import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.atom
import IMP.algebra
import IMP.core
import IMP.pmi.tools
import IMP.pmi.output

# Define particle properties
ids = ["A", "B", "C"]       # Particle types
radii = [20.0, 4.0, 10.0]    # Radii for each type
counts = [8, 8, 16]         # Number of particles for each type
color_map = {
    "A": "red",
    "B": "green",
    "C": "blue"
}

# Create the IMP model and PMI system
mdl = IMP.Model()
system = IMP.pmi.topology.System(mdl)
state = system.create_state()

# Add spherical particles for each type
for p_type, radius, count in zip(ids, radii, counts):
    for i in range(count):
        name = f"{p_type}_{i+1}"
        particle = state.create_molecule(name=name, sequence="X", chain_id=p_type)
        color = color_map.get(p_type, "gray")  # Default color if not found
        particle.add_representation(resolutions=[1], color=color)

# Build the hierarchy
hierarchy = system.build()
print(hierarchy.get_child(0).get_children())

# Assign random positions/radii to each molecule
bbox = IMP.algebra.BoundingBox3D(
    IMP.algebra.Vector3D(-100, -100, -100),
    IMP.algebra.Vector3D(100, 100, 100)
)
for molecule in IMP.atom.get_by_type(hierarchy, IMP.atom.MOLECULE_TYPE):
    p_type = molecule.get_name().split('_')[0]
    radius = radii[ids.index(p_type)] 

    # Get the child particle representation
    particle = molecule.get_child(0).get_particle()
    xyzr = IMP.core.XYZR.setup_particle(particle)
    xyzr.set_radius(radius)
    xyzr.set_coordinates(IMP.algebra.get_random_vector_in(bbox))

# Shuffle configuration to reduce overlaps
IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=20.0)

# Set up RMF output
output = IMP.pmi.output.Output()
output.init_rmf("spherical_particles.rmf3", [hierarchy])  # Initialize the RMF file
output.write_rmf("spherical_particles.rmf3")  # Write the RMF file

# Optional: Print the particle hierarchy
#IMP.atom.show_with_representations(hierarchy)