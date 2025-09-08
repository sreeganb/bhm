import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.atom
import IMP.algebra
import IMP.core
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints
import IMP.pmi.dof
import IMP.pmi.macros
import numpy as np
from simple_rex import SimpleReplicaExchange

#--------------------------------------------------------------------------
# Simple Distance Restraint Class that works with PMI infrastructure
#--------------------------------------------------------------------------
class SimpleDistanceRestraint(IMP.pmi.restraints.RestraintBase):
    """A simple distance restraint between two particles"""

    def __init__(self, model, particle1, particle2, 
                 distancemin=0, distancemax=100, kappa=1.0,
                 label=None, weight=1.):
        """Setup simple distance restraint between two particles.
        @param model The IMP model
        @param particle1 First particle (IMP.Particle)
        @param particle2 Second particle (IMP.Particle)
        @param distancemin The minimum distance
        @param distancemax The maximum distance
        @param kappa The harmonic parameter
        @param label A unique label for outputs and restraint names
        @param weight Weight of restraint
        """
        super().__init__(model, label=label, weight=weight)
        
        # Create upper and lower bound score functions
        ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
        ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)

        # Create and add the distance restraints
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts1,
                                     particle1, particle2))
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts2,
                                     particle1, particle2))
        
        # Print info
        particle1_name = particle1.get_name() if particle1.get_name() else "Particle1"
        particle2_name = particle2.get_name() if particle2.get_name() else "Particle2"
        print("Created simple distance restraint between %s and %s" % 
              (particle1_name, particle2_name))

#--------------------------------------------------------------------------
# Modified System Builder using simple spheres
#--------------------------------------------------------------------------
class SimpleParticleSystemBuilder:
    def __init__(self, ntype, copy_numbers, radii, colors):
        self.ntype = ntype
        self.copy_numbers = copy_numbers
        self.radii = radii  # radius for each particle type
        self.colors = colors

    def build_system(self, box_size=400.0):
        """Build system with simple spherical particles"""
        mdl = IMP.Model()
        
        # Create bounding box for initial placement
        bb = IMP.algebra.BoundingBox3D(
            IMP.algebra.Vector3D(-box_size, -box_size, -box_size),
            IMP.algebra.Vector3D(box_size, box_size, box_size)
        )
        
        particles = []
        particle_types = []
        movers = []
        
        for ptype in range(self.ntype):
            for copy in range(self.copy_numbers[ptype]):
                # Create simple particle
                p = IMP.Particle(mdl)
                p.set_name(f"type{ptype}_copy{copy}")
                
                # Set up as XYZR (coordinates + radius)
                random_pos = IMP.algebra.get_random_vector_in(bb)
                sphere = IMP.algebra.Sphere3D(random_pos, self.radii[ptype])
                d = IMP.core.XYZR.setup_particle(p, sphere)
                d.set_coordinates_are_optimized(True)
                movers.append(IMP.core.BallMover(mdl, p, self.radii[ptype]*2))
                movers[-1].set_was_used(True)
                IMP.display.Colored.setup_particle(
                    p, IMP.display.Color(1,0,0) if ptype==0 else
                       (IMP.display.Color(0,1,0) if ptype==1 else IMP.display.Color(0,0,1))
                )

                # Add to hierarchy so PMI tools work
                h = IMP.atom.Hierarchy.setup_particle(mdl, p)
                
                # Add mass (optional, for dynamics)
                mass = 4.0/3.0 * np.pi * (self.radii[ptype]**3)
                IMP.atom.Mass.setup_particle(mdl, p, mass)
                
                particles.append(p)
                particle_types.append(ptype)
        
        # Create a simple hierarchy root
        root = IMP.atom.Hierarchy.setup_particle(mdl, mdl.add_particle("root"))
        for p in particles:
            root.add_child(IMP.atom.Hierarchy(mdl, p))


        return mdl, root, particles, particle_types, movers

#--------------------------------------------------------------------------
# Scoring function class using PMI infrastructure
#--------------------------------------------------------------------------
class ScoringFunction:
    def __init__(self, model, particles, particle_types, ntype, copy_numbers):
        self.model = model
        self.particles = particles
        self.particle_types = particle_types
        self.ntype = ntype
        self.copy_numbers = copy_numbers
        self.output_objects = []
        
    def add_excluded_volume_restraint(self):
        """Use PMI's excluded volume restraint"""
        evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects=self.particles,
            resolution=1)
        evr.add_to_model()
        self.output_objects.append(evr)
        print(f"Added excluded volume restraint")
        return evr
    
    def add_pair_distance_restraints(self):
        """
        Add distance restraints between particle pairs
        @param interaction_params: dict with keys as (type1, type2) tuples
                                  and values as (dist_min, dist_max, kappa)
        """
        restraints_added = []
        
        # create a single distance restraint between a type 0
        # particle and a type 1 particle
        dr = SimpleDistanceRestraint(
            self.model, self.particles[0], self.particles[1],
            distancemin=0.0, distancemax=49.0, kappa=1.0,
            label="dist_type0_0_type1_0"
        )
        dr.add_to_model()
        self.output_objects.append(dr)
        
#        # Get particles by type
#        def get_particles_by_type(ptype):
#            return [p for i, p in enumerate(self.particles) 
#                   if self.particle_types[i] == ptype]
#        
#        for (type1, type2), (dist_min, dist_max, kappa) in interaction_params.items():
#            particles_type1 = get_particles_by_type(type1)
#            particles_type2 = get_particles_by_type(type2)
#            
#            if type1 == type2:
#                # Same type interactions (avoid double counting)
#                for i in range(len(particles_type1)):
#                    for j in range(i+1, len(particles_type1)):
#                        r = SimplePairDistanceRestraint(
#                            self.model, particles_type1[i], particles_type1[j],
#                            dist_min, dist_max, kappa,
#                            label=f"dist_type{type1}_{i}_{j}")
#                        r.add_to_model()
#                        self.output_objects.append(r)
#                        restraints_added.append(r)
#            else:
#                # Different type interactions
#                for i, p1 in enumerate(particles_type1):
#                    for j, p2 in enumerate(particles_type2):
#                        r = SimplePairDistanceRestraint(
#                            self.model, p1, p2,
#                            dist_min, dist_max, kappa,
#                            label=f"dist_type{type1}_{i}_type{type2}_{j}")
#                        r.add_to_model()
#                        self.output_objects.append(r)
#                        restraints_added.append(r)
#        
#        print(f"Added {len(restraints_added)} distance restraints")
        return dr

#--------------------------------------------------------------------------
# Main execution
#--------------------------------------------------------------------------
if __name__ == "__main__":
    # Your original system parameters
    ntype = 3
    copy_numbers = [8, 8, 16]
    radii = [24.0, 14.0, 16.0]  # Different radii for each type
    colors = ['red', 'blue', 'green']
    
    # Build system with simple particles
    builder = SimpleParticleSystemBuilder(ntype, copy_numbers, radii, colors)
    mdl, hierarchy, particles, particle_types, movers = builder.build_system(box_size=400.0)
    
    print(f"Created {len(particles)} particles")
    
    # Shuffle configuration
    IMP.pmi.tools.shuffle_configuration(hierarchy, max_translation=20.0)
    
    # Save initial configuration
    output = IMP.pmi.output.Output()
    output.init_rmf("initial_particles.rmf3", [hierarchy])
    output.write_rmf("initial_particles.rmf3")
    
    # Set up scoring function
    sf = ScoringFunction(mdl, particles, particle_types, ntype, copy_numbers)
    
    # Add excluded volume using PMI's implementation
    sf.add_excluded_volume_restraint()
    sf.add_pair_distance_restraints()
#    
#    # Define interaction parameters: (type1, type2): (min_dist, max_dist, kappa)
#    # Type 0 = A (8 copies), Type 1 = B (8 copies), Type 2 = C (16 copies)
#    interactions = {
#        (0, 0): (8.0, 15.0, 1.0),   # A-A interactions
#        (0, 1): (10.0, 20.0, 0.5),  # A-B interactions  
#        (1, 2): (6.0, 12.0, 2.0),   # B-C interactions
#        # Add more as needed
#    }
#    
#    sf.add_pair_distance_restraints(interactions)
#    
#    # Set up degrees of freedom using PMI
    dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
    for p in particles:
        dof.create_rigid_body(p)
    #sm = IMP.core.SerialMover(movers)
    print(f"Degrees of freedom: {dof.get_movers()}")
#
#    # Create flexible beads or rigid bodies
#    for particle in particles:
#        # Make each particle flexible (it needs to be a hierarchy)
#        h = IMP.atom.Hierarchy(mdl, particle)
#        dof.create_flexible_bead(h, max_trans=5.0)
#        
#    print(f"Degrees of freedom: {dof.get_movers()}")
#    
#    # Save shuffled configuration
#    output_shuffled = IMP.pmi.output.Output()
#    output_shuffled.init_rmf("shuffled_particles.rmf3", [hierarchy])
#    output_shuffled.write_rmf("shuffled_particles.rmf3")
#    
#    # Run replica exchange using PMI's built-in functionality
#    rex = IMP.pmi.macros.ReplicaExchange(
#        mdl, 
#        root_hier=hierarchy,
#        monte_carlo_sample_objects=dof.get_movers(),
#        output_objects=sf.output_objects,
#        monte_carlo_steps=10,
#        number_of_frames=100,
#        global_output_directory="output/",
#        atomistic=False)
    
    rex = SimpleReplicaExchange(
        model=mdl,
        root_hier=hierarchy,
        monte_carlo_sample_objects=dof.get_movers(),
        output_objects=sf.output_objects,
        monte_carlo_steps=10,
        number_of_frames=100,
        global_output_directory="output/"
    )

    rex.execute_macro()
