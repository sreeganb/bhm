#!/usr/bin/env python
#------------------------------------------
# test_spherical_particles.py 
# create a simple particle system with four particles and run an MCMC sampler
# to optimize their positions. Furthermore, make this simple MC sampler into a
# two level MCMC sampler with two conformationally distinct states
#------------------------------------------
import IMP
import IMP
import IMP.isd
import IMP.core
import IMP.container
import IMP.atom
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry

# Define a class for the particle system
class ParticleSystem:
    def __init__(self):
        self.model = IMP.Model()
        self.particles = []
        self.rs = IMP.RestraintSet(self.model, "PairRestraintSet")
    
    def create_particles(self, num_particles):
        for _ in range(num_particles):
            p = IMP.Particle(self.model)
            self.particles.append(p)
            IMP.core.XYZ.setup_particle(p, IMP.algebra.Vector3D(1.0, 1.0, 0.0))
            IMP.atom.Mass.setup_particle(p, 1.0)
        print("particles = : ", self.particles)
        return self.particles

    def create_pairwise_restraints(self):
        self.particle_pairs = []
        #print("model = : ", self.model())
        for i in range(len(self.particles) - 1):
            hu = IMP.core.HarmonicUpperBound(1.5, 0.1)
            dps = IMP.core.DistancePairScore(hu)
            self.particle_pairs.append((self.particles[i], self.particles[i+1]))
            r = IMP.core.PairRestraint(self.model, dps, (self.particles[i].get_index(), self.particles[i+1].get_index()))
            print("restraints = : ", r)
            self.rs.add_restraint(r)
        return self.rs, dps

    def run_mcmc_sampler(self, num_steps, temperature, mvr):
        #mc.set_temperature(temperature)
        restraints, dps = self.create_pairwise_restraints()
        mc = IMP.core.MonteCarlo(self.model)
        mc.add_movers(mvr)
        mc.set_scoring_function(dps)
        
        #for r in restraints:
        #    mc.add_restraint(r)
        mc.optimize(num_steps)
        return mc

    def get_optimized_coordinates(self):
        return [IMP.core.XYZ(p).get_coordinates() for p in self.particles]
        
# Create an instance of the ParticleSystem class
particle_system = ParticleSystem()

# Create four particles
particle_system.create_particles(4)

# Define movers for the particles
movers = []
radius = 1.0
for p in particle_system.particles:
    print("particle = : ", p.get_index())
    mover = IMP.core.BallMover(particle_system.model, p.get_index(), radius)
    movers.append(mover)

# Run the MCMC sampler
sampler = particle_system.run_mcmc_sampler(1000, 1.0, movers)

# Get the optimized coordinates of the particles
optimized_coordinates = particle_system.get_optimized_coordinates()

# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")