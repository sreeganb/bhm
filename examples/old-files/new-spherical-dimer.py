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
import IMP.pmi.output
import RMF
import IMP.rmf
import random
#import two_level_mcmc as tlmcmc

# Define a class for the particle system
class ParticleSystem:
    def __init__(self):
        self.model = IMP.Model()
        self.hier = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model))
        self.particles = []
        self.rs = IMP.RestraintSet(self.model, "PairRestraintSet")

    def create_particles(self, num_particles):
        # create 4 particles with coordinates (0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3)
        i = 0
        for _ in range(num_particles):
            j = random.uniform(-11.0,10.0)
            k = random.uniform(-11.0,10.0)
            l = random.uniform(-11.0,10.0)
            p = IMP.Particle(self.model) # create a particle
            self.particles.append(p) # add the particle to the list of particles
            dr = IMP.core.XYZR.setup_particle(p) # create a decorator for the particle
            dr.set_coordinates(IMP.algebra.Vector3D(j, k, l))
            if i == 0 or i == num_particles // 2 - 1 or i == num_particles // 2 or i == num_particles - 1:
                dr.set_radius(1.8)
            else:
                dr.set_radius(1.0)
            dr.set_coordinates_are_optimized(True)
            IMP.atom.Mass.setup_particle(p, 0.5)
            # these are newly added, if the code breaks, comment them out
            hr = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "molecule" + str(i)))
            hr.add_child(IMP.atom.Hierarchy.setup_particle(dr))
            self.hier.add_child(hr)
            i = i + 1
        return self.particles

    def create_bayesian_restraints(self):
        # create a Bayesian restraint for the particle system
        # create a Jeffreys prior for the distance parameter between two adjacent particles
        hubound = IMP.core.HarmonicUpperBound(2.5, 20.0)
        hlbound = IMP.core.HarmonicLowerBound(1.1, 20.0)
        for i in range(len(self.particles) // 2 -1):
            r1 = IMP.core.DistanceRestraint(self.model, hubound, self.particles[i], self.particles[i+1])
            r2 = IMP.core.DistanceRestraint(self.model, hlbound, self.particles[i], self.particles[i+1])
            ambres = IMP.isd.AmbiguousRestraint(self.model, 20.0, r1, r2)
            self.rs.add_restraint(ambres)
            r1 = IMP.core.DistanceRestraint(self.model, hubound, self.particles[len(self.particles) - i -1], self.particles[len(self.particles) - i - 2])
            r2 = IMP.core.DistanceRestraint(self.model, hlbound, self.particles[len(self.particles) - i -1], self.particles[len(self.particles) - i - 2])
            ambres = IMP.isd.AmbiguousRestraint(self.model, 20.0, r1, r2)
            self.rs.add_restraint(ambres)
        # adding excluded volume restraint so that the particles do not overlap
        ev = IMP.core.ExcludedVolumeRestraint(self.particles, 20.0, 1, "EV")
        self.rs.add_restraint(ev)
        # Repulsive restraint between particles belonging to the two different strings
        hubound = IMP.core.HarmonicUpperBound(11.0, 10.0)
        hlbound = IMP.core.HarmonicLowerBound(6.0, 10.0)
        for i in range(len(self.particles) // 2):
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hubound, self.particles[i], self.particles[len(self.particles) // 2 + i]))
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hlbound, self.particles[i], self.particles[len(self.particles) // 2 + i]))
            print("connection between particles: ", self.particles[i], self.particles[len(self.particles) // 2 + i])
        # Create a nuisance particle to denote the end to end distance of the two strings
        # Give it a non-informative prior (Jeffrey's prior maybe)
        self.chi_is_sampled = True
        chiminnuis = 1.0
        chimaxnuis = 30.0
        chiinit = 15.0
        chimin = 0.01
        chimax = 100.0
        #chitrans = 0.5
        chi = IMP.pmi.tools.SetupNuisance(self.model, chiinit, chiminnuis, chimaxnuis, self.chi_is_sampled).get_particle()
        self.rs.add_restraint(IMP.isd.UniformPrior(self.model, chi, 10000.0, chimax, chimin))

        return self.rs

    def create_pairwise_restraints(self):
        # create pairwise restraints between the particles
        #self.particle_pairs = []
        #cpc = IMP.container.ConsecutivePairContainer(self.model, self.particles)
        #print("close pairs: ", cpc)
        #hps = IMP.core.HarmonicDistancePairScore(10.0, 100.0)
        #print("harmonic distance pair score: ", hps)    
        #self.rs = IMP.container.PairsRestraint(hps, cpc)
        # add a simple distance restraint between consecutive particles
        hubound = IMP.core.HarmonicUpperBound(4.0, 10.0)
        hlbound = IMP.core.HarmonicLowerBound(2.0, 10.0)
        for i in range(len(self.particles) - 1):
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hubound, self.particles[i], self.particles[i+1]))
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hlbound, self.particles[i], self.particles[i+1]))
        # adding excluded volume restraint so that the particles do not overlap
        ev = IMP.core.ExcludedVolumeRestraint(self.particles, 10.0, 1, "EV")
        self.rs.add_restraint(ev)
        return self.rs

    def run_mcmc_sampler(self, num_steps, temperature):
        # run an MCMC sampler to optimize the positions of the particles 
        # Extend this to a two-level MCMC sampler with two conformationally distinct states
        #restraints = self.create_pairwise_restraints()
        restraints = self.create_bayesian_restraints()
        mc = IMP.core.MonteCarlo(self.model)
        mc.set_log_level(IMP.SILENT)
        mc.set_scoring_function([restraints])
        mc.set_kt(temperature) # this is how temperature is set in MonteCarlo module
        bmvr = [IMP.core.BallMover(self.model, x, 0.5) for x in self.particles]
        #mv = IMP.core.SerialMover(bmvr)
        mc.add_movers(bmvr) # add the movers to the MonteCarlo object
        f = RMF.create_rmf_file("strings-dimeric-beads.rmf3")
        IMP.rmf.add_hierarchy(f, self.hier)
        IMP.rmf.add_restraints(f, [restraints])
        IMP.rmf.add_particles(f, self.particles)
        os = IMP.rmf.SaveOptimizerState(self.model, f)
        for i in range(10):
            os.update_always("initial conformation")
            os.set_log_level(IMP.SILENT)
            #mc.add_optimizer_state(os)
            mc.optimize(num_steps)
            print("number of accepted MCMC steps: ", mc.get_number_of_accepted_steps())
        
        return mc

    def get_optimized_coordinates(self):
        return [IMP.core.XYZR(p).get_coordinates() for p in self.particles]
    
    def define_global_parameter(self):
        # The end to end distance of the 10 particle setup is the global parameter
        # Calculate the end to end distance
        dfirst = IMP.core.XYZR(self.model, self.particles[0]).get_coordinates()
        dlast = IMP.core.XYZR(self.model, self.particles[(len(self.particles)) // 2 - 1]).get_coordinates()
        end_to_end1 = IMP.algebra.get_distance(dfirst, dlast) # calculate the distance between the first and last particles
        print(f"end_to_end1 = {end_to_end1} between particles {self.particles[0]} and {self.particles[(len(self.particles)) // 2 -1]}")
        dfirst = IMP.core.XYZR(self.model, self.particles[(len(self.particles)) // 2 ]).get_coordinates()
        dlast = IMP.core.XYZR(self.model, self.particles[-1]).get_coordinates()
        end_to_end2 = IMP.algebra.get_distance(dfirst, dlast) # calculate the distance between the first and last particles
        print(f"end_to_end2 = {end_to_end2} between particles {self.particles[(len(self.particles)) // 2 ]} and {self.particles[-1]}")
        return end_to_end1, end_to_end2

# Create an instance of the ParticleSystem class
particle_system = ParticleSystem()

# Create four particles
particle_system.create_particles(20)

# Get the optimized coordinates of the particles
optimized_coordinates = particle_system.get_optimized_coordinates()

# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")

# Run the MCMC sampler
sampler = particle_system.run_mcmc_sampler(2000, 2.0)

# Get the optimized coordinates of the particles
optimized_coordinates = particle_system.get_optimized_coordinates()
# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")
eedis1, eedis2 = particle_system.define_global_parameter()