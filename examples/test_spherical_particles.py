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
            j = random.uniform(-10.0,10.0)
            k = random.uniform(-10.0,10.0)
            l = random.uniform(-10.0,10.0)
            p = IMP.Particle(self.model) # create a particle
            self.particles.append(p) # add the particle to the list of particles
            dr = IMP.core.XYZR.setup_particle(p) # create a decorator for the particle
            dr.set_coordinates(IMP.algebra.Vector3D(j, k, l))
            dr.set_radius(1.0)
            dr.set_coordinates_are_optimized(True)
            IMP.atom.Mass.setup_particle(p, 0.5)
            # these are newly added, if the code breaks, comment them out
            hr = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "molecule" + str(i)))
            hr.add_child(IMP.atom.Hierarchy.setup_particle(dr))
            self.hier.add_child(hr)
            i = i + 1
        print("particles = : ", self.particles)
        return self.particles

    def bayesian_restraints(self):
        # create a Bayesian restraint for the particle system
        # create a Jeffreys prior for the distance parameter between two adjacent particles
        

        return 0


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
        restraints = self.create_pairwise_restraints()
        mc = IMP.core.MonteCarlo(self.model)
        mc.set_log_level(IMP.SILENT)
        mc.set_scoring_function([restraints])
        mc.set_kt(temperature) # this is how temperature is set in MonteCarlo module
        bmvr = [IMP.core.BallMover(self.model, x, 0.5) for x in self.particles]
        #mv = IMP.core.SerialMover(bmvr)
        mc.add_movers(bmvr) # add the movers to the MonteCarlo object
        f = RMF.create_rmf_file("balls-move.rmf3")
        IMP.rmf.add_hierarchy(f, self.hier)
        IMP.rmf.add_restraints(f, [restraints])
        IMP.rmf.add_particles(f, self.particles)
        #o = IMP.pmi.output.Output()
        os = IMP.rmf.SaveOptimizerState(self.model, f)
        for i in range(10):
            os.update_always("initial conformation")
            os.set_log_level(IMP.SILENT)
            #os.set_simulator(mc)
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
        dlast = IMP.core.XYZR(self.model, self.particles[-1]).get_coordinates()
        end_to_end = IMP.algebra.get_distance(dfirst, dlast) # calculate the distance between the first and last particles
        print("end to end distance: ", end_to_end)
        return end_to_end

# Create an instance of the ParticleSystem class
particle_system = ParticleSystem()

# Create four particles
particle_system.create_particles(10)

# Get the optimized coordinates of the particles
optimized_coordinates = particle_system.get_optimized_coordinates()
print("Initial optimized_coordinates = : ", optimized_coordinates)

# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")

# Run the MCMC sampler
sampler = particle_system.run_mcmc_sampler(200, 2.0)

# Get the optimized coordinates of the particles
optimized_coordinates = particle_system.get_optimized_coordinates()
print("optimized_coordinates = : ", optimized_coordinates)

# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")
eedis = particle_system.define_global_parameter()
print("end to end distance as calculated: ", eedis)
