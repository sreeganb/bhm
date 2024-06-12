#!/usr/bin/env python
#------------------------------------------
# spherical_dimer.py 
# create a simple particle system with a dimer of ten particles 
# and run an MCMC sampler
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
        self.particles1 = []
        self.particles2 = []
        self.rs = IMP.RestraintSet(self.model, "PairRestraintSet")

    def string_of_beads(self, num_particles, iter):
        # create a string of beads with 10 particles
        pts = []
        i = 0
        for _ in range(num_particles):
            j = random.uniform(-10.0,10.0)
            k = random.uniform(-10.0,10.0)
            l = random.uniform(-10.0,10.0)
            p = IMP.Particle(self.model)
            dr = IMP.core.XYZR.setup_particle(p)
            dr.set_coordinates(IMP.algebra.Vector3D(j, k, l))
            if i == 0 or i == num_particles - 1:
                dr.set_radius(2.0)
            else:
                dr.set_radius(1.0)
            dr.set_coordinates_are_optimized(True)
            IMP.atom.Mass.setup_particle(p, 0.5)
            pts.append(p)
            hr = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "molecule" + str(i) + "_iter" + str(iter)))
            hr.add_child(IMP.atom.Hierarchy.setup_particle(dr))
            self.hier.add_child(hr)
            i = i + 1

        return pts

    def create_particles(self, num_particles):
        # create two sets of particle systems
        self.particles1 = self.string_of_beads(num_particles, 1)
        self.particles2 = self.string_of_beads(num_particles, 2)        
        # create hierarchies for set 1 and set 2
        #hier_set1 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "set1"))
        #hier_set2 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "set2"))
        
        #for i in range(num_particles):
        #    hr1 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "molecule" + str(i) + "_set1"))
        #    hr1.add_child(IMP.atom.Hierarchy.setup_particle(dr1))
        #    hier_set1.add_child(hr1)
            
        #    hr2 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "molecule" + str(i) + "_set2"))
        #    hr2.add_child(IMP.atom.Hierarchy.setup_particle(dr2))
        #    hier_set2.add_child(hr2)
        
        #self.hier.add_child(hier_set1)
        #self.hier.add_child(hier_set2)
        #print("hierarchy: ", hier_set1)
        
        #print("particles set 1 = : ", self.particles1)
        #print("particles set 2 = : ", self.particles2)
        
        return self.particles1, self.particles2
    
    def create_bayesian_restraints(self):
        # create a Bayesian restraint for the particle system
        # create a Jeffreys prior for the distance parameter between two adjacent particles
        hubound = IMP.core.HarmonicUpperBound(4.8, 10.0)
        hlbound = IMP.core.HarmonicLowerBound(1.8, 10.0)
        for i in range(len(self.particles1) -1):
            r1 = IMP.core.DistanceRestraint(self.model, hubound, self.particles1[i], self.particles1[i+1])
            r2 = IMP.core.DistanceRestraint(self.model, hlbound, self.particles1[i], self.particles1[i+1])
            r3 = IMP.core.DistanceRestraint(self.model, hubound, self.particles2[i], self.particles2[i+1])
            r4 = IMP.core.DistanceRestraint(self.model, hlbound, self.particles2[i], self.particles2[i+1])
            ambres1 = IMP.isd.AmbiguousRestraint(self.model, 10.0, r1, r2)
            self.rs.add_restraint(ambres1)
            ambres2 = IMP.isd.AmbiguousRestraint(self.model, 10.0, r3, r4)
            self.rs.add_restraint(ambres2)
            #print("ambiguous restraint: ", ambres)
        # adding excluded volume restraint so that the particles do not overlap
        ev = IMP.core.ExcludedVolumeRestraint(self.particles1, 10.0, 1, "EV")
        self.rs.add_restraint(ev)
        ev = IMP.core.ExcludedVolumeRestraint(self.particles2, 10.0, 1, "EV")
        self.rs.add_restraint(ev)
        # add a distance restraint between the middle particles of the two strings
        self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hubound, self.particles1[4], self.particles2[4]))
        self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hlbound, self.particles1[4], self.particles2[4]))

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
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hubound, self.particles1[i], self.particles1[i+1]))
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hlbound, self.particles1[i], self.particles1[i+1]))
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hubound, self.particles2[i], self.particles2[i+1]))
            self.rs.add_restraint(IMP.core.DistanceRestraint(self.model, hlbound, self.particles2[i], self.particles2[i+1]))
        # adding excluded volume restraint so that the particles do not overlap
        ev = IMP.core.ExcludedVolumeRestraint(self.particles1, 10.0, 1, "EV")
        self.rs.add_restraint(ev)
        return self.rs

    def run_mcmc_sampler(self, num_steps, temperature):
        # run an MCMC sampler to optimize the positions of the particles 
        # Extend this to a two-level MCMC sampler with two conformationally distinct states
        #restraints = self.create_pairwise_restraints()
        restraints = self.create_bayesian_restraints()
        #print("ambiguous restraint: ", ambres)
        mc = IMP.core.MonteCarlo(self.model)
        mc.set_log_level(IMP.SILENT)
        mc.set_scoring_function([restraints])
        mc.set_kt(temperature) # this is how temperature is set in MonteCarlo module
        # add ball movers for two sets of strings of particles
        bmvr = [IMP.core.BallMover(self.model, x, 0.5) for x in self.particles1 + self.particles2]
        mc.add_movers(bmvr) # add the movers to the MonteCarlo object
        f = RMF.create_rmf_file("balls-move.rmf3")
        IMP.rmf.add_restraints(f, [restraints])
        IMP.rmf.add_particles(f, self.particles1 + self.particles2)
        #IMP.rmf.add_particles(f, self.particles2)
        IMP.rmf.add_hierarchy(f, self.hier)
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
        return [IMP.core.XYZR(p).get_coordinates() for p in self.particles1], [IMP.core.XYZR(p).get_coordinates() for p in self.particles2]
    
    def define_global_parameter(self):
        # The end to end distance of the 10 particle setup is the global parameter
        # Calculate the end to end distance
        dfirst = IMP.core.XYZR(self.model, self.particles[0]).get_coordinates()
        dlast = IMP.core.XYZR(self.model, self.particles[-1]).get_coordinates()
        end_to_end = IMP.algebra.get_distance(dfirst, dlast) # calculate the distance between the first and last particles
        return end_to_end

# Create an instance of the ParticleSystem class
particle_system = ParticleSystem()

# Create four particles
particle_system.create_particles(10)

# Get the optimized coordinates of the particles
optimized_coordinates = particle_system.get_optimized_coordinates()
#print("Initial optimized_coordinates = : ", optimized_coordinates)

# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")

# Run the MCMC sampler
sampler = particle_system.run_mcmc_sampler(200, 2.0)

# Get the optimized coordinates of the particles
#optimized_coordinates1, optimized_coordinates2 = particle_system.get_optimized_coordinates()
#print("optimized_coordinates = : ", optimized_coordinates)
# Print the optimized coordinates
for i, coord in enumerate(optimized_coordinates):
    print(f"Particle {i+1}: {coord}")
eedis = particle_system.define_global_parameter()
print("end to end distance as calculated: ", eedis)