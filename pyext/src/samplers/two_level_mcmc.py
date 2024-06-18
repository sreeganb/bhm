#!/usr/bin/env python
"""MCMC protocol for the Bayesian hierarchical modeling
   This involves two stages of parameter estimation.
   The first stage is a simple MCMC step, whereas the 
   second stage involves the joint posterior estimation of the 
   local and global parameters. Use the core IMP library for 
   most of the MCMC steps."""
import IMP
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.restraints
import RMF
import IMP.rmf
import IMP.pmi.tools

class MCMCsampler():
    # MCMC sampler for the hierarchical system
    def __init__(self, root_hier, temperature, num_steps):
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        #****************************************************************************************
        # IMP.pmi.tools is where most of the extraction of information from the hierarchy is done
        # First try to just add movers for the first system and check if the MCMC will work
        #****************************************************************************************
        self.rs = IMP.pmi.tools.get_restraint_set(self.m)
        num_strings = root_hier.get_number_of_children()
        num_beads = root_hier.get_child(0).get_number_of_children()
        self.particles = []
        for i in range(num_strings):
            string = root_hier.get_child(i)
            for j in range(num_beads):
                self.particles.append(string.get_child(j).get_particle())
            #self.particles.append([string.get_child(j).get_particle() for j in range(num_beads)])
        # Setup the MCMC parameters
        IMP.pmi.tools.shuffle_configuration(self.root_hier, max_translation=10)
        mc = IMP.core.MonteCarlo(self.m)
        mc.set_kt(temperature)
        sf = IMP.core.RestraintsScoringFunction(self.rs, "SF")
        mc.set_scoring_function(sf)
        bmvr = [IMP.core.BallMover(self.m, x, 0.5) for x in self.particles]
        mc.add_movers(bmvr)
        # Saving the frames to RMF file
        f = RMF.create_rmf_file("string_of_beads.rmf3")
        IMP.rmf.add_hierarchy(f, self.root_hier)
        IMP.rmf.add_particles(f, self.particles)
        IMP.rmf.add_restraints(f, [self.rs])
        #o = IMP.pmi.output.Output()
        #os = IMP.rmf.SaveOptimizerState(self.m, f)
        #os.update_always("initial conformation")
        #os.set_log_level(IMP.SILENT)
        #os.set_simulator(mc)
        # update the decorator (average) 
        #mc.add_optimizer_state(os)
        mc.optimize(num_steps)

class TwoLevelMCMC:
    def __init__(self, root_hier, temperature, num_steps):
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()

    def _create_global_particle(self):
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
        
        self.global_particle = IMP.Particle(self.m)
        self.global_particle.set_name('global')
        self.global_particle.add_attribute(IMP.FloatKey('mean'), 0)
        self.global_particle.add_attribute(IMP.FloatKey('precision'), 1)
        
