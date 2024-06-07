#!/usr/bin/env python
#----------------------------------------------------------------------
# test_string_beads.py 
# Create a system of string of beads and run hierarchical MCMC sampler
# Use the IMP hierarchy type to represent the system 
# Use a class for each separate restraint 
#----------------------------------------------------------------------

import IMP
import IMP.isd
import IMP.core
import IMP.container
import IMP.atom
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import numpy as np
import RMF
import IMP.rmf
import IMP.display
import IMP.pmi.output

class system_of_beads():
    # system representation:
    # create a system of beads and strings
    def __init__(self):
        self.model = IMP.Model()
#    
    def _create_beads(self, num_strings, num_beads):
        # Top level of the hierarchy is the root_hier
        root_hier = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "root")) # create a hierarchy
        for i in range(num_strings):
            # Second level of the hierarchy is the string
            pc = IMP.Particle(self.model, "string" + str(i))
            pcd = IMP.atom.Hierarchy.setup_particle(pc)            
            root_hier.add_child(pcd) # first child which is the string as a whole
            for j in range(num_beads):
                #dis = np.random.uniform(0, 10)
                p = IMP.Particle(self.model, "bead" + str(i) + str(j))
                IMP.display.Colored.setup_particle(p, IMP.display.get_display_color(i))
                #pd = IMP.atom.Hierarchy.setup_particle(p)
                dr = IMP.core.XYZR.setup_particle(p) # create a decorator for the particle
                dr.set_coordinates(IMP.algebra.get_random_vector_in(IMP.algebra.Sphere3D(IMP.algebra.Vector3D(0, 0, 0), 20.0)))
                if j == 0 or j == num_beads - 1:
                    dr.set_radius(2.0)
                else:
                    dr.set_radius(1.0)
                dr.set_coordinates_are_optimized(True)
                IMP.atom.Mass.setup_particle(p, 1.0)
                pcd.add_child(p) # first grandchild which is a bead
        return root_hier

class ConnectBeadsRestraint(IMP.pmi.restraints.RestraintBase):
    # A simple distance restraint between two beads connecting the 
    # string of beads    
    def __init__(self, root_hier, num_strings, num_beads,
                 distancemin=0, distancemax=100, kappa=1.0,
                 label=None, weight=1.):
        """Setup distance restraint.
        @param root_hier The hierarchy to select from
        @param distancemin The minimum dist
        @param distancemax The maximum dist
        @param kappa The harmonic parameter
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        """
        self.model = root_hier.get_model()
        super(ConnectBeadsRestraint, self).__init__(self.model, label=label,
            weight=weight)

        for j in range(num_strings):
            string = root_hier.get_child(j)
            particles = [string.get_child(i).get_particle() for i in range(num_beads)]

            ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
            ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)

            for i in range(0, num_beads-1):
                particle1 = particles[i]
                particle2 = particles[i + 1]        
                print("Created distance restraint between "
                      "%s and %s" % (particle1.get_name(),
                                     particle2.get_name()))

                self.rs.add_restraint(
                    IMP.core.DistanceRestraint(self.model, ts1,
                                           particle1,
                                           particle2))
                self.rs.add_restraint(
                    IMP.core.DistanceRestraint(self.model, ts2,
                                           particle1,
                                           particle2))
        if (num_strings > 1):
            ts1 = IMP.core.HarmonicUpperBound(distancemax + 2.0, kappa)
            ts2 = IMP.core.HarmonicLowerBound(distancemin + 2.0, kappa)
            for i in range(0, num_beads - 1):
                for j in range(num_strings):
                    particle1 = root_hier.get_child(j).get_child(i).get_particle()
                    particle2 = root_hier.get_child(j).get_child(i+1).get_particle()
                    self.rs.add_restraint(
                        IMP.core.DistanceRestraint(self.model, ts1,
                                               particle1,
                                               particle2))
                    self.rs.add_restraint(
                        IMP.core.DistanceRestraint(self.model, ts2,
                                               particle1,
                                               particle2))
        def get_score(self):
            print(self.rs.get_score())
            return self.rs.get_score()

class MCMCsampler():
    # MCMC sampler for the hierarchical system
    def __init__(self, root_hier, num_strings, num_beads, temperature, num_steps):
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        #****************************************************************************************
        # IMP.pmi.tools is where most of the extraction of information from the hierarchy is done
        #****************************************************************************************
        self.rs = IMP.pmi.tools.get_restraint_set(self.m)
        self.particles = []
        for i in range(num_strings):
            string = root_hier.get_child(i)
            for j in range(num_beads):
                self.particles.append(string.get_child(j).get_particle())
            #self.particles.append([string.get_child(j).get_particle() for j in range(num_beads)])
        # Setup the MCMC parameters
        IMP.pmi.tools.shuffle_configuration(self.root_hier, max_translation=10)
        mc = IMP.core.MonteCarlo(self.m)
        mc.set_log_level(IMP.SILENT)
        mc.set_kt(temperature)
        sf = IMP.core.RestraintsScoringFunction(self.rs, "SF")
        mc.set_scoring_function(sf)
        #print("particles = ", self.particles)
        bmvr = [IMP.core.BallMover(self.m, x, 0.5) for x in self.particles]
        mc.add_movers(bmvr)
        # Saving the frames to RMF file
        f = RMF.create_rmf_file("string_of_beads.rmf3")
        IMP.rmf.add_hierarchy(f, self.root_hier)
        #IMP.rmf.add_particles(f, self.particles)
        #print("restraints set :", self.rs)
        IMP.rmf.add_restraints(f, [self.rs])
        #o = IMP.pmi.output.Output()
        os = IMP.rmf.SaveOptimizerState(self.m, f)
        os.update_always("initial conformation")
        os.set_log_level(IMP.SILENT)
        #os.set_simulator(mc)
        mc.add_optimizer_state(os)
        mc.optimize(num_steps)
        print("number of accepted MCMC steps: ", mc.get_number_of_accepted_steps())
        
if __name__ == "__main__":
    num_beads = 10
    num_strings = 2
    system = system_of_beads()
    root_hier = system._create_beads(num_strings, num_beads)   
    cbr = ConnectBeadsRestraint(root_hier, num_strings, num_beads, 1.0, 4.0, kappa = 10.0, label = "disres")
    cbr.add_to_model()
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(root_hier)
    evr.add_to_model()
    MCMCsampler(root_hier, num_strings, num_beads, 3.0, 1000)
