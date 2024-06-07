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

class system_of_beads():
    # system representation:
    # create a system of beads and strings
    def __init__(self):
        self.model = IMP.Model()
#    
    def _create_beads(self, num_beads):
        ps = [IMP.core.XYZR.setup_particle(IMP.Particle(self.model))
        for i in range(0, num_beads)]
        root_hier = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "root"))
        i = 0
        for p in ps:
            p.set_radius(1)
            p.set_coordinates_are_optimized(True)
            IMP.atom.Mass.setup_particle(p, 30)
            p.set_coordinates(IMP.algebra.get_random_vector_in(IMP.algebra.Sphere3D((0.0,0.0,0.0), 50.0)))
            IMP.display.Colored.setup_particle(p, IMP.display.get_display_color(i))
            hr = IMP.atom.Hierarchy.setup_particle(
                IMP.Particle(self.model, "molecule " + str(i)))
            hr.add_child(IMP.atom.Hierarchy.setup_particle(p))
            root_hier.add_child(hr)
            i = i + 1
        print(root_hier.get_children())
        #pl = []
        #for i in range(0, len(ps) // 2):
        #    pp = (ps[2 * i], ps[2 * i + 1])
        #    pl.append(pp)
        #    # create OK staring position
        #    pp[1].set_coordinates(IMP.algebra.get_random_vector_on(pp[0].get_sphere()))
        #    IMP.display.Colored.setup_particle(pp[0], IMP.display.get_display_color(i))
        #    IMP.display.Colored.setup_particle(pp[1], IMP.display.get_display_color(i))
        #    hr = IMP.atom.Hierarchy.setup_particle(
        #        IMP.Particle(self.model, "molecule " + str(i)))
        #    hr.add_child(IMP.atom.Hierarchy.setup_particle(pp[0]))
        #    hr.add_child(IMP.atom.Hierarchy.setup_particle(pp[1]))
        #    root_hier.add_child(hr)
        #print("particle list: ", pl)
        return root_hier
        # create a system of beads
#        pp = IMP.Particle(self.model) # parent particle
#        tr = IMP.core.HierarchyTraits("String of beads") 
#        root_hier = IMP.core.Hierarchy.setup_particle(pp, tr) # create a hierarchy
#        pp.set_name("root")
#        for i in range(num_strings):
#            pc = IMP.Particle(self.model)
#            pcd = IMP.core.Hierarchy.setup_particle(pc, tr)
#            root_hier.add_child(pcd) # first child which is the string as a whole
#            pc.set_name("string" + str(i))
#            for j in range(num_beads):
#                dis = np.random.uniform(0, 10)
#                p = IMP.Particle(self.model)
#                pd = IMP.core.Hierarchy.setup_particle(p, tr)
#                p.set_name("bead" + str(j))
#                dr = IMP.core.XYZR.setup_particle(p) # create a decorator for the particle
#                dr.set_coordinates(IMP.algebra.Vector3D(dis, -1.0*dis , dis-1.0))
#                if i == 0 or i == num_beads - 1:
#                    dr.set_radius(2.0)
#                else:
#                    dr.set_radius(1.0)
#                dr.set_coordinates_are_optimized(True)
#                IMP.atom.Mass.setup_particle(p, 1.0)
#                pcd.add_child(pd) # first grandchild which is a bead
#        return root_hier

class ConnectBeadsRestraint(IMP.pmi.restraints.RestraintBase):
    # A simple distance restraint between two beads connecting the 
    # string of beads    
    def __init__(self, root_hier, particle1, particle2,
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
        ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
        ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)

        model = root_hier.get_model()
        super(ConnectBeadsRestraint, self).__init__(model, label=label,
                                                weight=weight)        
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

class MCMCsampler():
    # MCMC sampler for the hierarchical system
    def __init__(self, root_hier, particles, temperature, num_steps):
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        #****************************************************************************************
        # IMP.pmi.tools is where most of the extraction of information from the hierarchy is done
        #****************************************************************************************
        self.rs = IMP.pmi.tools.get_restraint_set(self.m)
        self.particles = particles
        # Setup the MCMC parameters
        mc = IMP.core.MonteCarlo(self.m)
        mc.set_log_level(IMP.SILENT)
        mc.set_kt(temperature)
        sf = IMP.core.RestraintsScoringFunction(self.rs, "SF")
        mc.set_scoring_function(sf)
        bmvr = [IMP.core.BallMover(self.m, x, 0.5) for x in self.particles]
        #mv = IMP.core.SerialMover(bmvr)
        mc.add_movers(bmvr)
        #mc.add_movers(bmvr) # add the movers to the MonteCarlo object
        #f = RMF.create_rmf_file("string_of_beads.rmf3")
        #IMP.rmf.add_hierarchy(f, self.root_hier)
        #IMP.rmf.add_restraints(f, [self.rs])
        #IMP.rmf.add_particles(f, self.particles)
        #o = IMP.pmi.output.Output()
        #os = IMP.rmf.SaveOptimizerState(self.m, f)
        #earray = []
        #for i in range(400):
        #os.update_always("initial conformation")
        #os.set_log_level(IMP.SILENT)
        #os.set_simulator(mc)
        #mc.add_optimizer_state(os)
        IMP.pmi.tools.shuffle_configuration(self.root_hier, max_translation=5)
        mc.optimize(num_steps)
        #print("number of accepted MCMC steps: ", mc.get_number_of_accepted_steps())
        #eted = self.define_global_parameter()
        #earray.append(eted)
        
if __name__ == "__main__":
    num_beads = 10
    num_strings = 1
    system = system_of_beads()
    for i in range(num_strings):
        root_hier = system._create_beads(num_beads)
    particles = [root_hier.get_child(i).get_particle() for i in range(num_beads)]
    #------------------------------------------------------------------
    #Some of the allowed usage of root_hierarchy type
    #------------------------------------------------------------------
    for i in range(num_beads - 1):
        n1 = particles[i]
        n2 = particles[i + 1]
        cbr = ConnectBeadsRestraint(root_hier, n1, n2, 1.0, 4.0, kappa = 10.0, label = "disres")
        cbr.add_to_model()
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(root_hier)
    evr.add_to_model()
#    rs = []
#    rs.append(cbr.get_restraint())
#    rs.append(evr.get_restraint())
    MCMCsampler(root_hier, particles, 1.0, 1000)
