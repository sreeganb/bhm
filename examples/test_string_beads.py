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
import IMP.bhm
import IMP.bhm.restraints.strings

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

class EstimateChi():

    def __init__ (self, root_hier, num_strings, num_beads):
        end_to_end = []
        self.model = root_hier.get_model()
        for i in range(num_strings):
            dfirst = IMP.core.XYZR(self.model, root_hier.get_child(i).get_child(0).get_particle()).get_coordinates()
            dlast = IMP.core.XYZR(self.model, root_hier.get_child(i).get_child(num_beads - 1).get_particle()).get_coordinates()
            end_to_end.append(IMP.algebra.get_distance(dfirst, dlast))
        print("end to end distances: ", end_to_end)
        
        # Save end to end distances to a file
        #with open('./derived_data/end_to_end_distances.txt', 'w') as f:
        #    for distance in end_to_end:
        #        f.write(str(distance) + '\n')

class MCMCsampler():
    # MCMC sampler for the hierarchical system
    def __init__(self, root_hier, temperature, num_steps):
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        #****************************************************************************************
        # IMP.pmi.tools is where most of the extraction of information from the hierarchy is done
        #****************************************************************************************
        self.rs = IMP.pmi.tools.get_restraint_set(self.m)
        num_beads = root_hier.get_child(0).get_number_of_children()
        num_strings = root_hier.get_number_of_children()
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
        EstimateChi(self.root_hier, num_strings, num_beads)
        print("number of accepted MCMC steps: ", mc.get_number_of_accepted_steps())
        
if __name__ == "__main__":
    num_beads = 10
    num_strings = 1
    system = system_of_beads()
    root_hier = system._create_beads(num_strings, num_beads)   
    cbr = IMP.bhm.restraints.strings.ConnectBeadsRestraint(root_hier, 1.0, 4.0, kappa = 10.0, label = "disres")
    cbr.add_to_model()
    #****************************************************************************************
    # Read in the data from the file end_to_end_data.txt
    # pass this on to the class EndToEndRestraint and create the restraint
    #****************************************************************************************
    etedata = np.loadtxt('./derived_data/end_to_end_data.txt')
    etr = IMP.bhm.restraints.strings.EndToEndRestraint(root_hier, etedata, label = "endtoend", weight = 1.0)
    etr.add_to_model()  # add restraint to model  
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(root_hier)
    evr.add_to_model()
    MCMCsampler(root_hier, 3.0, 1000)