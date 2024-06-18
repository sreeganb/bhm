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
import IMP.pmi.macros
import IMP.pmi.dof
import IMP.rmf
import IMP.display
import IMP.pmi.output
import IMP.bhm
import IMP.bhm.restraints.strings
import IMP.bhm.samplers.two_level_mcmc
import IMP.bhm.system_representation.build

#class EstimateChi():

#    def __init__ (self, root_hier, num_strings, num_beads):
#        end_to_end = []
#        self.model = root_hier.get_model()
#        for i in range(num_strings):
#            dfirst = IMP.core.XYZR(self.model, root_hier.get_child(i).get_child(0).get_particle()).get_coordinates()
#            dlast = IMP.core.XYZR(self.model, root_hier.get_child(i).get_child(num_beads - 1).get_particle()).get_coordinates()
#            end_to_end.append(IMP.algebra.get_distance(dfirst, dlast))
#        print("end to end distances: ", end_to_end)
        
        # Save end to end distances to a file
        #with open('./derived_data/end_to_end_distances.txt', 'w') as f:
        #    for distance in end_to_end:
        #        f.write(str(distance) + '\n')
        
if __name__ == "__main__":
    num_systems = 2
    num_strings = np.array([1, 2])
    #num_strings = np.array([1])
    num_beads = 10
    
    system_rep = IMP.bhm.system_representation.build.model()
    root_hier, build_sys = system_rep._create_beads(num_systems, num_strings, num_beads)
    #print(root_hier[1].get_child(1).get_children())
    etedata = np.loadtxt('./derived_data/end_to_end_data.txt')
    for i in range(num_systems):
        cbr = IMP.bhm.restraints.strings.ConnectBeadsRestraint(root_hier[i], 1.0, 4.0, kappa = 10.0, label = "disres")
        cbr.add_to_model()  # add restraint to model
        etr = IMP.bhm.restraints.strings.EndToEndRestraint(root_hier[i], etedata, label = "endtoend", weight = 1.0)
        etr.add_to_model()  # add restraint to model
        evr = "evr_" + str(i)
        evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(root_hier)
        evr.add_to_model()
    print("degrees of freedom: ", build_sys[0].execute_macro()[1])
    IMP.bhm.samplers.two_level_mcmc.MCMCsampler(root_hier[1], 3.0, 100)