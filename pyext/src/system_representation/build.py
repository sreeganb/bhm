#!/usr/bin/env python
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
import IMP.rmf
import IMP.display
import IMP.pmi.output
import IMP.bhm
import IMP.pmi.macros
import IMP.pmi.dof

class model():
    # system representation:
    # create a system of beads and strings
    def __init__(self):
        self.model = IMP.Model()
    def _create_beads(self, num_systems, num_strings, num_beads):
        """
        num_systems: number of systems (to treat the monomer and dimer simultaneously)
        num_strings: number of strings, array type with entries corresponding to the number of strings (2 for a dimer, 1 for a monomer)
        num_beads: number of beads in each string
        """
        # Top level of the hierarchy is the root_hier
        hier = []
        bs = []
        for k in range(num_systems):
            h = "root_hier" + str(k)
            h = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "root"+str(k))) # create a hierarchy
            bsys = IMP.pmi.macros.BuildSystem(self.model)
            for i in range(int(num_strings[k])):
                # Second level of the hierarchy is the string
                pc = IMP.Particle(self.model, "string" + str(i))
                pcd = IMP.atom.Hierarchy.setup_particle(pc)            
                h.add_child(pcd)
                for j in range(num_beads):
                    p = IMP.Particle(self.model, "bead" + str(i) + str(j))
                    IMP.display.Colored.setup_particle(p, IMP.display.get_display_color(i))
                    dr = IMP.core.XYZR.setup_particle(p)
                    dr.set_coordinates(IMP.algebra.get_random_vector_in(IMP.algebra.Sphere3D(IMP.algebra.Vector3D(0, 0, 0), 20.0)))
                    if j == 0 or j == num_beads - 1:
                        dr.set_radius(2.0)
                    else:
                        dr.set_radius(1.0)
                    dr.set_coordinates_are_optimized(True)
                    IMP.atom.Mass.setup_particle(p, 1.0)
                    pcd.add_child(p)
            #temp1, temp2 = bsys.execute_macro()
            #print(temp2)
            bsys.add_state(h)
            hier, dof = bsys.execute_macro()
            print("degrees of freedom: ", dof)
            print("leaves: ", IMP.atom.get_leaves(h))
            bs.append(bsys)
            hier.append(h)

        return hier, bs