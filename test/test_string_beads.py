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
import numpy as np

class system_of_beads():
    def __init__(self):
        self.model = IMP.Model()
#    
    def _create_beads(self, num_beads, num_strings):
        # create a system of beads
        pp = IMP.Particle(self.model) # parent particle 
        root_hier = IMP.core.Hierarchy.setup_particle(pp)
        pp.set_name("root")
        for i in range(num_strings):
            pc = IMP.Particle(self.model)
            pcd = IMP.core.Hierarchy.setup_particle(pc)
            root_hier.add_child(pcd) # first child which is the string as a whole
            pc.set_name("string" + str(i))
            for j in range(num_beads):
                dis = np.random.uniform(0, 10)
                p = IMP.Particle(self.model)
                pd = IMP.core.Hierarchy.setup_particle(p)
                p.set_name("bead" + str(j))
                dr = IMP.core.XYZR.setup_particle(p) # create a decorator for the particle
                dr.set_coordinates(IMP.algebra.Vector3D(dis, -1.0*dis , dis-1.0))
                if i == 0 or i == num_beads - 1:
                    dr.set_radius(2.0)
                else:
                    dr.set_radius(1.0)
                dr.set_coordinates_are_optimized(True)
                IMP.atom.Mass.setup_particle(p, 1.0)
                pcd.add_child(pd) # first grandchild which is a bead
        return root_hier

class ConnectBeadsRestraint(IMP.pmi.restraints.RestraintBase):
        
    def __init__(self, root_hier, tuple_selection1, tuple_selection2,
                 distancemin=0, distancemax=100, resolution=1.0, kappa=1.0,
                 label=None, weight=1.):
        """Setup distance restraint.
        @param root_hier The hierarchy to select from
        @param tuple_selection1 (resnum, resnum, molecule name, copy
               number (=0))
        @param tuple_selection2 (resnum, resnum, molecule name, copy
               number (=0))
        @param distancemin The minimum dist
        @param distancemax The maximum dist
        @param resolution For selecting particles
        @param kappa The harmonic parameter
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number.
        """
        ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
        ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)

        model = root_hier.get_model()
        copy_num1 = 0
        if len(tuple_selection1) > 3:
            copy_num1 = tuple_selection1[3]
        copy_num2 = 0
        if len(tuple_selection2) > 3:
            copy_num2 = tuple_selection2[3]

        sel1 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection1[2],
                                  residue_index=tuple_selection1[0],
                                  copy_index=copy_num1)
        particles1 = sel1.get_selected_particles()
        sel2 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection2[2],
                                  residue_index=tuple_selection2[0],
                                  copy_index=copy_num2)
        particles2 = sel2.get_selected_particles()

        super(ConnectBeadsRestraint, self).__init__(model, label=label,
                                                weight=weight)
        print(self.name)

        print("Created distance restraint between "
              "%s and %s" % (particles1[0].get_name(),
                             particles2[0].get_name()))

        if len(particles1) > 1 or len(particles2) > 1:
            raise ValueError("more than one particle selected")

        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts1,
                                       particles1[0],
                                       particles2[0]))
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts2,
                                       particles1[0],
                                       particles2[0]))

if __name__ == "__main__":
    num_beads = 10
    num_strings = 2
    system = system_of_beads()
    root_hier = system._create_beads(num_beads, num_strings)
    rest = ConnectBeadsRestraint()
    