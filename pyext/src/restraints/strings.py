"""@namespace IMP.pmi.restraints.basic
End to end restraint for the model problem
"""

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.tools
import IMP.pmi.restraints
import numpy as np

class ConnectBeadsRestraint(IMP.pmi.restraints.RestraintBase):
    # A simple distance restraint between two beads connecting the 
    # string of beads    
    def __init__(self, root_hier,
                 distancemin=0, distancemax=100, kappa=1.0,
                 label=None, weight=1.):
        self.model = root_hier.get_model()
        num_strings = root_hier.get_number_of_children()
        num_beads = root_hier.get_child(0).get_number_of_children()

        super(ConnectBeadsRestraint, self).__init__(self.model, label=label,
            weight=weight)
        print(self.name)

        for j in range(num_strings):
            self.rset = IMP.RestraintSet(self.model, "connectbeads"+str(j))
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

                self.rset.add_restraint(IMP.core.DistanceRestraint(self.model, ts1,
                                           particle1,
                                           particle2))

                self.rset.add_restraint(
                    IMP.core.DistanceRestraint(self.model, ts2,
                                           particle1,
                                           particle2))
            self.rs.add_restraint(self.rset)
    
class EndToEndRestraint(IMP.pmi.restraints.RestraintBase):
    '''Distance restraint for the end-to-end distance of a string of beads
    '''
    def __init__(self, root_hier, etedata, label = None, weight = 1.0):
        '''
        input two particles, read in the experimental end to end distance value, 
        construct a forward model for the data point and add the noise model which
        is a gaussian function with a free parameter sigma where sigma^2 is the variance
        '''
        self.model = root_hier.get_model()
        super().__init__(self.model, label=label, weight=weight)
        print(self.name)
        
        # create nuisance particles
        self.sigma = IMP.pmi.tools.SetupNuisance(
            self.model, 1., 0.01, 100., isoptimized=True).get_particle()
        num_strings = root_hier.get_number_of_children()
        num_beads = root_hier.get_child(0).get_number_of_children()
        sel_tuple = []
        distances = []
        # create a restraint set for each system separately
        for i in range(num_strings):
            self.rset = IMP.RestraintSet(self.model, "endtoend"+str(i))
            self.d1 = root_hier.get_child(i).get_child(0).get_particle() # select the first bead
            self.d2 = root_hier.get_child(i).get_child(num_beads - 1).get_particle() # select the last bead
            pair = (self.d1, self.d2)
            for j in range(len(etedata)):
                distances.append(etedata[j])
                sel_tuple.append(pair)
                self.rset.add_restraint(IMP.bhm.EndtoendRestraint(self.model, sel_tuple[0][0], sel_tuple[0][1], self.sigma, etedata[j]))
            self.rs.add_restraint(self.rset)
        #rs_priors.add_restraint(IMP.isd.JeffreysRestraint(self.model,
        #                                                       self.sigma))
        
    def get_output(self):
        """Get the output of the restraint to be used by the IMP.pmi.output
        object"""
        output = super().get_output()

        output["EndtoEndRestraint between particles_"] = str(self.d1)

        return output

    def get_likelihood(self):
        """Get the unweighted likelihood of the restraint"""
        likelihood = 1
        for restraint in self.rs:
            likelihood *= restraint.get_probability()

        return likelihood