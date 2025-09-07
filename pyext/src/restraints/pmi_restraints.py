#!/usr/bin/env python

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
from math import log

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
    _include_in_rmf = True
    def __init__(self, root_hier, etedata, label = None, weight = 1.0):
        '''
        input two particles, read in the experimental end to end distance value, 
        construct a forward model for the data point and add the noise model which
        is a gaussian function with a free parameter sigma where sigma^2 is the variance
        '''
        self.model = root_hier.get_model()
        super(EndToEndRestraint, self).__init__(self.model, label=label,
            weight=weight)
        print(self.name)
        self.sigma_dictionary = {}
        
        # create nuisance particles and add a uniform prior to it
        self.sigma_is_sampled = True
        sigmaname = "SIGMA"
        self.rssig = self._create_restraint_set("PriorSig")
        self.sigma = self.create_sigma(sigmaname)
        
        # get all parameters from the root hierarchy
        num_strings = root_hier.get_child(0).get_number_of_children()
        sel_tuple = []
        distances = []
        # create a restraint set for each system separately
        for i in range(num_strings):
            #self.rset = IMP.RestraintSet(self.model, "endtoend"+str(i))
            # get the first and last bead of the string
            self.d1 = root_hier.get_child(0).get_child(i).get_child(0).get_child(0).get_child(0).get_particle() # select the first bead
            self.d2 = root_hier.get_child(0).get_child(i).get_child(0).get_child(0).get_children()[-1].get_particle() # select the last bead
            pair = (self.d1, self.d2)
            for j in range(len(etedata)):
                distances.append(etedata[j])
                sel_tuple.append(pair)
                #self.rset.add_restraint(IMP.bhm.EndtoendRestraint(self.model, sel_tuple[0][0], sel_tuple[0][1], self.sigma, etedata[j]))
                self.rs.add_restraint(IMP.bhm.EndtoendRestraint(self.model, sel_tuple[0][0], sel_tuple[0][1], self.sigma, etedata[j]))
            #self.rs.add_restraint(self.rset)
        
    def create_sigma(self, name):
        """ This is called internally. Creates a nuisance
        on the structural uncertainty """
        if name in self.sigma_dictionary:
            return self.sigma_dictionary[name][0]

        sigmainit = 5.0
        sigmaminnuis = 3.0
        sigmamaxnuis = 30.0
        sigmamin = 3.0
        sigmamax = 30.0
        sigmatrans = 0.5
        sigma = IMP.pmi.tools.SetupNuisance(
            self.model, sigmainit, sigmaminnuis, sigmamaxnuis,
            self.sigma_is_sampled, name=name).get_particle()
        self.sigma_dictionary[name] = (
            sigma,
            sigmatrans,
            self.sigma_is_sampled)
        self.rssig.add_restraint(
            IMP.isd.UniformPrior(
                self.model,
                sigma,
                100.0,
                sigmamax,
                sigmamin))
        return sigma
    
    def get_restraint_for_rmf(self):
        """ get the dummy restraints to be displayed in the rmf file """
        return self.rs, self.rssig
    
    def get_restraint_sets(self):
        """ get the restraint set """
        return self.rs

    def get_output(self):
        """Get the output of the restraint to be used by the IMP.pmi.output
           object"""
        output = super().get_output()
        output["EndtoEndRestraint_Score_"
                ] = str(-log(self.rs.unprotected_evaluate(None)))
        dis0 = IMP.core.XYZ(self.d1)
        dis1 = IMP.core.XYZ(self.d2)
        output["EndtoEndRestraint_Distance_"
                ] = str(IMP.core.get_distance(dis0, dis1))
        
        for sigmaname in self.sigma_dictionary:
            output["EndtoEndRestraint_Sigma_" +
                   str(sigmaname) + self._label_suffix] = str(
                    self.sigma_dictionary[sigmaname][0].get_scale())
        
        return output

    def get_likelihood(self):
        """Get the unweighted likelihood of the restraint"""
        likelihood = 1
        for restraint in self.rs:
            likelihood *= restraint.get_probability()

        return likelihood
    
    def get_movers(self):
        """ Get all need data to construct a mover in IMP.pmi.dof class"""
        movers = []
        if self.sigma_is_sampled:
            for sigmaname in self.sigma_dictionary:
                mover_name = \
                    "Nuisances_EndtoEndRestraint_Sigma_" \
                    + str(sigmaname) + "_" + self.label
                particle = self.sigma_dictionary[sigmaname][0]
                maxstep = (self.sigma_dictionary[sigmaname][1])
                mv = IMP.core.NormalMover(
                    [particle], IMP.FloatKeys([IMP.FloatKey("nuisance")]),
                    maxstep)
                mv.set_name(mover_name)
                movers.append(mv)

        return movers