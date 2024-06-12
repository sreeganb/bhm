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
#import IMP.bhm.restraint

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
    def __init__(self, root_hier,
                 distancemin=0, distancemax=100, kappa=1.0,
                 label=None, weight=1.):
        self.model = root_hier.get_model()
        num_beads = root_hier.get_child(0).get_number_of_children()
        num_strings = root_hier.get_number_of_children()
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

class EndToEndRestraint(object):
    '''Distance restraint for the end-to-end distance of a string of beads
    '''
    def __init__(self, root_hier, etedata):
        '''
        input two particles, read in the experimental end to end distance value, 
        construct a forward model for the data point and add the noise model which
        is a gaussian function with a free parameter sigma where sigma^2 is the variance
        '''
        self.model = root_hier.get_model()
        IMP.Restraint.__init__(self, self.model, "EndToEndRestraint %1%")
        # create restraint sets to join restraints
        self.rs = IMP.RestraintSet(self.model, "likelihood")
        self.rs_priors = IMP.RestraintSet(self.model, "priors")

        # create nuisance particles
        self.sigma = IMP.pmi.tools.SetupNuisance(
            self.model, 1., 0.01, 100., isoptimized=True).get_particle()

        num_beads = root_hier.get_child(0).get_number_of_children()
        num_strings = root_hier.get_number_of_children()
        sel_tuple = []
        distances = []
        for i in range(num_strings):
            self.d1 = root_hier.get_child(i).get_child(0).get_particle() # select the first bead
            self.d2 = root_hier.get_child(i).get_child(num_beads - 1).get_particle() # select the last bead
            pair = (self.d1, self.d2)
            for j in range(len(etedata)):
                distances.append(etedata[j])
                sel_tuple.append(pair)
        print("the tuple list is: ", sel_tuple)
        rend = IMP.bhm.EndtoendRestraint(self.model, sel_tuple, self.sigma, distances)

    def add_to_model(self):
        """Add the restraints to the model."""
        for rs in [self.rs, self.rs_priors]:
            IMP.pmi.tools.add_restraint_to_model(self.model, rs)

    def get_restraint(self):
        return self.rs

    def get_restraint_for_rmf(self):
        """Get the restraint for visualization in an RMF file."""
        return self.rs

    def get_particles_to_sample(self):
        """Get any created particles which should be sampled."""
        out = {}
        out["Nuisances_Sigma"] = ([self.sigma], .1)
        return out

    def get_output(self):
        """Get outputs to write to stat files."""
        output = {}
        self.m.update()
        likelihood_score = self.rs.unprotected_evaluate(None)
        prior_score = self.rs_priors.unprotected_evaluate(None)
        output["_TotalScore"] = likelihood_score + prior_score
        output["EndtoendLikelihood_Score"] = likelihood_score
        output["EndtoendPrior_Score"] = prior_score
        output["Endtoend_Sigma"] = self.sigma.get_scale()
        return output

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
        with open('./end_to_end_distances.txt', 'w') as f:
            for distance in end_to_end:
                f.write(str(distance) + '\n')

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
        IMP.rmf.add_restraints(f, [self.rs])
        #o = IMP.pmi.output.Output()
        os = IMP.rmf.SaveOptimizerState(self.m, f)
        os.update_always("initial conformation")
        os.set_log_level(IMP.SILENT)
        #os.set_simulator(mc)
        # update the decorator (average) 
        mc.add_optimizer_state(os)
        mc.optimize(num_steps)
        EstimateChi(self.root_hier, num_strings, num_beads)
        print("number of accepted MCMC steps: ", mc.get_number_of_accepted_steps())
        
if __name__ == "__main__":
    num_beads = 10
    num_strings = 1
    system = system_of_beads()
    root_hier = system._create_beads(num_strings, num_beads)   
    cbr = ConnectBeadsRestraint(root_hier, 1.0, 4.0, kappa = 10.0, label = "disres")
    cbr.add_to_model()
    #****************************************************************************************
    # Read in the data from the file end_to_end_data.txt
    # pass this on to the class EndToEndRestraint and create the restraint
    #****************************************************************************************
    etedata = np.loadtxt('end_to_end_data.txt')
    etr = EndToEndRestraint(root_hier, etedata)
    
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(root_hier)
    evr.add_to_model()
    MCMCsampler(root_hier, num_strings, num_beads, 3.0, 1000)
