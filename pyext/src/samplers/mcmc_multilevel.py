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
import IMP.atom

class _RMFRestraints:
    """All restraints that are written out to the RMF file"""
    def __init__(self, model, user_restraints):
        self._rmf_rs = IMP.pmi.tools.get_restraint_set(model, rmf=True)
        self._user_restraints = user_restraints if user_restraints else []

    def __len__(self):
        return (len(self._user_restraints)
                + self._rmf_rs.get_number_of_restraints())

    def __bool__(self):
        return len(self) > 0

    def __getitem__(self, i):
        class FakePMIWrapper:
            def __init__(self, r):
                self.r = IMP.RestraintSet.get_from(r)

            def get_restraint(self):
                return self.r

        lenuser = len(self._user_restraints)
        if 0 <= i < lenuser:
            return self._user_restraints[i]
        elif 0 <= i - lenuser < self._rmf_rs.get_number_of_restraints():
            r = self._rmf_rs.get_restraint(i - lenuser)
            return FakePMIWrapper(r)
        else:
            raise IndexError("Out of range")

class MCMCsampler():
    # MCMC sampler for the hierarchical system
    def __init__(self, root_hier, dof, temperature, num_steps, file_name, 
                 output_objects,
                 rmf_output_objects, 
                 stat_file_name_suffix, global_output_directory,
                 initial_rmf_name_suffix="initial",
                 nframes_write_coordinates=1
                 ):
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        self.output_objects = output_objects
        self.vars = {}
        self.vars["global_output_directory"] = global_output_directory
        self.vars["initial_rmf_name_suffix"] = initial_rmf_name_suffix
        self.vars["number_of_mcmc_steps"] = num_steps
        self.vars["nframes_write_coordinates"] = nframes_write_coordinates
        self._rmf_restraints = _RMFRestraints(self.m, None)
        self.rmf_output_objects = rmf_output_objects
        samplers = []
        sampler_mc = None
        #****************************************************************************************
        # IMP.pmi.tools is where most of the extraction of information from the hierarchy is done
        # First try to just add movers for the first system and check if the MCMC will work
        #****************************************************************************************
        self.rs = IMP.pmi.tools.get_restraint_set(self.m)
        f = RMF.create_rmf_file(file_name)
        IMP.rmf.add_hierarchy(f, self.root_hier)
        IMP.rmf.add_restraints(f, [self.rs])
        #--------------------------------------------------
        # Setup the stat file and print the output onto it
        #--------------------------------------------------
        atomistic = False
        self.vars["atomistic"] = atomistic
        output = IMP.pmi.output.Output(atomistic=self.vars["atomistic"])
        if (isinstance(root_hier, IMP.atom.Hierarchy)
                and not root_hier.get_parent()):
            if self.output_objects is not None:
                print("isinstance is true man")
                self.output_objects.append(
                    IMP.pmi.io.TotalScoreOutput(self.m))
        # Set up the stat file
        print("Setting up stat file")
        globaldir = self.vars["global_output_directory"] + "/"
        stat_file = globaldir + stat_file_name_suffix + "_bhm.out"
        self.vars["stat_file"] = stat_file
        self.m.update() # update the model before writing the stat file
        if self.output_objects is not None:
                output.init_stat2(stat_file,
                                  self.output_objects,
                                  extralabels=["rmf_file", "rmf_frame_index"])
        #--------------------------------------------------
        # setup the RMF files for saving coordinates
        #--------------------------------------------------
        output_hierarchies = [self.root_hier]
        print("Setting up and writing initial rmf coordinate file")
        init_suffix = globaldir + self.vars["initial_rmf_name_suffix"]
        output.init_rmf(init_suffix + ".rmf3",
                            output_hierarchies,
                            listofobjects=self.rmf_output_objects)
        if self._rmf_restraints:
            output.add_restraints_to_rmf(
                init_suffix + ".rmf3",
                self._rmf_restraints)
        output.write_rmf(init_suffix  + ".rmf3")
        output.close_rmf(init_suffix + ".rmf3")

        print("Setting up production rmf files")
        rmfname = globaldir +"monomer"+ ".rmf3"
        output.init_rmf(rmfname, output_hierarchies,
                        listofobjects=self.rmf_output_objects)

        if self._rmf_restraints:
            output.add_restraints_to_rmf(rmfname, self._rmf_restraints)
        #----------------------------------------------------------------------
        # Using the IMP.pmi.samplers.MonteCarlo class to perform the MCMC steps
        #----------------------------------------------------------------------
        sampler_mc = IMP.pmi.samplers.MonteCarlo(self.m, dof.get_movers(), temperature)
        if self.output_objects is not None:
                self.output_objects.append(sampler_mc)
                if self.rmf_output_objects is not None:
                    self.rmf_output_objects.append(sampler_mc)
                samplers.append(sampler_mc)
        #----------------------------------------------------------------------
        # Run the MCMC for nframes number of steps
        #----------------------------------------------------------------------
        nframes = 100 # number of frames to run MCMC
        for i in range(nframes):
            sampler_mc.optimize(self.vars["number_of_mcmc_steps"])
            score = IMP.pmi.tools.get_restraint_set(
                    self.m).evaluate(False)
            output.set_output_entry("score", score)
            if i % self.vars["nframes_write_coordinates"] == 0:
                print('--- writing coordinates')
                output.write_rmf(rmfname)
                output.set_output_entry("rmf_file", rmfname)
                #output.set_output_entry("rmf_frame_index", ntimes_at_low_temp)
            else:
                output.set_output_entry("rmf_file", rmfname)
                output.set_output_entry("rmf_frame_index", '-1')
            if self.output_objects is not None:
                output.write_stat2(stat_file)
        #-----------------------------------------------------------------------
        #mc.set_kt(temperature)
        #sf = IMP.core.RestraintsScoringFunction(self.rs, "SF")
        #mc.set_scoring_function([self.rs])
        #o = IMP.pmi.output.Output()
        #os = IMP.rmf.SaveOptimizerState(self.m, f)
        #os.update_always("initial conformation")
        #os.set_log_level(IMP.SILENT)
        #sampler_mc.optimize(num_steps)
        #print("frame number : ", sampler_mc.get_frame_number())
        #sampler_mc.get_output()   

#        mc = IMP.core.MonteCarlo(self.m)
#        mc.set_kt(temperature)
#        sf = IMP.core.RestraintsScoringFunction(self.rs, "SF")
#        mc.set_scoring_function(sf)
#        sm = dof.get_movers()
#        IMP.set_log_level(IMP.SILENT)
#        mc.add_movers(sm)
#        # Saving the frames to RMF file
#        o = IMP.pmi.output.Output()
#        os = IMP.rmf.SaveOptimizerState(self.m, f)
#        os.update_always("initial conformation")
#        os.set_log_level(IMP.SILENT)
#        mc.add_optimizer_state(os)
#        mc.optimize(num_steps)
#        print("number of accepted steps: ", mc.get_number_of_accepted_steps())


class TwoLevelMCMC:
    def __init__(self, root_hier, dof, temperature, num_steps):
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
        
