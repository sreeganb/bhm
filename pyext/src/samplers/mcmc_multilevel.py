#!/usr/bin/env python
"""
MCMC protocol for the Bayesian hierarchical modeling.
This involves two stages of parameter estimation.
The first stage is a simple MCMC step, whereas the 
second stage involves the joint posterior estimation of the 
local and global parameters. Use the core IMP library for 
most of the MCMC steps.
"""

import IMP
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.restraints
import RMF
import IMP.rmf
import IMP.atom
import os

class RMFRestraints:
    """All restraints that are written out to the RMF file."""
    def __init__(self, model, user_restraints):
        print("Initializing RMFRestraints")
        self._rmf_rs = IMP.pmi.tools.get_restraint_set(model, rmf=True)
        self._user_restraints = user_restraints if user_restraints else []
        print(f"Initialized RMFRestraints with {len(self._user_restraints)} user restraints")

    def __len__(self):
        return len(self._user_restraints) + self._rmf_rs.get_number_of_restraints()

    def __bool__(self):
        return len(self) > 0

    def __getitem__(self, i):
        class FakePMIWrapper:
            def __init__(self, r):
                self.r = IMP.RestraintSet.get_from(r)

            def get_restraint(self):
                return self.r

        len_user = len(self._user_restraints)
        if 0 <= i < len_user:
            return self._user_restraints[i]
        elif 0 <= i - len_user < self._rmf_rs.get_number_of_restraints():
            r = self._rmf_rs.get_restraint(i - len_user)
            return FakePMIWrapper(r)
        else:
            raise IndexError("Out of range")

class MCMCSampler:
    """MCMC sampler for the hierarchical system."""
    def __init__(self, root_hier, dof, temperature, num_steps,  
                 output_objects, rmf_output_objects, stat_file_name_suffix, 
                 file_name="bhm0",
                 global_output_directory="output_dir", 
                 initial_rmf_name_suffix="initial",
                 nframes_write_coordinates=1, 
                 do_clean_first=True, 
                 do_create_directories=True, 
                 atomistic=False, 
                 test_mode=False, 
                 use_nestor=False):
        print("Initializing MCMCSampler")
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        self.output_objects = output_objects
        self.vars = {
            "global_output_directory": global_output_directory,
            "do_clean_first": do_clean_first,
            "do_create_directories": do_create_directories,
            "initial_rmf_name_suffix": initial_rmf_name_suffix,
            "number_of_mcmc_steps": num_steps,
            "nframes_write_coordinates": nframes_write_coordinates,
            "atomistic": atomistic,
            "production_rmf_file_name": file_name,
        }
        self.vars["geometries"] = None
        self._rmf_restraints = RMFRestraints(self.m, None)
        self.rmf_output_objects = rmf_output_objects
        self.nest = use_nestor
        self.test_mode = test_mode
        self.samplers = []
        self.sampler_mc = None
        self.stat_file_name_suffix = stat_file_name_suffix
        self.file_name = file_name

        self.run_mcmc(dof, temperature)
        print("MCMCSampler initialized")
        
    def add_geometries(self, geometries):
        if self.vars["geometries"] is None:
            self.vars["geometries"] = list(geometries)
        else:
            self.vars["geometries"].extend(geometries)

    def setup_output(self):
        """Setup the output directories and stat file."""
        globaldir = self.vars["global_output_directory"] + "/"
        if not self.test_mode and not self.nest:
            if self.vars["do_clean_first"]:
                print("Cleaning output directory")
                pass
            if self.vars["do_create_directories"]:
                print(f"Creating output directory: {globaldir}")
                os.makedirs(globaldir, exist_ok=True)

        stat_file = globaldir + self.stat_file_name_suffix + "_bhm.out"
        self.vars["stat_file"] = stat_file
        self.m.update()  # Update the model before writing the stat file
        output = IMP.pmi.output.Output(atomistic=self.vars["atomistic"])
        if isinstance(self.root_hier, IMP.atom.Hierarchy) and not self.root_hier.get_parent():
            if self.output_objects is not None:
                print("Appending TotalScoreOutput to output_objects")
                self.output_objects.append(IMP.pmi.io.TotalScoreOutput(self.m))
        if self.output_objects is not None:
            print(f"Initializing stat file: {stat_file}")
            output.init_stat2(stat_file, self.output_objects, extralabels=["rmf_file", "rmf_frame_index"])

    def run_mcmc(self, dof, temperature):
        """Setup the RMF files and run the MCMC simulation."""
        print("Setting up RMF files")
        output_hierarchies = [self.root_hier]
        output = IMP.pmi.output.Output(atomistic=self.vars["atomistic"])
        globaldir = self.vars["global_output_directory"] + "/"
        init_suffix = globaldir + self.vars["initial_rmf_name_suffix"]

        print("Setting up and writing initial RMF coordinate file")
        output.init_rmf(init_suffix + ".rmf3", 
                        output_hierarchies,
                        geometries=self.vars["geometries"], 
                        listofobjects=self.rmf_output_objects)
        if self._rmf_restraints:
            print("Adding restraints to initial RMF file")
            output.add_restraints_to_rmf(init_suffix + ".rmf3", self._rmf_restraints)
        output.write_rmf(init_suffix + ".rmf3")
        output.close_rmf(init_suffix + ".rmf3")

        prod_suffix = globaldir + self.vars["production_rmf_file_name"]        
        print("Setting up production RMF coordinate file")
        rmfname = prod_suffix + ".rmf3"
        output.init_rmf(rmfname, 
                        output_hierarchies, 
                        geometries=self.vars["geometries"],
                        listofobjects=self.rmf_output_objects)
        if self._rmf_restraints:
            print("Adding restraints to production RMF file")
            output.add_restraints_to_rmf(rmfname, self._rmf_restraints)

        output.write_rmf(rmfname)

        self.vars["rmfname"] = rmfname  # Store the production RMF file name in vars

        print("Running setup_output()")
        self.setup_output()

        print("movers are: ", dof.get_movers())
        self.sampler_mc = IMP.pmi.samplers.MonteCarlo(self.m, dof.get_movers(), temperature)

        if self.output_objects is not None:
            print("Appending sampler_mc to output_objects")
            self.output_objects.append(self.sampler_mc)

        if self.rmf_output_objects is not None:
            print("Appending sampler_mc to rmf_output_objects")
            self.rmf_output_objects.append(self.sampler_mc)
        self.samplers.append(self.sampler_mc)

        nframes = 100  # Number of frames to run MCMC
        print(f"RMF file name: {rmfname}")
        #print("rmf dictionary", output.dictionary_rmfs)

        stat_file = self.vars["stat_file"]

        print(f"Starting MCMC loop for {nframes} frames")
        for i in range(nframes):
            self.sampler_mc.optimize(self.vars["number_of_mcmc_steps"])
            score = IMP.pmi.tools.get_restraint_set(self.m).evaluate(False)
            print(f"Frame {i}: score is {score}")
            output.set_output_entry("score", score)
            if i % self.vars["nframes_write_coordinates"] == 0:
                print(f"Frame {i}: Writing coordinates to {rmfname}")
                if rmfname not in output.dictionary_rmfs:
                    print(f"ERROR: RMF file {rmfname} is NOT registered in dictionary_rmfs")
                print("elements of the dictionary:", output.dictionary_rmfs[rmfname][3])
                try:
                    output.write_rmf(rmfname)
                except KeyError as e:
                    print(f"KeyError: {e}")
                    print("Available keys:", output.dictionary_rmfs[rmfname][2])
                output.set_output_entry("rmf_file", rmfname)
            else:
                output.set_output_entry("rmf_file", rmfname)
                output.set_output_entry("rmf_frame_index", '-1')
            if self.output_objects is not None:
                print(f"Frame {i}: Writing stat file {stat_file}")
                output.write_stat2(stat_file)
        output.close_rmf(rmfname)

class TwoLevelMCMC:
    """Class for two-level MCMC sampling."""
    def __init__(self, root_hier, dof, temperature, num_steps):
        print("Initializing TwoLevelMCMC")
        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        print("TwoLevelMCMC initialized")

    def _create_global_particle(self):
        """Create a nuisance particle to denote the end-to-end distance of the two strings."""
        print("Creating global particle")
        self.chi_is_sampled = True
        chiminnuis = 1.0
        chimaxnuis = 30.0
        chiinit = 15.0
        chimin = 0.01
        chimax = 100.0
        chi = IMP.pmi.tools.SetupNuisance(self.model, chiinit, chiminnuis, chimaxnuis, self.chi_is_sampled).get_particle()
        self.rs.add_restraint(IMP.isd.UniformPrior(self.model, chi, 10000.0, chimax, chimin))

        self.global_particle = IMP.Particle(self.m)
        self.global_particle.set_name('global')
        self.global_particle.add_attribute(IMP.FloatKey('mean'), 0)
        self.global_particle.add_attribute(IMP.FloatKey('precision'), 1)
        print("Global particle created")