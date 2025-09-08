import IMP
import IMP.core
import IMP.atom
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.tools
import os
import math

class SimpleReplicaExchange:
    """
    Replica exchange for simple spherical particles.
    Closely follows the original ReplicaExchange class structure but without protein-specific features.
    """
    
    def __init__(self, 
                 model, 
                 root_hier,
                 monte_carlo_sample_objects,
                 output_objects=None,
                 rmf_output_objects=None,
                 monte_carlo_temperature=1.0,
                 replica_exchange_minimum_temperature=1.0,
                 replica_exchange_maximum_temperature=2.5,
                 replica_exchange_swap=True,
                 num_sample_rounds=1,
                 monte_carlo_steps=10,
                 number_of_frames=1000,
                 nframes_write_coordinates=1,
                 global_output_directory="./",
                 rmf_dir="rmfs/",
                 score_moved=False):
        """
        Constructor following original ReplicaExchange pattern.
        
        @param model The IMP model
        @param root_hier Root hierarchy containing particles
        @param monte_carlo_sample_objects Movers for MC sampling
        @param output_objects Objects to include in output
        @param rmf_output_objects Objects to include in RMF
        @param monte_carlo_temperature MC temperature
        @param replica_exchange_minimum_temperature Min temp for REX
        @param replica_exchange_maximum_temperature Max temp for REX
        @param replica_exchange_swap Enable temperature swap
        @param num_sample_rounds Number of MC rounds per cycle
        @param monte_carlo_steps Number of MC steps per round
        @param number_of_frames Number of frames to run
        @param nframes_write_coordinates How often to write coordinates
        @param global_output_directory Output directory
        @param rmf_dir RMF subdirectory
        @param score_moved Cache scoring for unmoved particles
        """
        
        self.model = model
        self.root_hier = root_hier
        self.monte_carlo_sample_objects = monte_carlo_sample_objects
        self.score_moved = score_moved
        
        # Output objects
        if output_objects is None:
            self.output_objects = []
        else:
            self.output_objects = output_objects
            
        if rmf_output_objects is None:
            self.rmf_output_objects = []
        else:
            self.rmf_output_objects = rmf_output_objects
        
        # Add total score output
        self.output_objects.append(IMP.pmi.io.TotalScoreOutput(self.model))
        self.rmf_output_objects.append(IMP.pmi.io.TotalScoreOutput(self.model))
        
        # Parameters
        self.vars = {}
        self.vars["monte_carlo_temperature"] = monte_carlo_temperature
        self.vars["replica_exchange_minimum_temperature"] = replica_exchange_minimum_temperature
        self.vars["replica_exchange_maximum_temperature"] = replica_exchange_maximum_temperature
        self.vars["replica_exchange_swap"] = replica_exchange_swap
        self.vars["num_sample_rounds"] = num_sample_rounds
        self.vars["monte_carlo_steps"] = monte_carlo_steps
        self.vars["number_of_frames"] = number_of_frames
        self.vars["nframes_write_coordinates"] = nframes_write_coordinates
        self.vars["global_output_directory"] = global_output_directory
        self.vars["rmf_dir"] = rmf_dir
        
    def execute_macro(self):
        """Execute the replica exchange simulation following original pattern"""
        
        temp_index_factor = 100000.0
        
        # Setup Monte Carlo sampler
        print("Setting up MonteCarlo")
        sampler_mc = IMP.pmi.samplers.MonteCarlo(
            self.model, 
            self.monte_carlo_sample_objects,
            self.vars["monte_carlo_temperature"],
            score_moved=self.score_moved
        )
        
        self.output_objects.append(sampler_mc)
        self.rmf_output_objects.append(sampler_mc)
        
        # Setup Replica Exchange
        print("Setting up ReplicaExchange")
        rex = IMP.pmi.samplers.ReplicaExchange(
            self.model,
            self.vars["replica_exchange_minimum_temperature"],
            self.vars["replica_exchange_maximum_temperature"],
            [sampler_mc],
            replica_exchange_object=None
        )
        
        self.replica_exchange_object = rex.rem
        
        myindex = rex.get_my_index()
        self.output_objects.append(rex)
        self.rmf_output_objects.append(rex)
        
        # Reset minimum temperature due to binary precision
        min_temp_index = int(min(rex.get_temperatures()) * temp_index_factor)
        
        # Setup directories
        globaldir = self.vars["global_output_directory"] + "/"
        rmf_dir = globaldir + self.vars["rmf_dir"]
        
        # Create directories
        try:
            os.makedirs(globaldir)
        except:
            pass
        try:
            os.makedirs(rmf_dir)
        except:
            pass
        
        # Add stopwatch
        sw = IMP.pmi.tools.Stopwatch()
        self.output_objects.append(sw)
        self.rmf_output_objects.append(sw)
        
        # Setup output
        output = IMP.pmi.output.Output(atomistic=False)
        
        # Update model before initial save
        self.model.update()
        
        # Setup initial RMF
        print("Setting up and writing initial rmf coordinate file")
        init_rmf_name = globaldir + "initial." + str(myindex) + ".rmf3"
        output.init_rmf(init_rmf_name, [self.root_hier], 
                       listofobjects=self.rmf_output_objects)
        output.write_rmf(init_rmf_name)
        output.close_rmf(init_rmf_name)
        
        # Setup production RMF
        print("Setting up production rmf files")
        rmfname = rmf_dir + "/" + str(myindex) + ".rmf3"
        output.init_rmf(rmfname, [self.root_hier],
                       listofobjects=self.rmf_output_objects)
        
        # Setup MPI values (or mock if not using MPI)
        try:
            mpivs = IMP.pmi.samplers.MPI_values(self.replica_exchange_object)
        except:
            # Simple mock class if MPI not available
            class MockMPIValues:
                def set_value(self, key, value): pass
                def get_percentile(self, key): return 0.0
            mpivs = MockMPIValues()
        
        ntimes_at_low_temp = 0
        
        # Show info
        if myindex == 0:
            self.show_info()
            
        self.replica_exchange_object.set_was_used(True)
        
        # Main simulation loop
        nframes = self.vars["number_of_frames"]
        print(f"Running {nframes} frames")
        
        for i in range(nframes):
            # Sampling rounds
            for nr in range(self.vars["num_sample_rounds"]):
                sampler_mc.optimize(self.vars["monte_carlo_steps"])
            
            # Evaluate score using IMP.pmi.tools method (as in original)
            score = IMP.pmi.tools.get_restraint_set(self.model).evaluate(False)
            mpivs.set_value("score", score)
            
            # Check if at minimum temperature
            my_temp_index = int(rex.get_my_temp() * temp_index_factor)
            save_frame = (min_temp_index == my_temp_index)
            
            # Update model before saving
            if save_frame:
                self.model.update()
                
            # Save coordinates if at minimum temperature
            if save_frame:
                print(f"--- frame {i} score {score}")
                
                if i % self.vars["nframes_write_coordinates"] == 0:
                    print('--- writing coordinates')
                    output.write_rmf(rmfname)
                    
                ntimes_at_low_temp += 1
            
            # Attempt temperature swap
            if self.vars["replica_exchange_swap"]:
                rex.swap_temp(i, score)
        
        # Close RMF file
        print("closing production rmf files")
        output.close_rmf(rmfname)
        
        print(f"Replica {myindex} complete")
        print(f"Saved {ntimes_at_low_temp} frames at minimum temperature")
        
    def show_info(self):
        """Display configuration info"""
        print("SimpleReplicaExchange: generates initial.*.rmf3 and rmfs/*.rmf3 for each replica")
        print("--- coordinates saved only at lowest temperature")
        print("--- variables:")
        keys = list(self.vars.keys())
        keys.sort()
        for v in keys:
            print("------", v.ljust(30), self.vars[v])