## \example bhm/helix.py
# In this example we will build a helix and define the parameters associated
# with this and define a prior 
# The data likelihood is defined by the end to end distance for the structure
# The error parameter is the simple sigma with an associated Jeffrey's prior
from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.crosslinking
import IMP.isd
import IMP.pmi.dof
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.tools
import IMP.pmi.macros
import IMP.atom
import sys
import IMP.rmf
from pmi_restraints import ConnectAtomsRestraint, DihedralHelixRestraint

if len(sys.argv) > 1:
    out_dir = sys.argv[1]
else:
    out_dir = "output"

#------------------------------------------------------------------------------
# Read in input information
#------------------------------------------------------------------------------
pdb_dir = "./data/pdb/"
fasta_dir = "./data/fasta/"
saxs_data = "./derived_data/saxs/ala19lys_bb.pdb.dat"
saxs_weight = 1.0

sequences = IMP.pmi.topology.Sequences(fasta_dir + "ala19lys1.fasta.txt")

mdl = IMP.Model() # Define the model
sys = IMP.pmi.topology.System(mdl) # define the system
st = sys.create_state() # define the state

alalys = st.create_molecule("A", sequence=sequences["ala19lys"])
a1 = alalys.add_structure(pdb_dir + "ala19lys_bb.pdb", chain_id = "A")
#------------------------------------------------------------------------------
# Add representation and test it
# Going with a pdb file that contains CA, O, N, C atoms, so 4 beads per
# residue
#------------------------------------------------------------------------------
alalys.add_representation(a1, resolutions=0) #, ideal_helix=True) 
hier = sys.build()
beads = IMP.atom.Selection(hier).get_selected_particles()

IMP.atom.show_with_representations(hier)
#------------------------------------------------------------------------------
# Define the degrees of freedom and create movers 
#------------------------------------------------------------------------------
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
lys_rb = alalys[0:3]
fx_bead = alalys[3:20]
rb1 = dof.create_rigid_body(lys_rb, max_trans=1.0, max_rot=0.5, nonrigid_parts = lys_rb & alalys.get_non_atomic_residues())

ala1 = dof.create_flexible_beads(fx_bead, max_trans=1.0)
#------------------------------------------------------------------------------
# Add restraints: in this case the data corresponds to some end to end distances
# as well as some SAXS profiles presumably. The distance restraint is added 
# as the mean of the distribution obtained from MD simulations. SAXS profile
# is obtained from the FOXS webserver. Adding contraints on the phi, psi angles
# would result in formation of a helix, but how do we contrain dihedrals in IMP?
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Add restraints for the monomer 
#------------------------------------------------------------------------------
# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT (modified as connectivity atoms restraint)
#
# Restrains residues/particles that are collected in sequence
# This should be used for any system without an atomic force field
# (e.g. CHARMM)
# Here, we pass root_hier to apply this restraint to the entire system
#cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(alalys)
#cr.add_to_model()           # add restraint to the model
#output_objects.append(cr)   # add restraint to the output
cr = ConnectAtomsRestraint(alalys)
cr.add_to_model()
output_objects.append(cr)

# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of both molecule chains to included_objects to
# apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest
# resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[alalys],
                                            resolution=1000)
output_objects.append(evr)
# -------------------------
# %%%%%% Dihedral restraint
# Just get the dihedrals and print as a first step
# -------------------------
# print("residues: ", IMP.atom.Residue(alalys, 1))
#for r in alalys.residues:
#    print(r.get_residue_type())
res = IMP.atom.get_by_type(hier, IMP.atom.RESIDUE_TYPE)
#print("residues: ", res)
for h in res[1:-2]:
    r = IMP.atom.Residue(h)
    phi = IMP.atom.get_phi_dihedral_atoms(r)
    #print("atoms: ", phi[0])
    d = IMP.core.get_dihedral(*[IMP.core.XYZ(x) for x in phi]) 
    #print(d)
    pos = []
    for i in range(4):
        pos.append(IMP.core.XYZ(phi[i]))
    #print(ParticleIndexAdaptor phi[0])
#    dih_rest = IMP.core.DihedralRestraint(mdl, ParticleIndexAdaptor phi[0], phi[1], phi[2], phi[3])
#    ouput_objects.append(dih_rest)
    #print(IMP.atom.get_phi_dihedral_atoms(r))
    #print(IMP.atom.Residue(r).get_residue_type())
#print("residues: ", mdl.residues)
#r = IMP.atom.Residue()
#phi = IMP.atom.get_phi_dihedral_atoms()


# -------------------------
# %%%%% SAXS RESTRAINT
#
# Scores the SAXS data against the predicted SAXS for the input_objects
# of the model
#
# The form factor type (ff_type) for calculating the model SAXS profile
# can be either:
#   IMP.saxs.RESIDUES - Residue-level calculation - requires a model
#                       at resolution=1
#   IMP.saxs.CA_ATOMS - Residue-level calculation - requires a model
#                       at atomic resolution
#   IMP.saxs.HEAVY_ATOMS - No explicit hydrogens - requires a model
#                          at atomic resolution
#   IMP.saxs.ALL_ATOMS - Uses explicit hydrogens - requires a model
#                        at atomic resolution with hydrogens
sr = IMP.pmi.restraints.saxs.SAXSRestraint(
    input_objects=[alalys],
    saxs_datafile=saxs_data,
    weight=saxs_weight,         # scaling factor for the SAXS score
    ff_type=IMP.saxs.CA_ATOMS,
    # Maximum q at which to compare SAXS curves. Defaults are:
    #   ~0.15 for residue-level calculations
    #   0.4 for HEAVY_ATOMS
    #   0.5 for ALL_ATOMS (up to 1.0)
    maxq=0.15)
#output_objects.append(sr)

#####################################################
#                      SAMPLING                     #
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(hier,
                                    max_translation=10)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)

evr.add_to_model()
#emr.add_to_model()
#slx.add_to_model()
#sr.add_to_model()

# Run replica exchange Monte Carlo sampling
rex = IMP.pmi.macros.ReplicaExchange(
    mdl,
    # pass the root hierarchy
    root_hier=hier,
    # pass all objects to be moved ( almost always dof.get_movers() )
    monte_carlo_sample_objects=dof.get_movers(),
    # The output directory for this sampling run.
    global_output_directory='run_manual1/output/',
    # Items in output_objects write information to the stat file.
    output_objects=output_objects,
    # Number of MC steps between writing frames
    monte_carlo_steps=1,
    # set >0 to store best PDB files (but this is slow)
    number_of_best_scoring_models=2,
    # Total number of frames to run / write to the RMF file.
    number_of_frames=200)

# Ok, now we finally do the sampling!
rex.execute_macro()

# Outputs are then analyzed using a separate analysis script.

# hrest = IMP.pmi.restraints.stereochemistry.HelixRestraint(hier, (0,20,alalys)) this didnt work, of course



# Write a single frame of RMF file to view it
#out = IMP.pmi.output.Output()
#out.init_rmf("ala19lys1-helix.rmf3", hierarchies=[hier])
#out.write_rmf("ala19lys1-helix.rmf3")
