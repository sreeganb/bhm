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
import IMP.pmi.restraints.crosslinking
import ihm.cross_linkers
import IMP.isd
import IMP.pmi.dof
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.tools
import IMP.pmi.macros
import IMP.atom
import sys
import IMP.rmf

if len(sys.argv) > 1:
    out_dir = sys.argv[1]
else:
    out_dir = "output"

#------------------------------------------------------------------------------
# Read in input information
#------------------------------------------------------------------------------
pdb_dir = "./data/pdb/"
fasta_dir = "./data/fasta/"
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
# Add restraints: in this case the data corresponds to some end to end distances
# as well as some SAXS profiles presumably. The distance restraint is added 
# as the mean of the distribution obtained from MD simulations. SAXS profile
# is obtained from the FOXS webserver. Adding contraints on the phi, psi angles
# would result in formation of a helix, but how do we contrain dihedrals in IMP?
# 
#------------------------------------------------------------------------------
# Use a ideal helix restraint class and check what happens
# hrest = IMP.pmi.restraints.stereochemistry.HelixRestraint(hier, (0,20,alalys)) this didnt work, of course



# Write a single frame of RMF file to view it
# out = IMP.pmi.output.Output()
# out.init_rmf("ala19lys1-helix.rmf3", hierarchies=[hier])
# out.write_rmf("ala19lys1-helix.rmf3")
