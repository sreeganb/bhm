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
import IMP.algebra
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.topology
import numpy as np
import os
import sys
import IMP.pmi.macros
import IMP.pmi.dof
import IMP.rmf
import IMP.display
import IMP.pmi.output
import IMP.bhm
import IMP.bhm.restraints.strings
import IMP.bhm.restraints.pmi_restraints
import IMP.bhm.samplers.two_level_mcmc
import IMP.bhm.samplers.mcmc_multilevel
import IMP.bhm.system_representation.build
#----------------------------------------------------------------------
# New build system for the monomer 
#----------------------------------------------------------------------
ids = ["A", "B", "C"]
mdl = IMP.Model()
sequence = 'A'*10
#--------------------------------------------------
# New code
# State 1
#--------------------------------------------------
s1 = IMP.pmi.topology.System(mdl)
st1 = s1.create_state()
num_mols = 1
mols1 = []
k = 0
for j in range(num_mols):
    m1 = st1.create_molecule("prot%s"%(ids[k]), sequence, chain_id ='%s'%(ids[j]))
    m1.add_representation(m1, resolutions = [1])
    mols1.append(m1)
    k += 1
r1_hier = s1.build()
dof_s1 = IMP.pmi.dof.DegreesOfFreedom(mdl)
IMP.pmi.tools.shuffle_configuration(r1_hier, max_translation=30)
for k in mols1:
     dof_s1.create_flexible_beads(k, max_trans = 2.0)
#IMP.atom.show_with_representations(r1_hier)
#------------------------------------
# state 2
#------------------------------------
s2 = IMP.pmi.topology.System(mdl)
st2 = s2.create_state()
num_mols = 2
mols2 = []
k = 1
for j in range(num_mols):
    m2 = st2.create_molecule("prot%s"%(ids[k]), sequence, chain_id ='%s'%(ids[k]))
    m2.add_representation(m2, resolutions = [1])
    mols2.append(m2)
    k += 1
r2_hier = s2.build()
dof_s2 = IMP.pmi.dof.DegreesOfFreedom(mdl)
IMP.pmi.tools.shuffle_configuration(r2_hier, max_translation=30)
for k in mols2:
     dof_s2.create_flexible_beads(k, max_trans = 2.0)
#IMP.atom.show_with_representations(r2_hier)
#--------------------------------------------------
#print("test ", r2_hier.get_child(0).get_child(1).get_child(0).get_child(0).get_children()[-1].get_particle())
for i in range(len(mols1)):
    cr1 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mols1[i])
    cr1.add_to_model()
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(r1_hier)
evr1.add_to_model()
for i in range(len(mols2)):
    cr2 = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mols2[i])
    cr2.add_to_model()
evr2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(r2_hier)
evr2.add_to_model()
# Add the end to end restraint to the model by reading the data from a file
etedata = np.loadtxt('./derived_data/synthetic_data_monomer.txt')
etr = IMP.bhm.restraints.pmi_restraints.EndToEndRestraint(r1_hier, etedata, label = "endtoend", weight = 1.0)
#etr.add_to_model()  # add restraint to model
#dof_s1.get_nuisances_from_restraint(etr)

etedata = np.loadtxt('./derived_data/synthetic_data_dimer.txt')
etr2 = IMP.bhm.restraints.pmi_restraints.EndToEndRestraint(r2_hier, etedata, label = "endtoend", weight = 1.0)
etr2.add_to_model()  # add restraint to model
dof_s2.get_nuisances_from_restraint(etr2)
#IMP.bhm.samplers.mcmc_multilevel.MCMCsampler(r1_hier, dof_s1, 2.0, 400, "monomer.rmf3")
IMP.bhm.samplers.mcmc_multilevel.MCMCsampler(r2_hier, dof_s2, 2.0, 1000, "dimer.rmf3")