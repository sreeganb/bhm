"""PMI-style restraints."""
import math
import IMP
import IMP.isd
import IMP.core
import IMP.container
import IMP.atom
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.restraints
from operator import itemgetter
import math
import regex

class ConnectAtomsRestraint(IMP.pmi.restraints.RestraintBase):
    """
    Restraint to keep the sequence of atoms connected because this is necessary to
    keep the protein structure intact. For pairs of atoms in the backbone representation
    add connectivity restraints between them. Not all pairs are connected, so only the 
    neighboring atoms in the list needs to be attached, except for the first and last 
    residues. IN this case, 
          O     
         //     
    N-CA-C-N-CA-C-N-CA-C
                \\
                 O 
    This set of atoms is useful for the definition of a dihedral restraint as well
    """
    def get_bond_length(self,first, last):
        """
        Function to pass the bond distances to be used as constraints based on the atom pair
        """
        #print(IMP.atom.Atom(first).get_atom_type(), IMP.atom.Atom(last).get_atom_type())
        atom1 = IMP.atom.Atom(first).get_atom_type()
        atom2 = IMP.atom.Atom(last).get_atom_type()
        ca_at = IMP.atom.AtomType("CA")
        c_at = IMP.atom.AtomType("C")
        o_at = IMP.atom.AtomType("O")
        n_at = IMP.atom.AtomType("N") 
        if (atom1 == n_at or atom1 == ca_at) and (atom2 == n_at or atom2 ==ca_at):
            dist = 1.46
        elif (atom1 == c_at or atom1 == ca_at) and (atom2 == c_at or atom2 ==ca_at):
            dist = 1.53
        elif (atom1 == c_at or atom1 == o_at) and (atom2 == c_at or atom2 ==o_at):
            dist = 1.43
        elif (atom1 == c_at or atom1 == n_at) and (atom2 == n_at or atom2 ==c_at):
            dist = 1.48
        else:
            print("wrong combination of atoms, check the input")
        return dist

    def __init__(self,objects,scale = 1.0,disorderedlength=False,upperharmonic=True,resolution=0,label=None):
        hiers = IMP.pmi.tools.input_adaptor(objects, resolution)
        if len(hiers) > 1:
            raise Exception("ConnectivityRestraint: only pass stuff from "
                            "one Molecule, please")
        hiers = hiers[0]
        m = list(hiers)[0].get_model()
        super(ConnectAtomsRestraint, self).__init__(m, label=label)

        self.kappa = 10  # spring constant used for the harmonic restraints changed from 10 to 20 since we have atoms now
        SortedSegments = []
        for h in hiers:
            try:
                start = IMP.atom.Hierarchy(h).get_children()[0]
            except:  # noqa: E722
                start = IMP.atom.Hierarchy(h)

            try:
                end = IMP.atom.Hierarchy(h).get_children()[-1]
            except:  # noqa: E722
                end = IMP.atom.Hierarchy(h)

            startres = IMP.pmi.tools.get_residue_indexes(start)[0]
            SortedSegments.append((start, end, startres))
        SortedSegments = sorted(SortedSegments, key=itemgetter(2))
        #print("length of the sorted segments is: ", len(SortedSegments))
        #print("these are the sorted segments: ", SortedSegments[2])
        # connect the particles
        self.particle_pairs = []
        nres = max(IMP.pmi.tools.get_residue_indexes(end))
        j = 0 # counter for the number of particles to be connected
        #--------------------------------------------------------------------
        # Loop over residues and atoms within residues and add connectivities
        #--------------------------------------------------------------------
        for x in range(nres): # this should loop over the number of residues
            for y in range(4): # hardcoded 4 here, but need to be properly changed based on the number of atoms in the coarse grained representation
                if (j!=nres*4-2):
                    first = SortedSegments[j+1][0]
                else:
                    break
                last = SortedSegments[j][1]
                if (y==3):
                    last = SortedSegments[j - 1][1]
                    first = SortedSegments[j + 1][0]
                    if (x == nres-1):
                        break
                apply_restraint = True
                j=j+1
                # Apply connectivity runless ALL of the following are true:
                #   - first and last both have RigidBodyMember decorators
                #   - first and last are both RigidMembers
                #   - first and last are part of the same RigidBody object

                # Check for both in a rigid body
                if IMP.core.RigidBodyMember.get_is_setup(first) and \
                        IMP.core.RigidBodyMember.get_is_setup(last) and \
                        IMP.core.RigidMember.get_is_setup(first) and \
                        IMP.core.RigidMember.get_is_setup(last):
                    # Check if the rigid body objects for each particle are
                    # the same object.
                    #  if so, skip connectivity restraint
                    if IMP.core.RigidBodyMember(first).get_rigid_body() \
                            == IMP.core.RigidBodyMember(last).get_rigid_body():
                        apply_restraint = False

                if apply_restraint:

                    nreslast = len(IMP.pmi.tools.get_residue_indexes(last))
                    lastresn = IMP.pmi.tools.get_residue_indexes(last)[-1]
                    nresfirst = len(IMP.pmi.tools.get_residue_indexes(first))
                    firstresn = IMP.pmi.tools.get_residue_indexes(first)[0]
                    residuegap = firstresn - lastresn - 1
                    if disorderedlength and (nreslast / 2 + nresfirst / 2
                                             + residuegap) > 20.0:
                        # calculate the distance between the sphere centers
                        # using Kohn PNAS 2004
                        optdist = math.sqrt(5 / 3) * 1.93 * \
                            (nreslast / 2 + nresfirst / 2 + residuegap) ** 0.6
                        if upperharmonic:
                            hu = IMP.core.HarmonicUpperBound(optdist, self.kappa)
                        else:
                            hu = IMP.core.Harmonic(optdist, self.kappa)
                        dps = IMP.core.DistancePairScore(hu)
                    else:  # default
                        #optdist = (0.0 + (float(residuegap) + 1.0) * 3.6) * scale
                        #print("first,last: ", first, last)
                        optdist = self.get_bond_length(first, last)
                        #print("optdist value: ", optdist)
                        if upperharmonic:  # default
                            hu = IMP.core.HarmonicUpperBound(optdist, self.kappa)
                        else:
                            hu = IMP.core.Harmonic(optdist, self.kappa)
                        dps = IMP.core.SphereDistancePairScore(hu)
                        
                    pt0 = last.get_particle()
                    pt1 = first.get_particle()
                    self.particle_pairs.append((pt0, pt1))
                    #print("particles list: ", pt0,pt1)
                    r = IMP.core.PairRestraint(
                        self.model, dps, (pt0.get_index(), pt1.get_index()))
                    
                    print("Adding sequence connectivity restraint between",
                        pt0.get_name(), " and ", pt1.get_name(), 'of distance',
                        optdist)
                    self.rs.add_restraint(r)

    def get_num_restraints(self):
        """ Returns number of connectivity restraints """
        return len(self.rs.get_restraints())

    def get_particle_pairs(self):
        """ Returns the list of connected particles pairs """
        return self.particle_pairs
 
# Restraint the dihedrals to form a helix
class DihedralHelixRestraint(object):
    """ Restrain a protein using the ideal helix dihedrals 
        from the Ramachandran plot. In this case the phi angle is 
    """
    def __init__(self, tuple_selection_quads, dihedrals, root_hier):
        self.m = root_hier.get_model()
        # create restraint sets to join restraints
        self.rs = IMP.RestraintSet(self.m, "likelihood")
        self.rs_priors = IMP.RestraintSet(self.m, "priors")

        # create nuisance particles
        self.kappa = IMP.pmi.tools.SetupNuisance(
            self.m, 1., 1e-5, 1e5, isoptimized=True).get_particle()

        # create NOE likelihood restraints
        for (tuple1, tuple2, tuple3, tuple4), angle in zip(
                tuple_selection_quads, dihedrals):
            # get atom 1
            p1 = self.get_atom_from_selection(root_hier, tuple1)
            p2 = self.get_atom_from_selection(root_hier, tuple2)
            p3 = self.get_atom_from_selection(root_hier, tuple3)
            p4 = self.get_atom_from_selection(root_hier, tuple4)

            # create restraint and add to set
            angle = angle * math.pi / 180.  # convert to radians
            r = IMP.isd.TALOSRestraint(self.m, p1, p2, p3, p4, [angle],
                                       self.kappa)
            r.set_name(
                "DihedralRestraint_{0}:{1}_{2}:{3}_{4}:{5}_{6}:{7}".format(
                    tuple1[0], tuple1[1], tuple2[0], tuple2[1],
                    tuple3[0], tuple3[1], tuple4[0], tuple4[1]))
            self.rs.add_restraint(r)

        # create prior restraints
        self.rs_priors.add_restraint(IMP.isd.vonMisesKappaJeffreysRestraint(
            self.m, self.kappa))

    def get_atom_from_selection(self, root_hier, sel_tuple):
        res_id, atom_name = sel_tuple
        sel = IMP.atom.Selection(root_hier,
                                 resolution=0,
                                 residue_index=res_id,
                                 atom_type=IMP.atom.AtomType(atom_name))
        return sel.get_selected_particles()[0]

    def add_to_model(self):
        """Add the restraints to the model."""
        for rs in [self.rs, self.rs_priors]:
            IMP.pmi.tools.add_restraint_to_model(self.m, rs)

    def get_restraint(self):
        return self.rs

    def get_restraint_for_rmf(self):
        """Get the restraint for visualization in an RMF file."""
        return self.rs

    def get_particles_to_sample(self):
        """Get any created particles which should be sampled."""
        out = {}
        out["Nuisances_Kappa"] = ([self.kappa], .1)
        return out

    def get_output(self):
        """Get outputs to write to stat files."""
        output = {}
        self.m.update()
        likelihood_score = self.rs.unprotected_evaluate(None)
        prior_score = self.rs_priors.unprotected_evaluate(None)
        output["_TotalScore"] = likelihood_score + prior_score
        output["DihedralLikelihood_Score"] = likelihood_score
        output["DihedralPrior_Score"] = prior_score
        output["Dihedral_Kappa"] = self.kappa.get_scale()
        return output
        

