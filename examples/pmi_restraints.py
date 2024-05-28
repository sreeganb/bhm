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
        atom1 = IMP.atom.Atom(first).get_atom_type()
        atom2 = IMP.atom.Atom(last).get_atom_type()
        ca_at = IMP.atom.AtomType("CA") # define IMP objects with same atom types to compare 
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

        self.kappa = 25  # spring constant used for the harmonic restraints changed from 10 to 20 since we have atoms now
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
                if (y==3): # include the C-CA particle pair
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
                        #optdist = self.get_bond_length(first, last)
                        bl = self.get_bond_length(first, last)
                        optdist = bl *1.0
                        if upperharmonic:  # default
                            hu = IMP.core.HarmonicUpperBound(optdist, self.kappa)
                        else:
                            hu = IMP.core.Harmonic(optdist, self.kappa)
                        #dps = IMP.core.SphereDistancePairScore(hu)
                        dps = IMP.core.DistancePairScore(hu)
                    pt0 = last.get_particle()
                    pt1 = first.get_particle()
                    self.particle_pairs.append((pt0, pt1))
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
class DihedralHelixRestraint(IMP.pmi.restraints.RestraintBase):
    """ Restrain a protein using the ideal helix dihedrals 
        from the Ramachandran plot. In this case the phi angle is 
    """
    def __init__(self, root_hier):
        self.m = root_hier.get_model()
        # create restraint sets to join restraints
        self.rs = IMP.RestraintSet(self.m, "likelihood")
        self.rs_priors = IMP.RestraintSet(self.m, "priors")

        # create nuisance particles
        self.kappa = IMP.pmi.tools.SetupNuisance(
            self.m, 1., 1e-5, 1e5, isoptimized=True).get_particle()
        res = IMP.atom.get_by_type(root_hier,IMP.atom.RESIDUE_TYPE)
            # Takes 4 particles as input and constraints the 
            # dihedrals based on a HarmonicWell scoring function
        k = 0
        for h in res[1:-2]:    
            anglemin = -60.0
            anglemax = -52.0
            strength = 25.0
            ts = IMP.core.HarmonicWell((math.pi * anglemin /180.0, math.pi * anglemax / 180.0), strength)
            # add the restraint
            ra = IMP.atom.Residue(h)
            phi = IMP.atom.get_phi_dihedral_atoms(ra)
            #phi = self.get_phi_psi_atom_tuple(root_hier, h)
            dresphi = IMP.core.DihedralRestraint(self.m, ts, phi[0], phi[1], phi[2], phi[3]) 
            self.rs.add_restraint(dresphi)
            # Same procedure for the psi dihedral angle
            anglemin = -53.0
            anglemax = -41.0
            strength = 25.0
            ts = IMP.core.HarmonicWell((math.pi * anglemin /180.0, math.pi * anglemax / 180.0), strength)
            # add the restraint
            ra = IMP.atom.Residue(h)
            psi = IMP.atom.get_psi_dihedral_atoms(ra)
            #phi = self.get_phi_psi_atom_tuple(root_hier, h)
            print("psi dihedral: ", psi)
            drespsi = IMP.core.DihedralRestraint(self.m, ts, psi[0], psi[1], psi[2], psi[3]) 
            self.rs.add_restraint(drespsi)

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

class DistanceHelixRestraint(IMP.pmi.restraints.RestraintBase):
    """A simple distance restraint with a nuisance particle to include structural uncertainty."""
    def create_sigma(self):
        """Create a nuisance on the length of the helix."""
        lengthinit = 20.0
        self.lengthissampled = True
        lengthminnuis = 4.0
        lengthmaxnuis = 40.0
        length = IMP.pmi.tools.SetupNuisance(self.model, lengthinit,
                                             lengthminnuis, lengthmaxnuis,
                                             self.lengthissampled
                                             ).get_particle()
        self.rs.add_restraint(
            IMP.isd.LognormalRestraint(length, 15.0, 6.0))
        
    def __init__(self, root_hier, tuple_selection1, tuple_selection2,
                 distancemin=0, distancemax=100, resolution=1.0, kappa=1.0,
                 label=None, weight=1.0):
        """Setup distance restraint.
        @param root_hier The hierarchy to select from
        @param tuple_selection1 (resnum, resnum, molecule name, copy
               number (=0))
        @param tuple_selection2 (resnum, resnum, molecule name, copy
               number (=0))
        @param distancemin The minimum dist
        @param distancemax The maximum dist
        @param resolution For selecting particles
        @param kappa The harmonic parameter
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number.
        @dof the degrees of freedom object
        By default this will take the CA (c-alpha) atoms between two chosen residues 
        and adds a distance restraint between them.
        """
        self.model = root_hier.get_model()
        self.rs = IMP.RestraintSet(self.model, "distance_helix_restraint")
        #print("model: ", self.model)

        ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
        ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)
        
        copy_num1 = 0
        if len(tuple_selection1) > 3:
            copy_num1 = tuple_selection1[3]
        copy_num2 = 0
        if len(tuple_selection2) > 3:
            copy_num2 = tuple_selection2[3]

        sel1 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection1[2],
                                  residue_index=tuple_selection1[0],
                                  copy_index=copy_num1)
        particles1 = sel1.get_selected_particles()
        sel2 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection2[2],
                                  residue_index=tuple_selection2[0],
                                  copy_index=copy_num2)
        particles2 = sel2.get_selected_particles()

        super(DistanceHelixRestraint, self).__init__(self.model, label=label,
                                                weight=weight)
        print(f"Created distance restraint between "
            f"{particles1[1].get_name()} and {particles2[1].get_name()}")
        '''
                Add the nuisance particle to the model        '''
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts1,
                                       particles1[1],
                                       particles2[1]))
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts2,
                                       particles1[1],
                                       particles2[1]))
        self.create_sigma()
