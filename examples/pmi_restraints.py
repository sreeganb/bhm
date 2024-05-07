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

class ConnectAtomsRestraint(IMP.pmi.restraints.RestraintBase):
    """
    Restraint to keep the sequence of atoms connected because this is necessary to
    keep the protein structure intact. For pairs of atoms in the backbone representation
    add connectivity restraints between them.
    """
    def __init__(self,objects,scale = 1.0,disorderedlength=False,upperharmonic=True,resolution=0,label=None):
        hiers = IMP.pmi.tools.input_adaptor(objects, resolution)
        if len(hiers) > 1:
            raise Exception("ConnectivityRestraint: only pass stuff from "
                            "one Molecule, please")
        hiers = hiers[0]
        m = list(hiers)[0].get_model()
        super(ConnectAtomsRestraint, self).__init__(m, label=label)

        self.kappa = 10  # spring constant used for the harmonic restraints
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
        for x in range(len(SortedSegments) - 1):

            last = SortedSegments[x][1]
            first = SortedSegments[x + 1][0]

            apply_restraint = True

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
                    optdist = (0.0 + (float(residuegap) + 1.0) * 3.6) * scale
                    if upperharmonic:  # default
                        hu = IMP.core.HarmonicUpperBound(optdist, self.kappa)
                    else:
                        hu = IMP.core.Harmonic(optdist, self.kappa)
                    dps = IMP.core.SphereDistancePairScore(hu)

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
        

