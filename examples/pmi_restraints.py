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

class ConnectAtomsRestraint(IMP.pmi.restraints.RestraintBase):
    """
    Restraint to keep the sequence of atoms connected because this is necessary to
    keep the protein structure intact
    """
    def __init__(self, 
                 objects,
                 scale = 1.0,
                 disorderedlength=False,
                 upperharmonic=True,
                 resolution=0,
                 label=None):
        self.m = root_hier.get_model()
        

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
        

