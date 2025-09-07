/**
 *  \file bhm/EndtoendRestraint.cpp
 *  \brief Restrain a list of particle pairs with a lognormal restraint.
 *  NOTE: for now, the derivatives are written to all variables.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bhm/EndtoendRestraint.h>
#include <IMP/isd/FNormal.h>
#include <IMP/isd/Scale.h>
#include <IMP/core/XYZ.h>
#include <IMP/UnaryFunction.h>
#include <cmath>

IMPBHM_BEGIN_NAMESPACE

EndtoendRestraint::EndtoendRestraint(Model *m, Particle *p0, Particle *p1,
                           Particle *sigma, double Vexp)
    : Restraint(m, "EndtoendRestraint%1%"),
      p0_(p0),
      p1_(p1),
      sigma_(sigma),
      Vexp_(Vexp) {}

/* Apply the restraint to two atoms, two Scales, one experimental value.
 */
double EndtoendRestraint::unprotected_evaluate(DerivativeAccumulator *accum) const {
  core::XYZ d0(p0_);
  core::XYZ d1(p1_);
  isd::Scale sigma_nuis(sigma_);
  /* compute Icalc */
  algebra::Vector3D c0 = d0.get_coordinates();
  algebra::Vector3D c1 = d1.get_coordinates();
  double diff = (c0 - c1).get_magnitude();
  double sigma_val = sigma_nuis.get_scale();
  double Icalc = 1.0 * pow(diff, 1);
  /* compute all arguments to FNormal */
  double FA = log(Vexp_);
  double FM = log(Icalc);
  double JA = 1.0 / Vexp_;
  IMP_NEW(isd::FNormal, lognormal, (FA, JA, FM, sigma_val));
  lognormal->set_was_used(true);
  /* get score */
  double score = lognormal->evaluate();
  const_cast<EndtoendRestraint *>(this)->set_chi(FA - FM);

  if (accum) {
    /* derivative for coordinates */
    double DFM = lognormal->evaluate_derivative_FM();
    double factor = 1 / diff; /* d(log(1.0*pow(diff, 1)))/d(diff) */
    algebra::Vector3D deriv = DFM * factor * (c0 - c1) / diff;
    d0.add_to_derivatives(deriv, *accum);
    d1.add_to_derivatives(-deriv, *accum);
    /* derivative for sigma */
    sigma_nuis.add_to_scale_derivative(lognormal->evaluate_derivative_sigma(),
                                       *accum);
    /* derivative for gamma */
  /*  gamma_nuis.add_to_scale_derivative(DFM / gamma_val, *accum); */
  }
  return score;
}

ModelObjectsTemp EndtoendRestraint::do_get_inputs() const {
  ParticlesTemp ret;
  ret.push_back(p0_);
  ret.push_back(p1_);
  ret.push_back(sigma_);
  return ret;
}

IMPBHM_END_NAMESPACE
