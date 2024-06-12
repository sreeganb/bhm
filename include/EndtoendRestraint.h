/**
 *  \file IMP/bhm/EndtoendRestraint.h
 *  \brief A lognormal restraint for the end to end distance of a string
 *  of beads.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBHM_Endtoend_RESTRAINT_H
#define IMPBHM_Endtoend_RESTRAINT_H

#include <IMP/bhm/bhm_config.h>
#include <IMP/SingletonScore.h>
#include <IMP/core/XYZ.h>
#include <IMP/Restraint.h>
#include <IMP/PairContainer.h>
#include <IMP/PairScore.h>

IMPBHM_BEGIN_NAMESPACE

//! Apply an end to end distance restraint between two particles.
class IMPBHMEXPORT EndtoendRestraint : public Restraint {
  Pointer<Particle> p0_;
  Pointer<Particle> p1_;
  Pointer<Particle> sigma_;
  double Vexp_;
  double chi_;
  void set_chi(double chi) { chi_ = chi; }

 public:
  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */
  EndtoendRestraint(Model *m, Particle *p0, Particle *p1, Particle *sigma,
                double Iexp);

  /* call for probability */
  double get_probability() const { return exp(-unprotected_evaluate(nullptr)); }

  double get_chi() const { return chi_; }

  virtual double unprotected_evaluate(IMP::DerivativeAccumulator *accum)
      const override;
  virtual IMP::ModelObjectsTemp do_get_inputs() const override;
  IMP_OBJECT_METHODS(EndtoendRestraint);
};

IMPBHM_END_NAMESPACE

#endif /* IMPBHM_Endtoend_RESTRAINT_H */

