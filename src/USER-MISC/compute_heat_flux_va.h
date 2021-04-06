/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/*------------------------------------------------------------------------
  Contributing Authors : Stephen Sanderson (UQ)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS

ComputeStyle(heat/flux/va,ComputeHeatFluxVA)

#else

#ifndef LMP_COMPUTE_HEAT_FLUX_VA_H
#define LMP_COMPUTE_HEAT_FLUX_VA_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeHeatFluxVA : public Compute {
  public:
    ComputeHeatFluxVA(class LAMMPS *, int, char **);
    virtual ~ComputeHeatFluxVA();
    void init();
    void init_list(int, class NeighList *);
    void compute_array();

  private:

    void compute_pairs();

    int me,nvalues,dir;
    int *which;

    double *values_local,*values_global;
    double pos,pos1,dt,nktv2p,ftm2v;
    double area;
    class NeighList *list;

  };

}

#endif
#endif

/* ERROR/WARNING messages:

 E: Illegal ... command

 Self-explanatory.  Check the input script syntax and compare to the
 documentation for the command.  You can use -echo screen as a
 command-line option when running LAMMPS to see the offending line.

 E: No pair style is defined for compute heat/flux/va

 Self-explanatory. Compute heat/flux/va requires the definition of a pair style.

 E: Pair style does not support compute heat/flux/va

 The pair style does not have a single() function, so it can
 not be invoked by compute heat/flux/va.

 W: compute heat/flux/va does not account for bond potentials

 W: compute heat/flux/va does not account for angle potentials

 W: compute heat/flux/va does not account for dihedral potentials

 W: compute heat/flux/va does not account for improper potentials

 W: compute heat/flux/va does not account for kspace contributions

 Compute heat/flux/va only accounts for pairwise additive interactions for
 the computation of local stress tensor components.

 */

