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

ComputeStyle(heat/flux/va/chunk,ComputeHeatFluxVAChunk)

#else

#ifndef LMP_COMPUTE_HEAT_FLUX_VA_H
#define LMP_COMPUTE_HEAT_FLUX_VA_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHeatFluxVAChunk : public Compute {
public:
  ComputeHeatFluxVAChunk(class LAMMPS *, int, char **);
  virtual ~ComputeHeatFluxVAChunk();
  void init();
  void setup();
  void init_list(int, class NeighList *);
  void compute_array();

private:

  void compute_flux();
  int find_crossing(int, int, double *, double *, int *, double *);
  int crossing_bin1d(int, int, double *, double *, int *, double *);
  int crossing_bin2d(int, int, double *, double *, int *, double *);
  int crossing_bin3d(int, int, double *, double *, int *, double *);
  int crossing_binsphere(int, int, double *, double *, int *, double *);
  int crossing_bincylinder(int, int, double *, double *, int *, double *);
  void allocate();

  int me;

  int biasflag;
  char *id_temp;
  class Compute *c_temp;

  char *id_ke, *id_pe;
  class Compute *c_ke, *c_pe;

  int nchunk,maxchunk;
  char *id_chunk;
  class ComputeChunkAtom *c_chunk;

  double **cvalues_local,**cvalues;
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

E: Could not find compute heat/flux/va compute ID

Self-explanatory.

E: Compute heat/flux/va compute ID does not compute ke/atom

Self-explanatory.

E: Compute heat/flux/va compute ID does not compute pe/atom

Self-explanatory.

E: Compute heat/flux/va requires chunk/atom compute

Self-explanatory.

E: Bias compute does not calculate a velocity bias

Self-explanatory.

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
