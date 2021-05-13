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
  int find_crossing(double *, double *);
  template<int>
  int crossing_bin1d(double *, double *, int *, double *, int **);
  int crossing_bin2d(double *, double *, int *, double *, int **);
  int crossing_bin3d(double *, double *, int *, double *, int **);
  // int crossing_binsphere(double *, double *, int *, double *);
  // int crossing_bincylinder(double *, double *, int *, double *);
  int x2chunk(double *);
  void allocate();

  int me;
  int deform_vremap;

  int biasflag;
  char *id_temp_c, *id_temp_k;
  class ComputeTempChunk *c_temp_c;
  class Compute *c_temp_k;

  char *id_pe;
  class Compute *c_pe;

  int nchunk,maxchunk,linechunk,maxlinechunk;
  char *id_chunk;
  class ComputeChunkAtom *c_chunk;

  double **cvalues_local,**cvalues;
  double *cfactor, *cfactor2d, *cfactor3d;
  int *c_ids, *c_ids2d, *c_ids3d;
  int **c_wrap, **c_wrap2d, **c_wrap3d;
  double c_hi[3];
  bool novoid[3];
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

E: Could not find compute heat/flux/va/chunk compute ID

Self-explanatory.

E: Compute heat/flux/va/chunk compute ID does not compute ke/atom

Self-explanatory.

E: Compute heat/flux/va/chunk compute ID does not compute pe/atom

Self-explanatory.

E: Compute heat/flux/va/chunk requires chunk/atom compute

Self-explanatory.

E: Bias compute does not calculate a velocity bias

Self-explanatory.

E: No pair style is defined for compute heat/flux/va/chunk

Self-explanatory. Compute heat/flux/va/chunk requires the
definition of a pair style.

E: Pair style does not support compute heat/flux/va/chunk

The pair style does not have a single() function, so it can
not be invoked by compute heat/flux/va/chunk.

W: compute heat/flux/va/chunk does not account for bond potentials

W: compute heat/flux/va/chunk does not account for angle potentials

W: compute heat/flux/va/chunk does not account for dihedral potentials

W: compute heat/flux/va/chunk does not account for improper potentials

W: compute heat/flux/va/chunk does not account for kspace contributions

Compute heat/flux/va/chunk only accounts for pairwise additive interactions for
the computation of local stress tensor components.

*/
