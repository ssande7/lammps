/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(temp/nhmesh/profile,ComputeTempProfileNHMesh)

#else

#ifndef LMP_COMPUTE_TEMP_PROFILE_NHMESH_H
#define LMP_COMPUTE_TEMP_PROFILE_NHMESH_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempProfileNHMesh : public Compute {
 public:
  ComputeTempProfileNHMesh(class LAMMPS *, int, char **);
  ~ComputeTempProfileNHMesh();
  void init();
  void setup();
  double compute_scalar();
  void compute_vector();
  void compute_array();

  void remove_bias(int, double *);
  void remove_bias_thr(int, double *, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_thr(int, double *, double *);
  void restore_bias_all();
  double memory_usage();

 protected:
  char *idcoupling;
  class ComputeCouplingNHMesh *coupling;
  int n_thermostats;

  int xflag,yflag,zflag,ncount,outflag;
  int nbinx,nbiny,nbinz,nbins;
  int ivx,ivy,ivz;
  double tfactor;

  int box_change,triclinic;
  int *periodicity;
  double *boxlo,*boxhi,*prd;
  double invdelta[3];

  int maxatom;
  int *bin;
  double **vbin,**binave;
  double *tbin,*tbinall;

  void (ComputeTempProfileNHMesh::*array_compute_fn)();
  virtual void compute_array_bin();
  virtual void compute_array_ke();

  virtual void dof_compute();
  virtual void bin_average();
  virtual void bin_setup();
  virtual void bin_assign();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for temp/nhmesh/profile does not exist

Self-explanatory.

E: Invalid coupling compute for temp/nhmesh/profile

Coupling compute must be derived from ComputeCouplingNHMesh.

E: Compute temp/nhmesh/profile cannot use vz for 2d systemx

Self-explanatory.

E: Compute temp/nhmesh/profile cannot bin z for 2d systems

Self-explanatory.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

W: Compute temp/nhmesh/profile can't be used as the temperature compute for fix
   temp/nhmesh with out bin.

The out bin flag change the compute_array function, making it return results
that are incompatible with fix temp/nhmesh. This could cause undefined
behaviour if used as the temperature compute for fix temp/nhmesh.

*/
