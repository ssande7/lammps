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

#include "compute_temp_profile.h"

namespace LAMMPS_NS {

class ComputeTempProfileNHMesh : public ComputeTempProfile {
 public:
  ComputeTempProfileNHMesh(class LAMMPS *, int, char **);
  ~ComputeTempProfileNHMesh();
  void compute_array();

 protected:
  char *idcoupling;
  class ComputeCouplingNHMesh *coupling;
  int n_thermostats;

  void (ComputeTempProfileNHMesh::*array_compute_fn)();
  virtual void compute_array_bin();
  virtual void compute_array_ke();

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
