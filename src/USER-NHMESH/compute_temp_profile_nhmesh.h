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

ComputeStyle(temp/profile/nhmesh,ComputeTempProfileNHMesh)

#else

#ifndef LMP_COMPUTE_TEMP_PROFILE_NHMESH_H
#define LMP_COMPUTE_TEMP_PROFILE_NHMESH_H

#include "compute_temp_profile.h"

namespace LAMMPS_NS {

class ComputeTempProfileNHMesh : public ComputeTempProfile {
 public:
  ComputeTempProfileNHMesh(class LAMMPS *, int, char **);
  ~ComputeTempProfileNHMesh();
  virtual void compute_array();

 protected:
  char *idcoupling;
  class ComputeNHMeshCouplingAtom *coupling;
  int n_thermostats;

  void (ComputeTempProfileNHMesh::*array_compute_fn)();
  virtual void compute_array_bin();
  virtual void compute_array_ke();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Compute ID for temp/profile/nhmesh does not exist

Self-explanatory.

E: Invalid coupling compute for temp/profile/nhmesh

Coupling compute must be derived from nhmesh/coupling/atom.

W: Compute temp/profile/nhmesh can't be used as the temperature compute for fix
   nvt/nhmesh with out bin.

The out bin flag changes the compute_array function, making it return results
that are incompatible with fix temp/nhmesh. This could cause undefined
behaviour if used as the temperature compute for fix temp/nhmesh.

*/
