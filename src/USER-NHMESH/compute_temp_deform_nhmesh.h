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

ComputeStyle(temp/deform/nhmesh,ComputeTempDeformNHMesh)

#else

#ifndef LMP_COMPUTE_TEMP_DEFORM_NHMESH_H
#define LMP_COMPUTE_TEMP_DEFORM_NHMESH_H

#include "compute_temp_deform.h"

namespace LAMMPS_NS {

class ComputeTempDeformNHMesh : public ComputeTempDeform {
 public:
  ComputeTempDeformNHMesh(class LAMMPS *, int, char **);
  virtual ~ComputeTempDeformNHMesh();
  virtual void compute_array();
  virtual double memory_usage();

 protected:
  char *idcoupling;
  class ComputeNHMeshCouplingAtom *coupling;
  int n_thermostats;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Compute ID for temp/deform/nhmesh does not exist

Self-explanatory.

E: Invalid coupling compute for temp/deform/nhmesh

Coupling compute must be derived from nhmesh/coupling/atom

*/
