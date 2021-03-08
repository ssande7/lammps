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

// TODO:
// * compute class to return kinetic temperature relevant to each thermostat
// * dof_compute function for DoFs per atom (stored in dof). Thermostat fix then
//   multiplies this by the sum of each row of the coupling matrix to get atom
//   DoFs per thermostat.
// * add protected variable to store total DoFs for scalar output??
// * add memory_usage() function

// NOTE:
// * takes compute_coupling_nhmesh or compute_coupling_conserve_nhmesh as input
// * returns (for compatability with tempflag):
//    + scalar: average kinetic temperature
//    + vector: vector of length 6 with total KE tensor
//    + array:  n_thermostats x 6 array. Each row is the KE tensor for that
//              thermostat. Summing the row and multiplying by tfactor gives the
//              temperature for each thermostat

#ifdef COMPUTE_CLASS

ComputeStyle(temp/nhmesh,ComputeTempNHMesh)

#else

#ifndef LMP_COMPUTE_TEMP_NHMESH_H
#define LMP_COMPUTE_TEMP_NHMESH_H

#include "compute_temp.h"

namespace LAMMPS_NS {

class ComputeTempNHMesh : public ComputeTemp {
 public:
  ComputeTempNHMesh(class LAMMPS *, int, char **);
  virtual ~ComputeTempNHMesh();
  virtual void compute_array();

 protected:
  double tfactor;

  char *idcoupling;

  // Compute that calculates atom-thermostat couplings
  class ComputeCouplingNHMesh *coupling;

  int n_thermostats;    // Number of thermostats controlling the atoms

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Compute ID for temp/nhmesh does not exist

Self-explanatory.  temp/nhmesh requires a compute that provides
particle-thermostat couplings

E: Invalid coupling compute for temp/nhmesh/profile

Coupling compute must be derived from nhmesh/coupling.

*/
