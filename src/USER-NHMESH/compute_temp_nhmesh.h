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

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempNHMesh : public Compute {
 public:
  ComputeTempNHMesh(class LAMMPS *, int, char **);
  virtual ~ComputeTempNHMesh();
  void init();
  void setup();
  virtual double compute_scalar();
  virtual void compute_vector();
  virtual void compute_array();

 protected:
  double tfactor;

  char *idcoupling;
  Compute *coupling;    // Compute that calculates atom-thermostat couplings

  int n_thermostats;    // Number of thermostats controlling the atoms

  virtual void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.

E: Number of thermostats must be > 0

Self-explanatory.

E: Compute ID for temp/nhmesh does not exist

Self-explanatory.  temp/nhmesh requires a compute that provides
particle-thermostat couplings

E: Compute doesn't output the correct array size for temp/nhmesh

Self-explanatory.

E: Temperature compute degrees of freedom < 0

Self-explanatory.

*/
