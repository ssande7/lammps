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
// * add protected variable to store total DoFs for scalar output

// NOTE:
// * takes compute_coupling_nhmesh or compute_coupling_conserve_nhmesh as input
// * returns:
//    + scalar: average kinetic temperature
//    + vector: vector of length n_thermostats with temperature of each
//    + array:  n_thermostats x 6 array of KE tensors. Summing these would give
//              the KE tensor of the entire group.

#ifdef COMPUTE_CLASS

ComputeStyle(temp,ComputeTemp)

#else

#ifndef LMP_COMPUTE_TEMP_NHMESH_H
#define LMP_COMPUTE_TEMP_NHMESH_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempNHMesh : public Compute {
 public:
  ComputeTempNHMesh(class LAMMPS *, int, char **);
  virtual ~ComputeTempNHMesh();
  void init() {}
  void setup();
  virtual double compute_scalar();
  virtual void compute_vector();

 protected:
  double tfactor;

  virtual void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
