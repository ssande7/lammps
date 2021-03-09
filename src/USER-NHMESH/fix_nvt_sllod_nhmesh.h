/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nvt/sllod/nhmesh,FixNVTSllodNHMesh)

#else

#ifndef LMP_FIX_NVT_SLLOD_NHMESH_H
#define LMP_FIX_NVT_SLLOD_NHMESH_H

#include "fix_nhmesh.h"

namespace LAMMPS_NS {

class FixNVTSllodNHMesh : public FixNHMesh {
 public:
  FixNVTSllodNHMesh(class LAMMPS *, int, char **);
  ~FixNVTSllodNHMesh() {}
  void init();

 private:
  int nondeformbias;

  void nhmesh_v_temp();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature for fix nvt/sllod/nhmesh does not have a bias

The specified compute must compute temperature with a bias.

E: Using fix nvt/sllod/nhmesh with inconsistent fix deform remap option

Fix nvt/sllod/nhmesh requires that deforming atoms have a velocity profile
provided by "remap v" as a fix deform option.

E: Using fix nvt/sllod/nhmesh with no fix deform defined

Self-explanatory.

*/
