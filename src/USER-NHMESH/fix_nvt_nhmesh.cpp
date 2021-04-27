/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nvt_nhmesh.h"
#include <cstring>

#include "group.h"
#include "modify.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVTNHMesh::FixNVTNHMesh(LAMMPS *lmp, int narg, char **arg) :
  FixNHMesh(lmp, narg, arg)
{
  // Create new temperature compute
  // id = fix-ID + _temp

  std::string tcmd = id + std::string("_temp");
  id_temp = new char[tcmd.size()+1];
  strcpy(id_temp, tcmd.c_str());

  tcmd += fmt::format(" {} temp/nhmesh {}", group->names[igroup], id_coupling);
  modify->add_compute(tcmd);
  tcomputeflag = 1;
}
