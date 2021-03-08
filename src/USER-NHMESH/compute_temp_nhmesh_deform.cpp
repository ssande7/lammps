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

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "compute_temp_nhmesh_deform.h"
#include "compute_coupling_nhmesh_atom.h"

#include <cstring>
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempDeformNHMesh::ComputeTempDeformNHMesh(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempDeform(lmp, narg-1, arg)
{
  idcoupling = utils::strdup(arg[narg-1]);

  int icoupling = modify->find_compute(idcoupling);
  if (icoupling < 0)
    error->all(FLERR,"Compute ID for temp/nhmesh/profile does not exist");
  Compute *comp = modify->compute[icoupling];
  coupling = dynamic_cast<ComputeCouplingNHMesh *>(comp);
  if (coupling == nullptr)
    error->all(FLERR,"Invalid coupling compute for temp/nhmesh/profile");
  n_thermostats = coupling->get_n_thermostats();

  array_flag = 1;
  size_array_rows = n_thermostats;
  size_array_cols = 6;
  extarray = 1;
  memory->create(array,size_array_rows,size_array_cols,
      "temp/nhmesh/deform:array");
}

/* ---------------------------------------------------------------------- */

ComputeTempDeformNHMesh::~ComputeTempDeformNHMesh()
{
  if (!copymode) {
    delete [] idcoupling;
    memory->destroy(array);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformNHMesh::compute_array()
{
  double lamda[3],vstream[3],vthermal[3];

  invoked_vector = update->ntimestep;

  coupling->compute_peratom();
  double **&couple_mat = coupling->array_atom;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  double massone,t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      for (int j = 0; j < n_thermostats; j++) {
        t[0] += couple_mat[i][j] * massone * vthermal[0]*vthermal[0];
        t[1] += couple_mat[i][j] * massone * vthermal[1]*vthermal[1];
        t[2] += couple_mat[i][j] * massone * vthermal[2]*vthermal[2];
        t[3] += couple_mat[i][j] * massone * vthermal[0]*vthermal[1];
        t[4] += couple_mat[i][j] * massone * vthermal[0]*vthermal[2];
        t[5] += couple_mat[i][j] * massone * vthermal[1]*vthermal[2];
      }
    }

  // array[i] = &data[i*n2], n2 = 6, so this reduces the whole array in one call
  MPI_Allreduce(t,array[0],6*n_thermostats,MPI_DOUBLE,MPI_SUM,world);
  for (int j = 0; j < n_thermostats; j++)
    for (int i = 0; i < 6; i++) array[j][i] *= force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeformNHMesh::memory_usage()
{
  double bytes = size_array_rows * size_array_cols * sizeof(double);
  return bytes + ComputeTempDeform::memory_usage();
}
