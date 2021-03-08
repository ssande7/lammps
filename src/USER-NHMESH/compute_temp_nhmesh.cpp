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

#include "compute_temp_nhmesh.h"
#include "compute_coupling_nhmesh_atom.h"

#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempNHMesh::ComputeTempNHMesh(LAMMPS *lmp, int narg, char **arg) :
  ComputeTemp(lmp, narg-1, arg)
{
  idcoupling = utils::strdup(arg[narg-1]);

  int icoupling = modify->find_compute(idcoupling);
  if (icoupling < 0)
    error->all(FLERR,"Compute ID for temp/nhmesh does not exist");
  Compute *comp = modify->compute[icoupling];
  coupling = dynamic_cast<ComputeCouplingNHMesh *>(comp);
  if (coupling == nullptr)
    error->all(FLERR,"Invalid coupling compute for temp/nhmesh");
  n_thermostats = coupling->get_n_thermostats();

  array_flag = 1;
  size_array_cols = 6;
  size_array_rows = n_thermostats;
  extarray = 1;

  memory->create(array,size_array_rows,size_array_cols,"temp/nhmesh:array");
}

/* ---------------------------------------------------------------------- */

ComputeTempNHMesh::~ComputeTempNHMesh()
{
  if (!copymode) {
    delete [] idcoupling;
    memory->destroy(array);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempNHMesh::compute_array()
{
  int i, j;

  invoked_vector = update->ntimestep;

  // TODO: check that this only recalculates once per time step
  coupling->compute_peratom();
  double **&couple_mat = coupling->array_atom;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[n_thermostats*6];
  for (j = 0; j < n_thermostats; j++)
    for (i = 0; i < 6; i++) t[j*6+i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      for (j = 0; j < n_thermostats; j++) {
        t[j*6+0] += couple_mat[i][j] * massone * v[i][0]*v[i][0];
        t[j*6+1] += couple_mat[i][j] * massone * v[i][1]*v[i][1];
        t[j*6+2] += couple_mat[i][j] * massone * v[i][2]*v[i][2];
        t[j*6+3] += couple_mat[i][j] * massone * v[i][0]*v[i][1];
        t[j*6+4] += couple_mat[i][j] * massone * v[i][0]*v[i][2];
        t[j*6+5] += couple_mat[i][j] * massone * v[i][1]*v[i][2];
      }
    }

  // array[i] = &data[i*n2], n2 = 6, so this reduces the whole array in one call
  MPI_Allreduce(t,array[0],6*n_thermostats,MPI_DOUBLE,MPI_SUM,world);
  for (j = 0; j < n_thermostats; j++)
    for (i = 0; i < 6; i++) array[j][i] *= force->mvv2e;
}
