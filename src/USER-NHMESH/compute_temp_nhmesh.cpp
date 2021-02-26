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
  Compute(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal compute temp/nhmesh command");

  n_thermostats = utils::inumeric(FLERR,arg[3],false,lmp);
  if (n_thermostats <= 0) error->all(FLERR,"Number of thermostats must be > 0");

  idcoupling = utils::strdup(arg[4]);

  scalar_flag = vector_flag = array_flag = 1;
  size_vector = 6;
  size_array_cols = n_thermostats;
  size_array_rows = 6;
  size_array_rows_variable = 1;
  extscalar = 0;
  extvector = 1;
  extarray = 1;
  tempflag = 1;

  vector = new double[size_vector];
  memory->create(array,size_array_rows,size_array_cols,"temp/nhmesh:array");
}

/* ---------------------------------------------------------------------- */

ComputeTempNHMesh::~ComputeTempNHMesh()
{
  if (!copymode) {
    delete [] vector;
    delete [] idcoupling;
    memory->destroy(array);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempNHMesh::init()
{
  // get pointer to coupling compute

  int icoupling = modify->find_compute(idcoupling);
  if (icoupling == -1)
    error->all(FLERR,"Compute ID for temp/nhmesh does not exist");
  coupling = modify->compute[icoupling];
  if (coupling->size_peratom_cols != n_thermostats)
    error->all(FLERR,
        "Compute doesn't output the correct array size for temp/nhmesh");
}

/* ---------------------------------------------------------------------- */

void ComputeTempNHMesh::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempNHMesh::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempNHMesh::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  // assumes sum of couplings for each atom is 1
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempNHMesh::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  // assumes sum of couplings for each atom is 1
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
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
