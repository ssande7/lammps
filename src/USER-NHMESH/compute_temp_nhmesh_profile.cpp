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

#include "compute_temp_nhmesh_profile.h"
#include "compute_coupling_nhmesh_atom.h"

#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "domain.h"
#include "memory.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{TENSOR,BIN};

/* ---------------------------------------------------------------------- */

ComputeTempProfileNHMesh::ComputeTempProfileNHMesh(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempProfile(lmp, narg-1, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal compute temp/nhmesh/profile command");

  scalar_flag = 1;
  extscalar = 0;
  tempflag = 1;
  tempbias = 1;

  idcoupling = utils::strdup(arg[narg-1]);

  int icoupling = modify->find_compute(idcoupling);
  if (icoupling < 0)
    error->all(FLERR,"Compute ID for temp/nhmesh/profile does not exist");
  Compute *comp = modify->compute[icoupling];
  coupling = dynamic_cast<ComputeCouplingNHMesh *>(comp);
  if (coupling == nullptr)
    error->all(FLERR,"Invalid coupling compute for temp/nhmesh/profile");
  n_thermostats = coupling->get_n_thermostats();

  if (outflag == TENSOR) {
    array_flag = 1;
    size_array_rows = n_thermostats;
    size_array_cols = 6;
    extarray = 1;
    memory->create(array,size_array_rows,size_array_cols,
        "temp/nhmesh/profile:array");
    array_compute_fn = &ComputeTempProfileNHMesh::compute_array_ke;
  } else {
    array_compute_fn = &ComputeTempProfileNHMesh::compute_array_bin;
    error->warning(FLERR,"Compute temp/nhmesh/profile can't be used as the "
                      "temperature compute for fix temp/nhmesh with out bin.");
  }
}

/* ---------------------------------------------------------------------- */

ComputeTempProfileNHMesh::~ComputeTempProfileNHMesh()
{
  if (!copymode) {
    delete [] idcoupling;
    if (outflag == TENSOR) memory->destroy(array);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempProfileNHMesh::compute_array()
{
  (this->*array_compute_fn)();
}

/* ---------------------------------------------------------------------- */

void ComputeTempProfileNHMesh::compute_array_ke() {
  int i, j, ibin;

  double vthermal[3];

  invoked_vector = update->ntimestep;

  bin_average();

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
      ibin = bin[i];
      if (xflag) vthermal[0] = v[i][0] - binave[ibin][ivx];
      else vthermal[0] = v[i][0];
      if (yflag) vthermal[1] = v[i][1] - binave[ibin][ivy];
      else vthermal[1] = v[i][1];
      if (zflag) vthermal[2] = v[i][2] - binave[ibin][ivz];
      else vthermal[2] = v[i][2];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      for (j = 0; j < n_thermostats; j++) {
        t[j*6+0] += couple_mat[i][j] * massone * vthermal[0]*vthermal[0];
        t[j*6+1] += couple_mat[i][j] * massone * vthermal[1]*vthermal[1];
        t[j*6+2] += couple_mat[i][j] * massone * vthermal[2]*vthermal[2];
        t[j*6+3] += couple_mat[i][j] * massone * vthermal[0]*vthermal[1];
        t[j*6+4] += couple_mat[i][j] * massone * vthermal[0]*vthermal[2];
        t[j*6+5] += couple_mat[i][j] * massone * vthermal[1]*vthermal[2];
      }
    }

  // array[i] = &data[i*n2], n2 = 6, so this reduces the whole array in one call
  MPI_Allreduce(t,array[0],6*n_thermostats,MPI_DOUBLE,MPI_SUM,world);
  for (j = 0; j < n_thermostats; j++)
    for (i = 0; i < 6; i++) array[j][i] *= force->mvv2e;
}

/* ---------------------------------------------------------------------- */

void ComputeTempProfileNHMesh::compute_array_bin() {
  ComputeTempProfile::compute_array();
}

