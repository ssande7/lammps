/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "fix_nvt_sllod_nhmesh.h"
#include "compute_coupling_nhmesh_atom.h"

#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "compute.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVTSllodNHMesh::FixNVTSllodNHMesh(LAMMPS *lmp, int narg, char **arg) :
  FixNHMesh(lmp, narg, arg)
{

  // create a new compute temp style
  // id = fix-ID + temp

  std::string cmd = id + std::string("_temp");
  id_temp = new char[cmd.size()+1];
  strcpy(id_temp,cmd.c_str());

  cmd+=fmt::format(" {} temp/nhmesh/deform {}",group->names[igroup],id_coupling);
  modify->add_compute(cmd);
  tcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixNVTSllodNHMesh::init()
{
  FixNHMesh::init();

  if (!temperature->tempbias)
    error->all(FLERR,"Temperature for fix nhmesh/nvt/sllod does not have a bias");

  nondeformbias = 0;
  if (strcmp(temperature->style,"temp/nhmesh/deform") != 0) nondeformbias = 1;

  // check fix deform remap settings

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"deform",6) == 0) {
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using fix nhmesh/nvt/sllod with inconsistent fix "
                   "deform remap option");
      break;
    }
  if (i == modify->nfix)
    error->all(FLERR,"Using fix nhmesh/nvt/sllod with no fix deform defined");
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

void FixNVTSllodNHMesh::nh_v_temp()
{
  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  if (nondeformbias) temperature->compute_scalar();

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ke_local[n_thermostats];
  double fac, massone;
  for (int i = 0; i < n_thermostats; i++) ke_local[i] = 0.0;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double h_two[6],vdelu[3];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fac = 0;
      for (int j = 0; j < n_thermostats; j++)
        fac += coupling->array_atom[i][j] * factor_eta[j];
      fac = exp(fac);
      vdelu[0] = h_two[0]*v[i][0] + h_two[5]*v[i][1] + h_two[4]*v[i][2];
      vdelu[1] = h_two[1]*v[i][1] + h_two[3]*v[i][2];
      vdelu[2] = h_two[2]*v[i][2];
      temperature->remove_bias(i,v[i]);
      v[i][0] = v[i][0]*fac - dthalf*vdelu[0];
      v[i][1] = v[i][1]*fac - dthalf*vdelu[1];
      v[i][2] = v[i][2]*fac - dthalf*vdelu[2];
      for (int j = 0; j < n_thermostats; j++)
        ke_local[j] += coupling->array_atom[i][j] * massone *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      temperature->restore_bias(i,v[i]);
    }
  }
  MPI_Allreduce(ke_local,ke_current,n_thermostats,MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < n_thermostats; i++) {
    ke_current[i] *= force->mvv2e;
    if (mesh_coupling_flag)
      for (int j = 0; j < n_thermostats; j++)
        ke_current[i] += mesh_coupling[i][j]*eta_mass[j]*eta_dot[j]*eta_dot[j];
  }
}
