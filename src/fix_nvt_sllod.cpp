// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
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

#include "fix_nvt_sllod.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// from FixNH:
enum{NOBIAS,BIAS};

// from FixDeform:
enum{NONE=0,FINAL,DELTA,SCALE,VEL,ERATE,TRATE,VOLUME,WIGGLE,VARIABLE};
/* ---------------------------------------------------------------------- */

FixNVTSllod::FixNVTSllod(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/sllod");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt/sllod");

  // default values
  p_sllod = 1;
  peculiar = 0;

  if (mtchain_default_flag) mtchain = 1;

  for (int i = 0; i < narg; ++i) {
    if (strcmp(arg[i], "p_sllod")==0) {
      if (++i >= narg) error->all(FLERR, "Illegal fix nvt/sllod command");
      p_sllod = utils::logical(FLERR, arg[i], false, lmp);
    } else if (strcmp(arg[i], "peculiar")==0) {
      if (++i >= narg) error->all(FLERR, "Illegal fix nvt/sllod command");
      peculiar = utils::logical(FLERR, arg[i], false, lmp);
    }
  }

  // create a new compute temp style
  // id = fix-ID + temp

  id_temp = utils::strdup(std::string(id) + "_temp");
  if (peculiar) modify->add_compute(fmt::format("{} {} temp",
                                  id_temp,group->names[igroup]));
  else modify->add_compute(fmt::format("{} {} temp/deform",
                                  id_temp,group->names[igroup]));

  tcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixNVTSllod::init()
{
  FixNH::init();

  if (!peculiar && !temperature->tempbias)
    error->all(FLERR,"Temperature for fix nvt/sllod does not have a bias");

  nondeformbias = 0;
  if (strcmp(temperature->style,"temp/deform") != 0) nondeformbias = 1;

  // check fix deform remap settings

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"deform",6) == 0) {
      auto def = dynamic_cast<FixDeform *>(modify->fix[i]);
      if (!peculiar && def->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using fix nvt/sllod with inconsistent fix deform "
                   "remap option");
      if (peculiar && def->remapflag != Domain::NO_REMAP)
        error->all(FLERR,"Using fix nvt/sllod with inconsistent fix deform "
                   "remap option");
      bool elongation = false;
      for (int j = 0; j < 3; ++j) {
        if (def->set[i].style) {
          elongation = true;
          if (def->set[j].style != TRATE)
            error->all(FLERR,"fix nvt/sllod requires the trate style for x/y/z deformation");
        }
      }
      for (int j = 3; j < 6; ++j) {
        if (def->set[j].style && def->set[j].style != ERATE) {
          if (elongation)
            error->all(FLERR,"fix nvt/sllod requires the erate style for xy/xz/yz deformation under mixed shear/extensional flow");
          else
            error->warning(FLERR, "Using non-constant shear rate with fix nvt/sllod");
        }
      }
      break;
    }
  if (i == modify->nfix)
    error->all(FLERR,"Using fix nvt/sllod with no fix deform defined");

  // Apply initial kick if we can
  if (!peculiar && !nondeformbias) {
    temperature->remove_bias_all();
    temperature->restore_bias_all();
    temperature->restore_bias_all();
  }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

void FixNVTSllod::nh_v_temp()
{
  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  if (which == BIAS) temperature->compute_scalar();

  double **v = atom->v;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double grad_u[6],vdelu[3];
  double* h_rate = domain->h_rate;
  double* h = domain->h;
  grad_u[0] = h_rate[0]/h[0];
  grad_u[1] = h_rate[1]/h[1];
  grad_u[2] = h_rate[2]/h[2];
  grad_u[3] = (h_rate[3] - grad_u[1]*h[3])/h[2];
  grad_u[4] = (h_rate[4] - grad_u[0]*h[4])/h[2];
  grad_u[5] = (h_rate[5] - grad_u[0]*h[5])/h[1];

  if (peculiar) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (which == BIAS) temperature->remove_bias(i,v[i]);
        v[i][0] = v[i][0]*factor_eta;
        v[i][1] = v[i][1]*factor_eta;
        v[i][2] = v[i][2]*factor_eta;
        if (which == BIAS) temperature->restore_bias(i,v[i]);
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (!p_sllod) temperature->remove_bias(i,v[i]);
        vdelu[0] = grad_u[0]*v[i][0] + grad_u[5]*v[i][1] + grad_u[4]*v[i][2];
        vdelu[1] = grad_u[1]*v[i][1] + grad_u[3]*v[i][2];
        vdelu[2] = grad_u[2]*v[i][2];
        if (p_sllod) temperature->remove_bias(i,v[i]);
        v[i][0] = v[i][0]*factor_eta - dthalf*vdelu[0];
        v[i][1] = v[i][1]*factor_eta - dthalf*vdelu[1];
        v[i][2] = v[i][2]*factor_eta - dthalf*vdelu[2];
        temperature->restore_bias(i,v[i]);
      }
    }
  }
}

void FixNVTSllod::nve_v()
{
  double dtfm, dtf2;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double grad_u[6], vfac[3];
  double* h_rate = domain->h_rate;
  double* h = domain->h;
  grad_u[0] = h_rate[0]/h[0];
  grad_u[1] = h_rate[1]/h[1];
  grad_u[2] = h_rate[2]/h[2];
  grad_u[3] = (h_rate[3] - grad_u[1]*h[3])/h[2];
  grad_u[4] = (h_rate[4] - grad_u[0]*h[4])/h[2];
  grad_u[5] = (h_rate[5] - grad_u[0]*h[5])/h[1];

  if (peculiar) {
    dtf2 = 0.5*dtf;
    vfac[0] = exp(-grad_u[0]*dtf2);
    vfac[1] = exp(-grad_u[1]*dtf2);
    vfac[2] = exp(-grad_u[2]*dtf2);
  }
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) dtfm = dtf / rmass[i];
      else dtfm = dtf / mass[type[i]];

      if (peculiar) {
        if (which == BIAS) temperature->remove_bias(i,v[i]);
        if (p_sllod) {
          // Add dtf2*p-SLLOD force separately so that pure shear is identical
          // between SLLOD and p-SLLOD. Using dtf2*(SLLOD_force + p-SLLOD_force)
          // causes numerical divergence.
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
          v[i][2] -= dtf2*grad_u[2]*grad_u[2]*x[i][2];
          v[i][1] -= dtf2*grad_u[3]*v[i][2] + dtf2*grad_u[1]*grad_u[1]*x[i][1];
          v[i][0] -= dtf2*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                     + dtf2*grad_u[0]*grad_u[0]*x[i][0];
          v[i][0] += dtfm*f[i][0];
          v[i][1] += dtfm*f[i][1];
          v[i][2] += dtfm*f[i][2];
          v[i][0] -= dtf2*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                     + dtf2*grad_u[0]*grad_u[0]*x[i][0];
          v[i][1] -= dtf2*grad_u[3]*v[i][2] + dtf2*grad_u[1]*grad_u[1]*x[i][1];
          v[i][2] -= dtf2*grad_u[2]*grad_u[2]*x[i][2];
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
        } else {
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
          v[i][1] -= dtf2*grad_u[3]*v[i][2];
          v[i][0] -= dtf2*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
          v[i][0] += dtfm*f[i][0];
          v[i][1] += dtfm*f[i][1];
          v[i][2] += dtfm*f[i][2];
          v[i][0] -= dtf2*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
          v[i][1] -= dtf2*grad_u[3]*v[i][2];
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];

          // Exact ODE solution for mixed flow, if e_aa != e_bb, e_aa != 0 and e_bb != 0.
          // v[i][0] = v[i][2]*grad_u[5]*grad_u[3]/(grad_u[2]-grad_u[1])
          //         * ((vfac[2]-vfac[0])/(grad_u[2]-grad_u[0]) - (vfac[1]-vfac[0])/(grad_u[1]-grad_u[0]))
          //         + v[i][1]*(vfac[1]-vfac[0])*grad_u[5]/(grad_u[1]-grad_u[0])
          //         + v[i][2]*(vfac[2]-vfac[0])*grad_u[4]/(grad_u[2]-grad_u[0])
          //         + (v[i][0]+f[i][0]*dtfm)*vfac[0];
          // v[i][1] = v[i][2]*(vfac[2]-vfac[1])*grad_u[3]/(grad_u[2]-grad_u[1])
          //         + (v[i][1]+f[i][1]*dtfm)*vfac[1];
          // v[i][2] = (v[i][2]+f[i][2]*dtfm)*vfac[2];
        }
        if (which == BIAS) temperature->restore_bias(i,v[i]);
      } else {
        v[i][0] += dtfm*f[i][0];
        v[i][1] += dtfm*f[i][1];
        v[i][2] += dtfm*f[i][2];
      }
    }
  }
}


/* ----------------------------------------------------------------------
   perform full-step update of positions
-----------------------------------------------------------------------*/

void FixNVTSllod::nve_x()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  double grad_u[6], xfac[3];
  double dtv2 = dtv*0.5;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x update by full step only for atoms in group
  
  if (peculiar) {
    double* h_rate = domain->h_rate;
    double* h = domain->h;
    grad_u[0] = h_rate[0]/h[0];
    grad_u[1] = h_rate[1]/h[1];
    grad_u[2] = h_rate[2]/h[2];
    grad_u[3] = (h_rate[3] - grad_u[1]*h[3])/h[2];
    grad_u[4] = (h_rate[4] - grad_u[0]*h[4])/h[2];
    grad_u[5] = (h_rate[5] - grad_u[0]*h[5])/h[1];
    xfac[0] = exp(grad_u[0]*dtv2);
    xfac[1] = exp(grad_u[1]*dtv2);
    xfac[2] = exp(grad_u[2]*dtv2);
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (peculiar) {
        x[i][0] *= xfac[0];
        x[i][1] *= xfac[1];
        x[i][2] *= xfac[2];
        x[i][1] += dtv2 * grad_u[3]*x[i][2];
        x[i][0] += dtv2 * (grad_u[5]*x[i][1] + grad_u[4]*x[i][2]);
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        x[i][0] += dtv2 * (grad_u[5]*x[i][1] + grad_u[4]*x[i][2]);
        x[i][1] += dtv2 * grad_u[3]*x[i][2];
        x[i][0] *= xfac[0];
        x[i][1] *= xfac[1];
        x[i][2] *= xfac[2];

        // Exact ODE solution for mixed flow, if e_aa != e_bb, e_aa != 0 and e_bb != 0.
        // x[i][0] = x[i][2]*grad_u[5]*grad_u[3]/(grad_u[2]-grad_u[1])
        //         * ((xfac[2]-xfac[0])/(grad_u[2]-grad_u[0]) - (xfac[1]-xfac[0])/(grad_u[1]-grad_u[0]))
        //         + x[i][1]*(xfac[1]-xfac[0])*grad_u[5]/(grad_u[1]-grad_u[0])
        //         + x[i][2]*(xfac[2]-xfac[0])*grad_u[4]/(grad_u[2]-grad_u[0])
        //         + (x[i][0]+v[i][0]*dtv)*xfac[0];
        // x[i][1] = x[i][2]*(xfac[2]-xfac[1])*grad_u[3]/(grad_u[2]-grad_u[1])
        //         + (x[i][1]+v[i][1]*dtv)*xfac[1];
        // x[i][2] = (x[i][2]+v[i][2]*dtv)*xfac[2];

      } else {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }
}

