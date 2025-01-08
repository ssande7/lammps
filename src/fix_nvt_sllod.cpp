// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   LAMMPS development team: developers@lammps.org, Sandia National Laboratories

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
#include "comm.h"
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

/* ---------------------------------------------------------------------- */

FixNVTSllod::FixNVTSllod(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/sllod");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt/sllod");

  // default values

  psllod_flag = 0;
  peculiar_flag = 1;
  if (mtchain_default_flag) mtchain = 1;

  // select SLLOD/p-SLLOD/g-SLLOD variant and velocity frame

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"psllod") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod psllod", error);
      psllod_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"peculiar") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod peculiar", error);
      peculiar_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else iarg++;
  }

  // create a new compute temp style
  // id = fix-ID + temp

  id_temp = utils::strdup(std::string(id) + "_temp");
  if (peculiar_flag) modify->add_compute(fmt::format("{} {} temp",id_temp,group->names[igroup]));
  else modify->add_compute(fmt::format("{} {} temp/deform",id_temp,group->names[igroup]));
  tcomputeflag = 1;
  nondeformbias = peculiar_flag;
}

/* ---------------------------------------------------------------------- */

void FixNVTSllod::init()
{
  FixNH::init();

  if (!peculiar_flag && !temperature->tempbias)
    error->all(FLERR,"Temperature for fix {} does not have a bias", style);

  nondeformbias = 0;
  if (strcmp(temperature->style,"temp/deform") != 0) nondeformbias = 1;

  // check fix deform remap settings

  auto deform = modify->get_fix_by_style("^deform");
  if (deform.size() < 1) error->all(FLERR,"Using fix {} with no fix deform defined", style);

  for (auto &ifix : deform) {
    auto f = dynamic_cast<FixDeform *>(ifix);
    if (f == nullptr) continue;
    if ((peculiar_flag && f->remapflag != Domain::NO_REMAP) ||
        (!peculiar_flag && f->remapflag != Domain::V_REMAP))
      error->all(FLERR,"Using fix {} with inconsistent fix deform remap option", style);

    // error on unsupported mixed flows
    bool elongation = false;
    for (int j = 0; j < 3; ++j) {
      if (f->set[j].style) {
        elongation = true;
        if (f->set[j].style != FixDeform::TRATE)
          error->all(FLERR,"fix {} requires the trate style for x,y,z deformation", style);
      }
    }
    for (int j = 3; j < 6; ++j) {
      if (f->set[j].style && f->set[j].style != FixDeform::ERATE) {
        if (elongation) error->all(FLERR,"fix {} requires the erate style for "
            "xy/xz/yz deformation under mixed shear/extensional flow", style);
        else if (comm->me == 0)
          error->warning(FLERR,"Using non-constant shear rate with fix nvt/sllod");
      }
    }
    if (comm->me == 0) {
      // Warn about fix deform settings that do not produce a constant flow tensor
      // No need to warn about xy + yz shear since this is handled in fix deform
      if (f->set[5].style && f->set[5].rate != 0.0 &&
          (f->set[3].style || domain->yz != 0.0) &&
          (f->set[4].style != FixDeform::ERATE ||
           f->set[5].style != FixDeform::ERATE ||
           (f->set[3].style && f->set[3].style != FixDeform::ERATE)))
        error->warning(FLERR,"Shearing xy with a yz tilt is only handled correctly "
            "if fix deform uses the erate style for xy, xz and yz");
      if (f->end_flag)
        error->warning(FLERR,"fix {} requires box deformation to occur with "
            "position updates to be strictly correct. Set the N parameter of "
            "fix deform to 0 to enable this.", style);

      if (!peculiar_flag)
        error->warning(FLERR,"fix {} will produce incorrect energy dissipation "
            "if the peculiar flag is not set", style);
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

void FixNVTSllod::nh_v_temp()
{
  // thermostat thermal velocity only
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  double **v = atom->v;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // Flow tensor = h_rate * h_inv
  // layout: xx, yy, zz, yz, xz, xy
  double grad_u[6];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,grad_u);

  if (peculiar_flag) {
    // update velocities in a time-reversible manner
    // atom velocities stored in peculiar frame w.r.t. streaming velocity
    // Remove/restore bias only for thermostat step if needed.
    //   Loops need to be split in this case so that timestep stays reversible,
    //   since some temperature computes require v hasn't changed after compute_scalar().
    //   Main use for this is something like a molecular thermostat.
    // Add dt4*p-SLLOD force as separate term so that pure shear is
    // identical between SLLOD and p-SLLOD. dt4*(SLLOD force + p-SLLOD force)
    // causes numerical divergence.

    double dt4 = 0.5*dthalf;
    double vfac[3];
    vfac[0] = exp(-grad_u[0]*dt4);
    vfac[1] = exp(-grad_u[1]*dt4);
    vfac[2] = exp(-grad_u[2]*dt4);
    if (which == FixNH::BIAS) {
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
          if (psllod_flag) {
            v[i][2] -= dt4*grad_u[2]*grad_u[2]*x[i][2];
            v[i][1] -= dt4*grad_u[3]*v[i][2] + dt4*grad_u[1]*grad_u[1]*x[i][1];
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                       + dt4*grad_u[0]*grad_u[0]*x[i][0];
          } else {
            v[i][1] -= dt4*grad_u[3]*v[i][2];
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
          }
        }
      }
      temperature->compute_scalar();
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          temperature->remove_bias(i,v[i]);
          v[i][0] *= factor_eta;
          v[i][1] *= factor_eta;
          v[i][2] *= factor_eta;
          temperature->restore_bias(i,v[i]);
        }
      }
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          if (psllod_flag) {
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                       + dt4*grad_u[0]*grad_u[0]*x[i][0];
            v[i][1] -= dt4*grad_u[3]*v[i][2] + dt4*grad_u[1]*grad_u[1]*x[i][1];
            v[i][2] -= dt4*grad_u[2]*grad_u[2]*x[i][2];
          } else {
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
            v[i][1] -= dt4*grad_u[3]*v[i][2];
          }
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
        }
      }
    } else {
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
          if (psllod_flag) {
            v[i][2] -= dt4*grad_u[2]*grad_u[2]*x[i][2];
            v[i][1] -= dt4*grad_u[3]*v[i][2] + dt4*grad_u[1]*grad_u[1]*x[i][1];
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                       + dt4*grad_u[0]*grad_u[0]*x[i][0];
            v[i][0] *= factor_eta;
            v[i][1] *= factor_eta;
            v[i][2] *= factor_eta;
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                       + dt4*grad_u[0]*grad_u[0]*x[i][0];
            v[i][1] -= dt4*grad_u[3]*v[i][2] + dt4*grad_u[1]*grad_u[1]*x[i][1];
            v[i][2] -= dt4*grad_u[2]*grad_u[2]*x[i][2];
          } else {
            v[i][1] -= dt4*grad_u[3]*v[i][2];
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
            v[i][0] *= factor_eta;
            v[i][1] *= factor_eta;
            v[i][2] *= factor_eta;
            v[i][0] -= dt4*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
            v[i][1] -= dt4*grad_u[3]*v[i][2];
          }
          v[i][0] *= vfac[0];
          v[i][1] *= vfac[1];
          v[i][2] *= vfac[2];
        }
      }
    }
  } else {
    // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
    // vdelu = SLLOD correction   = -vthermal . grad_u
    //    OR = p-SLLOD correction = -vthermal . grad_u - x . grad_u . grad_u
    //                            = -v . grad_u
    // Integrating as below is not reversible!
    if (nondeformbias) temperature->compute_scalar();
    double vdelu[3];
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (!psllod_flag) temperature->remove_bias(i,v[i]);
        vdelu[0] = grad_u[0]*v[i][0] + grad_u[5]*v[i][1] + grad_u[4]*v[i][2];
        vdelu[1] = grad_u[1]*v[i][1] + grad_u[3]*v[i][2];
        vdelu[2] = grad_u[2]*v[i][2];
        if (psllod_flag) temperature->remove_bias(i,v[i]);
        v[i][0] = v[i][0]*factor_eta - dthalf*vdelu[0];
        v[i][1] = v[i][1]*factor_eta - dthalf*vdelu[1];
        v[i][2] = v[i][2]*factor_eta - dthalf*vdelu[2];
        temperature->restore_bias(i,v[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions with streaming velocity
-----------------------------------------------------------------------*/

void FixNVTSllod::nve_x()
{
  if (!peculiar_flag) return FixNH::nve_x();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x update by full step only for atoms in group
  // identical for SLLOD and p-SLLOD
  // velocity stored as peculiar, so need to add in streaming velocity here

  double dtv2 = dtv*0.5;
  double grad_u[6], xfac[3];
  MathExtra::multiply_shape_shape(domain->h_rate, domain->h_inv, grad_u);
  xfac[0] = exp(grad_u[0]*dtv2);
  xfac[1] = exp(grad_u[1]*dtv2);
  xfac[2] = exp(grad_u[2]*dtv2);

  // Fix deform keeps the box centre fixed under elongation,
  // and the lower corner fixed under shear, so adjust for that
  // to avoid an apparent drift relative to the box and prevent
  // extra atom exchanges between MPI ranks
  double xmid[3];
  for (int i = 0; i < 3; ++i) {
    xmid[i] = (domain->boxhi[i] + domain->boxlo[i])/2.;
  }
  double *xlo = domain->boxlo;

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      x[i][0] = xmid[0] + (x[i][0] - xmid[0])*xfac[0];
      x[i][1] = xmid[1] + (x[i][1] - xmid[1])*xfac[1];
      x[i][2] = xmid[2] + (x[i][2] - xmid[2])*xfac[2];
      x[i][1] += dtv2 * grad_u[3]*(x[i][2] - xlo[2]);
      x[i][0] += dtv2 * (grad_u[5]*(x[i][1] - xlo[1]) + grad_u[4]*(x[i][2] - xlo[2]));
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
      x[i][0] += dtv2 * (grad_u[5]*(x[i][1] - xlo[1]) + grad_u[4]*(x[i][2] - xlo[2]));
      x[i][1] += dtv2 * grad_u[3]*(x[i][2] - xlo[2]);
      x[i][0] = xmid[0] + (x[i][0] - xmid[0])*xfac[0];
      x[i][1] = xmid[1] + (x[i][1] - xmid[1])*xfac[1];
      x[i][2] = xmid[2] + (x[i][2] - xmid[2])*xfac[2];
    }
  }
}
