// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "fix_deform.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "irregular.h"
#include "kspace.h"
#include "lattice.h"
#include "math_const.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDeform::FixDeform(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
irregular(nullptr), set(nullptr)
{
  const std::string thiscmd = fmt::format("fix {}", style);
  if (narg < 4) utils::missing_cmd_args(FLERR, thiscmd, error);

  no_change_box = 1;
  restart_global = 1;
  pre_exchange_migrate = 1;
  end_flag = 1;
  need_flip_change = 0;
  allow_flip_change = 1;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) error->all(FLERR, "Fix {} Nevery must be >= 0", style);
  else if (nevery == 0) end_flag = 0;

  // arguments for child classes

  std::unordered_set<std::string> child_parameters;
  std::unordered_map<std::string, int> child_styles;
  int nskip;
  if (utils::strmatch(style, "^deform/pressure")) {
    child_parameters.insert("box");
    child_styles.insert({{"pressure", 4}, {"pressure/mean", 4}, {"erate/rescale", 3}, {"volume", 2}});
  }

  // set defaults

  set = new Set[6];
  memset(set, 0, 6 * sizeof(Set));

  // parse all parameter/style arguments for this parent and also child classes
  // for child classes, simply store them in leftover_iarg and skip over them

  triclinic = domain->triclinic;

  int index;
  int iarg = 4;

  while (iarg < narg) {
    if ((strcmp(arg[iarg], "x") == 0)
        || (strcmp(arg[iarg], "y") == 0)
        || (strcmp(arg[iarg], "z") == 0)) {

      if (strcmp(arg[iarg], "x") == 0) index = 0;
      else if (strcmp(arg[iarg], "y") == 0) index = 1;
      else if (strcmp(arg[iarg], "z") == 0) index = 2;

      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, thiscmd, error);
      if (strcmp(arg[iarg + 1], "final") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, thiscmd + " final", error);
        set[index].style = FINAL;
        set[index].flo = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].fhi = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "delta") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, thiscmd + " delta", error);
        set[index].style = DELTA;
        set[index].dlo = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].dhi = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "scale") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " scale", error);
        set[index].style = SCALE;
        set[index].scale = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "vel") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " vel", error);
        set[index].style = VEL;
        set[index].vel = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "erate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " erate", error);
        set[index].style = ERATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "trate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " trate", error);
        set[index].style = TRATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "volume") == 0) {
        set[index].style = VOLUME;
        iarg += 2;
      } else if (strcmp(arg[iarg + 1], "wiggle") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, thiscmd + " wiggle", error);
        set[index].style = WIGGLE;
        set[index].amplitude = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].tperiod = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].tperiod <= 0.0)
          error->all(FLERR, "Illegal fix {} wiggle period, must be positive", style);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "variable") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, thiscmd + " variable", error);
        set[index].style = VARIABLE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2])
          error->all(FLERR, "Illegal fix {} variable name {}", style, arg[iarg + 2]);
        if (strstr(arg[iarg + 3], "v_") != arg[iarg + 3])
          error->all(FLERR, "Illegal fix {} variable name {}", style, arg[iarg + 3]);
        delete[] set[index].hstr;
        delete[] set[index].hratestr;
        set[index].hstr = utils::strdup(&arg[iarg + 2][2]);
        set[index].hratestr = utils::strdup(&arg[iarg + 3][2]);
        iarg += 4;
      } else if (child_styles.find(arg[iarg + 1]) != child_styles.end()) {
        nskip = child_styles[arg[iarg + 1]];
        if (iarg + nskip > narg)
          utils::missing_cmd_args(FLERR, fmt::format("fix {} {}", style, arg[iarg + 1]), error);
        for (int i = 0; i < nskip; i++) leftover_iarg.push_back(iarg + i);
        iarg += nskip;
      } else error->all(FLERR, "Illegal fix {} command argument: {}", style, arg[iarg + 1]);

    } else if ((strcmp(arg[iarg], "xy") == 0)
               || (strcmp(arg[iarg], "xz") == 0)
               || (strcmp(arg[iarg], "yz") == 0)) {

      if (triclinic == 0) error->all(FLERR,"Fix {} tilt factors require triclinic box", style);
      if (strcmp(arg[iarg], "xy") == 0) index = 5;
      else if (strcmp(arg[iarg], "xz") == 0) index = 4;
      else if (strcmp(arg[iarg], "yz") == 0) index = 3;

      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, thiscmd, error);
      if (strcmp(arg[iarg + 1], "final") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " final", error);
        set[index].style = FINAL;
        set[index].ftilt = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "delta") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " delta", error);
        set[index].style = DELTA;
        set[index].dtilt = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "vel") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " vel", error);
        set[index].style = VEL;
        set[index].vel = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "erate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " erate", error);
        set[index].style = ERATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "trate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, thiscmd + " trate", error);
        set[index].style = TRATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "wiggle") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, thiscmd + " wiggle", error);
        set[index].style = WIGGLE;
        set[index].amplitude = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].tperiod = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].tperiod <= 0.0)
          error->all(FLERR, "Illegal fix {} wiggle period, must be positive", style);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "variable") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, thiscmd + " variable", error);
        set[index].style = VARIABLE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2])
          error->all(FLERR, "Illegal fix {} variable name {}", style, arg[iarg + 2]);
        if (strstr(arg[iarg + 3], "v_") != arg[iarg + 3])
          error->all(FLERR, "Illegal fix {} variable name {}", style, arg[iarg + 3]);
        delete[] set[index].hstr;
        delete[] set[index].hratestr;
        set[index].hstr = utils::strdup(&arg[iarg + 2][2]);
        set[index].hratestr = utils::strdup(&arg[iarg + 3][2]);
        iarg += 4;
      } else if (child_styles.find(arg[iarg + 1]) != child_styles.end()) {
        nskip = child_styles[arg[iarg + 1]];
        if (iarg + nskip > narg)
         utils::missing_cmd_args(FLERR, fmt::format("fix {} {}", style, arg[iarg + 1]), error);
        for (int i = 0; i < nskip; i++) leftover_iarg.push_back(iarg + i);
        iarg += nskip;
      } else error->all(FLERR, "Illegal fix {} command argument: {}", style, arg[iarg + 1]);
    } else if (child_parameters.find(arg[iarg]) != child_parameters.end()) {
      if (child_styles.find(arg[iarg + 1]) != child_styles.end()) {
        nskip = child_styles[arg[iarg + 1]];
        if (iarg + nskip > narg)
         utils::missing_cmd_args(FLERR, fmt::format("fix {} {}", style, arg[iarg + 1]), error);
        for (int i = 0; i < nskip; i++) leftover_iarg.push_back(iarg + i);
        iarg += nskip;
      } else error->all(FLERR, "Illegal fix {} command argument: {}", style, arg[iarg + 1]);
    } else break;
  }

  // read options from end of input line

  iarg_options_start = iarg;
  options(narg - iarg, &arg[iarg]);

  // no x remap effectively moves atoms within box, so set restart_pbc

  if (remapflag != Domain::X_REMAP) restart_pbc = 1;

  // setup dimflags used by other classes to check for volume-change conflicts

  for (int i = 0; i < 6; i++)
    if (set[i].style == NONE) dimflag[i] = 0;
    else dimflag[i] = 1;

  if (dimflag[0]) box_change |= BOX_CHANGE_X;
  if (dimflag[1]) box_change |= BOX_CHANGE_Y;
  if (dimflag[2]) box_change |= BOX_CHANGE_Z;
  if (dimflag[3]) box_change |= BOX_CHANGE_YZ;
  if (dimflag[4]) box_change |= BOX_CHANGE_XZ;
  if (dimflag[5]) box_change |= BOX_CHANGE_XY;

  // no tensile deformation on shrink-wrapped dims
  // b/c shrink wrap will change box-length

  for (int i = 0; i < 3; i++)
    if (set[i].style && (domain->boundary[i][0] >= 2 || domain->boundary[i][1] >= 2))
      error->all(FLERR, "Cannot use fix {} on a shrink-wrapped boundary", style);

  // no tilt deformation on shrink-wrapped 2nd dim
  // b/c shrink wrap will change tilt factor in domain::reset_box()

  if (set[3].style && (domain->boundary[2][0] >= 2 || domain->boundary[2][1] >= 2))
    error->all(FLERR, "Cannot use fix {} tilt on a shrink-wrapped 2nd dim", style);
  if (set[4].style && (domain->boundary[2][0] >= 2 || domain->boundary[2][1] >= 2))
    error->all(FLERR, "Cannot use fix {} tilt on a shrink-wrapped 2nd dim", style);
  if (set[5].style && (domain->boundary[1][0] >= 2 || domain->boundary[1][1] >= 2))
    error->all(FLERR, "Cannot use fix {} tilt on a shrink-wrapped 2nd dim", style);

  // apply scaling to FINAL,DELTA,VEL,WIGGLE since they have dist/vel units

  int flag = 0;
  for (int i = 0; i < 6; i++)
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == VEL || set[i].style == WIGGLE) flag = 1;

  double xscale, yscale, zscale;
  if (flag && scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // for 3,4,5: scaling is in 1st dimension, e.g. x for xz

  double map[6];
  map[0] = xscale; map[1] = yscale; map[2] = zscale;
  map[3] = yscale; map[4] = xscale; map[5] = xscale;

  for (int i = 0; i < 3; i++) {
    if (set[i].style == FINAL) {
      set[i].flo *= map[i];
      set[i].fhi *= map[i];
    } else if (set[i].style == DELTA) {
      set[i].dlo *= map[i];
      set[i].dhi *= map[i];
    } else if (set[i].style == VEL) {
      set[i].vel *= map[i];
    } else if (set[i].style == WIGGLE) {
      set[i].amplitude *= map[i];
    }
  }

  for (int i = 3; i < 6; i++) {
    if (set[i].style == FINAL) set[i].ftilt *= map[i];
    else if (set[i].style == DELTA) set[i].dtilt *= map[i];
    else if (set[i].style == VEL) set[i].vel *= map[i];
    else if (set[i].style == WIGGLE) set[i].amplitude *= map[i];
  }

  // for VOLUME, setup links to other dims
  // fixed, dynamic1, dynamic2
  // only check for parent, otherwise child will check

  if (strcmp(style, "deform") == 0) {
    for (int i = 0; i < 3; i++) {
      if (set[i].style != VOLUME) continue;
      int other1 = (i + 1) % 3;
      int other2 = (i + 2) % 3;

      // Cannot use VOLUME option without at least one deformed dimension
      if (set[other1].style == NONE || set[other1].style == VOLUME)
        if (set[other2].style == NONE || set[other2].style == VOLUME)
          error->all(FLERR, "Fix {} volume setting is invalid", style);

      if (set[other1].style == NONE) {
        set[i].substyle = ONE_FROM_ONE;
        set[i].fixed = other1;
        set[i].dynamic1 = other2;
      } else if (set[other2].style == NONE) {
        set[i].substyle = ONE_FROM_ONE;
        set[i].fixed = other2;
        set[i].dynamic1 = other1;
      } else if (set[other1].style == VOLUME) {
        set[i].substyle = TWO_FROM_ONE;
        set[i].fixed = other1;
        set[i].dynamic1 = other2;
      } else if (set[other2].style == VOLUME) {
        set[i].substyle = TWO_FROM_ONE;
        set[i].fixed = other2;
        set[i].dynamic1 = other1;
      } else {
        set[i].substyle = ONE_FROM_TWO;
        set[i].dynamic1 = other1;
        set[i].dynamic2 = other2;
      }
    }
  }

  // set varflag

  varflag = 0;
  for (int i = 0; i < 6; i++)
    if (set[i].style == VARIABLE) varflag = 1;

  // set initial values at time fix deform is issued

  for (int i = 0; i < 3; i++) {
    set[i].lo_initial = domain->boxlo[i];
    set[i].hi_initial = domain->boxhi[i];
    set[i].vol_initial = domain->xprd * domain->yprd * domain->zprd;
  }
  set[3].tilt_initial = domain->yz;
  set[4].tilt_initial = domain->xz;
  set[5].tilt_initial = domain->xy;

  // reneighboring only forced if flips can occur due to shape changes

  if (flipflag && (set[3].style || set[4].style || set[5].style))
    force_reneighbor = 1;
  next_reneighbor = -1;

  flip = 0;

  if (force_reneighbor) irregular = new Irregular(lmp);
  else irregular = nullptr;
}

/* ---------------------------------------------------------------------- */

FixDeform::~FixDeform()
{
  if (set) {
    for (int i = 0; i < 6; i++) {
      delete[] set[i].hstr;
      delete[] set[i].hratestr;
    }
  }
  delete[] set;

  delete irregular;

  // reset domain's h_rate = 0.0, since this fix may have made it non-zero

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  h_rate[0] = h_rate[1] = h_rate[2] =
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixDeform::setmask()
{
  int mask = 0;
  if (force_reneighbor) mask |= PRE_EXCHANGE;
  if (end_flag) mask |= END_OF_STEP;
  else mask |= POST_INTEGRATE | POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeform::init()
{
  // error if more than one fix deform
  // domain, fix nvt/sllod, compute temp/deform only work on single h_rate

  if (modify->get_fix_by_style("deform").size() > 1)
    error->all(FLERR, "More than one fix deform");

  // Kspace setting

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // elapsed time for entire simulation, including multiple runs if defined

  double delt = (update->endstep - update->beginstep) * update->dt;

  // check variables for VARIABLE style

  for (int i = 0; i < 6; i++) {
    if (set[i].style != VARIABLE) continue;
    set[i].hvar = input->variable->find(set[i].hstr);
    if (set[i].hvar < 0)
      error->all(FLERR, "Variable name {} for fix {} does not exist", set[i].hstr, style);
    if (!input->variable->equalstyle(set[i].hvar))
      error->all(FLERR, "Variable {} for fix {} is invalid style", set[i].hstr, style);
    set[i].hratevar = input->variable->find(set[i].hratestr);
    if (set[i].hratevar < 0)
      error->all(FLERR, "Variable name {} for fix {} does not exist", set[i].hratestr, style);
    if (!input->variable->equalstyle(set[i].hratevar))
      error->all(FLERR, "Variable {} for fix {} is invalid style", set[i].hratestr, style);
  }

  // set start/stop values for box size and shape
  // if single run, start is current values
  // if multiple runs enabled via run start/stop settings,
  //   start is value when fix deform was issued
  // if VARIABLE or NONE, no need to set

  for (int i = 0; i < 3; i++) {
    if (update->firststep == update->beginstep) {
      set[i].lo_start = domain->boxlo[i];
      set[i].hi_start = domain->boxhi[i];
      set[i].vol_start = domain->xprd * domain->yprd * domain->zprd;
    } else {
      set[i].lo_start = set[i].lo_initial;
      set[i].hi_start = set[i].hi_initial;
      set[i].vol_start = set[i].vol_initial;
    }

    if (set[i].style == FINAL) {
      set[i].lo_stop = set[i].flo;
      set[i].hi_stop = set[i].fhi;
    } else if (set[i].style == DELTA) {
      set[i].lo_stop = set[i].lo_start + set[i].dlo;
      set[i].hi_stop = set[i].hi_start + set[i].dhi;
    } else if (set[i].style == SCALE) {
      double shift = 0.5 * set[i].scale * (set[i].hi_start - set[i].lo_start);
      set[i].lo_stop = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_stop = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
    } else if (set[i].style == VEL) {
      set[i].lo_stop = set[i].lo_start - 0.5 * delt * set[i].vel;
      set[i].hi_stop = set[i].hi_start + 0.5 * delt * set[i].vel;
    } else if (set[i].style == ERATE) {
      double shift = 0.5 * delt * set[i].rate * (set[i].hi_start - set[i].lo_start);
      set[i].lo_stop = set[i].lo_start - shift;
      set[i].hi_stop = set[i].hi_start + shift;
      if (set[i].hi_stop <= set[i].lo_stop)
        error->all(FLERR, "Final box dimension due to fix {} is < 0.0", style);
    } else if (set[i].style == TRATE) {
      double shift = 0.5 * ((set[i].hi_start - set[i].lo_start) * exp(set[i].rate * delt));
      set[i].lo_stop = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_stop = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
    } else if (set[i].style == WIGGLE) {
      double shift = 0.5 * set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
      set[i].lo_stop = set[i].lo_start - shift;
      set[i].hi_stop = set[i].hi_start + shift;
    }
  }

  for (int i = 3; i < 6; i++) {
    if (update->firststep == update->beginstep) {
      if (i == 5) set[i].tilt_start = domain->xy;
      else if (i == 4) set[i].tilt_start = domain->xz;
      else if (i == 3) set[i].tilt_start = domain->yz;
    } else set[i].tilt_start = set[i].tilt_initial;

    if (set[i].style == FINAL) {
      set[i].tilt_stop = set[i].ftilt;
    } else if (set[i].style == DELTA) {
      set[i].tilt_stop = set[i].tilt_start + set[i].dtilt;
    } else if (set[i].style == VEL) {
      set[i].tilt_stop = set[i].tilt_start + delt * set[i].vel;
    } else if (set[i].style == ERATE) {
      // Account for effect of TRATE elongation on shearing
      double arate = 0.0, brate = 0.0, h_bb;
      if (i == 3) {
        if (set[1].style == TRATE) arate = set[1].rate;
        if (set[2].style == TRATE) brate = set[2].rate;
        h_bb = set[2].hi_start - set[2].lo_start;
      }
      if (i == 4) {
        if (set[0].style == TRATE) arate = set[0].rate;
        if (set[2].style == TRATE) brate = set[2].rate;
        h_bb = set[2].hi_start - set[2].lo_start;
      }
      if (i == 5) {
        if (set[0].style == TRATE) arate = set[0].rate;
        if (set[1].style == TRATE) brate = set[1].rate;
        h_bb = set[1].hi_start - set[1].lo_start;
      }
      // TODO: set up/use a nearly_equal function for
      //       floating point comparisons?
      if (arate == 0.0) {
        if (brate == 0.0)
          set[i].tilt_stop = set[i].rate * h_bb * delt;
        else
          set[i].tilt_stop = set[i].rate * h_bb / brate * (exp(brate * delt)-1.0);
      } else {
        if (brate == 0.0)
          set[i].tilt_stop = set[i].rate * h_bb / arate * (exp(arate * delt)-1.0);
        else if (arate == brate)
          set[i].tilt_stop = set[i].rate * h_bb * delt * exp(arate * delt);
        else
          set[i].tilt_stop = set[i].rate * h_bb / (brate - arate) * (exp(brate * delt) - exp(arate * delt));
      }
      set[i].tilt_stop += set[i].tilt_start * exp(arate * delt);
    } else if (set[i].style == TRATE) {
      set[i].tilt_stop = set[i].tilt_start * exp(set[i].rate * delt);
    } else if (set[i].style == WIGGLE) {
      double shift = set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
      set[i].tilt_stop = set[i].tilt_start + shift;

      // compute min/max for WIGGLE = extrema tilt factor will ever reach

      if (set[i].amplitude >= 0.0) {
        if (delt < 0.25 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start;
          set[i].tilt_max = set[i].tilt_start + shift;
        } else if (delt < 0.5 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        } else if (delt < 0.75 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - shift;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        } else {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        }
      } else {
        if (delt < 0.25 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - shift;
          set[i].tilt_max = set[i].tilt_start;
        } else if (delt < 0.5 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start;
        } else if (delt < 0.75 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start + shift;
        } else {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        }
      }
    }
  }

  // if using tilt TRATE, then initial tilt must be non-zero

  for (int i = 3; i < 6; i++)
    if (set[i].style == TRATE && set[i].tilt_start == 0.0)
      error->all(FLERR, "Cannot use fix {} trate on a box with zero tilt", style);

  // set domain->h_rate values for use by domain and other fixes/computes
  // initialize all rates to 0.0
  // cannot set here for VOLUME,WIGGLE,VARIABLE since not constant
  // TRATE also not constant, but is known at init

  h_rate = domain->h_rate;
  h_ratelo = domain->h_ratelo;

  for (int i = 0; i < 3; i++) {
    h_rate[i] = h_ratelo[i] = 0.0;
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == SCALE || set[i].style == VEL ||
        set[i].style == ERATE) {
      double dlo_dt, dhi_dt;
      if (delt != 0.0) {
        dlo_dt = (set[i].lo_stop - set[i].lo_start) / delt;
        dhi_dt = (set[i].hi_stop - set[i].hi_start) / delt;
      } else dlo_dt = dhi_dt = 0.0;
      h_rate[i] = dhi_dt - dlo_dt;
      h_ratelo[i] = dlo_dt;
    } else if (set[i].style == TRATE) {
      h_rate[i] = set[i].rate * (set[i].hi_start - set[i].lo_start);
      h_ratelo[i] = -0.5 * h_rate[i];
    }
  }

  for (int i = 3; i < 6; i++) {
    h_rate[i] = 0.0;
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == VEL) {
      if (delt != 0.0)
        h_rate[i] = (set[i].tilt_stop - set[i].tilt_start) / delt;
      else h_rate[i] = 0.0;
    } else if (set[i].style == ERATE) {
        double arate = 0.0, h_bb;
        if (i == 3) {
          if (set[1].style == TRATE) arate = set[1].rate;
          h_bb = set[2].hi_start - set[2].lo_start;
        }
        if (i == 4) {
          if (set[0].style == TRATE) arate = set[0].rate;
          h_bb = set[2].hi_start - set[2].lo_start;
        }
        if (i == 5) {
          if (set[0].style == TRATE) arate = set[0].rate;
          h_bb = set[1].hi_start - set[1].lo_start;
        }
        h_rate[i] = set[i].rate*h_bb + set[i].tilt_start*arate;
    } else if (set[i].style == TRATE) {
      h_rate[i] = set[i].rate*set[i].tilt_start;
    }
  }

  // Account for xy shear with yz tilt

  if (set[5].style == ERATE && set[5].rate != 0.0 && set[4].style == ERATE) {
    h_rate[4] += set[5].rate*set[3].tilt_start;
    set[4].tilt_stop += calc_xz_correction(delt);
  }

  // if yz changes and will cause box flip, then xy cannot be changing
  // yz = [3], xy = [5]
  // this is b/c the flips would induce continuous changes in xz
  //   in order to keep the edge vectors of the flipped shape matrix
  //   an integer combination of the edge vectors of the unflipped shape matrix
  // VARIABLE for yz is error, since no way to calculate if box flip occurs
  // WIGGLE lo/hi flip test is on min/max oscillation limit, not tilt_stop
  // only trigger actual errors if flipflag is set
  // ERATE for xy and xz is accounted for, so allow in that case.

  if (set[3].style && set[5].style && set[3].style != ERATE && set[4].style != ERATE) {
    int flag = 0;
    double lo,hi;
    if (flipflag && set[3].style == VARIABLE)
      error->all(FLERR, "Fix {} cannot use yz variable with xy", style);
    if (set[3].style == WIGGLE) {
      lo = set[3].tilt_min;
      hi = set[3].tilt_max;
    } else lo = hi = set[3].tilt_stop;
    if (flipflag) {
      if (lo / (set[1].hi_start - set[1].lo_start) < -0.5 ||
          hi / (set[1].hi_start - set[1].lo_start) > 0.5) flag = 1;
      if (set[1].style) {
        if (lo / (set[1].hi_stop - set[1].lo_stop) < -0.5 ||
            hi / (set[1].hi_stop - set[1].lo_stop) > 0.5) flag = 1;
      }
      if (flag)
        error->all(FLERR, "Fix {} is changing yz too much with xy", style);
    }
  }

  // Variables/computes may not be valid during post_integrate,
  // so require end_flag for VARIABLE style
  if (!end_flag) {
    for (int i = 0; i < 6; ++i) {
      if (set[i].style == VARIABLE)
        error->all(FLERR,"Fix deform cannot use the variable style with nevery = 0");
    }
  }

  // detect if any rigid fixes exist so rigid bodies can be rescaled
  // rfix[] = vector with pointers to each fix rigid

  rfix.clear();

  for (auto &ifix : modify->get_fix_list())
    if (ifix->rigid_flag) rfix.push_back(ifix);

  if (!end_flag && utils::strmatch(update->integrate_style,"^respa")) {
    auto respa = dynamic_cast<Respa*>(update->integrate);
    nlevels_respa = respa->nlevels;
    step_respa = respa->step;
    nloop0_respa = respa->loop[nlevels_respa-1];
    kspace_level_respa = respa->level_kspace;
    for (int iloop = nlevels_respa-2; iloop >= 0; --iloop)
      nloop0_respa *= respa->loop[iloop];
    allow_flip_change = 0;
  } else allow_flip_change = 1;
}

/* ----------------------------------------------------------------------
  box flipped on previous step
  reset box tilts for flipped config and create new box in domain
  image_flip() adjusts image flags due to box shape change induced by flip
  remap() puts atoms outside the new box back into the new box
  perform irregular on atoms in lamda coords to migrate atoms to new procs
  important that image_flip comes before remap, since remap may change
    image flags to new values, making eqs in doc of Domain:image_flip incorrect
------------------------------------------------------------------------- */

void FixDeform::pre_exchange()
{
  if (flip == 0) return;

  domain->yz = set[3].tilt_target = set[3].tilt_flip;
  domain->xz = set[4].tilt_target = set[4].tilt_flip;
  domain->xy = set[5].tilt_target = set[5].tilt_flip;
  domain->set_global_box();
  domain->set_local_box();

  domain->image_flip(flipxy, flipxz, flipyz);
  domain->remap_all();

  domain->x2lamda(atom->nlocal);
  irregular->migrate_atoms();
  domain->lamda2x(atom->nlocal);

  flip = 0;
}

/* ---------------------------------------------------------------------- */

void FixDeform::end_of_step()
{
  nsteps = update->ntimestep - update->beginstep;
  nsteps_total = update->endstep - update->beginstep;
  dt = update->dt;
  update_box();
}

/* ---------------------------------------------------------------------- */

void FixDeform::post_integrate()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixDeform::post_integrate_respa(int ilevel, int iloop) {
  // Update box on RESPA level 0 since that is where x updates occur
  if (ilevel == 0) {
    // Set correct number of timesteps and dt based on RESPA level.
    // NOTE: This reduces the limit on the maximum number of timesteps by
    // a factor of nloop0_respa (the total number of inner timesteps).
    nsteps = (update->ntimestep - 1 - update->beginstep)*nloop0_respa + iloop + 1;
    nsteps_total = (update->endstep - update->beginstep)*nloop0_respa;
    dt = step_respa[ilevel];

    // Turn off k-space flag when calling update_box() to avoid
    // unnecessary kspace->setup() calls
    int old_kspace_flag = kspace_flag;
    kspace_flag = 0;
    update_box();
    kspace_flag = old_kspace_flag;

  } else if (need_flip_change && ilevel == nlevels_respa-1) {
    // Last step needed a box flip, so allow that now.
    // Need to flip on outer level so that reneighbouring occurs.
    // Don't change nsteps, nsteps_total, dt since positions haven't been
    // integrated yet.
    allow_flip_change = 1;
    update_box();
    allow_flip_change = 0;
    need_flip_change = 0;
  }

  // Box is changing every inner step, so kspace->setup() must be called
  // on kspace level
  if (kspace_flag && ilevel == kspace_level_respa) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   Update the box
 ---------------------------------------------------------------------- */
void FixDeform::update_box()
{
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  // set new box size for strain-based dims

  apply_strain();

  // set new box size for VOLUME dims that are linked to other dims
  // NOTE: still need to set h_rate for these dims

  apply_volume();

  if (varflag) modify->addstep_compute(update->ntimestep + nevery);

  update_domain();

  // redo KSpace coeffs since box has changed

  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   apply strain controls
------------------------------------------------------------------------- */

void FixDeform::apply_strain()
{
  // for NONE, target is current box size
  // for TRATE, set target directly based on current time, also set h_rate
  // for WIGGLE, set target directly based on current time, also set h_rate
  // for VARIABLE, set target directly via variable eval, also set h_rate
  // for others except VOLUME, target is linear value between start and stop

  double delta = nsteps;
  if (delta != 0.0) delta /= nsteps_total;

  for (int i = 0; i < 3; i++) {
    if (set[i].style == NONE) {
      set[i].lo_target = domain->boxlo[i];
      set[i].hi_target = domain->boxhi[i];
    } else if (set[i].style == TRATE) {
      double delt = nsteps * dt;
      double shift = 0.5 * ((set[i].hi_start - set[i].lo_start) * exp(set[i].rate * delt));
      set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
      h_rate[i] = set[i].rate * domain->h[i];
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (set[i].style == WIGGLE) {
      double delt = nsteps * dt;
      double shift = 0.5 * set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
      set[i].lo_target = set[i].lo_start - shift;
      set[i].hi_target = set[i].hi_start + shift;
      h_rate[i] = MY_2PI / set[i].tperiod * set[i].amplitude *
        cos(MY_2PI * delt / set[i].tperiod);
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (set[i].style == VARIABLE) {
      double del = input->variable->compute_equal(set[i].hvar);
      set[i].lo_target = set[i].lo_start - 0.5 * del;
      set[i].hi_target = set[i].hi_start + 0.5 * del;
      h_rate[i] = input->variable->compute_equal(set[i].hratevar);
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (set[i].style == FINAL || set[i].style == DELTA || set[i].style == SCALE ||
               set[i].style == VEL || set[i].style == ERATE) {
      set[i].lo_target = set[i].lo_start + delta * (set[i].lo_stop - set[i].lo_start);
      set[i].hi_target = set[i].hi_start + delta * (set[i].hi_stop - set[i].hi_start);
    }
  }

  // for triclinic, set new box shape
  // for NONE, target is current tilt
  // for TRATE, set target directly based on current time. set h_rate later after box flips
  // for ERATE, set target accounting for other deformation rates. set h_rate later after box flips
  // for WIGGLE, set target directly based on current time. also set h_rate
  // for VARIABLE, set target directly via variable eval. also set h_rate
  // for other styles, target is linear value between start and stop values

  if (triclinic) {
    for (int i = 3; i < 6; i++) {
      if (set[i].style == NONE) {
        if (i == 5) set[i].tilt_target = domain->xy;
        else if (i == 4) set[i].tilt_target = domain->xz;
        else if (i == 3) set[i].tilt_target = domain->yz;
      } else if (set[i].style == TRATE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        set[i].tilt_target = set[i].tilt_start * exp(set[i].rate * delt);
        h_rate[i] = set[i].rate * domain->h[i];
      } else if (set[i].style == ERATE) {
        // Solve ODE for a,b,c box vectors accounting for elongation caused by TRATE.
        // This is needed for SLLOD to be correct under mixed flow.
        double delt = nsteps * dt;
        double arate = 0.0, brate = 0.0, h_bb;
        if (i == 3) {
          if (set[1].style == TRATE) arate = set[1].rate;
          if (set[2].style == TRATE) brate = set[2].rate;
          h_bb = set[2].hi_start - set[2].lo_start;
        }
        if (i == 4) {
          if (set[0].style == TRATE) arate = set[0].rate;
          if (set[2].style == TRATE) brate = set[2].rate;
          h_bb = set[2].hi_start - set[2].lo_start;
        }
        if (i == 5) {
          if (set[0].style == TRATE) arate = set[0].rate;
          if (set[1].style == TRATE) brate = set[1].rate;
          h_bb = set[1].hi_start - set[1].lo_start;
        }
        // TODO: nearly_equal function for floating point comparisons?
        if (arate == 0.0) {
          if (brate == 0.0)
            set[i].tilt_target = set[i].rate*h_bb*delt;
          else
            set[i].tilt_target = set[i].rate*h_bb/brate * (exp(brate*delt)-1.0);
        } else {
          if (brate == 0.0)
            set[i].tilt_target = set[i].rate*h_bb/arate * (exp(arate*delt)-1.0);
          else if (arate == brate)
            set[i].tilt_target = set[i].rate*h_bb*delt * exp(arate*delt);
          else
            set[i].tilt_target = set[i].rate*h_bb/(brate-arate) * (exp(brate*delt)-exp(arate*delt));
        }
        set[i].tilt_target += set[i].tilt_start*exp(arate*delt);
      } else if (set[i].style == WIGGLE) {
        double delt = nsteps * dt;
        set[i].tilt_target = set[i].tilt_start +
          set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
        h_rate[i] = MY_2PI / set[i].tperiod * set[i].amplitude *
          cos(MY_2PI * delt / set[i].tperiod);
      } else if (set[i].style == VARIABLE) {
        double delta_tilt = input->variable->compute_equal(set[i].hvar);
        set[i].tilt_target = set[i].tilt_start + delta_tilt;
        h_rate[i] = input->variable->compute_equal(set[i].hratevar);
      } else {
        set[i].tilt_target = set[i].tilt_start + delta * (set[i].tilt_stop - set[i].tilt_start);
      }
    }

    // Correct for effects of deformation on xz tilt
    if (set[5].style == ERATE && set[5].rate != 0.0 && set[4].style == ERATE)
      set[4].tilt_target += calc_xz_correction(nsteps * dt);
  }
}

/* ----------------------------------------------------------------------
   apply volume controls
------------------------------------------------------------------------- */

void FixDeform::apply_volume()
{
  for (int i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;

    int dynamic1 = set[i].dynamic1;
    int dynamic2 = set[i].dynamic2;
    int fixed = set[i].fixed;
    double v0 = set[i].vol_start;
    double shift = 0.0;

    if (set[i].substyle == ONE_FROM_ONE) {
      shift = 0.5 * (v0 / (set[dynamic1].hi_target - set[dynamic1].lo_target) /
             (set[fixed].hi_start - set[fixed].lo_start));
    } else if (set[i].substyle == ONE_FROM_TWO) {
      shift = 0.5 * (v0 / (set[dynamic1].hi_target - set[dynamic1].lo_target) /
             (set[dynamic2].hi_target - set[dynamic2].lo_target));
    } else if (set[i].substyle == TWO_FROM_ONE) {
      shift = 0.5 * sqrt(v0 * (set[i].hi_start - set[i].lo_start) /
                 (set[dynamic1].hi_target - set[dynamic1].lo_target) /
                 (set[fixed].hi_start - set[fixed].lo_start));
    }

    h_rate[i] = (2.0 * shift / (domain->boxhi[i] - domain->boxlo[i]) - 1.0) / update->dt;
    h_ratelo[i] = -0.5 * h_rate[i];

    set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
    set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
  }
}

/* ----------------------------------------------------------------------
   Update box domain
------------------------------------------------------------------------- */

void FixDeform::update_domain()
{
  // tilt_target can be large positive or large negative value
  // add/subtract box lengths until tilt_target is closest to current value
  // need to know final xy tilt first since yz adjusts c vector by multiple of b vector
  // adjust xz last to account for adjustments made by yz


  if (triclinic) {
    double *h = domain->h;
    for (int i : {5, 3, 4}) {
      int idenom = 0;
      if (i == 3) idenom = 1; // yz
      else idenom = 0;        // xz || xy
      double denom = set[idenom].hi_target - set[idenom].lo_target;
      double denom_inv = 1.0 / denom;

      double current = h[i] / h[idenom];

      while (set[i].tilt_target * denom_inv - current > 0.0) {
        set[i].tilt_target -= denom;
        if (i == 3) set[4].tilt_target -= set[5].tilt_target;
      }
      while (set[i].tilt_target * denom_inv - current < 0.0) {
        set[i].tilt_target += denom;
        if (i == 3) set[4].tilt_target += set[5].tilt_target;
      }
      if (fabs(set[i].tilt_target * denom_inv - 1.0 - current) <
          fabs(set[i].tilt_target * denom_inv - current)) {
        set[i].tilt_target -= denom;
        if (i == 3) set[4].tilt_target -= set[5].tilt_target;
      }
    }
  }

  // if any tilt ratios exceed 0.5, set flip = 1 and compute new tilt values
  // do not flip in x or y if non-periodic (can tilt but not flip)
  //   this is b/c the box length would be changed (dramatically) by flip
  // if yz tilt exceeded, adjust C vector by one B vector
  // if xz tilt exceeded, adjust C vector by one A vector
  // if xy tilt exceeded, adjust B vector by one A vector
  // check yz first since it may change xz, then xz check comes after
  // if end_flag = 1, flip is performed on current timestep, before reneighboring in pre_exchange()
  // if end_flag = 0, flip is performed on next timestep

  if (triclinic && flipflag) {
    double xprd = set[0].hi_target - set[0].lo_target;
    double yprd = set[1].hi_target - set[1].lo_target;
    double xprdinv = 1.0 / xprd;
    double yprdinv = 1.0 / yprd;
    if (set[3].tilt_target * yprdinv < -0.5 ||
        set[3].tilt_target * yprdinv > 0.5 ||
        set[4].tilt_target * xprdinv < -0.5 ||
        set[4].tilt_target * xprdinv > 0.5 ||
        set[5].tilt_target * xprdinv < -0.5 ||
        set[5].tilt_target * xprdinv > 0.5) {
      if (!allow_flip_change) {
        // For rRESPA, can only flip in outer timestep, but could be integrating
        // in inner timestep. Just flag to calculate new box in outer rRESPA
        // level of next timestep.
        need_flip_change = 1;

      } else {
        set[3].tilt_flip = set[3].tilt_target;
        set[4].tilt_flip = set[4].tilt_target;
        set[5].tilt_flip = set[5].tilt_target;

        flipxy = flipxz = flipyz = 0;

        if (domain->yperiodic) {
          if (set[3].tilt_flip * yprdinv < -0.5) {
            set[3].tilt_flip += yprd;
            set[4].tilt_flip += set[5].tilt_flip;
            flipyz = 1;
          } else if (set[3].tilt_flip * yprdinv > 0.5) {
            set[3].tilt_flip -= yprd;
            set[4].tilt_flip -= set[5].tilt_flip;
            flipyz = -1;
          }
        }
        if (domain->xperiodic) {
          if (set[4].tilt_flip * xprdinv < -0.5) {
            set[4].tilt_flip += xprd;
            flipxz = 1;
          }
          if (set[4].tilt_flip * xprdinv > 0.5) {
            set[4].tilt_flip -= xprd;
            flipxz = -1;
          }
          if (set[5].tilt_flip * xprdinv < -0.5) {
            set[5].tilt_flip += xprd;
            flipxy = 1;
          }
          if (set[5].tilt_flip * xprdinv > 0.5) {
            set[5].tilt_flip -= xprd;
            flipxy = -1;
          }
        }

        flip = 0;
        if (flipxy || flipxz || flipyz) flip = 1;
        if (flip) next_reneighbor = update->ntimestep + (end_flag ? 1 : 0);
      }
    }
  }

  // For styles with h_rate dependence on xy/xz/yz, need to set h_rate after
  // box flips so that the correct streaming velocity can be recovered by
  // fix nvt/sllod
  if (triclinic) {
    double *h = domain->h;

    for (int i = 3; i < 6; i++) {
      if (set[i].style == TRATE) {
        h_rate[i] = set[i].rate * domain->h[i];
      } else if (set[i].style == ERATE) {
        // Solve ODE for a,b,c vectors accounting for elongation caused by TRATE.
        // This is needed for SLLOD to be correct under mixed flow.
        // TODO: do other elongation styles need to be accounted for where possible?
        double delt = nsteps * dt;
        double arate = 0.0;
        if (i == 3) {
          if (set[1].style == TRATE) arate = set[1].rate;
          h_rate[i] = set[i].rate * (set[2].hi_target - set[2].lo_target) + arate * set[i].tilt_target;
        }
        if (i == 4) {
          if (set[0].style == TRATE) arate = set[0].rate;
          h_rate[i] = set[i].rate * (set[2].hi_target - set[2].lo_target) + arate * set[i].tilt_target;
        }
        if (i == 5) {
          if (set[0].style == TRATE) arate = set[0].rate;
          h_rate[i] = set[i].rate * (set[1].hi_target - set[1].lo_target) + arate * set[i].tilt_target;
        }
      }
    }
    // TODO: use nearly_equal for check on set[5].rate?
    if (set[5].style == ERATE && set[5].rate != 0.0 && set[4].style == ERATE)
      h_rate[4] += set[5].rate*set[3].tilt_target;
  }

  // convert atoms and rigid bodies to lamda coords

  if (remapflag == Domain::X_REMAP) {
    domain->x2lamda(atom->nlocal, groupbit);

    for (auto &ifix : rfix)
      ifix->deform(0);
  }

  // reset global and local box to new size/shape
  // only if deform fix is controlling the dimension

  if (dimflag[0]) {
    domain->boxlo[0] = set[0].lo_target;
    domain->boxhi[0] = set[0].hi_target;
  }
  if (dimflag[1]) {
    domain->boxlo[1] = set[1].lo_target;
    domain->boxhi[1] = set[1].hi_target;
  }
  if (dimflag[2]) {
    domain->boxlo[2] = set[2].lo_target;
    domain->boxhi[2] = set[2].hi_target;
  }
  if (triclinic) {
    if (dimflag[3]) domain->yz = set[3].tilt_target;
    if (dimflag[4]) domain->xz = set[4].tilt_target;
    if (dimflag[5]) domain->xy = set[5].tilt_target;
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert atoms and rigid bodies back to box coords

  if (remapflag == Domain::X_REMAP) {
    domain->lamda2x(atom->nlocal, groupbit);

    for (auto &ifix : rfix)
      ifix->deform(1);
  }
}

/* ----------------------------------------------------------------------
   Calculate correction to xz tilt due to xy shear with yz tilt.
   NOTE: only considers xx, yy, zz deformation with TRATE
         and xy, xz, yz with ERATE. (i.e. constant flow tensor)
   Non-zero xy rate and xz.style == ERATE is assumed.
   Requires xz style of ERATE. If xz style is NONE then changes aren't
   tracked properly (e.g. if there is pressure control on the xz tilt).
------------------------------------------------------------------------- */
double FixDeform::calc_xz_correction(double delt) {
  // Solve ODE for xy component of xz tilt factor
  double g_xy = set[5].rate;
  double g_yz = set[3].rate;
  double h_yz0 = set[3].tilt_start;
  if (set[3].style == ERATE && g_yz != 0.0) {
    double e_xx = 0.0, e_yy = 0.0, e_zz = 0.0;
    if (set[0].style == TRATE) e_xx = set[0].rate;
    if (set[1].style == TRATE) e_yy = set[1].rate;
    if (set[2].style == TRATE) e_zz = set[2].rate;
    double h_zz0 = set[2].hi_start - set[2].lo_start;
    if (e_xx == e_zz) {
      if (e_xx == 0.0) {
        if (e_yy == 0.0) {
          // e_xx = e_yy = e_zz = 0
          // (Pure shear)
          return g_xy * (h_yz0 * delt + 0.5 * g_yz * h_zz0 * delt * delt);
        } else {
          // e_xx = e_zz = 0, e_yy != 0
          // (Shear + y extension - non vol. preserving)
          double yyfac = (exp(e_yy * delt) - 1.0) / e_yy;
          return g_xy * g_yz * h_zz0 / e_yy * (yyfac - delt) + g_xy * h_yz0 * yyfac;
        }
      } else {
        if(e_yy == 0.0) {
          // e_xx = e_zz != 0, e_yy = 0
          // (Shear + xz extension - non vol. preserving)
          double xfac = exp(e_xx*delt);
          return g_xy * g_yz * h_zz0 / e_zz * (delt * xfac - (xfac - 1.0) / e_xx)
                 + g_xy * h_yz0 * ((xfac - 1.0) / e_xx);
        } else if (e_yy == e_zz) {
          // e_xx = e_yy = e_zz != 0
          // (Shear + xyz extension - non vol. preserving)
          return g_xy * (h_yz0 * delt + 0.5 * g_xy * g_yz * h_zz0 * delt * delt) * exp(e_xx * delt);
        } else {
          // e_xx = e_zz != 0, e_yy != e_xx, e_yy != 0
          // (Shear + xyz extension - possibly vol. preserving)
          double xfac = exp(e_xx * delt);
          double yfac = exp(e_yy * delt);
          double xyfac = (yfac - xfac) / (e_yy - e_xx);
          return g_xy * g_yz * h_zz0 / (e_zz - e_yy) * (delt * xfac - xyfac)
                 + g_xy * h_yz0 * xyfac;
        }
      }
    } else if (e_xx == 0.0) {
      if (e_yy == 0.0) {
        // e_xx = e_yy = 0, e_zz != 0
        // (Shear + z extension - non vol. preserving)
        return g_xy * g_yz * h_zz0 / e_zz * ((exp(e_zz * delt) - 1.0) / e_zz - delt)
               + g_xy * h_yz0 * delt;
      } else if (e_yy == e_zz) {
        // e_xx = 0, e_yy = e_zz != 0
        // (Shear + yz extension - non vol. preserving)
        double yfac = exp(e_yy*delt);
        return g_xy * g_yz * h_zz0 / e_yy * (delt * yfac - (yfac - 1.0) / e_yy)
               + g_xy * h_yz0 * ((yfac - 1.0) / e_yy);
      } else {
        // e_xx = 0, e_yy != 0, e_zz != 0, e_yy != e_zz
        // (Shear + yz extension - possibly vol. preserving)
        double yfac = (exp(e_yy * delt) - 1.0) / e_yy;
        return g_xy * g_yz * h_zz0 / (e_zz - e_yy) * ((exp(e_zz) - 1.0) / e_zz - yfac)
               + g_xy*h_yz0*yfac;
      }
    } else if (e_zz == 0.0) {
      if (e_yy == 0.0) {
        // e_xx != 0, e_yy = e_zz = 0
        // (Shear + x extension - non vol. preserving)
        double xfac = (exp(e_xx * delt) - 1.0) / e_xx;
        return g_xy * g_yz * h_zz0 / e_xx * (xfac - delt) + g_xy * h_yz0 * xfac;
      } else if (e_xx == e_yy) {
        // e_xx = e_yy != 0, e_zz = 0
        // (Shear + xy extension - non vol. preserving)
        double xfac = exp(e_xx * delt);
        return g_xy * g_yz / e_xx * ((1.0 - xfac) / e_xx + delt * xfac)
               + g_xy * h_yz0 * (delt * xfac);
      } else {
        // e_xx != 0, e_yy != 0, e_xx != e_yy, e_zz = 0
        // (Shear + xy extension - possibly vol. preserving)
        double xfac = exp(e_xx * delt);
        double yfac = exp(e_yy * delt);
        double xyfac = (yfac - xfac) / (e_yy - e_xx);
        return g_xy * g_yz * h_zz0 / e_yy * (xyfac + (1.0 - xfac) / e_xx)
               + g_xy * h_yz0 * xyfac;
      }
    } else {
      if (e_yy == 0.0) {
        // e_xx != 0, e_zz != 0, e_xx != e_zz, e_yy = 0
        // (Shear + xz extension - possibly vol. preserving)
        double xfac = exp(e_xx * delt);
        double zfac = exp(e_zz * delt);
        return g_xy * g_yz * h_zz0 / e_zz * ((zfac - xfac) / (e_zz - e_xx)
               + (1.0 - xfac) / e_xx) + g_xy * h_yz0 * ((1.0 - xfac) / e_xx);
      } else if (e_yy == e_zz) {
        // e_xx != 0, e_yy != 0, e_zz = e_yy, e_xx != e_zz
        // (Shear + xyz extension - possibly vol. preserving)
        double xfac = exp(e_xx*delt);
        double yfac = exp(e_yy*delt);
        double xyfac = (yfac - xfac) / (e_yy - e_xx);
        return g_xy * g_yz * h_zz0 / (e_yy - e_xx) * (delt * yfac - xyfac)
               + g_xy * h_yz0 * xyfac;
      } else if (e_yy == e_xx) {
        // e_xx != 0, e_yy = e_xx, e_zz != 0, e_xx != e_zz
        // (Shear + xyz extension - possibly vol. preserving)
        double xfac = exp(e_xx * delt);
        double zfac = exp(e_zz * delt);
        return g_xy * g_yz * h_zz0 / (e_zz - e_yy) * ((zfac - xfac) / (e_zz - e_xx)
               - delt * xfac) + g_xy * h_yz0 * (delt * xfac);
      } else {
        // e_xx != 0, e_yy != 0, e_zz != 0, e_xx != e_zz, e_xx != e_yy, e_yy != e_zz
        // (Shear + xyz extension - possibly vol. preserving)
        double xfac = exp(e_xx * delt);
        double yfac = exp(e_yy * delt);
        double zfac = exp(e_zz * delt);
        double xzfac = (zfac - xfac) / (e_zz - e_xx);
        double xyfac = (yfac - xfac) / (e_yy - e_xx);
        return g_xy * g_yz * h_zz0 / (e_zz - e_yy) * (xzfac - xyfac) + g_xy * h_yz0 * xyfac;
      }
    }
  } else {
    // h_yz is constant
    return g_xy * h_yz0 * delt;
  }
}

/* ----------------------------------------------------------------------
   write Set data to restart file
------------------------------------------------------------------------- */

void FixDeform::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = 6 * sizeof(Set);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(set, sizeof(Set), 6, fp);
  }
}

/* ----------------------------------------------------------------------
   use selected state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixDeform::restart(char *buf)
{
  int samestyle = 1;
  Set *set_restart = (Set *) buf;
  for (int i = 0; i < 6; ++i) {
    // restore data from initial state
    set[i].lo_initial = set_restart[i].lo_initial;
    set[i].hi_initial = set_restart[i].hi_initial;
    set[i].vol_initial = set_restart[i].vol_initial;
    set[i].tilt_initial = set_restart[i].tilt_initial;
    // check if style settings are consistent (should do the whole set?)
    if (set[i].style != set_restart[i].style)
      samestyle = 0;
    if (set[i].substyle != set_restart[i].substyle)
      samestyle = 0;
  }
  if (!samestyle)
    error->all(FLERR, "Fix {} settings not consistent with restart", style);
}

/* ---------------------------------------------------------------------- */

void FixDeform::options(int narg, char **arg)
{
  const std::string thiscmd = fmt::format("fix {}", style);
  if (narg < 0) utils::missing_cmd_args(FLERR, thiscmd, error);

  remapflag = Domain::X_REMAP;
  scaleflag = 1;
  flipflag = 1;

  // arguments for child classes

  std::unordered_map<std::string, int> child_options;
  if (utils::strmatch(style, "^deform/pressure")) {
    child_options.insert({{"couple", 2}, {"max/rate", 2}, {"normalize/pressure", 2},
                          {"vol/balance/p", 2}});
  }

  // parse all optional arguments for this parent and also child classes
  // for child classes, simply store them in leftover_iarg and skip over them

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "remap") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, thiscmd + " remap", error);
      if (strcmp(arg[iarg + 1], "x") == 0) remapflag = Domain::X_REMAP;
      else if (strcmp(arg[iarg + 1], "v") == 0) remapflag = Domain::V_REMAP;
      else if (strcmp(arg[iarg + 1], "none") == 0) remapflag = Domain::NO_REMAP;
      else error->all(FLERR, "Illegal fix {} remap command: {}", style, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, thiscmd + " units", error);
      if (strcmp(arg[iarg + 1], "box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0) scaleflag = 1;
      else error->all(FLERR, "Illegal fix {} units command: {}", style, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "flip") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, thiscmd + " flip", error);
      flipflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (child_options.find(arg[iarg]) != child_options.end()) {
      auto nskip = child_options[arg[iarg]];
      if (iarg + nskip > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} {}", style, arg[iarg]), error);
      for (int i = 0; i < nskip; i++) leftover_iarg.push_back(iarg + i);
      iarg += nskip;
    } else error->all(FLERR, "Unknown fix {} keyword: {}", style, arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of Irregular
------------------------------------------------------------------------- */

double FixDeform::memory_usage()
{
  double bytes = 0.0;
  if (irregular) bytes += irregular->memory_usage();
  return bytes;
}
