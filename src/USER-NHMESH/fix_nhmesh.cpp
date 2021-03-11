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

#include "fix_nhmesh.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "comm.h"
#include "neighbor.h"
#include "modify.h"
#include "compute.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "compute_nhmesh_coupling_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};

/* ----------------------------------------------------------------------
   NHMesh integrator for coupled and/or linear combinations of
   Nose-Hoover thermostats
 ---------------------------------------------------------------------- */

FixNHMesh::FixNHMesh(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_temp(nullptr), eta(nullptr), eta_dot(nullptr),
  eta_dotdot(nullptr), eta_mass(nullptr), t_start(nullptr), t_stop(nullptr),
  t_freq(nullptr), t_current(nullptr), t_target(nullptr), ke_current(nullptr),
  ke_target(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix temp/nhmesh command");

  restart_global = 1;
  dynamic_group_allow = 1;
  time_integrate = 1;
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 0;
  ecouple_flag = 1;

  // default values

  //drag = 0.0;
  nc_tchain = 1;
  eta_mass_flag = 1;
  mesh_coupling_flag = 0;

  tcomputeflag = 0;
  id_temp = nullptr;


  // used by FixNVTSllod to preserve non-default value

  // mtchain_default_flag = 1;

  id_coupling = utils::strdup(arg[3]);
  int icoupling = modify->find_compute(id_coupling);
  if (icoupling < 0)
    error->all(FLERR,"Coupling ID of fix temp/nhmesh doesn't exist");
  coupling = dynamic_cast<ComputeNHMeshCouplingAtom*>(modify->compute[icoupling]);
  if (coupling == nullptr)
    error->all(FLERR,"Invalid coupling compute for temp/nhmesh");

  n_thermostats = coupling->n_thermostats;
  if (n_thermostats <= 0)
    error->all(FLERR,"Number of thermostats for temp/nhmesh must be > 0");

  int iarg = 4, i;


  double t_period[n_thermostats];
  t_start = new double[n_thermostats];
  t_stop = new double[n_thermostats];
  t_freq = new double[n_thermostats];
  t_current = new double[n_thermostats];
  t_target = new double[n_thermostats];
  is_real = new int[n_thermostats];
  ke_current = new double[n_thermostats];
  ke_target = new double[n_thermostats];
  mesh_dof = new double[n_thermostats];
  tdof = new double[n_thermostats];
  factor_eta = new double[n_thermostats];
  for (i = 0; i < n_thermostats; i++) {
    mesh_dof[i] = 0.0;
    is_real[i] = 1;
  }

  // TODO: accept vectors here
  for (i = 0; i < n_thermostats; i++) {
    if (iarg+1 > narg)
      error->all(FLERR,"Illegal fix temp/nhmesh command");
    if (i == 0 && strcmp(arg[iarg], "all") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix temp/nhmesh command");
      for (; i < n_thermostats; i++) {
        t_start[i] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        t_target[i] = t_start[i];
      }
      iarg+=2;
    } else {
      t_start[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
      t_target[i] = t_start[i];
    }
  }
  for (i = 0; i < n_thermostats; i++) {
    if (iarg+1 > narg)
      error->all(FLERR,"Illegal fix temp/nhmesh command");
    if (i == 0 && strcmp(arg[iarg], "all") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix temp/nhmesh command");
      for (; i < n_thermostats; i++) {
        t_stop[i] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (t_start[i] <= 0.0 || t_stop[i] <= 0.0)
          error->all(FLERR,
                     "Target temperature for fix temp/nhmesh cannot be 0.0");
      }
      iarg+=2;
    } else {
      t_stop[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
      if (t_start[i] <= 0.0 || t_stop[i] <= 0.0)
        error->all(FLERR,
                   "Target temperature for fix temp/nhmesh cannot be 0.0");
    }
  }
  for (i = 0; i < n_thermostats; i++) {
    if (iarg+1 > narg)
      error->all(FLERR,"Illegal fix temp/nhmesh command");
    if (i == 0 && strcmp(arg[iarg], "all") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix temp/nhmesh command");
      for (; i < n_thermostats; i++)
        t_period[i] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else t_period[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
  }

  // process keywords

  while (iarg < narg) {

    if (strcmp(arg[iarg],"tloop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix temp/nhmesh command");
      nc_tchain = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nc_tchain < 0) error->all(FLERR,"Illegal fix temp/nhmesh command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"couple") == 0) {
      if (iarg + n_thermostats*n_thermostats > narg)
        error->all(FLERR,"Illegal fix temp/nhmesh command");
      mesh_coupling_flag = 1;
      // TODO: Accept vector/array input
      memory->create(mesh_coupling, n_thermostats,
          n_thermostats, "temp/nhmesh:mesh_coupling");
      for (i = 0; i < n_thermostats; i++)
        for (int j = 0; j < n_thermostats; j++) {
          mesh_coupling[i][j] = utils::numeric(FLERR,arg[iarg++],false,lmp);
          if (mesh_coupling[i][j] < 0 || mesh_coupling[i][j] > 1)
            error->all(FLERR,"Thermostat couplings must be between 0 and 1");
          mesh_dof[i] += mesh_coupling[i][j];
        }
    } else if (strcmp(arg[iarg],"virtual") == 0) {
      iarg++;
      while (iarg < narg && utils::is_integer(arg[iarg])) {
        int ival = utils::inumeric(FLERR,arg[iarg++],false,lmp);
        if (ival < 1 || ival > n_thermostats)
          error->all(FLERR,"Virtual thermostats must be > 0 and <= N");
        is_real[ival-1] = 0;
      }
    } else error->all(FLERR,"Illegal fix temp/nhmesh command");
  }

  // convert input periods to frequencies

  for (i = 0; i < n_thermostats; i++) t_freq[i] = 1.0 / t_period[i];

  // Nose/Hoover temp init

  eta = new double[n_thermostats];
  eta_dot = new double[n_thermostats];
  eta_dotdot = new double[n_thermostats];
  for (i = 0; i < n_thermostats; i++) {
    eta[i] = eta_dot[i] = eta_dotdot[i] = 0.0;
  }
  eta_mass = new double[n_thermostats];

  size_vector = 4*n_thermostats;
}

/* ---------------------------------------------------------------------- */

FixNHMesh::~FixNHMesh()
{
  if (copymode) return;

  // delete temperature compute if fix created it

  if (tcomputeflag) modify->delete_compute(id_temp);
  delete [] id_temp;

  if (mesh_coupling_flag) memory->destroy(mesh_coupling);
  delete [] mesh_dof;

  delete [] t_start;
  delete [] t_stop;
  delete [] t_current;
  delete [] t_target;
  delete [] ke_current;
  delete [] ke_target;
  delete [] t_freq;

  delete [] eta;
  delete [] eta_dot;
  delete [] eta_dotdot;
  delete [] eta_mass;

}

/* ---------------------------------------------------------------------- */

int FixNHMesh::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNHMesh::init()
{
  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/nhmesh does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  // set timesteps and frequencies

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  dto = dthalf;

  // tdrag_factor = 1.0 - (update->dt * t_freq * drag / nc_tchain);

  boltz = force->boltz;
  nktv2p = force->nktv2p;

  if (strstr(update->integrate_style,"respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
    dto = 0.5*step_respa[0];
  }

}

/* ----------------------------------------------------------------------
   compute T before integrator starts
------------------------------------------------------------------------- */

void FixNHMesh::setup(int /*vflag*/)
{
  int i;

  compute_temp_current();

  // masses and initial forces on thermostat variables

  for (i = 0; i < n_thermostats; i++)
    eta_mass[i] = tdof[i] * boltz * t_target[i] / (t_freq[i]*t_freq[i]);

  for (i = 0; i < n_thermostats; i++) {
    eta_dotdot[i] = ke_current[i] - tdof[i] * boltz * t_target[i];
    eta_dotdot[i] /= eta_mass[i];
  }
}

/* ----------------------------------------------------------------------
   Compute temperature of each thermostat from KE of particles and other
   thermostats
   Also updates tdof
------------------------------------------------------------------------- */

void FixNHMesh::compute_temp_current() {

  // NOTE: this also updates the coupling matrix for dof calculation. If another
  // temperature compute can be used then code here should make sure
  // coupling->compute_peratom() is called
  temperature->compute_array();
  double **&ke = temperature->array;

  // This is also called by temperature compute - could make friend class to
  // avoid recalculation, but that would probably limit modularity
  double natoms_temp = group->count(igroup);

  for (int i = 0; i < n_thermostats; i++) {
    tdof[i] = is_real[i]*coupling->therm_sum[i]*temperature->dof/natoms_temp
              + mesh_dof[i];
    // update ke_target since dof might have changed
    ke_target[i] = tdof[i] * boltz * t_target[i];
    ke_current[i] = is_real[i] * (ke[i][0] + ke[i][1] + ke[i][2]);
    ke_current[i] *= force->mvv2e;
    if (mesh_coupling_flag)
      for (int j = 0; j < n_thermostats; j++)
        ke_current[i] += mesh_coupling[i][j] * eta_mass[j]*eta_dot[j]*eta_dot[j];
    t_current[i] = tdof[i] == 0 ? t_target[i] : ke_current[i] / (boltz * tdof[i]);
  }
}

/* ----------------------------------------------------------------------
   1st half of Verlet update
------------------------------------------------------------------------- */

void FixNHMesh::initial_integrate(int /*vflag*/)
{
  // update eta_dot

  compute_temp_current();
  compute_temp_target();
  nhmesh_temp_integrate();

  // update particles

  nve_v();

  nve_x();

}

/* ----------------------------------------------------------------------
   2nd half of Verlet update
------------------------------------------------------------------------- */

void FixNHMesh::final_integrate()
{
  nve_v();

  // Update temperature
  // This also recomputes the coupling matrix based on new x values, and
  // accounting for reneighbouring, which is needed for temp_integrate()
  compute_temp_current();

  // update eta_dot

  nhmesh_temp_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNHMesh::initial_integrate_respa(int /*vflag*/, int ilevel, int /*iloop*/)
{
  // set timesteps by level

  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  // outermost level - update eta_dot and omega_dot, apply to v
  // all other levels - NVE update of v
  // x,v updates only performed for atoms in group

  if (ilevel == nlevels_respa-1) {

    // update eta_dot

    compute_temp_current();
    compute_temp_target();
    nhmesh_temp_integrate();

    // recompute pressure to account for change in KE
    // t_current is up-to-date, but compute_temperature is not
    // compute appropriately coupled elements of mvv_current

    // NOTE: as above, deleted code here, leaving comment in case it matters

    nve_v();

  } else nve_v();

  // innermost level - also update x only for atoms in group
  // if barostat, perform 1/2 step remap before and after

  if (ilevel == 0) nve_x();

}

/* ---------------------------------------------------------------------- */

void FixNHMesh::final_integrate_respa(int ilevel, int /*iloop*/)
{
  // set timesteps by level

  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  // outermost level - update eta_dot and omega_dot, apply via final_integrate
  // all other levels - NVE update of v

  if (ilevel == nlevels_respa-1) final_integrate();
  else nve_v();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixNHMesh::write_restart(FILE *fp)
{
  int nsize = size_restart_global();

  double *list;
  memory->create(list,nsize,"nh:list");

  pack_restart_data(list);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
  }

  memory->destroy(list);
}

/* ----------------------------------------------------------------------
    calculate the number of data to be packed
------------------------------------------------------------------------- */

int FixNHMesh::size_restart_global()
{
  int nsize = 1;              // n_thermostats
  nsize += 2*n_thermostats;   // eta, eta_dot

  return nsize;
}

/* ----------------------------------------------------------------------
   pack restart data
------------------------------------------------------------------------- */

int FixNHMesh::pack_restart_data(double *list)
{
  int n = 0;

  list[n++] = n_thermostats;
  for (int i = 0; i < n_thermostats; i++)
    list[n++] = eta[i];
  for (int i = 0; i < n_thermostats; i++)
    list[n++] = eta_dot[i];

  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixNHMesh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  int m = static_cast<int> (list[n++]);
  if (m != n_thermostats)
    error->all(FLERR,
        "fix temp/nhmesh restarted with incompatible coupling compute");
  for (int i = 0; i < n_thermostats; i++)
    eta[i] = list[n++];
  for (int i = 0; i < n_thermostats; i++)
    eta_dot[i] = list[n++];
}

/* ---------------------------------------------------------------------- */

int FixNHMesh::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tcomputeflag) {
      modify->delete_compute(id_temp);
      tcomputeflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->array_flag == 0 ||
        temperature->size_array_rows != n_thermostats ||
        temperature->size_array_cols < 3)
      error->all(FLERR,
                 "Fix_modify temperature ID is not compatible with nhmesh");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for fix modify is not for group all");

    return 2;

  }

  return 0;
}

/* ---------------------------------------------------------------------- */

double FixNHMesh::compute_scalar()
{
  int i;
  double energy = 0.0;

  // Thermostat energy is kinetic + potential:
  // KE = Sum(0.5 * eta_mass * eta_dot^2)
  // PE = Sum(dof * k * t_target * eta)
  for (i = 0; i < n_thermostats; i++)
    energy += ke_target[i] * eta[i] + 0.5*eta_mass[i]*eta_dot[i]*eta_dot[i];

  return energy;
}

/* ----------------------------------------------------------------------
   return a single element of the following vectors, in this order:
      eta[n_thermostats], eta_dot[n_thermostats], PE_eta[n_thermostats],
      KE_eta_dot[n_thermostats]
------------------------------------------------------------------------- */

double FixNHMesh::compute_vector(int n)
{
  int ilen;

  ilen = n_thermostats;
  if (n < ilen) return eta[n];
  n -= ilen;

  ilen = n_thermostats;
  if (n < ilen) return eta_dot[n];
  n -= ilen;

  ilen = n_thermostats;
  if (n < ilen) return ke_target[n] * eta[n];
  n -= ilen;

  ilen = n_thermostats;
  if (n < ilen) return 0.5*eta_mass[n]*eta_dot[n]*eta_dot[n];
  n -= ilen;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

void FixNHMesh::reset_target(double t_new)
{
  // Called by fix temper, so needs to take scalar. Maybe look at treating it as
  // a scaling factor instead of the actual new temperature? In that case, would
  // need a separate array to store config temperatures, and multiply that by
  // scale_fac to give t_target, etc.
  for (int i = 0; i < n_thermostats; i++)
    t_target[i] = t_start[i] = t_stop[i] = t_new;
}

/* ---------------------------------------------------------------------- */

void FixNHMesh::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  dt16 = 0.0625 * update->dt;
  dto = dthalf;

  // If using respa, then remap is performed in innermost level

  if (strstr(update->integrate_style,"respa"))
    dto = 0.5*step_respa[0];
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixNHMesh::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"n_thermostats") == 0) {
    return &n_thermostats;
  }
  dim=1;
  if (strcmp(str,"t_target") == 0) {
    return &t_target;
  } else if (strcmp(str,"t_start") == 0) {
    return &t_start;
  } else if (strcmp(str,"t_stop") == 0) {
    return &t_stop;
  } else if (strcmp(str,"eta") == 0) {
    return &eta;
  } else if (strcmp(str,"ke_current")) {
    return &ke_current;
  }
  dim=2;
  if (mesh_coupling_flag && strcmp(str,"mesh_coupling")) {
    return &mesh_coupling;
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   perform half-step update of chain thermostat variables
------------------------------------------------------------------------- */

void FixNHMesh::nhmesh_temp_integrate()
{
  int i;
  double expfac[n_thermostats];
  // ke_current should already be up to date from compute_temp_current();
  // for (i = 0; i < n_thermostats; i++)
  //   ke_current[i] = tdof[i] * boltz * t_current[i];

  // Update masses, to preserve initial freq, if flag set

  if (eta_mass_flag)
    for (i = 0; i < n_thermostats; i++)
      eta_mass[i] = tdof[i] * boltz * t_target[i] / (t_freq[i]*t_freq[i]);

  for (i = 0; i < n_thermostats; i++) {
    if (eta_mass[i] > 0.0)
      eta_dotdot[i] = (ke_current[i] - ke_target[i])/eta_mass[i];
    else eta_dotdot[i] = 0.0;
  }

  double ncfac = 1.0/nc_tchain;
  double eta_dot_step[n_thermostats];
  for (int iloop = 0; iloop < nc_tchain; iloop++) {

    // Update thermostat velocities
    for (i = 0; i < n_thermostats; i++) {
      if (mesh_coupling_flag) {
        expfac[i] = 0.0;
        for (int j = 0; j < n_thermostats; j++)
          expfac[i] += mesh_coupling[j][i]*eta_dot[j];
        expfac[i] = exp(-ncfac*dt16*expfac[i]);
      } else expfac[i] = 1.0;
      eta_dot_step[i] = eta_dot[i] * expfac[i];
      eta_dot_step[i] += eta_dotdot[i] * ncfac*dt8;
      eta_dot_step[i] *= expfac[i];
    }
    for (i = 0; i < n_thermostats; i++) {
      if (mesh_coupling_flag) {
        expfac[i] = 0.0;
        for (int j = 0; j < n_thermostats; j++)
          expfac[i] += mesh_coupling[j][i]*eta_dot_step[j];
        expfac[i] = exp(-ncfac*dt16*expfac[i]);
      } // else expfac[i] = 1.0; // should already be set from prev. loop
      eta_dot[i] = eta_dot_step[i] * expfac[i];
      eta_dot[i] += eta_dotdot[i] * ncfac*dt8;
      eta_dot[i] *= expfac[i];
    }

    for (i = 0; i < n_thermostats; i++)
      factor_eta[i] = -ncfac*dthalf*eta_dot[i];
    nhmesh_v_temp();

    // Recalculate current temperatures since nhmesh_v_temp() updates ke_current
    // only.
    // TODO: Look into just rescaling ke_current instead of recalculating.
    //       Might not be possible.
    for (i = 0; i < n_thermostats; i++)
      t_current[i] = tdof[i] == 0 ? t_target[i]
                                  : ke_current[i] / (boltz * tdof[i]);

    for (i = 0; i < n_thermostats; i++) {
      if (eta_mass[i] > 0.0)
        eta_dotdot[i] = (ke_current[i] - ke_target[i])/eta_mass[i];
      else eta_dotdot[i] = 0.0;
    }

    // Update thermostat 'positions'
    for (i = 0; i < n_thermostats; i++)
      eta[i] += ncfac*dthalf*eta_dot[i];

    // Update thermostat velocities
    for (i = 0; i < n_thermostats; i++) {
      if (mesh_coupling_flag) {
        expfac[i] = 0.0;
        for (int j = 0; j < n_thermostats; j++) {
          expfac[i] += mesh_coupling[j][i]*eta_dot[j];
        }
        expfac[i] = exp(-ncfac*dt16*expfac[i]);
      } // else expfac[i] = 1.0; // should already be set
      eta_dot_step[i] = eta_dot[i] * expfac[i];
    }
    for (i = 0; i < n_thermostats; i++) {
      if (expfac[i] != 1.0)
        for (int j = 0; j < n_thermostats; j++)
          eta_dotdot[i] += mesh_coupling[i][j]*eta_mass[j]/eta_mass[i]*
            (eta_dot_step[j]*eta_dot_step[j] - eta_dot[j]*eta_dot[j]);
      eta_dot_step[i] += eta_dotdot[i] * ncfac*dt8;
      eta_dot_step[i] *= expfac[i];
    }
    for (i = 0; i < n_thermostats; i++) {
      if (mesh_coupling_flag) {
        expfac[i] = 0.0;
        for (int j = 0; j < n_thermostats; j++)
          expfac[i] += mesh_coupling[j][i]*eta_dot_step[j];
        expfac[i] = exp(-ncfac*dt16*expfac[i]);
      } // else expfac[i] = 1.0;
      eta_dot[i] = eta_dot_step[i] * expfac[i];
    }
    for (i = 0; i < n_thermostats; i++) {
      if (expfac[i] != 1.0)
        for (int j = 0; j < n_thermostats; j++)
          eta_dotdot[i] += mesh_coupling[i][j]*eta_mass[j]/eta_mass[i]*
            (eta_dot_step[j]*eta_dot_step[j] - eta_dot[j]*eta_dot[j]);
      eta_dot[i] += eta_dotdot[i] * ncfac*dt8;
      eta_dot[i] *= expfac[i];
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities
-----------------------------------------------------------------------*/

void FixNHMesh::nve_v()
{
  double dtfm;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm*f[i][0];
        v[i][1] += dtfm*f[i][1];
        v[i][2] += dtfm*f[i][2];
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
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

void FixNHMesh::nve_x()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x update by full step only for atoms in group

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of velocities
   also updates ke_current
-----------------------------------------------------------------------*/

void FixNHMesh::nhmesh_v_temp()
{
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

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        fac = 0;
        for (int j = 0; j < n_thermostats; j++)
          fac += is_real[j] * coupling->array_atom[i][j] * factor_eta[j];
        fac = exp(fac);
        v[i][0] *= fac;
        v[i][1] *= fac;
        v[i][2] *= fac;
        for (int j = 0; j < n_thermostats; j++)
          ke_local[j] += coupling->array_atom[i][j] * massone *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      }
    }
  } else if (which == BIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        temperature->remove_bias(i,v[i]);
        fac = 0;
        for (int j = 0; j < n_thermostats; j++)
          fac += is_real[j] * coupling->array_atom[i][j] * factor_eta[j];
        fac = exp(fac);
        v[i][0] *= fac;
        v[i][1] *= fac;
        v[i][2] *= fac;
        for (int j = 0; j < n_thermostats; j++)
          ke_local[j] += coupling->array_atom[i][j] * massone *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
        temperature->restore_bias(i,v[i]);
      }
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

/* ----------------------------------------------------------------------
   compute target temperature and kinetic energy
-----------------------------------------------------------------------*/

void FixNHMesh::compute_temp_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  for (int i = 0; i < n_thermostats; i++) {
    t_target[i] = t_start[i] + delta * (t_stop[i] - t_start[i]);
    ke_target[i] = tdof[i] * boltz * t_target[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixNHMesh::memory_usage()
{
  double bytes = 0.0;
  if (mesh_coupling_flag)
    bytes += n_thermostats*n_thermostats * sizeof(double);
  // Should stack-allocated arrays be included here too?
  return bytes;
}
