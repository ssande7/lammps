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

#include "compute_coupling_nhmesh_atom.h"

#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCouplingNHMesh::ComputeCouplingNHMesh(LAMMPS *lmp, int narg, char **arg)
  : Compute(lmp, narg, arg), coupling(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal compute nhmesh/coupling command");

  int i;

  // defaults
  for (i = 0; i < 3; i++) grid_decay[i] = 1;

  n_thermostats = utils::inumeric(FLERR,arg[3],false,lmp);
  if (n_thermostats <= 0) error->all(FLERR,"Number of thermostats must be > 0");

  fill_remainder = 0;
  int iarg = 5;
  if (strcmp(arg[4], "grid") == 0) {
    heuristic = GRID;
    for (i = 0; i < 3; i++) {
      if (narg < iarg+1)
        error->all(FLERR,"Illegal compute nhmesh/coupling command");
      if (strcmp(arg[iarg], "span") == 0) {
        grid_lo[i] = domain->boxlo[i];
        grid_hi[i] = domain->boxhi[i];
        grid_span[i] = 1;
        iarg++;
      } else {
        if (narg < iarg+2)
          error->all(FLERR,"Illegal compute nhmesh/coupling command");
        int ivar = input->variable->find(arg[iarg]);
        if (ivar >= 0) {
          if (!input->variable->equalstyle(ivar))
            error->all(FLERR,
                "Compute nhmesh/coupling grid variables must be equal style");
          // TODO: evaluate during run to allow variable volume etc?
          grid_lo[i] = input->variable->compute_equal(ivar);
        } else grid_lo[i] = utils::numeric(FLERR,arg[iarg],false,lmp);
        ivar = input->variable->find(arg[iarg+1]);
        if (ivar >= 0) {
          if (!input->variable->equalstyle(ivar))
            error->all(FLERR,
                "Compute nhmesh/coupling grid variables must be equal style");
          // TODO: evaluate during run to allow variable volume etc?
          grid_hi[i] = input->variable->compute_equal(ivar);
        } grid_hi[i] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (grid_lo[i] > grid_hi[i])
          error->all(FLERR,"Illegal grid boundaries for nhmesh/coupling");
        grid_span[i] = 0;
        iarg += 2;
      }
    }
    int n_check = 1;
    if (narg < iarg+3)
      error->all(FLERR,"Illegal compute nhmesh/coupling command");
    for (i=0; i<3; i++) {
      int ivar = input->variable->find(arg[iarg]);
      if (ivar >= 0) {
        if (!input->variable->equalstyle(ivar))
          error->all(FLERR,
              "Compute nhmesh/coupling grid variables must be equal style");
        // TODO: evaluate during run to allow variable volume etc?
        grid_n[i] = input->variable->compute_equal(ivar);
        iarg++;
      } else grid_n[i] = utils::inumeric(FLERR,arg[iarg++],false,lmp);
      n_check *= grid_n[i];
      // Make sure final grid plane isn't a periodic
      // duiplicate of the first if spanning
      if (domain->periodicity[i] && grid_span[i] && grid_n[i] > 1)
        grid_hi[i] = (grid_n[i]-1) * grid_hi[i]/grid_n[i];
    }
    if (n_thermostats != n_check)
      error->all(FLERR,"Illegal grid dimensions for nhmesh/coupling - "
                       "nx*ny*nz must equal N");
    if (narg > iarg+1) {
      if (narg != iarg+3)
        error->all(FLERR,"Illegal compute nhmesh/coupling command");
      for (i=0; i<3; i++) {
        int ivar = input->variable->find(arg[iarg]);
        if (ivar >= 0) {
          if (!input->variable->equalstyle(ivar))
            error->all(FLERR,
                "Compute nhmesh/coupling grid variables must be equal style");
          // TODO: evaluate during run to allow variable volume etc?
          grid_decay[i] = input->variable->compute_equal(ivar);
          iarg++;
        } else grid_decay[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
        if (grid_decay[i] <= 0 || grid_decay[i] > grid_n[i])
          error->all(FLERR,"Illegal grid decay length for nhmesh/coupling");
      }
    }
    memory->create(grid_pts,n_thermostats,3,"nhmesh/coupling:grid_pts");
  } else if (strcmp(arg[4], "points") == 0) {
    heuristic = POINTS;
    memory->create(points,n_thermostats-1,4,"nhmesh/coupling:points");
    points_str = new char**[n_thermostats-1];
    points_varflag = new int*[n_thermostats-1];
    points_anyvar = 0;
    int j,ivar;
    fill_remainder = 1;
    if (iarg < narg && strcmp(arg[iarg], "nofill") == 0) {
      fill_remainder = 0;
      iarg++;
    }
    for (i=0; i<n_thermostats-1; i++) {
      if (narg < iarg+4)
        error->all(FLERR,"Illegal compute nhmesh/coupling command");
      points_varflag[i] = new int[4];
      points_str[i] = new char*[4];
      for (int j = 0; j < 4; j++) {
        ivar = input->variable->find(arg[iarg]);
        if (ivar >= 0) {
          points_varflag[i][j] = 1;
          points_anyvar = 1;
          if (!input->variable->equalstyle(ivar))
            error->all(FLERR,
                "Compute nhmesh/coupling points variables must be equal style");
          points_str[i][j] = utils::strdup(arg[iarg++]);
        } else {
          points_str[i][j] = nullptr;
          points_varflag[i][j] = 0;
          if (strcmp(arg[iarg], "NULL") == 0) {
            points[i][j] = NAN;
            iarg++;
          } else points[i][j] = utils::numeric(FLERR,arg[iarg++],false,lmp);
        }
      }
    }

    if (narg == iarg) points_decay = LINEAR;
    else if (narg != iarg+1)
      error->all(FLERR,"Illegal compute nhmesh/coupling command");
    else if (strcmp(arg[iarg],"linear")==0) points_decay = LINEAR;
    else if (strcmp(arg[iarg],"gaussian")==0) points_decay = GAUSSIAN;
    else if (strcmp(arg[iarg],"exp")==0) points_decay = EXP;
    else error->all(FLERR,"Unknown decay style for nhmesh/coupling points command");
    iarg++;
  } else error->all(FLERR,"Unknown nhmesh/coupling heuristic");

  peratom_flag = 1;
  size_peratom_cols = n_thermostats;

  nmax = 0;
  memory->create(therm_sum, n_thermostats, "nhmesh/coupling:therm_sum");
}

/* ---------------------------------------------------------------------- */

ComputeCouplingNHMesh::~ComputeCouplingNHMesh()
{
  if (!copymode) {
    memory->destroy(coupling);
    memory->destroy(therm_sum);

    if (heuristic == GRID) {
      memory->destroy(grid_pts);
    } else if (heuristic == POINTS) {
      memory->destroy(points);
      for (int i = 0; i < n_thermostats-1; i++) {
        for (int j = 0; j < 4; j++)
          if (points_str[i][j]) delete [] points_str[i][j];
        delete [] points_str[i];
        delete [] points_varflag[i];
      }
      delete [] points_str;
      delete [] points_varflag;
    }
    memory->destroy(coupling);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeCouplingNHMesh::init() {}

/* ---------------------------------------------------------------------- */

void ComputeCouplingNHMesh::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  update_heuristics();
}

/* ---------------------------------------------------------------------- */

void ComputeCouplingNHMesh::update_heuristics() {
  switch (heuristic) {
    case POINTS:
      if (points_anyvar) {
        double pt_var;
        int ivar;
        for (int i=0; i<n_thermostats-1; i++)
          for (int j = 0; j < 4; j++) {
            if (points_varflag[i][j]) {
              ivar = input->variable->find(points_str[i][j]);
              points[i][j] = input->variable->compute_equal(ivar);
            }
          }
      }
      break;
    case GRID:
      for (int i = 0; i < 3; i++) {
        if (grid_n[i] > 1)
          grid_dlength[i] = grid_decay[i]*(grid_hi[i]-grid_lo[i])/(grid_n[i]-1);
        else grid_dlength[i] = 0;
      }
      int gid[3];
      for (int j = 0; j < n_thermostats; j++) {
        gid[0] = j / (grid_n[1]*grid_n[2]);
        gid[1] = (j - gid[0]*grid_n[1]*grid_n[2]) / grid_n[2];
        gid[2] = j - gid[0]*grid_n[1]*grid_n[2] - gid[1]*grid_n[2];
        for (int i = 0; i < 3; i++) {
          if (grid_n[i] > 1)
            grid_pts[j][i] = grid_lo[i] + gid[i]*
              (grid_hi[i]-grid_lo[i])/(grid_n[i]-1);
          else
            grid_pts[j][i] = grid_lo[i];
        }
      }
      break;
  }
}

/* ---------------------------------------------------------------------- */

double ComputeCouplingNHMesh::calc_weight(double *x, int &j) {
  static const double r2pi = sqrt(2*MathConst::MY_PI);
  switch (heuristic) {
    case POINTS:
      {
        if (fill_remainder && j == n_thermostats-1) return 0;
        double pi[3];
        for (int i = 0; i < 3; i++) {
          if (std::isnan(points[j][i])) pi[i] = x[i];
          else pi[i] = points[j][i];
        }
        domain->remap_near(pi, x);
        double r2,dx,dy,dz;
        dx = fabs(x[0]-pi[0]);
        dy = fabs(x[1]-pi[1]);
        dz = fabs(x[2]-pi[2]);
        r2 = dx*dx + dy*dy + dz*dz;

        switch (points_decay) {
          case LINEAR:
            return (r2 < points[j][3]*points[j][3])
              ? 1-sqrt(r2)/points[j][3]
              : 0.0;
          case GAUSSIAN:
            return 1/(points[j][3]*r2pi)*exp(-r2/(2*points[j][3]*points[j][3]));
          case EXP:
            return exp(-points[j][3]*sqrt(r2));
        }
      }
      break;
    case GRID:
      {
        int i;
        double pi[3], dx;
        for (i = 0; i < 3; i++) {
          if (grid_dlength[i] != 0)
            pi[i] = grid_pts[j][i];
          else pi[i] = x[i];
        }
        domain->remap_near(pi, x);
        double wt = 1.0;
        for (i = 0; i < 3; i++) {
          dx = fabs(x[i] - pi[i]);
          if (grid_dlength[i] > 0 && dx > grid_dlength[i]) return 0.0;
          else if (grid_dlength[i] != 0.0) wt *= 1.0 - dx/grid_dlength[i];
        }
        return wt;
      }
      break;
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void ComputeCouplingNHMesh::compute_peratom()
{
  // grow coupling array if needed

  if (atom->nmax > nmax) {
    memory->destroy(coupling);
    nmax = atom->nmax;
    memory->create(coupling,nmax,n_thermostats,"nhmesh/coupling:atom");
    array_atom = coupling;
  }

  int i, j;

  invoked_peratom = update->ntimestep;

  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  double wt_sum;
  double local_sum[n_thermostats];
  for (j = 0; j < n_thermostats; j++) local_sum[j] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      wt_sum = 0;
      double * __restrict ci = coupling[i];
      for (j = 0; j < n_thermostats; j++) {
        ci[j] = calc_weight(x[i], j);
        wt_sum += ci[j];
      }
      if (wt_sum > 1) {
        for (j = 0; j < n_thermostats; j++) coupling[i][j] /= wt_sum;
      } else if (fill_remainder && wt_sum < 1)
        coupling[i][n_thermostats-1] = 1.0-wt_sum;
      for (j = 0; j < n_thermostats; j++) local_sum[j] += coupling[i][j];
    } else for (j = 0; j < n_thermostats; j++) coupling[i][j] = 0.0;
  // Calculating sum for each thermostat here to save an extra loop over atoms
  MPI_Allreduce(local_sum,therm_sum,n_thermostats,MPI_DOUBLE,MPI_SUM,world);
}
