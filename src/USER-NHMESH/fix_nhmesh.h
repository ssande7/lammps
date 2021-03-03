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

/* -------------------------------------------------------------------------
TODO:
 * strip out barostat
 * change thermostat implementation
 * create compute for temperature in class constructor
 * add thermostat-thermostat coupling + input handling for it
 * handle velocity bias for SLLOD etc.
 * conservation of row sums in coupling matrix handled by conserve compute
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
NOTE:
 * fix_npt inherits FixNH directly, so no need to worry about things done for
   that
 * temperature of each thermostat comes from compute_temp_nhmesh
 * R matrix (for particle-thermostat couplings) from compute_coupling_nhmesh
 * Needs to perform time integration like fix_nh for proper thermostat
   "particle" dynamics
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
INPUT:
fix fix-ID grp-ID temp/nhmesh c_couple [t_start] [t_stop] [t_period] kwargs
  * c_couple     = compute that returns the coupling matrix mapping
                   particles to thermostats.
  * [t_start]    = vector of initial temperatures of length N
  * [t_stop]     = vector of final temperatures of length N
  * [t_period]   = vector of coupling time periods of length N
  * kwargs (optional):
     + couple args:
        args     = thermostat-thermostat coupling coefficients, in the format
                   M_11, M_12, ..., M_1N, M_21, ..., M_2N, ..., M_N1, ..., M_NN
                   these are the extent to which each thermostat controls each
                   other thermostat
     + tloop     = number of sub-cycles to perform in thermostat integration
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(temp/nhmesh,FixNHMesh)

#else

#ifndef LMP_FIX_NHMESH_H
#define LMP_FIX_NHMESH_H


#include "fix.h"                // IWYU pragma: export

namespace LAMMPS_NS {

class FixNHMesh : public Fix {
 public:
  FixNHMesh(class LAMMPS *, int, char **);
  virtual ~FixNHMesh();
  int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  double compute_scalar();            // Should return thermostat total energy
  virtual double compute_vector(int); // Returns vector of thermostat states
  void write_restart(FILE *);
  virtual int pack_restart_data(double *); // pack restart data
  virtual void restart(char *);
  int modify_param(int, char **);
  void reset_target(double);
  void reset_dt();
  virtual void *extract(const char*,int &);
  double memory_usage();

 protected:
  int which;
  double dtv,dtf,dthalf,dt4,dt8,dt16,dto;
  double boltz,nktv2p;
  double *tdof;         // dof controlled by each thermostat

  int n_thermostats;
  double *t_start, *t_stop;
  double *t_current, *t_target;
  double *ke_current, *ke_target;
  double *t_freq;

  int nlevels_respa;
  double *step_respa;

  char *id_coupling;
  class ComputeCouplingNHMesh *coupling;

  char *id_temp;
  class Compute *temperature;
  int tcomputeflag;                // 1 = compute was created by fix
                                   // 0 = created externally

  double **mesh_coupling;          // Thermostat-thermostat couplings
  int mesh_coupling_flag;          // 0/1 if couplings aren't/are defined
  double *mesh_dof;

  double *eta,*eta_dot;            // chain thermostat for particles
  double *eta_dotdot;
  double *eta_mass;

  int nc_tchain;                   // Number of sub-cycles
  double *factor_eta;

  int eta_mass_flag;               // 1 if eta_mass updated, 0 if not.

  void nhmesh_temp_integrate();

  virtual void nve_x();            // may be overwritten by child classes
  virtual void nve_v();
  virtual void nhmesh_v_temp();
  virtual void compute_temp_current();
  virtual void compute_temp_target();
  virtual int size_restart_global();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix temp/nhmesh command

Self-explanatory.

E: Number of thermostats for temp/nhmesh must be > 0

Self-explanatory.

E: Target temperature for fix temp/nhmesh cannot be 0.0

Self-explanatory.

E: Temperature ID for fix nhmesh does not exist

Self-explanatory.

E: Coupling ID of fix temp/nhmesh doesn't exist

Self-explanatory.

E: Invalid coupling compute for temp/nhmesh

Self-explanatory.

E: Thermostat couplings must be between 0 and 1

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

E: Fix_modify temperature ID is not compatible with nhmesh

Temperature compute must produce a global array of the kinetic energy tensors
of each thermostat, or at least the xx, yy, and zz components

*/
