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
fix fix-ID grp-ID temp/nhmesh N [t_start] [t_stop] [t_period] c_couple kwargs
  * N            = Number of thermostats
  * [t_start]    = vector of initial temperatures of length N
  * [t_stop]     = vector of final temperatures of length N
  * [t_period]   = vector of coupling time periods of length N
  * c_couple     = compute that returns the coupling matrix mapping
                   particles to thermostats.
  * kwargs (optional):
     + couple args:
        args     = N vectors of coupling coefficients, corresponding to
                   rows of the thermostat-thermostat coupling matrix (ie.
                   degree to which each thermostat controls each other
                   thermostat)
     + conserve arg:
        arg      = on or off (default on)
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
  double dtv,dtf,dthalf,dt4,dt8,dto;
  double boltz,nktv2p,tdof;

  int n_thermostats;
  double t_start,t_stop;
  double t_current,t_target,ke_target;
  double t_freq;

  double drag,tdrag_factor;        // drag factor on particle thermostat

  int nlevels_respa;
  double *step_respa;

  char *id_temp;
  class Compute *temperature;
  int tcomputeflag;                // 1 = compute was created by fix
                                   // 0 = created externally

  double *eta,*eta_dot;            // chain thermostat for particles
  double *eta_dotdot;
  double *eta_mass;
  int mtchain;                     // length of chain
  int mtchain_default_flag;        // 1 = mtchain is default

  int nc_tchain;
  double factor_eta;

  int eta_mass_flag;               // 1 if eta_mass updated, 0 if not.

  void nhc_temp_integrate();

  virtual void nve_x();            // may be overwritten by child classes
  virtual void nve_v();
  virtual void nh_v_temp();
  virtual void compute_temp_target();
  virtual int size_restart_global();
};

}

#endif
#endif

/* ERROR/WARNING messages:

TODO

E: Number of thermostats for temp/nhmesh must be > 0

Self-explanatory.

E: Temperature ID for fix nhmesh does not exist

Self-explanatory.

*/
