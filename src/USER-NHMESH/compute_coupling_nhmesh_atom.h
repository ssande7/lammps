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
INFO:
 * compute class to give fraction of influence of each thermostat on each
   atom
 * can use grid or point-based heuristics, with potential to plug in others
TODO:
 * give correct return from memory_usage()
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
INPUT:
 compute ID grp-ID nhmesh/coupling N heuristic
 * N          = number of thermostats (should match fix temp/nhmesh)
 * heuristic  = grid args or points args
    + grid xlo xhi ylo yhi zlo zhi nx ny nz decayx decayy decayz
      - xlo, xhi, ylo, yhi, zlo, zhi
                    = extents of grid
      - nx, ny, nz  = number of thermostats in each direction (min. 1)
      - decayx, decayy, decayz
                    = number of grid points away at which point influence
                      reaches 0. Default is 1, max is n<dir> (in which case the
                      influence does not decay in that direction).  Influence
                      is constant outside the region (equal to the border
                      value)
    + points [p_1] ... [p_N-1] (decay_type)
       - [p_x]      = vector of [x, y, z, rad] for each thermostat except the
                      last giving points and decay lengths for each. The last
                      thermostat will control any atoms not (fully) conrolled by
                      the other thermostats
       - decay_type = 'linear' or 'gaussian' or 'exp' or 'inv'
                      If the sum of thermostat influences on a given atom is > 1
                      it will be normalised.
          * linear    = (default) linear decay from 1 at the point to 0 at rad
          * gaussian  = gaussian distribution around point with std. of rad
          * exp       = exp(-rad*dist) where dist is distance from point
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(nhmesh/coupling,ComputeTemp)

#else

#ifndef LMP_COMPUTE_COUPLING_NHMESH_H
#define LMP_COMPUTE_COUPLING_NHMESH_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCouplingNHMesh : public Compute {
 public:
  ComputeCouplingNHMesh(class LAMMPS *, int, char **);
  virtual ~ComputeCouplingNHMesh();
  void init();
  void setup();
  virtual void compute_peratom();

 protected:
  int n_thermostats;    // Number of thermostats controlling the atoms
  int nmax;             // Max. number of atoms
  double *therm_sum;    // Sums of dofs controlled by each thermostat
  double **coupling;    // particle-thermostat couplings

  enum {
    GRID,
    POINTS
  } heuristic;

  int grid_n[3];
  int **grid_idtherm;
  double grid_decay[3];
  double grid_dlength[3];
  double grid_lo[3], grid_hi[3];

  double **points;
  char **points_str;
  int *points_varflag;
  int points_anyvar;
  enum {
    LINEAR,
    GAUSSIAN,
    EXP
  } points_decay;
  virtual void update_heuristics();
  virtual double calc_weight(double *, int&);

  // TODO: use dynamic_cast<ComputeCouplingNHMesh *> in fix to access therm_sum
  //       returns nullptr if fail
  friend class FixNHMesh;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.

E: Illegal grid dimensions for nhmesh/coupling - nx*ny*nz must equal N

Self-explanatory.

E: Illegal grid boundaries for nhmesh/coupling

Self-explanatory.

E: Illegal grid boundaries for nhmesh/coupling

Self-explanatory.

E: Number of thermostats must be > 0

Self-explanatory.

E: Unknown nhmesh/coupling heuristic

Self-explanatory.

E: Compute nhmesh/coupling points variables must be vector style

Self-explanatory.

E: Compute nhmesh/coupling point variables must return a vector of length 4

Self-explanatory.

*/
