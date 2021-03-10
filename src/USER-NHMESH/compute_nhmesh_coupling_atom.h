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
 compute ID grp-ID nhmesh/coupling/atom N heuristic
 * N          = number of thermostats (should match fix temp/nhmesh)
 * (units) arg = optional units:
                  box
                  lattice
                  reduced - unimplemented. Maybe not needed.
 * heuristic = thermostat partitioning heuristic: grid or points
    + grid xlo xhi ylo yhi zlo zhi nx ny nz decayx decayy decayz
      - xlo, xhi, ylo, yhi, zlo, zhi
                    = extents of grid. Replace xlo xhi with span to use
                    bounding box (at time of creation only, doesn't change as
                    system progresses.. yet..)
      - nx, ny, nz  = number of thermostats in each direction (min. 1)
      - decayx, decayy, decayz
                    = number of grid points away at which point influence
                      reaches 0. Default is 1, max is n<dir> (in which case the
                      influence does not decay in that direction).  Influence
                      is constant outside the region (equal to the border
                      value)
    + points (nofill) [p_1] ... [p_N-1] (decay_type)
       - nofill     = optional, set to exclude final thermostat that covers any
                      un- or partially-thermostatted particles
       - [p_x]      = vector of [x, y, z, rad] for each thermostat except the
                      last giving points and decay lengths for each. The last
                      thermostat will control any atoms not (fully) conrolled by
                      the other thermostats if nofill is not set. Use NULL  for
                      x or y or z to extend the point infinitely in that
                      direction (eg. as a line/plane)
       - decay_type = 'linear' or 'gaussian' or 'exp' or 'inv'
                      If the sum of thermostat influences on a given atom is > 1
                      it will be normalised.
          * linear    = (default) linear decay from 1 at the point to 0 at rad
          * gaussian  = gaussian distribution around point with std. of rad
          * exp       = exp(-rad*dist) where dist is distance from point
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(nhmesh/coupling/atom,ComputeNHMeshCouplingAtom)

#else

#ifndef LMP_COMPUTE_NHMESH_COUPLING_ATOM_H
#define LMP_COMPUTE_NHMESH_COUPLING_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeNHMeshCouplingAtom : public Compute {
 public:
  ComputeNHMeshCouplingAtom(class LAMMPS *, int, char **);
  virtual ~ComputeNHMeshCouplingAtom();
  void init();
  void setup();
  virtual void compute_peratom();

  virtual int get_n_thermostats() {return n_thermostats;}

 protected:
  int n_thermostats;    // Number of thermostats controlling the atoms
  int nmax;             // Max. number of atoms
  double *therm_sum;    // Sums of dofs controlled by each thermostat
  double **coupling;    // particle-thermostat couplings

  int fill_remainder;

  enum {
    GRID,
    POINTS
  } heuristic;

  
  enum{BOX,LATTICE,REDUCED} scaleflag;
  double scale[3];

  int grid_n[3];
  double **grid_pts;
  double grid_decay[3];
  double grid_dlength[3];
  double grid_lo[3], grid_hi[3];
  int grid_span[3];

  double **points;
  char ***points_str;
  int **points_varflag;
  int points_anyvar;
  enum {
    LINEAR,
    GAUSSIAN,
    EXP
  } points_decay;
  virtual void update_heuristics();
  virtual double calc_weight(double *, int&);

  // FixNHMesh needs access to therm_sum for efficiency reasons
  // could expose via public get method instead, but that seems less safe
  friend class FixNHMesh;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.

E: Illegal grid dimensions for nhmesh/coupling/atom - nx*ny*nz must equal N

Self-explanatory.

E: Illegal grid boundaries for nhmesh/coupling/atom

Self-explanatory.

E: Illegal grid boundaries for nhmesh/coupling/atom

Self-explanatory.

E: Number of thermostats must be > 0

Self-explanatory.

E: Unknown nhmesh/coupling/atom heuristic

Self-explanatory.

E: Compute nhmesh/coupling/atom points variables must be vector style

Self-explanatory.

E: Compute nhmesh/coupling/atom grid variables must be equal style

Self-explanatory.

E: Compute nhmesh/coupling/atom decay radii can not be NULL

Self-explanatory.

E: Unknown decay style for nhmesh/coupling/atom points command

Self-explanatory.

*/
