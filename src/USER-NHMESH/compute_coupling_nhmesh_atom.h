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
 * compute class to give fraction of influence of each thermostat on each
   atom
 * can use grid or point-based heuristics, with potential to plug in others
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
INPUT:
 compute ID grp-ID nhmesh/coupling N heuristic
 * N          = number of thermostats (should match fix temp/nhmesh)
 * heuristic  = grid args or points args
    + grid region nx ny nz kwargs
      - region      = block style region outlining the area over which the grid
                      of thermostats is placed (thermostats can still influence
                      particles outside this region)
      - nx, ny, nz  = number of thermostats in each direction
      - kwargs:
         * decay argx argy argz:
            + arg<dir> = number of grid points away at which point influence
                         reaches 0. Default is 1, max is n<dir> (in which case
                         the influence does not decay in that direction).
                         Influence is constant outside the region (equal to the
                         border value)
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
          * inv       = 1/(dist)^rad
------------------------------------------------------------------------- */

