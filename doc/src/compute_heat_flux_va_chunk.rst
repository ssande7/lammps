.. index:: compute heat/flux/va/chunk

compute heat/flux/va/chunk command
=========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID heat/flux/va/chunk pe-ID chunk-ID keyword args

* ID, group-ID are documented in :doc:`compute <compute>` command
* heat/flux/va/chunk = style name of this compute command
* pe-ID = ID of a compute that calculates per-atom potential energy
* chunk-ID = ID of a compute chunk/atom command
* zero or more keyword/arg pairs may be appended

   .. parsed-literal::

        *bias* arg = bias-ID
          bias-ID = ID of a temperature compute of style temp/chunk
        *kbias* arg = kbias-ID
          kbias-ID = ID of a temperature compute that calculates a velocity bias

Examples
""""""""

.. code-block:: LAMMPS

   compute myFlux all heat/flux/va/chunk myPE myChunks
   compute myFlux all heat/flux/va/chunk myPE myChunks bias myTempChunk
   compute myFlux all heat/flux/va/chunk myPE myChunks bias myTempChunk kbias myTempDeform

Description
"""""""""""

Define a computate that calculates the local heat flux vector in each chunk
using the volume averaged (VA) method as described by :ref:`(Smith) <Smith>`.
The heat flux in a chunk is calculated as

.. math::
   \mathbf{J}_c &= \frac{1}{V_c} \left[ \sum_{i \in c} e_i \mathbf{v}_i - \sum_{i} \sum_{j} \mathbf{r}_{ij} \mathbf{F}_{ij} \cdot \mathbf{v}_i l_{ij} \right] \\

The first term is the kinetic (convective) component, where :math:`e_i` is the
per-atom energy (potential and kinetic) as calculated by the computes *ke-ID*
and *pe-ID*, and the sum is over all atoms in chunk c. The second term is the
configurational component, and :math:`l_{ij}` is a selector function which gives
the fraction of the vector :math:`\mathbf{r}_{ij}` within the chunk.  Note that
this means atoms that are not inside a chunk can still contribute to that
chunk's heat flux, including atoms that would otherwise be excluded due to the
specifacation of a region in compute chunk/atom.  The configurational component
considers all atom pairs in which at least one atom is within the specified
compute group.

.. warning::

   Only pairwise interactions are currently supported in the calculation of the
   configurational component. All other forces, including bond, angle, dihedral,
   improper, and kspace, are ignored. Hence, this compute is not suitable for
   molecular systems.

.. note::

   The configurational component does not include any Lennard-Jones tail
   corrections added by the :doc:`pair_modify tail yes <pair_modify>`
   command, since those are contributions to the global system.

Removal of a velocity bias is supported via the *bias* keyword. It takes the
ID of a compute of type temp/chunk as an argument (see
:doc:`compute temp/chunk <compute_temp_chunk>`).  The temp/chunk compute should
use *chunk-ID* as its chunk/atom compute, and specify the "com yes" option.
With velocity bias removed, :math:`\mathbf{v}_i` in the the kinetic component of
the equation for :math:`\mathbf{J}_c` is replaced by the peculiar velocity, and
the velocity in the configurational component becomes
:math:`\left( \mathbf{v}_i - \mathbf{v}_c \right)`, where
:math:`\mathbf{v}_c` is the center of mass velocity of the chunk as calculated
by :doc:`compute temp/chunk <compute_temp_chunk>`.  An alternative bias can be
specified using the *kbias* keyword for the kinetic component only (for example,
it may be desirable in a SLLOD simulation to use
:doc:`compute temp/deform <compute_temp_deform>` to calculate peculiar
velocities).  If *kbias* is not specified but *bias* is, both components are
calculated using the same velocity bias.

.. note::

   The use of center of mass velocity is an approximation to the exact solution,
   as discussed in :ref:`(Smith) <Smith>`. The most accurate results will be
   obtained when there is minimal difference in the local streaming velocity
   within a chunk.

The compute takes two arguments which are IDs of other
:doc:`computes <compute>`.  One calculates per-atom  potential energy
(\ *pe-ID)*\ , and the other assigns chunk IDs to atoms (\ *chunk-ID*\ ). Unlike
:doc:`compute heat/flux <compute_heat_flux>`, this compute calculates kinetic
energy internally, accounting for a velocity bias if required.

.. note::

   These other computes, as well as temperature computes provided for bias
   calculation, should provide values for all the atoms in the group this
   compute specifies.  That means the other computes could use the same group as
   this compute, or they can just use group "all" (or any group whose atoms are
   a superset of the atoms in this compute's group).  LAMMPS does not check for
   this.

In LAMMPS, chunks are collections of atoms defined by a
:doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  Only chunks defined by spatial bins are
allowed (and not those defined by type or molecule), as the compute uses the
bin volume in the heat flux calculation.  See the :doc:`compute chunk/atom <compute_chunk_atom>`
doc page and the :doc:`Howto chunk <Howto_chunk>` doc page for details of how
chunks can be defined and examples of how they can be used to measure
properties of a system.

.. note::

   Unlike :doc:`compute heat/flux <compute_heat_flux>`, this compute includes
   volume in the heat flux calculation, and hence a 1/`:math:`{V}` scaling
   factor is not required in post-processing.

.. note::

   The :doc:`compute pe/atom <compute_pe_atom>` command has options for which
   terms to include in its calculation (pair, bond, etc).  The heat
   flux calculation will thus include exactly the same terms in the kinetic
   (convective) component. However, this has no influence on the calculation of
   the configurational component (see warning above).

See the documentation of :doc:`compute heat/flux <compute_heat_flux>` for a
discussion of how to calculate thermal conductivity.

Output info
"""""""""""

This compute calculates a global array of size *Nchunks* by 6.
The first 3 components of each row are the :math:`x`, :math:`y`, :math:`z`
components of the chunk's heat flux vector,
i.e. (:math:`J_x`, :math:`J_y`, :math:`J_z`).
The next 3 components are the :math:`x`, :math:`y`, :math:`z` components
of just the kinetic (convective) portion of the flux, i.e. the
first term in the equation for :math:`\mathbf{J}`.
Each component can be accessed by row indices 1-Nchunks and column indices 1-6.
These values can be used by any command that uses global array values from a
compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

The array values calculated by this compute are "intensive", meaning
they are independent of the number of atoms in the simulation.

The array values will be in energy/area/time :doc:`units <units>`.

Restrictions
""""""""""""

This command is part of the USER-MISC package. It is only enabled if LAMMPS is
built with that package.  See the :doc:`Build package <Build_package>` doc page
for more information.

Only two-body pair interactions are supported, as the pair->single() class
method is required.  All other interactions such as intra-molecular interactions
and long range (kspace) interactions are ignored.

2D simulations are not supported.

If a velocity bias is subtracted, the ID of the chunk/atom compute used by the
compute temp/chunk command (*bias-ID*) must match the ID supplied to this command as
*chunk-ID*.

Related commands
""""""""""""""""

:doc:`compute heat/flux <compute_heat_flux>`,
:doc:`compute chunk/atom <compute_chunk_atom>`,
:doc:`compute temp/chunk <compute_temp_chunk>`,
:doc:`fix ave/correlate <fix_ave_correlate>`,
:doc:`variable <variable>`

Default
"""""""

none

----------

.. _Smith:

**(Smith)**  Smith, Daivis, Todd, J. Chem. Phys. 150, 064103 (2019).
