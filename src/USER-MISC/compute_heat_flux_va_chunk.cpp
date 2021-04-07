/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/*------------------------------------------------------------------------
  Contributing Authors : Stephen Sanderson (UQ)
  --------------------------------------------------------------------------*/

#include "compute_heat_flux_va_chunk.h"

#include "arg_info.h"
#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};
enum{BOX,LATTICE,REDUCED};           // copied from compute_chunk_atom
enum{ONCE,NFREQ,EVERY};              // used in several files

/* ---------------------------------------------------------------------- */

ComputeHeatFluxVAChunk::ComputeHeatFluxVAChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cvalues(nullptr), cvalues_local(nullptr)
{
  if (narg < 6) error->all(FLERR,"Illegal compute stress/mop command");


  MPI_Comm_rank(world,&me);

  id_ke = utils::strdup(arg[3]);
  id_pe = utils::strdup(arg[4]);

  // KE and PE computes

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  if (ike < 0 || ipe < 0)
    error->all(FLERR,"Could not find compute heat/flux/va compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,"Compute heat/flux/va compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");

  // Chunk compute

  id_chunk = utils::strdup(arg[5]);
  int ichunk = modify->find_compute(id_chunk);
  if (ichunk < 0) error->all(FLERR,
      "Chunk/atom compute does not exist for compute heat/flux/va");

  // Velocity bias

  biasflag = NOBIAS;
  c_temp = nullptr;
  if (narg > 6)
    if (strcmp(arg[6],"bias")==0) {
      biasflag = BIAS;
      id_temp = utils::strdup(arg[7]);
      int itemp = modify->find_compute(id_temp);
      if (itemp < 0)
        error->all(FLERR,"Could not find compute heat/flux/va compute ID");
      c_temp = modify->compute[itemp];
      if (!c_temp->tempbias)
        error->all(FLERR,"Bias compute does not calculate a velocity bias");
    }

  c_chunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[ichunk]);
  if (c_chunk == nullptr)
    error->all(FLERR,"Compute heat/flux/va requires chunk/atom compute");

  // Need chunks with volume
  if (c_chunk->which == ArgInfo::TYPE || c_chunk->which == ArgInfo::MOLECULE ||
      c_chunk->which == ArgInfo::COMPUTE || c_chunk->which == ArgInfo::FIX ||
      c_chunk->which == ArgInfo::VARIABLE)
    error->all(FLERR,"Unsupported chunk/atom compute for compute heat/flux/va");

  // Can't discard empty chunks since they could still have some heat flux due
  // to the configurational component
  if (c_chunk->compress)
    error->all(FLERR,"Unsupported chunk/atom compute for compute heat/flux/va");

  // Warn about idsflag == ONCE - invalid if atoms can move into different
  // chunks
  if (c_chunk->idsflag == ONCE)
    error->warning(FLERR,"Compute chunk/atom using ids once. This can cause "
                   "incorrect heat flux results if atoms move between chunks.");

  // 3D only since calculations are based on chunk volumes
  if (domain->dimension < 3)
    error->all(FLERR,
        "Compute heat/flux/va incompatible with simulation dimension");

  c_ke = modify->compute[ike];
  c_pe = modify->compute[ipe];

  array_flag = 1;
  extarray = 0;
  size_array_cols = 6;
  size_array_rows = nchunk;
  size_array_rows_variable = 1;

  array = cvalues_local = nullptr;
  maxchunk = 0;
  nchunk = 1;
  allocate();
  cvalues = array;


}

/* ---------------------------------------------------------------------- */

ComputeHeatFluxVAChunk::~ComputeHeatFluxVAChunk()
{

  memory->destroy(array);
  memory->destroy(cvalues_local);

}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::allocate()
{
  size_array_rows = nchunk;
  if (nchunk > maxchunk) {
    maxchunk = nchunk;
    memory->grow(array,size_array_rows,size_array_cols,"heat/flux/va:array");
    memory->grow(cvalues_local,size_array_rows,size_array_cols,"heat/flux/va:cvalues_local");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::init()
{

  // Conversion constants

  // nktv2p = force->nktv2p;
  // ftm2v = force->ftm2v;

  // Error check

  // Compute stress/mop requires fixed simulation box
  // TODO: CHECK IF THIS IS A PROBLEM
  // if (domain->box_change_size || domain->box_change_shape || domain->deform_flag)
  //   error->all(FLERR, "Compute stress/mop requires a fixed simulation box");

  // This compute requires a pair style with pair_single method implemented

  if (force->pair == nullptr)
    error->all(FLERR,"No pair style is defined for compute heat/flux/va");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute heat/flux/va");

  // Warnings

  if (me==0) {

    //Compute stress/mop only accounts for pair interactions.
    // issue a warning if any intramolecular potential or Kspace is defined.

    if (force->bond!=nullptr)
      error->warning(FLERR,"compute heat/flux/va does not account for bond potentials");
    if (force->angle!=nullptr)
      error->warning(FLERR,"compute heat/flux/va does not account for angle potentials");
    if (force->dihedral!=nullptr)
      error->warning(FLERR,"compute heat/flux/va does not account for dihedral potentials");
    if (force->improper!=nullptr)
      error->warning(FLERR,"compute heat/flux/va does not account for improper potentials");
    if (force->kspace!=nullptr)
      error->warning(FLERR,"compute heat/flux/va does not account for kspace contributions");
  }

  // request an occasional half neighbor list
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::setup()
{
  nchunk = c_chunk->setup_chunks();
  allocate();
  c_chunk->bin_volumes();
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}


/* ----------------------------------------------------------------------
   compute output array
   ------------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::compute_array()
{
  invoked_array = update->ntimestep;

  // update chunk binning
  nchunk = c_chunk->setup_chunks();
  allocate();
  if (c_chunk->nchunkflag == EVERY) c_chunk->bin_volumes();

  compute_flux();

  // Pack output array
  int m,n;
  for (m=0; m<nchunk; m++) {
    for (n=0; n<3; n++) cvalues[m][n] += cvalues[m][n+3];
    for (n=0; n<6; n++) {
      if (  c_chunk->which == ArgInfo::BINSPHERE
         || c_chunk->which == ArgInfo::BINCYLINDER)
        cvalues[m][n] /= c_chunk->chunk_volume_vec[m];
      else cvalues[m][n] /= c_chunk->chunk_volume_scalar;
      // TODO: fix units
    }
  }

}


/*------------------------------------------------------------------------
  compute heat flux per chunk
  -------------------------------------------------------------------------*/

void ComputeHeatFluxVAChunk::compute_flux()
{
  int i,j,m,n,ii,jj,inum,jnum,itype,jtype;
  double rsq,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  for (m=0; m<nchunk; m++)
    for (n=0; n<6; n++)
      cvalues_local[m][n] = 0.0;

  neighbor->build_one(list);

  if (!(c_ke->invoked_flag & Compute::INVOKED_PERATOM)) {
    c_ke->compute_peratom();
    c_ke->invoked_flag |= Compute::INVOKED_PERATOM;
  }
  if (!(c_pe->invoked_flag & Compute::INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= Compute::INVOKED_PERATOM;
  }

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  c_chunk->compute_ichunk();
  int *ichunk = c_chunk->ichunk;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  double xi[3];
  double *xj;
  double rij[3];
  double fij[3];
  double confpair[3];
  double pairdot;
  double r;
  double *vi;
  double *vj;

  double cfactor[nchunk];
  int c_ids[nchunk];
  int n_intersect;

  if (biasflag) {
    if (c_temp->invoked_scalar != update->ntimestep) c_temp->compute_scalar();
    c_temp->remove_bias_all();
  }

  //Compute configurational contribution
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    vi = v[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Configurational component
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // skip if neither i nor j are in group
      if (!(mask[i] & groupbit || mask[j] & groupbit)) continue;

      // Can't skip if (ichunk[i] == 0 && ichunk[j] == 0)
      // in case rij still crosses a chunk.

      xj = x[j];
      rij[0] = xj[0];
      rij[1] = xj[1];
      rij[2] = xj[2];
      domain->remap_near(rij, xi);
      rij[0] -= xi[0];
      rij[1] -= xi[1];
      rij[2] -= xi[2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      r = sqrt(rsq);

      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
      fij[0] = rij[0]/r * fpair;
      fij[1] = rij[1]/r * fpair;
      fij[2] = rij[2]/r * fpair;

      if (newton_pair || j < nlocal) {
        // pair only stored once
        pairdot = (fij[0]*(vi[0]+vj[0]) +
                   fij[1]*(vi[1]+vj[1]) +
                   fij[2]*(vi[2]+vj[2]))/2;
      } else {
        // ghost atom, pair is repeated on another process
        pairdot = (fij[0]*vi[0] + fij[1]*vi[1] + fij[2]*vi[2])/2;
      }

      confpair[0] = rij[0] * pairdot;
      confpair[1] = rij[1] * pairdot;
      confpair[2] = rij[2] * pairdot;

      n_intersect = find_crossing(ichunk[i],ichunk[j],xi,rij,c_ids,cfactor);
      for (int ci=0, c=c_ids[ci]; c < n_intersect; c++) {
        cvalues_local[c][0] -= cfactor[ci] * confpair[0];
        cvalues_local[c][1] -= cfactor[ci] * confpair[1];
        cvalues_local[c][2] -= cfactor[ci] * confpair[2];
      }

    }
  }

  // Kinetic component
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      int ci = ichunk[i];
      if (ci > 0) {
        double ei = c_ke->vector_atom[i] + c_pe->vector_atom[i];
        cvalues_local[ci-1][0] += ei * v[i][0];
        cvalues_local[ci-1][1] += ei * v[i][1];
        cvalues_local[ci-1][2] += ei * v[i][2];
      }
    }


  if (biasflag) c_temp->restore_bias_all();

  MPI_Allreduce(cvalues_local[0],cvalues[0],nchunk*6, MPI_DOUBLE,MPI_SUM,world);
}

/*------------------------------------------------------------------------
  calculate fraction of vector rij in each chunk it crosses through
  -------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::find_crossing(
    int ci, int cj, double *xi, double *rij, int *c_ids, double *cfactor)
{
  if (c_chunk->scaleflag == REDUCED) {
    domain->x2lamda(xi,xi);
    domain->x2lamda(rij, rij);
  }
  // No need to convert back afterwards since these values are no longer needed

  if (c_chunk->which == ArgInfo::BINSPHERE)
    return crossing_binsphere(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BINCYLINDER)
    return crossing_bincylinder(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN1D)
    return crossing_bin1d(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN2D)
    return crossing_bin2d(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN3D)
    return crossing_bin3d(ci, cj, xi, rij, c_ids, cfactor);

  error->all(FLERR,
      "Unsupported compute chunk/atom style. This should never be reached.");

  return 0;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin1d(
    int ci, int cj, double *xi, double *rij, int *c_ids, double *cfactor)
{
  if (ci == cj && ci > 0) {
    cfactor[0] = 1.0;
    c_ids[0] = ci-1;
    return 1;
  }

  int n_intersect;
  int dim = c_chunk->dim[0];
  int dir = rij[dim] > 0 ? 1 : -1;
  double delta = c_chunk->delta[0];
  double offset = c_chunk->offset[0];
  double invdelta = c_chunk->invdelta[0];
  double ndelta = rij[dim]/delta;
  double x;
  double *boxlo,*boxhi,*prd;
  int periodicity = domain->periodicity[dim];
  if (periodicity) {
    if (c_chunk->scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  for (int i = 0; i <= ndelta; i++) {
    x = xi[dim] + i*delta*dir;
    if (periodicity) {
      if (x < boxlo[dim]) x += prd[dim];
      if (x >= boxhi[dim]) x -= prd[dim];
    }
    int ibin = static_cast<int>((x - offset)*invdelta);
    if (x < offset) ibin--;

    if (i == 0 && ibin >= 0) {
      cfactor[n_intersect] = fabs(offset+delta*ibin - x);
      if (dir < 0) cfactor[n_intersect] = delta - cfactor[n_intersect];
      cfactor[n_intersect] /= rij[dim];
      c_ids[n_intersect] = ibin;
      n_intersect++;
    } else if (i > ndelta-1 && ibin >= 0) {
      cfactor[n_intersect] = fabs(offset+delta*ibin - x);
      if (dir > 0) cfactor[n_intersect] = delta - cfactor[n_intersect];
      cfactor[n_intersect] /= rij[dim];
      c_ids[n_intersect] = ibin;
      n_intersect++;
    } else if (ibin >= 0) {
      cfactor[n_intersect] = fabs(delta/rij[dim]);
      c_ids[n_intersect] = ibin;
      n_intersect++;
    }
  }

  return n_intersect;
}
