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
#include "compute_temp_chunk.h"
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

ComputeHeatFluxVAChunk::ComputeHeatFluxVAChunk(
    LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cvalues(nullptr), cvalues_local(nullptr),
  cfactor(nullptr), c_ids(nullptr), c_wrap(nullptr)
{
  if (narg < 6) error->all(FLERR,"Illegal compute heat/flux/va/chunk command");


  MPI_Comm_rank(world,&me);

  id_pe = utils::strdup(arg[3]);

  // KE and PE computes

  int ipe = modify->find_compute(id_pe);
  if (ipe < 0)
    error->all(FLERR,"Could not find compute heat/flux/va/chunk compute ID");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");

  // Chunk compute

  id_chunk = utils::strdup(arg[4]);
  int ichunk = modify->find_compute(id_chunk);
  if (ichunk < 0) error->all(FLERR,
      "Chunk/atom compute does not exist for compute heat/flux/va/chunk");

  // Velocity bias

  biasflag = NOBIAS;
  c_temp_c = nullptr;
  c_temp_k = nullptr;
  int iarg = 5;
  while (narg > iarg) {
    if (strcmp(arg[iarg],"bias")==0) {
      biasflag = BIAS;
      id_temp_c = utils::strdup(arg[iarg+1]);
      int itemp = modify->find_compute(id_temp_c);
      if (itemp < 0)
        error->all(FLERR,"Could not find compute heat/flux/va/chunk compute ID");
      c_temp_c = dynamic_cast<ComputeTempChunk *>(modify->compute[itemp]);
      if (c_temp_c == nullptr)
        error->all(FLERR,"Bias must be a compute of type temp/chunk");
      if (!c_temp_c->comflag)
        error->all(FLERR,
            "temp/chunk compute ID does not calculate chunk velocity");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kbias")==0) {
      id_temp_k = utils::strdup(arg[iarg+1]);
      int itemp = modify->find_compute(id_temp_k);
      if (itemp < 0)
        error->all(FLERR,"Could not find compute heat/flux/va/chunk compute ID");
      c_temp_k = modify->compute[itemp];
      if (!c_temp_k->tempbias)
        error->all(FLERR,"kbias compute does not calculate a velocity bias");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute heat/flux/va/chunk command");
  }
  if (c_temp_k != nullptr && !biasflag)
    error->warning(FLERR,
        "Kinetic bias specified, without bias keyword. Bias will be ignored.");
  if (biasflag && c_temp_k == nullptr) c_temp_k = (Compute *)c_temp_c;

  c_chunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[ichunk]);
  if (c_chunk == nullptr)
    error->all(FLERR,"Compute heat/flux/va/chunk requires chunk/atom compute");

  // Need chunks with volume
  if (c_chunk->which == ArgInfo::TYPE || c_chunk->which == ArgInfo::MOLECULE ||
      c_chunk->which == ArgInfo::COMPUTE || c_chunk->which == ArgInfo::FIX ||
      c_chunk->which == ArgInfo::VARIABLE)
    error->all(FLERR,
        "Unsupported chunk/atom compute for compute heat/flux/va/chunk");

  // Can't discard empty chunks since they could still have some heat flux due
  // to the configurational component
  if (c_chunk->compress)
    error->all(FLERR,
        "Unsupported chunk/atom compute for compute heat/flux/va/chunk");

  // 3D only since calculations are based on chunk volumes
  if (domain->dimension < 3)
    error->all(FLERR,
        "Compute heat/flux/va/chunk incompatible with simulation dimension");

  // Need atom map for chunk IDs of neighbors
  if (!atom->map_style)
    error->all(FLERR, "Compute heat/flux/va/chunk requires an atom map");

  // Warn about idsflag == ONCE - invalid if atoms can move into different
  // chunks
  if (c_chunk->idsflag == ONCE)
    error->warning(FLERR,"Compute chunk/atom using ids once. This can cause "
                   "incorrect heat flux results if atoms move between chunks.");

  if (c_chunk->regionflag)
    error->warning(FLERR,"Compute heat/flux/va/chunk may produce incorrect "
        "results when compute chunk/atom defines a region if that region "
        "is not aligned to chunk boundaries.");

  c_pe = modify->compute[ipe];

  array_flag = 1;
  extarray = 0;
  size_array_cols = 6;
  size_array_rows = nchunk;
  size_array_rows_variable = 1;

  array = nullptr;
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

  memory->destroy(cfactor);
  memory->destroy(c_ids);
  if (domain->deform_vremap) memory->destroy(c_wrap);

}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::allocate()
{
  size_array_rows = nchunk;
  if (nchunk > maxchunk) {
    maxchunk = nchunk;
    memory->grow(array,size_array_rows,
        size_array_cols,"heat/flux/va/chunk:array");
    memory->grow(cvalues_local,size_array_rows,
        size_array_cols,"heat/flux/va/chunk:cvalues_local");
    cvalues = array;

    // cfactor stores fraction of rij in each segment. Allocate nchunk+2 in case
    // rij crosses the full region of chunks, with xi and xj outside the region
    memory->grow(cfactor, nchunk+2, "heat/flux/va/chunk:cfactor");
    memory->grow(c_ids, nchunk+2, "heat/flux/va/chunk:c_ids");
    if (domain->deform_vremap)
      memory->grow(c_wrap, nchunk+2, 3, "heat/flux/va/chunk:c_wrap");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::init()
{

  // Conversion constants

  // nktv2p = force->nktv2p;
  // ftm2v = force->ftm2v;

  // Error check

  // This compute requires a pair style with pair_single method implemented

  if (force->pair == nullptr)
    error->all(FLERR,"No pair style is defined for compute heat/flux/va/chunk");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute heat/flux/va/chunk");

  // Warnings

  if (me==0) {

    // Only accounts for pair interactions at this stage.
    // Issue a warning if any intramolecular potential or Kspace is defined.

    if (force->bond!=nullptr) error->warning(FLERR,
          "compute heat/flux/va/chunk does not account for bond potentials");
    if (force->angle!=nullptr) error->warning(FLERR,
        "compute heat/flux/va/chunk does not account for angle potentials");
    if (force->dihedral!=nullptr) error->warning(FLERR,
        "compute heat/flux/va/chunk does not account for dihedral potentials");
    if (force->improper!=nullptr) error->warning(FLERR,
        "compute heat/flux/va/chunk does not account for improper potentials");
    if (force->kspace!=nullptr) error->warning(FLERR,
        "compute heat/flux/va/chunk does not account for kspace contributions");
  }

  // Request an occasional full neighbor list
  // Need full list since half list would mean calculation requires
  // velocity and chunk id of neighbors, which breaks under MPI
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVAChunk::setup()
{
  // Check that chunk definitions match here since compute temp/chunk looks for
  // cchunk in init()
  if (biasflag && c_temp_c->cchunk != c_chunk)
    error->all(FLERR,"Chunk/atom compute used by bias must match chunk "
                     "compute of compute heat/flux/va/chunk");
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
  if (biasflag) {
    // update compute temp/chunk as well if it's used so vcm_compute() can be
    // called without the overhead of compute_scalar()
    c_temp_c->nchunk = nchunk;
    if (c_temp_c->nchunk > c_temp_c->maxchunk) c_temp_c->allocate();
  }
  if (c_chunk->nchunkflag == EVERY) c_chunk->bin_volumes();

  compute_flux();

  // Pack output array
  int m,n;
  for (m=0; m<nchunk; m++) {
    for (n=0; n<6; n++) {
      if (  c_chunk->which == ArgInfo::BINSPHERE
         || c_chunk->which == ArgInfo::BINCYLINDER)
        cvalues[m][n] /= c_chunk->chunk_volume_vec[m];
      else cvalues[m][n] /= c_chunk->chunk_volume_scalar;
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

  if (!(c_pe->invoked_flag & Compute::INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= Compute::INVOKED_PERATOM;
  }

  c_chunk->compute_ichunk();
  int *ichunk = c_chunk->ichunk;

  if (biasflag) {
    // Only need vcm, so no need to call compute_scalar, etc.
    if (c_temp_c->comstep != update->ntimestep) c_temp_c->vcm_compute();

    // compute temp/chunk just uses vcmall for bias removal, so no need to
    // recompute for c_temp_k with compute_scalar() if it's the same as c_temp_c
    if (c_temp_k != c_temp_c && c_temp_k->invoked_scalar != update->ntimestep)
      c_temp_k->compute_scalar();
  }

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  double rij[3];
  double fij[3];
  double confpair[3];
  double pairdot, pairdot_c = 0;
  double vi[3];
  double vj[3];

  int n_seg;

  // Compute kinetic component
  double massone;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double ei;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      int c = ichunk[i]-1;
      if (c >= 0) {
        if (biasflag) c_temp_k->remove_bias(i, v[i]);

        // Need peculiar kinetic energy
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        ei = massone*(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])
             + c_pe->vector_atom[i];

        cvalues_local[c][3] += ei * v[i][0];
        cvalues_local[c][4] += ei * v[i][1];
        cvalues_local[c][5] += ei * v[i][2];

        if (biasflag) c_temp_k->restore_bias(i, v[i]);
      }
    }

  //Compute configurational contribution
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    vi[0] = v[i][0];
    vi[1] = v[i][1];
    vi[2] = v[i][2];
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

      rij[0] = x[i][0] - x[j][0];
      rij[1] = x[i][1] - x[j][1];
      rij[2] = x[i][2] - x[j][2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
      fij[0] = rij[0] * fpair;
      fij[1] = rij[1] * fpair;
      fij[2] = rij[2] * fpair;

      pairdot = (fij[0]*vi[0] + fij[1]*vi[1] + fij[2]*vi[2])*0.5;

      confpair[0] = rij[0] * pairdot;
      confpair[1] = rij[1] * pairdot;
      confpair[2] = rij[2] * pairdot;

      n_seg = find_crossing(x[i],x[j]);
      for (int ci=0; ci < n_seg; ci++) {
        int c = c_ids[ci];
        if (c >= 0) {
          if (biasflag) {
            double vcm[3] = {c_temp_c->vcmall[c][0],
                             c_temp_c->vcmall[c][1],
                             c_temp_c->vcmall[c][2]};
            if (domain->deform_vremap) {
              double *h_rate = domain->h_rate;
              // TODO: check wrap_i and wrap_j are correct
              vcm[0] += h_rate[0]*c_wrap[ci][0] +
                        h_rate[5]*c_wrap[ci][1] +
                        h_rate[4]*c_wrap[ci][2];
              vcm[1] += h_rate[1]*c_wrap[ci][1] + h_rate[3]*c_wrap[ci][2];
              vcm[2] += h_rate[2]*c_wrap[ci][2];
            }
            pairdot_c = -(fij[0]*vcm[0] + fij[1]*vcm[1] + fij[2]*vcm[2])*0.5;
          }
          cvalues_local[c][0] += cfactor[ci]*(confpair[0] + rij[0]*pairdot_c);
          cvalues_local[c][1] += cfactor[ci]*(confpair[1] + rij[1]*pairdot_c);
          cvalues_local[c][2] += cfactor[ci]*(confpair[2] + rij[2]*pairdot_c);
        }
      }
    }
  }

  for (m = 0; m < nchunk; m++) {
    cvalues_local[m][0] += cvalues_local[m][3];
    cvalues_local[m][1] += cvalues_local[m][4];
    cvalues_local[m][2] += cvalues_local[m][5];
  }

  MPI_Allreduce(cvalues_local[0],cvalues[0],nchunk*6, MPI_DOUBLE,MPI_SUM,world);
}

/*------------------------------------------------------------------------
  Calculate fraction of vector rij in each chunk it crosses through, and
  whether reaching that chunk requires wrapping around a periodic boundary.
  Output stored in c_ids, cfactor and c_wrap
  -------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::find_crossing(double *xi_in, double *xj_in)
{
  double xi[3] = {xi_in[0], xi_in[1], xi_in[2]};
  double xj[3] = {xj_in[0], xj_in[1], xj_in[2]};
  if (c_chunk->scaleflag == REDUCED) {
    domain->x2lamda(xi, xi);
    domain->x2lamda(xj, xj);
  }

  double rij[3] = {xj[0]-xi[0], xj[1]-xi[1], xj[2]-xi[2]};

  if (c_chunk->which == ArgInfo::BINSPHERE)
    return crossing_binsphere(xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BINCYLINDER)
    return crossing_bincylinder(xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN1D)
    return crossing_bin1d(xi, rij, 0, c_ids, cfactor, c_wrap);

  if (c_chunk->which == ArgInfo::BIN2D)
    return crossing_bin2d(xi, rij, c_ids, cfactor, c_wrap);

  if (c_chunk->which == ArgInfo::BIN3D)
    return crossing_bin3d(xi, rij, c_ids, cfactor, c_wrap);

  error->all(FLERR,
      "Unsupported compute chunk/atom style. This should never be reached.");

  return 0;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin1d(
    double *xi, double *rij, int cdim, int *c_ids, double *cfactor, int **c_wrap)
{
  int dim = c_chunk->dim[cdim];
  int ci, cj;
  double xiremap = xi[dim];
  double xjremap = xiremap+rij[dim];
  int periodicity = domain->periodicity[dim];
  double *boxlo,*boxhi,*prd;
  int wrap_i = 0, wrap_j = 0;
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
    if (xiremap < boxlo[dim]) {xiremap += prd[dim]; wrap_i = -1;}
    if (xiremap >= boxhi[dim]) {xiremap -= prd[dim]; wrap_i = 1;}
    if (xjremap < boxlo[dim]) {xjremap += prd[dim]; wrap_j = -1;}
    if (xjremap >= boxhi[dim]) {xjremap -= prd[dim]; wrap_j = 1;}
  }
  double offset = c_chunk->offset[cdim];
  double invdelta = c_chunk->invdelta[cdim];
  int nlayers = c_chunk->nlayers[cdim];

  // Get chunks of end points
  // Can't rely on c_chunk->ichunk since atom may not have a chunk id if it's
  // from a different group or owned by a different MPI process
  ci = static_cast<int>((xiremap - offset) * invdelta) + 1;
  if (xiremap < offset) ci--;
  if (ci > nlayers || ci < 0) ci = 0;
  cj = static_cast<int>((xjremap - offset) * invdelta) + 1;
  if (xjremap < offset) cj--;
  if (cj > nlayers || cj < 0) cj = 0;

  // skip if in the same chunk, catching edge cases
  double delta = c_chunk->delta[cdim];
  double rijabs = fabs(rij[dim]);
  int dir = rij[dim] > 0 ? 1 : -1;
  double rijinv = 1/rijabs;
  int n_seg = 0;
  if (ci == cj) {
    if (ci == 0) {
      // Could still cross chunks while not being in one
      double upper, lower;
      if (dir > 0) lower = xiremap;
      else {
        lower = xiremap + rij[dim];
        // lower always inside simulation box
        if (periodicity && lower < boxlo[dim]) lower += prd[dim];
      }
      upper = lower + rijabs;
      bool wrapped =
        (periodicity && lower <= offset+prd[dim] && upper >= offset+prd[dim]);
      if (lower <= offset && upper > offset || wrapped) {
        // All chunks between xi and xj
        if (wrapped) offset += prd[dim];
        if (dir > 0) cfactor[n_seg] = (offset - lower)*rijinv;
        else cfactor[n_seg] = (upper - offset+nlayers*delta)*rijinv;
        c_ids[n_seg++] = -1;
        for (; n_seg <= nlayers; n_seg++) {
          cfactor[n_seg] = delta*rijinv;
          c_ids[n_seg] = n_seg-1;
          if (domain->deform_vremap) c_wrap[n_seg][dim] = dir;
        }
        if (dir > 0) cfactor[n_seg] = (upper - offset+nlayers*delta)*rijinv;
        else cfactor[n_seg] = (offset - lower)*rijinv;
        c_ids[n_seg++] = -1;
        return n_seg;
      }

    // nlayers == 1 is tricky to handle efficiently
    // in the loop below, so catch it here
    } else if (nlayers == 1 && periodicity && prd[dim] > nlayers*delta) {
      // Check if wrapping through unchunked zone
      if (dir > 0 && xiremap + rij[dim] >= boxhi[dim]) {
        cfactor[0] = (offset+nlayers*delta - xiremap)*rijinv;
        c_ids[0] = ci-1;
        cfactor[2] = (xjremap - offset)*rijinv;
        c_ids[2] = ci-1;
        if (domain->deform_vremap) {
          c_wrap[0][dim] = wrap_i;
          c_wrap[2][dim] = wrap_j;
        }
        cfactor[1] = 1 - cfactor[0] - cfactor[2];
        c_ids[1] = -1;
        return 3;
      } else if (dir < 0 && xiremap + rij[dim] < boxlo[dim]) {
        cfactor[0] = (xiremap - offset)*rijinv;
        c_ids[0] = ci-1;
        cfactor[2] = (offset+nlayers*delta - xjremap)*rijinv;
        c_ids[2] = ci-1;
        if (domain->deform_vremap) {
          c_wrap[0][dim] = wrap_i;
          c_wrap[2][dim] = wrap_j;
        }
        cfactor[1] = 1 - cfactor[0] - cfactor[2];
        c_ids[1] = -1;
        return 3;
      }
    }

    // Doesn't cross any other chunks
    cfactor[0] = 1.0;
    c_ids[0] = ci-1;
    if (domain->deform_vremap) c_wrap[0][dim] = wrap_i;
    return 1;
  }

  double xcovered = 0;
  double intersection;
  int i = ci - dir;
  int cur_wrap = wrap_i;

  do {
    i += dir;
    // No periodicity check for this since it needs to happen in non-periodic
    // systems when chunks don't cover full system
    if (i < 0) i = nlayers;
    if (i > nlayers) i = 0;

    if (i > 0) intersection = offset+delta*(i+(dir-1)/2);
    else if (dir < 0) intersection = offset+delta*nlayers;
    else /* (dir > 0) */ intersection = offset;

    // probably more optimisation to do here. cfactor is just delta * rijinv
    // for all chunks except first, last, and any section of unchunked space
    cfactor[n_seg] = fabs(xiremap-intersection);
    if (periodicity && cfactor[n_seg] > prd[dim]/2) {
      cfactor[n_seg] = fabs(prd[dim] - cfactor[n_seg]);
      cur_wrap += dir;
    }

    if (i == cj && cfactor[n_seg] > rijabs - xcovered) {
      cfactor[n_seg] = (rijabs - xcovered) * rijinv;
      c_ids[n_seg] = i-1;
      if (domain->deform_vremap) c_wrap[n_seg][dim] = cur_wrap;
      n_seg++;
      break;
    }

    if (cfactor[n_seg] > 0) {
      xcovered += cfactor[n_seg];
      cfactor[n_seg] *= rijinv;
      c_ids[n_seg] = i-1;
      if (domain->deform_vremap) c_wrap[n_seg][dim] = cur_wrap;
      n_seg++;
      xiremap = intersection;
    }

  } while (i != cj);

  return n_seg;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin2d(
    double *xi, double *rij, int *c_ids, double *cfactor, int **c_wrap)
{
  int *nlayers = c_chunk->nlayers;
  int *c_ids1 = new int[nlayers[0]];
  int **c_wrap1;
  if (domain->deform_vremap)
    memory->create(c_wrap1, nlayers[0], 3, "heat/flux/va/chunk:c_wrap2d");
  double *cfactor1 = new double[nlayers[0]];
  int n_seg1 = crossing_bin1d(xi,rij,0,c_ids1,cfactor1,c_wrap1);

  if (n_seg1 == 0) return 0;

  int n_seg = 0;
  int dim1 = c_chunk->dim[0];
  int dim2 = c_chunk->dim[1];
  double xi2[3] = {xi[0], xi[1], xi[2]};
  double rij2[3] = {0.0, 0.0, 0.0};
  for (int i1 = 0; i1 < n_seg1; i1++) {
    xi2[dim2] += rij2[dim2];
    rij2[dim2] = rij[dim2]*cfactor1[i1];
    int n_new=crossing_bin1d(xi2,rij2,1,&c_ids[n_seg],&cfactor[n_seg],&c_wrap[n_seg]);

    for (int n = n_seg; n < n_seg+n_new; n++) {
      cfactor[n] *= cfactor1[i1];
      if (domain->deform_vremap) c_wrap[n][dim1] = c_wrap1[i1][dim1];
      if (c_ids[n] < 0 || c_ids1[i1] < 0)
        c_ids[n] = -1;
      else
        c_ids[n] += c_ids1[i1]*nlayers[1];
    }
    n_seg += n_new;
  }

  delete [] c_ids1;
  delete [] cfactor1;
  if (domain->deform_vremap) memory->destroy(c_wrap1);
  return n_seg;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin3d(
    double *xi, double *rij, int *c_ids, double *cfactor, int **c_wrap)
{
  int *nlayers = c_chunk->nlayers;
  int *c_ids12 = new int[nlayers[0]*nlayers[1]];
  int **c_wrap12;
  if (domain->deform_vremap)
    memory->create(c_wrap12, nlayers[0]*nlayers[1], 3, "heat/flux/va/chunk:c_wrap3d");
  double *cfactor12 = new double[nlayers[0]*nlayers[1]];
  int n_seg12 = crossing_bin2d(xi,rij,c_ids12,cfactor12,c_wrap12);

  if (n_seg12 == 0) return 0;

  int n_seg = 0;
  int dim1 = c_chunk->dim[0];
  int dim2 = c_chunk->dim[1];
  int dim3 = c_chunk->dim[2];
  double xi3[3], rij3[3];
  xi3[0] = xi[0]; xi3[1] = xi[1]; xi3[2] = xi[2];
  rij3[0] = rij3[1] = rij3[2] = 0;
  for (int i12 = 0; i12 < n_seg12; i12++) {
    xi3[dim3] += rij3[dim3];
    rij3[dim3] = rij[dim3]*cfactor12[i12];

    int n_new=crossing_bin1d(xi3,rij3,2,&c_ids[n_seg],&cfactor[n_seg],&c_wrap[n_seg]);
    for (int n = n_seg; n < n_seg+n_new; n++) {
      cfactor[n] *= cfactor12[i12];
      if (domain->deform_vremap) {
        c_wrap[n][dim1] = c_wrap12[i12][dim1];
        c_wrap[n][dim2] = c_wrap12[i12][dim2];
      }
      if (c_ids[n] < 0 || c_ids12[i12] < 0)
        c_ids[n] = -1;
      else
        c_ids[n] += c_ids12[i12]*nlayers[2];
    }
    n_seg += n_new;
  }

  delete [] c_ids12;
  delete [] cfactor12;
  if (domain->deform_vremap) memory->destroy(c_wrap12);
  return n_seg;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_binsphere(
    double *xi, double *rij, int *c_ids, double *cfactor)
{
  error->all(FLERR, "Spherical binning is currently unimplemented in "
      "compute heat/flux/va/chunk");
  return 0;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bincylinder(
    double *xi, double *rij, int *c_ids, double *cfactor)
{
  error->all(FLERR, "Cylinder binning is currently unimplemented in "
      "compute heat/flux/va/chunk");
  return 0;
}
