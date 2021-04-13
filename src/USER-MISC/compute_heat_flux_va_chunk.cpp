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

ComputeHeatFluxVAChunk::ComputeHeatFluxVAChunk(
    LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cvalues(nullptr), cvalues_local(nullptr)
{
  if (narg < 6) error->all(FLERR,"Illegal compute heat/flux/va/chunk command");


  MPI_Comm_rank(world,&me);

  id_ke = utils::strdup(arg[3]);
  id_pe = utils::strdup(arg[4]);

  // KE and PE computes

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  if (ike < 0 || ipe < 0)
    error->all(FLERR,"Could not find compute heat/flux/va/chunk compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,
        "Compute heat/flux/va/chunk compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");

  // Chunk compute

  id_chunk = utils::strdup(arg[5]);
  int ichunk = modify->find_compute(id_chunk);
  if (ichunk < 0) error->all(FLERR,
      "Chunk/atom compute does not exist for compute heat/flux/va/chunk");

  // Velocity bias

  biasflag = NOBIAS;
  c_temp = nullptr;
  if (narg > 6) {
    if (strcmp(arg[6],"bias")==0) {
      biasflag = BIAS;
      id_temp = utils::strdup(arg[7]);
      int itemp = modify->find_compute(id_temp);
      if (itemp < 0)
        error->all(FLERR,"Could not find compute heat/flux/va/chunk compute ID");
      c_temp = modify->compute[itemp];
      if (!c_temp->tempbias)
        error->all(FLERR,"Bias compute does not calculate a velocity bias");
    } else error->all(FLERR,"Illegal compute heat/flux/va/chunk command");
  }

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

  // Warn about idsflag == ONCE - invalid if atoms can move into different
  // chunks
  if (c_chunk->idsflag == ONCE)
    error->warning(FLERR,"Compute chunk/atom using ids once. This can cause "
                   "incorrect heat flux results if atoms move between chunks.");

  // 3D only since calculations are based on chunk volumes
  if (domain->dimension < 3)
    error->all(FLERR,
        "Compute heat/flux/va/chunk incompatible with simulation dimension");

  // Need atom map for chunk IDs of neighbors
  if (!atom->map_style)
    error->all(FLERR, "Compute heat/flux/va/chunk requires an atom map");

  c_ke = modify->compute[ike];
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

    if (force->bond!=nullptr)
      error->warning(FLERR,"compute heat/flux/va/chunk does not account for bond potentials");
    if (force->angle!=nullptr)
      error->warning(FLERR,"compute heat/flux/va/chunk does not account for angle potentials");
    if (force->dihedral!=nullptr)
      error->warning(FLERR,"compute heat/flux/va/chunk does not account for dihedral potentials");
    if (force->improper!=nullptr)
      error->warning(FLERR,"compute heat/flux/va/chunk does not account for improper potentials");
    if (force->kspace!=nullptr)
      error->warning(FLERR,"compute heat/flux/va/chunk does not account for kspace contributions");
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
  int idj;
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

  double *xi;
  double rij[3];
  double fij[3];
  double confpair[3];
  double pairdot;
  double *vi;

  // cfactor stores fraction of rij in each segment. Allocate nchunk+2 in case
  // rij crosses the full region of chunks, with xi and xj outside the region
  double cfactor[nchunk+2];
  //
  // c_ids calculation might store one extra value, so overallocate by 1
  int c_ids[nchunk+3];

  // Number of rij segments
  int n_intersect;

  if (biasflag) {
    if (c_temp->invoked_scalar != update->ntimestep) c_temp->compute_scalar();
    c_temp->remove_bias_all();
  }

  //Compute configurational contribution
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xi = x[i];
    vi = v[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Configurational component
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      idj = atom->map(atom->tag[j]);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // skip if neither i nor j are in group
      if (!(mask[i] & groupbit || mask[j] & groupbit)) continue;

      // Can't skip if (ichunk[i] == 0 && ichunk[idj] == 0)
      // in case rij still crosses a chunk.

      rij[0] = xi[0] - x[j][0];
      rij[1] = xi[1] - x[j][1];
      rij[2] = xi[2] - x[j][2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
      fij[0] = rij[0] * fpair;
      fij[1] = rij[1] * fpair;
      fij[2] = rij[2] * fpair;

      if (newton_pair || j < nlocal)
        pairdot = (fij[0]*(vi[0]+v[idj][0]) +
                   fij[1]*(vi[1]+v[idj][1]) +
                   fij[2]*(vi[2]+v[idj][2]))*0.5;
      else
        pairdot = (fij[0]*vi[0] + fij[1]*vi[1] + fij[2]*vi[2])*0.5;

      confpair[0] = rij[0] * pairdot;
      confpair[1] = rij[1] * pairdot;
      confpair[2] = rij[2] * pairdot;

      n_intersect = find_crossing(ichunk[i],ichunk[idj],xi,x[j],c_ids,cfactor);
      for (int ci=0; ci < n_intersect; ci++) {
        int c = c_ids[ci];
        if (c >= 0) {
          cvalues_local[c][0] += cfactor[ci] * confpair[0];
          cvalues_local[c][1] += cfactor[ci] * confpair[1];
          cvalues_local[c][2] += cfactor[ci] * confpair[2];
        }
      }
    }
  }

  // Kinetic component
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      int c = ichunk[i]-1;
      if (c >= 0) {
        double ei = c_ke->vector_atom[i] + c_pe->vector_atom[i];
        cvalues_local[c][3] += ei * v[i][0];
        cvalues_local[c][4] += ei * v[i][1];
        cvalues_local[c][5] += ei * v[i][2];
      }
    }

  if (biasflag) c_temp->restore_bias_all();

  for (m = 0; m < nchunk; m++) {
    cvalues_local[m][0] += cvalues_local[m][3];
    cvalues_local[m][1] += cvalues_local[m][4];
    cvalues_local[m][2] += cvalues_local[m][5];
  }

  MPI_Allreduce(cvalues_local[0],cvalues[0],nchunk*6, MPI_DOUBLE,MPI_SUM,world);
}

/*------------------------------------------------------------------------
  calculate fraction of vector rij in each chunk it crosses through
  -------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::find_crossing(
    int ci, int cj, double *xi_in, double *xj_in, int *c_ids, double *cfactor)
{
  double xi[3] = {xi_in[0], xi_in[1], xi_in[2]};
  double xj[3] = {xj_in[0], xj_in[1], xj_in[2]};
  if (c_chunk->scaleflag == REDUCED) {
    domain->x2lamda(xi, xi);
    domain->x2lamda(xj, xj);
  }

  double rij[3] = {xj[0]-xi[0], xj[1]-xi[1], xj[2]-xi[2]};

  if (c_chunk->which == ArgInfo::BINSPHERE)
    return crossing_binsphere(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BINCYLINDER)
    return crossing_bincylinder(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN1D)
    return crossing_bin1d(ci, cj, xi, rij, 0, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN2D)
    return crossing_bin2d(ci, cj, xi, rij, c_ids, cfactor);

  if (c_chunk->which == ArgInfo::BIN3D)
    return crossing_bin3d(ci, cj, xi, rij, c_ids, cfactor);

  error->all(FLERR,
      "Unsupported compute chunk/atom style. This should never be reached.");

  return 0;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin1d(int ci, int cj, double *xi,
    double *rij, int cdim, int *c_ids, double *cfactor)
{
  if (ci == cj && ci > 0) {
    cfactor[0] = 1.0;
    c_ids[0] = ci-1;
    return 1;
  }
  int n_intersect = 0;
  int dim = c_chunk->dim[cdim];
  int dir = rij[dim] > 0 ? 1 : -1;
  double delta = c_chunk->delta[cdim];
  double offset = c_chunk->offset[cdim];
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

  double xlast = xi[dim];
  double rijinv = 1/fabs(rij[dim]);
  double intersection;
  int i = ci-dir;
  do {
    i+=dir;
    if (i < 0) i += c_chunk->nlayers[cdim]+1;
    if (i > c_chunk->nlayers[cdim]) {
      if (periodicity) i -= c_chunk->nlayers[cdim]+1;
      else break;   // This shouldn't happen
    }

    if (i == cj)
      intersection = (xi[dim]+rij[dim]);
    else {
      if (i > 0) intersection = offset+delta*(i+(dir-1)/2);
      else if (dir < 0) intersection = offset+delta*c_chunk->nlayers[cdim];
      else /* (dir > 0) */ intersection = offset;
    }

    if (periodicity) {
      if (xlast < boxlo[dim]) xlast += prd[dim];
      if (xlast >= boxhi[dim]) xlast -= prd[dim];
    }

    cfactor[n_intersect] = fabs(xlast-intersection);
    if (periodicity && cfactor[n_intersect] > prd[dim]/2)
      cfactor[n_intersect] = fabs(prd[dim] - cfactor[n_intersect]);

    if (cfactor[n_intersect] > 0) {
      cfactor[n_intersect] *= rijinv;
      c_ids[n_intersect] = i-1;
      n_intersect++;
      xlast = intersection;
    }

  } while (i != cj);
  // (i != cj) assumes rij <= half prd if direction is periodic, but this is
  // already assumed in many other places so shouldn't be a problem.

  return n_intersect;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin2d(
    int ci, int cj, double *xi, double *rij, int *c_ids, double *cfactor)
{
  if (ci == cj && ci > 0) {
    cfactor[0] = 1.0;
    c_ids[0] = ci-1;
    return 1;
  }
  int *nlayers = c_chunk->nlayers;
  int c_ids1[nlayers[0]];
  double cfactor1[nlayers[0]];
  int dim1 = c_chunk->dim[0];
  int ci1, cj1, ci2, cj2;
  ci1 = (ci-1)/nlayers[1];
  ci2 = (ci-1) - ci1*nlayers[1];
  cj1 = (cj-1)/nlayers[1];
  ci1++; ci2++; cj1++;
  int n_intersect1 = crossing_bin1d(ci1,cj1,xi,rij,0,c_ids1,cfactor1);

  if (n_intersect1 == 0) return 0;

  int n_intersect = 0;
  int dim2 = c_chunk->dim[1];
  int periodicity = domain->periodicity[dim2];
  double boxlo, boxhi, prd;
  if (periodicity) {
    if (c_chunk->scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda[dim2];
      boxhi = domain->boxhi_lamda[dim2];
      prd = domain->prd_lamda[dim2];
    } else {
      boxlo = domain->boxlo[dim2];
      boxhi = domain->boxhi[dim2];
      prd = domain->prd[dim2];
    }
  }
  double xbin;
  double invdelta = 1/c_chunk->delta[dim2];
  double xi2[3], rij2[3];
  xi2[0] = xi[0]; xi2[1] = xi[1]; xi2[2] = xi[2];
  rij2[0] = rij2[1] = rij2[2] = 0;
  for (int i1 = 0; i1 < n_intersect1; i1++) {
    xi2[dim2] += rij2[dim2];
    rij2[dim2] = rij[dim2]*cfactor[i1];
    xbin = xi2[dim2]+rij2[dim2];
    if (periodicity) {
      if (xbin < boxlo) xbin += prd;
      if (xbin >= boxhi) xbin -= prd;
    }
    cj2 = static_cast<int>((xbin - c_chunk->offset[dim2])*invdelta);
    if (xbin < c_chunk->offset[dim2]) cj2--;
    cj2++;

    int n_new = crossing_bin1d(ci2, cj2, xi2, rij2, 1,
        &c_ids[n_intersect], &cfactor[n_intersect]);
    for (int n = n_intersect; n < n_intersect+n_new; n++) {
      cfactor[n] *= cfactor1[i1];
      if (c_ids[n] < 0 || c_ids1[i1] < 0)
        c_ids[n] = -1;
      else
        c_ids[n] += c_ids1[i1]*nlayers[1];
    }
    n_intersect += n_new;
    ci2 = cj2;
  }

  return n_intersect;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bin3d(
    int ci, int cj, double *xi, double *rij, int *c_ids, double *cfactor)
{
  if (ci == cj && ci > 0) {
    cfactor[0] = 1.0;
    c_ids[0] = ci-1;
    return 1;
  }
  int *nlayers = c_chunk->nlayers;
  int c_ids12[nlayers[0]*nlayers[1]];
  double cfactor12[nlayers[0]];
  int dim1 = c_chunk->dim[0];
  int ci12, cj12, ci3, cj3;
  ci12 = (ci-1)/nlayers[2];
  ci3 = (ci-1) - ci12*nlayers[2];
  cj12 = (cj-1)/nlayers[2];
  ci12++; ci3++; cj12++;
  int n_intersect12 = crossing_bin2d(ci12,cj12,xi,rij,c_ids12,cfactor12);

  if (n_intersect12 == 0) return 0;

  int n_intersect = 0;
  int dim3 = c_chunk->dim[2];
  int periodicity = domain->periodicity[dim3];
  double boxlo = domain->boxlo[dim3];
  double boxhi = domain->boxhi[dim3];
  double prd = domain->prd[dim3];
  double invdelta = 1/c_chunk->delta[dim3];
  double xbin;
  double xi3[3], rij3[3];
  xi3[0] = xi[0]; xi3[1] = xi[1]; xi3[2] = xi[2];
  rij3[0] = rij3[1] = rij3[2] = 0;
  for (int i12 = 0; i12 < n_intersect12; i12++) {
    xi3[dim3] += rij3[dim3];
    rij3[dim3] = rij[dim3]*cfactor[i12];
    xbin = xi3[dim3]+rij3[dim3];
    if (periodicity) {
      if (xbin < boxlo) xbin += prd;
      if (xbin >= boxhi) xbin -= prd;
    }
    cj3 = static_cast<int>(
        (xbin - c_chunk->offset[dim3])*invdelta);
    if (xbin < c_chunk->offset[dim3]) cj3--;
    cj3++;

    int n_new = crossing_bin1d(ci3, cj3, xi3, rij3, 1,
        &c_ids[n_intersect], &cfactor[n_intersect]);
    for (int n = n_intersect; n < n_intersect+n_new; n++) {
      cfactor[n] *= cfactor12[i12];
      if (c_ids[n] < 0 || c_ids12[i12] < 0)
        c_ids[n] = -1;
      else
        c_ids[n] += c_ids12[i12]*nlayers[2];
    }
    n_intersect += n_new;
    ci3 = cj3;
  }

  return n_intersect;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_binsphere(
    int ci, int cj, double *xi, double *rij, int *c_ids, double *cfactor)
{
  error->all(FLERR, "Spherical binning is currently unimplemented in "
      "compute heat/flux/va/chunk");
  return 0;
}

/*------------------------------------------------------------------------*/

int ComputeHeatFluxVAChunk::crossing_bincylinder(
    int ci, int cj, double *xi, double *rij, int *c_ids, double *cfactor)
{
  error->all(FLERR, "Cylinder binning is currently unimplemented in "
      "compute heat/flux/va/chunk");
  return 0;
}
