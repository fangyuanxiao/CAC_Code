#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "comm_brick.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "math_extra.h"
#include "universe.h"
#include "error.h"
#include "memory.h"
#include "update.h"

using namespace CAC_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define BIG 1.0e20

enum{SINGLE,MULTI};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};

CommBrick::CommBrick(CAC *cac) : Comm(cac)
{
   style = 0;
   layout = LAYOUT_UNIFORM;
   init_buffers();
}

/*----------------------------------------------------------------------------*/

CommBrick::~CommBrick()
{
  free_swap();
  if (sendlist) for (int i = 0; i < maxswap; i++) memory->destroy(sendlist[i]);
  memory->sfree(sendlist);
  if (esendlist) for (int i = 0; i < maxswap; i++) memory->destroy(esendlist[i]);
  memory->sfree(esendlist);
  memory->destroy(maxsendlist);
  memory->destroy(emaxsendlist);
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ----------------------------------------------------------------------
 *    initialize comm buffers and other data structs local to CommBrick
 *    ------------------------------------------------------------------------- */
void CommBrick::init_buffers()
{
  multilo = multihi = NULL;
  cutghostmulti = NULL;
  // bufextra = max size of one exchanged atom
  //      = allowed overflow of sendbuf in exchange()
  //atomvec, fix reset these 2 maxexchange values if needed
  //only necessary if their size > BUFEXTRA

  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
  maxswap = 6;
  allocate_swap(maxswap);
  sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"comm:sendlist");
  esendlist= (int **) memory->smalloc(maxswap*sizeof(int *),"comm:esendlist");

  memory->create(maxsendlist,maxswap,"comm:maxsendlist");
  memory->create(emaxsendlist,maxswap,"comm:emaxsendlist");

  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
    emaxsendlist[i] = BUFMIN;
    memory->create(sendlist[i],BUFMIN,"comm:sendlist[i]");
    memory->create(esendlist[i],BUFMIN,"comm:sendlist[i]");
  }
}

/* ----------------------------------------------------------------------
 *    allocation of swap info
 *    ------------------------------------------------------------------------- */

void CommBrick::allocate_swap(int n)
{
  memory->create(sendnum,n,"comm:sendnum");
  memory->create(esendnum,n,"comm:esendnum");
  memory->create(recvnum,n,"comm:recvnum");
  memory->create(erecvnum,n,"comm:erecvnum");
  memory->create(sendproc,n,"comm:sendproc");
  memory->create(recvproc,n,"comm:recvproc");
  memory->create(size_forward_recv,n,"comm:size");
  memory->create(esize_forward_recv,n,"comm:esize");
  memory->create(size_reverse_send,n,"comm:size");
  memory->create(esize_reverse_send,n,"comm:esize");
  memory->create(size_reverse_recv,n,"comm:size");
  memory->create(esize_reverse_recv,n,"comm:esize");
  memory->create(slablo,n,"comm:slablo");
  memory->create(slabhi,n,"comm:slabhi");
  memory->create(firstrecv,n,"comm:firstrecv");
  memory->create(efirstrecv,n,"comm:efirstrecv");
  memory->create(pbc_flag,n,"comm:pbc_flag");
  memory->create(pbc,n,6,"comm:pbc");
}

/* ----------------------------------------------------------------------
 *    free memory for swaps
 *    ------------------------------------------------------------------------- */
void CommBrick::free_swap()
{
  memory->destroy(sendnum);
  memory->destroy(esendnum);
  memory->destroy(recvnum);
  memory->destroy(erecvnum);
  memory->destroy(sendproc);
  memory->destroy(recvproc);
  memory->destroy(size_forward_recv);
  memory->destroy(esize_forward_recv);
  memory->destroy(size_reverse_send);
  memory->destroy(esize_reverse_send);
  memory->destroy(size_reverse_recv);
  memory->destroy(esize_reverse_recv);
  memory->destroy(slablo);
  memory->destroy(slabhi);
  memory->destroy(firstrecv);
  memory->destroy(efirstrecv);
  memory->destroy(pbc_flag);
  memory->destroy(pbc);
}

/* ---------------------------------------------------------------------- */

void CommBrick::init()
{
  Comm::init();

  // memory for multi-style communication
  if (mode == MULTI && multilo == NULL) {
    // if (me==0) fprintf(screen, "testing mode MULTI= %d\n", atom->ntypes+element->ntypes+1);
     //allocate_multi(maxswap);
     //memory->create(cutghostmulti,element->ntypes+1,3,"comm:cutghostmulti");   
  }
  if (cac->element_flag)
  {
    allocate_multi(maxswap);
    memory->create(cutghostmulti, element->ntypes+1,3,"comm:cutghostmulti");
    //error->all(FLERR, "testing from comm init about element flag\n");
  } 
  //if (mode == SINGLE && multilo) {
    //if (me==0) fprintf(screen, "testing = %d\n", 1);
  //} 
}

/* ----------------------------------------------------------------------
   allocation of multi-type swap info
------------------------------------------------------------------------- */

void CommBrick::allocate_multi(int n)
{
  multilo = memory->create(multilo,n,element->ntypes+1,"comm:multilo");
  multihi = memory->create(multihi,n,element->ntypes+1,"comm:multihi");
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single mode sets slab boundaries (slablo,slabhi) based on max cutoff
   multi mode sets type-dependent slab boundaries (multilo,multihi)
------------------------------------------------------------------------- */

void CommBrick::setup()
{
  // cutghost[] = max distance at which ghost atoms need to be acquired
  // for orthogonal:
  //   cutghost is in box coords = neigh->cutghost in all 3 dims
  // for triclinic:
  //   neigh->cutghost = distance between tilted planes in box coords
  //   cutghost is in lamda coords = distance between those planes
  // for multi:
  //   cutghostmulti = same as cutghost, only for each atom type

  int i,j;
  int ntypes = atom->ntypes;
  int netypes = element->ntypes;
  double *prd,*sublo,*subhi;

  double cut = MAX(neighbor->cutneighmax,cutghostuser);
  double ele_size = element->max_size;
  double cut_all = neighbor->cut_all;
  //neighbor->cut_all = cut_all;
  //neighbor->cut_allsq = cut_all*cut_all; 

  if (triclinic == 0) {
    prd = domain->prd;
    sublo = domain->sublo;
    subhi = domain->subhi;
    cutghost[0] = cutghost[1] = cutghost[2] = cut_all;
    if (cac->element_flag)
    {
      double *size = element->size;
      for (i=1; i<=netypes; i++)
        cutghostmulti[i][0] = cutghostmulti[i][1] = cutghostmulti[i][2] = 
          cut+size[i];
    }
    
  }
  
  // recvneed[idim][0/1] = # of procs away I recv atoms from, within cutghost
  //   0 = from left, 1 = from right
  //   do not cross non-periodic boundaries, need[2] = 0 for 2d
  // sendneed[idim][0/1] = # of procs away I send atoms to
  //   0 = to left, 1 = to right
  //   set equal to recvneed[idim][1/0] of neighbor proc
  // maxneed[idim] = max procs away any proc recvs atoms in either direction
  // layout = UNIFORM = uniform sized sub-domains:
  //   maxneed is directly computable from sub-domain size
  //     limit to procgrid-1 for non-PBC
  //   recvneed = maxneed except for procs near non-PBC
  //   sendneed = recvneed of neighbor on each side
  // layout = NONUNIFORM = non-uniform sized sub-domains:
  //   compute recvneed via updown() which accounts for non-PBC
  //   sendneed = recvneed of neighbor on each side
  //   maxneed via Allreduce() of recvneed

  int *periodicity = domain->periodicity;
  int left,right;

  if (layout == LAYOUT_UNIFORM) {
    maxneed[0] = static_cast<int> (cutghost[0] * procgrid[0] / prd[0]) + 1;
    maxneed[1] = static_cast<int> (cutghost[1] * procgrid[1] / prd[1]) + 1;
    maxneed[2] = static_cast<int> (cutghost[2] * procgrid[2] / prd[2]) + 1;


    if (domain->dimension == 2) maxneed[2] = 0;
    if (!periodicity[0]) maxneed[0] = MIN(maxneed[0],procgrid[0]-1);
    if (!periodicity[1]) maxneed[1] = MIN(maxneed[1],procgrid[1]-1);
    if (!periodicity[2]) maxneed[2] = MIN(maxneed[2],procgrid[2]-1);

    if (!periodicity[0]) {
      recvneed[0][0] = MIN(maxneed[0],myloc[0]);
      recvneed[0][1] = MIN(maxneed[0],procgrid[0]-myloc[0]-1);
      left = myloc[0] - 1;
      if (left < 0) left = procgrid[0] - 1;
      sendneed[0][0] = MIN(maxneed[0],procgrid[0]-left-1);
      right = myloc[0] + 1;
      if (right == procgrid[0]) right = 0;
      sendneed[0][1] = MIN(maxneed[0],right);
    } else recvneed[0][0] = recvneed[0][1] =
             sendneed[0][0] = sendneed[0][1] = maxneed[0];

    if (!periodicity[1]) {
      recvneed[1][0] = MIN(maxneed[1],myloc[1]);
      recvneed[1][1] = MIN(maxneed[1],procgrid[1]-myloc[1]-1);
      left = myloc[1] - 1;
      if (left < 0) left = procgrid[1] - 1;
      sendneed[1][0] = MIN(maxneed[1],procgrid[1]-left-1);
      right = myloc[1] + 1;
      if (right == procgrid[1]) right = 0;
      sendneed[1][1] = MIN(maxneed[1],right);
    } else recvneed[1][0] = recvneed[1][1] =
             sendneed[1][0] = sendneed[1][1] = maxneed[1];

    if (!periodicity[2]) {
      recvneed[2][0] = MIN(maxneed[2],myloc[2]);
      recvneed[2][1] = MIN(maxneed[2],procgrid[2]-myloc[2]-1);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = MIN(maxneed[2],procgrid[2]-left-1);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = MIN(maxneed[2],right);
    } else recvneed[2][0] = recvneed[2][1] =
             sendneed[2][0] = sendneed[2][1] = maxneed[2];
  }

  // allocate comm memory
  nswap = 2 * (maxneed[0]+maxneed[1]+maxneed[2]);
  if (nswap > maxswap) grow_swap(nswap);

  // setup parameters for each exchange:
  // sendproc = proc to send to at each swap
  // recvproc = proc to recv from at each swap
  // for mode SINGLE:
  //   slablo/slabhi = boundaries for slab of atoms to send at each swap
  //   use -BIG/midpt/BIG to insure all atoms included even if round-off occurs
  //   if round-off, atoms recvd across PBC can be < or > than subbox boundary
  //   note that borders() only loops over subset of atoms during each swap
  //   treat all as PBC here, non-PBC is handled in borders() via r/s need[][]
  // for mode MULTI:
  //   multilo/multihi is same, with slablo/slabhi for each atom type
  // pbc_flag: 0 = nothing across a boundary, 1 = something across a boundary
  // pbc = -1/0/1 for PBC factor in each of 3/6 orthogonal/triclinic dirs
  // for triclinic, slablo/hi and pbc_border will be used in lamda (0-1) coords
  // 1st part of if statement is sending to the west/south/down
  // 2nd part of if statement is sending to the east/north/up

  int dim,ineed;

  int iswap = 0;
  for (dim = 0; dim < 3; dim++) {
    for (ineed = 0; ineed < 2*maxneed[dim]; ineed++) {
      pbc_flag[iswap] = 0;
      pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] =
        pbc[iswap][3] = pbc[iswap][4] = pbc[iswap][5] = 0;
      if (ineed % 2 == 0) {
        sendproc[iswap] = procneigh[dim][0];
        recvproc[iswap] = procneigh[dim][1];
        if (mode == SINGLE) {
          if (ineed < 2) slablo[iswap] = -BIG;
          else slablo[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
          slabhi[iswap] = sublo[dim] + cutghost[dim];
        } 
        if (cac->element_flag) {
          for (i = 1; i <= netypes; i++) {
            if (ineed < 2) multilo[iswap][i] = -BIG;
            else multilo[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
            multihi[iswap][i] = sublo[dim] + cutghostmulti[i][dim];
          }
        }
        if (myloc[dim] == 0) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = 1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = 1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = 1;
          }
        }
      } else {
        sendproc[iswap] = procneigh[dim][1];
        recvproc[iswap] = procneigh[dim][0];
        if (mode == SINGLE) {
          slablo[iswap] = subhi[dim] - cutghost[dim];
          if (ineed < 2) slabhi[iswap] = BIG;
          else slabhi[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
        } 
        if (cac->element_flag) {
          for (i = 1; i <= netypes; i++) {
            multilo[iswap][i] = subhi[dim] - cutghostmulti[i][dim];
            if (ineed < 2) multihi[iswap][i] = BIG;
            else multihi[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
          }
          //if (me==0) fprintf(screen, "testing from comm setup\n");
        }
        if (myloc[dim] == procgrid[dim]-1) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = -1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = -1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = -1;
          }
        }
      }
      iswap++;
    }
  } 
}

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
------------------------------------------------------------------------- */

void CommBrick::grow_swap(int n)
{
  free_swap();
  allocate_swap(n);
  if (cac->element_flag) {
     free_multi();
     allocate_multi(n);
  }
  sendlist = (int **)
    memory->srealloc(sendlist,n*sizeof(int *),"comm:sendlist");
  esendlist = (int **)
    memory->srealloc(esendlist,n*sizeof(int *),"comm:esendlist");
  memory->grow(maxsendlist,n,"comm:maxsendlist");
  memory->grow(emaxsendlist,n,"comm:emaxsendlist");
  for (int i = maxswap; i < n; i++) {
    maxsendlist[i] = BUFMIN;
    emaxsendlist[i] = BUFMIN;
    memory->create(sendlist[i],BUFMIN,"comm:sendlist[i]");
    memory->create(esendlist[i],BUFMIN,"comm:esendlist[i]");
  }
  maxswap = n; 
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside a stencil proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommBrick::exchange()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2,nlocal;
  double lo,hi,value;
  double **x;
  double *sublo,*subhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  
  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  //if (me==0) fprintf(screen, "map_style = %d\n", map_style);
  
  atom->nghost = 0;
  element->nghost = 0;
  
  atom->avec->clear_bonus();

  // insure send buf is large enough for single atom
  // bufextra = max size of one atom = allowed overflow of sendbuf
  // fixes can change per-atom size requirement on-the-fly
 
  int bufextra_old = bufextra;
  //if (comm->me==0) fprintf(screen, "bufextra = %d\n", bufextra);

  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA; 
  //if (comm->me==0) fprintf(screen, "bufextra = %d\n", bufextra);
  //error->all(FLERR, "testing from comm exchange\n");
  if (bufextra > bufextra_old)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");

  // subbox bounds for orthogonal or triclinic

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  // loop over dimensions

  int dimension = domain->dimension;

  for (int dim = 0; dim < dimension; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom
    
    x = atom->x;
    lo = sublo[dim];
    hi = subhi[dim];
    nlocal = atom->nlocal;
    i = nsend = 0;

    while (i < nlocal) {
      if (x[i][dim] < lo || x[i][dim] >= hi) {
        if (nsend > maxsend) grow_send(nsend,1);
        nsend += avec->pack_exchange(i,&buf_send[nsend]);
        avec->copy(nlocal-1,i,1);
        nlocal--;
      } else i++;
    }
    atom->nlocal = nlocal;
    
    // send/recv atoms in both directions
    // send size of message first so receiver can realloc buf_recv if needed
    // if 1 proc in dimension, no send/recv
    //   set nrecv = 0 so buf_send atoms will be lost
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    if (procgrid[dim] == 1) nrecv = 0;
    else {
      MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
                   &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,MPI_STATUS_IGNORE);
      nrecv = nrecv1;
      if (procgrid[dim] > 2) {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
                     &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,MPI_STATUS_IGNORE);
        nrecv += nrecv2;
      }
      if (nrecv > maxrecv) grow_recv(nrecv);

      MPI_Irecv(buf_recv,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
                world,&request);
      MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
      MPI_Wait(&request,MPI_STATUS_IGNORE);

      if (procgrid[dim] > 2) {
        MPI_Irecv(&buf_recv[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
                  world,&request);
        MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
      }
    }
   
    // check incoming atoms to see if they are in my box
    // if so, add to my list
    // box check is only for this dimension,
    //   atom may be passed to another proc in later dims

    m = 0;
    while (m < nrecv) {
      value = buf_recv[m+dim+1];
      if (value >= lo && value < hi) m += avec->unpack_exchange(&buf_recv[m]);
      else m += static_cast<int> (buf_recv[m]);
    } 
  }

  // added for element exchange 
  // the basic idea is the same as the exchange of the atoms
  //
  if (cac->element_flag)
  {
    //if (comm->me==0) fprintf(screen, "testing cac->element_flag\n");
    for (int dim = 0; dim < dimension; dim++) {
      // fill buffer with elements leaving my box, using < and >=
      // when element is deleted, fill it with last element

      x = element->x;
      lo = sublo[dim];
      hi = subhi[dim];
      nlocal = element->nlocal;
      i = nsend = 0;
      //if (comm->me==0) fprintf(screen, "maxsend = %d\n",maxsend);
      //error->all(FLERR, "testing the maxsend\n");
      while (i < nlocal) {
        if (x[i][dim]<lo || x[i][dim] >=hi) {
          if (nsend > maxsend) grow_send(nsend,1);
          nsend +=element->pack_exchange(i,&buf_send[nsend]);
          element->copy(nlocal-1,i,1);
          nlocal--; 

        } else i++;
      }
      element->nlocal = nlocal;

      // send/recv elements in both directions
      // send size of message first so receiver can realloc buf_recv if needed
      // if 1 proc in dimension, no send/recv
      //   set nrecv = 0 so buf_send elements will be lost
      // if 2 procs in dimension, single send/recv
      // if more than 2 procs in dimension, send/recv to both neighbors

      if (procgrid[dim] == 1) nrecv = 0;
      else {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
                   &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,MPI_STATUS_IGNORE);
        nrecv = nrecv1;
        if (procgrid[dim]>2) {
          MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
                     &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,MPI_STATUS_IGNORE);
          nrecv +=nrecv2;
        }
        if (nrecv > maxrecv) grow_recv(nrecv);

        MPI_Irecv(buf_recv,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
                world,&request);
        MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        if (procgrid[dim] > 2) {
          MPI_Irecv(&buf_recv[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
                    world,&request);
          MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
          MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
      }
    
    // check incoming elements to see if they are in my box
    // if so, add to my list
    // box check is only for this dimension,
    //   element may be passed to another proc in later dims

    m = 0;
    while (m < nrecv) {
      value = buf_recv[m+dim+1];
      if (value >=lo && value < hi) m += element->unpack_exchange(&buf_recv[m]);
      else m += static_cast<int> (buf_recv[m]);
    } 
    }
    //error->all(FLERR, "testing from comm exchange\n");
  }

  bigint n = element->nlocal;
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (sum != element->nelements)
    error->all(FLERR,"Element lost");
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommBrick::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrick::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommBrick::borders()
{
  int i,n,itype,iswap,dim,ineed,twoneed;
  int nsend,nrecv,sendflag,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;
  
  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {
      
      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in sendlist for use in future timesteps

      x = atom->x;
      if (mode == SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }
      
      nsend = 0;

      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requires data
      // e.g. due to non-PBC or non-uniform sub-domains

      if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
      else sendflag = 1;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus only atoms in bordergroup
      // can only limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (mode == SINGLE) {
            for (i = nfirst; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            for (i = nfirst; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }

        }else {
          if (mode == SINGLE) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            for (i = atom->nlocal; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        }
      }
      
      // pack up list of border atoms

      if (nsend*size_border > maxsend) grow_send(nsend*size_border,0);
      //if (comm->me==0) fprintf(screen, "ghost_velocity = %d\n", ghost_velocity);     
      if (ghost_velocity)
        n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,
                                  pbc_flag[iswap],pbc[iswap]);
      else
        n = avec->pack_border(nsend,sendlist[iswap],buf_send,
                              pbc_flag[iswap],pbc[iswap]);

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
        if (nrecv*size_border > maxrecv) grow_recv(nrecv*size_border);
        if (nrecv) MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
                             recvproc[iswap],0,world,&request);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

      // unpack buffer

      if (ghost_velocity)
        avec->unpack_border_vel(nrecv,atom->nlocal+atom->nghost,buf);
      else
        avec->unpack_border(nrecv,atom->nlocal+atom->nghost,buf);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }     
  }

  // insure send/recv buffers are long enough for all forward & reverse comm
  
  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);
  
  if (cac->element_flag) {
  // do swaps for elements

  iswap = 0;
  esmax = ermax = 0;

  for (dim = 0; dim<3;dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {
      
      // find elements within slab boundaries lo/hi using <= and >=
      // check elements between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent elements indices in sendlist for use in future timesteps

      x = element->x;
      type = element->type;
      mlo = multilo[iswap];
      mhi = multihi[iswap];
      if (ineed % 2 ==0) {
        nfirst = nlast;
        nlast = element->nlocal + element->nghost;
      } 

      nsend = 0;
      
      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requiers data
      // e.g. due to non-PBC or non-uniform sub-domains

      if (ineed/2 >= sendneed[dim][ineed%2]) sendflag = 0;
      else sendflag = 1;

      // find send elements according to SINGLE vs MULTI
      // all elements eligible versus only elements in bordergroup
      // can only limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned and ghost

      if (sendflag) {
        if (!bordergroup || ineed >=2) {
          for (i = nfirst; i<nlast; i++) {
            itype = type[i];
            if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
              // fprintf(screen, "me = %d, i = %d, dim =%d\n", me, i, dim);
              if (nsend == emaxsendlist[iswap]) egrow_list(iswap, nsend);
              esendlist[iswap][nsend++] = i;
            }
          }
        }
      }

     // pack up list of border elements

      //if (comm->me==0) fprintf(screen, "size_border = %d\n", size_border);
      //if (me==0) fprintf(screen, "maxsend = %d\n", maxsend);
      if (nsend*esize_border > maxsend) grow_send(nsend*esize_border,0);
      //if (me==0) fprintf(screen, "ghost_velocity = %d\n", ghost_velocity);
      //if (ghost_velocity)
      //  n = element->pack_border_vel(nsend,sendlist[iswap],buf_send,
     //                                pbc_flag[iswap],pbc[iswap]);
     n = element->pack_border(nsend, esendlist[iswap],buf_send,
                              pbc_flag[iswap],pbc[iswap]);

     //fprintf(screen, "nsend = %d, me = %d, pbc_flag = %d\n", nsend, me, pbc_flag[iswap]);
      
      // swap elements with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my element arrays
      // if swapping with self, simply copy, no messages
       
      if(sendproc[iswap] !=me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
        if (nrecv*esize_border > maxrecv) grow_recv(nrecv*esize_border);
        if (nrecv) MPI_Irecv(buf_recv,nrecv*esize_border,MPI_DOUBLE,
                             recvproc[iswap],0,world,&request);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

      // unpack buffer

      element->unpack_border(nrecv, element->nlocal+element->nghost,buf);
      
      esmax = MAX(esmax,nsend);
      ermax = MAX(ermax,nrecv);
      esendnum[iswap] = nsend;
      erecvnum[iswap] = nrecv;
      esize_forward_recv[iswap] = nrecv*esize_forward;
      esize_reverse_send[iswap] = nrecv*esize_reverse;
      esize_reverse_recv[iswap] = nsend*esize_reverse;
      efirstrecv[iswap] = element->nlocal+element->nghost;
      element->nghost += nrecv;
      iswap++;
    } 
  }
 } 

 // insure send/recv buffers are long enough for all forward & reverse comm
 int emax = MAX(emaxforward*esmax,emaxreverse*ermax);
 if (emax > maxsend) grow_send(emax,0);
 emax = MAX(emaxforward*ermax,emaxreverse*esmax);
 if (emax > maxrecv) grow_recv(emax);

 // After borders communication, check if the atoms in each element need to be interpolated to put them into atom x
  
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrick::grow_list(int iswap, int n)
{
  maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist[iswap],maxsendlist[iswap],"comm:sendlist[iswap]");
}

/*----------------------------------------------------------------------
  similar function for growing esendlist
------------------------------------------------------------------------*/
void CommBrick::egrow_list(int iswap, int n)
{
  emaxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(esendlist[iswap], emaxsendlist[iswap],"comm:esendlist[iswap]");
   
}

/*------------------------------------------------------------------------------
  forward communication invoked by a pair
  nsize used only to set recv buffer limit
--------------------------------------------------------------------------------*/

void CommBrick::forward_comm_pair(Pair *pair)
{
   int iswap,n;
   double *buf;
   MPI_Request request;

   int nsize = pair->comm_forward;

   for (iswap = 0; iswap<nswap; iswap++) {
     
     n = pair->pack_forward_comm(sendnum[iswap],sendlist[iswap],
		                 buf_send,pbc_flag[iswap],pbc[iswap]);

     // exchange with another proc
     // if self, set recv buffer to send buffer

     if (sendproc[iswap] != me) {
	if (recvnum[iswap])
          MPI_Irecv(buf_recv, nsize*recvnum[iswap],MPI_DOUBLE,
	            recvproc[iswap],0,world,&request);
	if (sendnum[iswap])
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (recvnum[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
	buf = buf_recv;
     } else buf = buf_send;

     // unpack buffer
     
     pair->unpack_forward_comm(recvnum[iswap],firstrecv[iswap],buf);
   }
}

/*------------------------------------------------------------------------
  forward communication of elements invoked by a pair
  nsize used only to set recv buffer limit
--------------------------------------------------------------------------*/

void CommBrick::eforward_comm_pair(Pair *pair)
{
  int iswap,n;
  double *buf;
  MPI_Request request;

  int nsize = pair->ecomm_forward;
  if (comm->me==0&&update->ntimestep>334) fprintf(screen, "nsize = %d\n",nsize);

 // error->all(FLERR,"test3 from eforward comm pair");
  for (int iswap = 0; iswap < nswap; iswap++) {
   if (sendproc[iswap] != me) {
     if (erecvnum[iswap])
       MPI_Irecv(buf_recv, erecvnum[iswap]*nsize, MPI_DOUBLE, 
		 recvproc[iswap],0,world,&request);
     if (comm->me==0 && update->ntimestep > 334) fprintf(screen, "test1 from eforward\n"); 
     MPI_Barrier(world);
     n = pair->epack_forward_comm(esendnum[iswap],esendlist[iswap],buf_send);
     MPI_Barrier(world);
     if (comm->me==0 && update->ntimestep >334) fprintf(screen,"test4 from eforward\n");
     if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
     MPI_Barrier(world);
     if (erecvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
     MPI_Barrier(world);
     if (comm->me==0 && update->ntimestep > 334) fprintf(screen,"test2 from eforward\n");
     pair->eunpack_forward_comm(erecvnum[iswap],efirstrecv[iswap],buf_recv);
     if (comm->me==0 && update->ntimestep > 334) fprintf(screen, "test3 from eforward\n");
  } else {
    pair->epack_forward_comm(esendnum[iswap],esendlist[iswap],buf_send);
    pair->eunpack_forward_comm(erecvnum[iswap],efirstrecv[iswap],buf_send);
  }
 }
}

/*----------------------------------------------------------------------
  forward communication of element coords every timestep 
  other per element attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------*/

void CommBrick::eforward_comm(int dummy)
{
   int n;
   MPI_Request request;
   double ***nodex = element->nodex;
   double *buf;

   //if (comm->me==0) fprintf(logfile, "comm_x_only: %d\n", comm_x_only);
   //error->all(FLERR, "test from eforward_comm");
   //now comm_x_only is 1

   for (int iswap = 0; iswap < nswap; iswap++) {
     if (sendproc[iswap] !=me) {
       if (esize_forward_recv[iswap]) 
	 MPI_Irecv(buf_recv, esize_forward_recv[iswap], MPI_DOUBLE,
		   recvproc[iswap],0,world,&request);
         n = element->pack_comm(esendnum[iswap],esendlist[iswap],
			         buf_send, pbc_flag[iswap],pbc[iswap]);
	 if (n) MPI_Send(buf_send,n,MPI_DOUBLE, sendproc[iswap],0,world);
	 if (esize_forward_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
         element->unpack_comm(erecvnum[iswap], efirstrecv[iswap],buf_recv);
     } else {
       element->pack_comm(esendnum[iswap], esendlist[iswap],
		      buf_send,pbc_flag[iswap],pbc[iswap]);
       element->unpack_comm(erecvnum[iswap],efirstrecv[iswap],buf_send);
     }
   }
   //error->all(FLERR, "test from eforward_comm");
}
/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommBrick::forward_comm(int dummy)
{
  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me) {
      if (comm_x_only) {
        if (size_forward_recv[iswap]) {
          if (size_forward_recv[iswap]) buf = x[firstrecv[iswap]];
          else buf = NULL;
          MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        }
        n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
      } else if (ghost_velocity) {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
        avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
      } else {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
        avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    } else {
      if (comm_x_only) {
        if (sendnum[iswap])
          avec->pack_comm(sendnum[iswap],sendlist[iswap],
                          x[firstrecv[iswap]],pbc_flag[iswap],pbc[iswap]);
      } else if (ghost_velocity) {
        avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
        avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_send);
      } else {
        avec->pack_comm(sendnum[iswap],sendlist[iswap],
                        buf_send,pbc_flag[iswap],pbc[iswap]);
        avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_send);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

bigint CommBrick::memory_usage()
{
  bigint bytes = 0;
  bytes += nprocs * sizeof(int);    // grid2proc
  for (int i = 0; i < nswap; i++)
    bytes += memory->usage(sendlist[i],maxsendlist[i]);
  bytes += memory->usage(buf_send,maxsend+bufextra);
  bytes += memory->usage(buf_recv,maxrecv);
  return bytes;
}

/*-----------------------------------------------------------------------------*/

void CommBrick::free_multi()
{
  memory->destroy(multilo);
  memory->destroy(multihi);
}
