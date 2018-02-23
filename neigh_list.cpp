#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define PGDELTA 1

enum{NSQ,BIN,MULTI};     // also in neighbor.cpp

/*-----------------------------------------------------------------------------------*/

NeighList::NeighList(CAC *cac) :
  Pointers(cac)
{
  maxatoms = 0;
  maxeles = 0;

  inum = einum = gnum = 0;
  ilist = NULL;
  numneigh = NULL;
  firstneigh = NULL;
  firstdouble = NULL;

  eilist = NULL;
  enumneigh = NULL;
  efirstneigh = NULL;

  dnum = 0;

  last_build = -1;

  iskip = NULL;
  ijskip = NULL;

  maxstencil = 0;
  stencil = NULL;
  stencilxyz = NULL;

  maxstencil_multi = 0;
  nstencil_multi = NULL;
  stencil_multi = NULL;
  distsq_multi = NULL;

  ipage = NULL;
  epage = NULL;
  dpage = NULL;
}

NeighList::~NeighList()
{
  if (!listcopy) {
    memory->destroy(ilist);
    memory->destroy(numneigh);
    memory->sfree(firstneigh);
    memory->sfree(firstdouble);

    memory->destroy(eilist);
    memory->destroy(enumneigh);
    memory->sfree(efirstneigh);

    delete [] ipage;
    delete [] epage;
    if (dnum) delete [] dpage;
  }

  delete [] iskip;
  memory->destroy(ijskip);

  if (maxstencil) memory->destroy(stencil);
  if (ghostflag) memory->destroy(stencilxyz);

  if (maxstencil_multi) {
 //   for (int i = 1; i <= atom->ntypes; i++) {
 //     memory->destroy(stencil_multi[i]);
 //     memory->destroy(distsq_multi[i]);
 //   }
    delete [] nstencil_multi;
    delete [] stencil_multi;
    delete [] distsq_multi;
  }
}

/* ----------------------------------------------------------------------
   copy skip info from request rq into list's iskip,ijskip
------------------------------------------------------------------------- */

void NeighList::copy_skip_info(int *rq_iskip, int **rq_ijskip)
{
  int ntypes = atom->ntypes;
  iskip = new int[ntypes+1];
  memory->create(ijskip,ntypes+1,ntypes+1,"neigh_list:ijskip");
  int i,j;
  for (i = 1; i <= ntypes; i++) iskip[i] = rq_iskip[i];
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      ijskip[i][j] = rq_ijskip[i][j];
}

/* ---------------------------------------------------------------------- */

void NeighList::setup_pages(int pgsize_caller, int oneatom_caller, 
                            int dnum_caller)
{
 // fprintf(screen, "me = %d\n",comm->me);
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;
  dnum = dnum_caller;

  int nmypage = comm->nthreads;
  ipage = new MyPage<int>[nmypage];
  for (int i = 0; i < nmypage; i++)
    ipage[i].init(oneatom,pgsize,PGDELTA);

  if (cac->element_flag) {
    epage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      epage[i].init(oneatom,pgsize,PGDELTA);
  }

  if (dnum) {
    dpage = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++)
      dpage[i].init(dnum*oneatom,dnum*pgsize,PGDELTA);
  }
  else dpage = NULL; 
}

/* ----------------------------------------------------------------------
   grow atom arrays to allow for nmax atoms
   triggered by more atoms on a processor
   caller knows if this list stores neighs of local atoms or local+ghost
------------------------------------------------------------------------- */

void NeighList::grow(int nmax)
{
  // skip if this list is already long enough to store nmax atoms

  if (nmax <= maxatoms) return;
  maxatoms = nmax;

  memory->destroy(ilist);
  memory->destroy(numneigh);
  memory->sfree(firstneigh);
  memory->sfree(firstdouble);

  memory->create(ilist,maxatoms,"neighlist:ilist");
  memory->create(numneigh,maxatoms,"neighlist:numneigh");
  firstneigh = (int **) memory->smalloc(maxatoms*sizeof(int *),
                                        "neighlist:firstneigh");

  if (dnum)
    firstdouble = (double **) memory->smalloc(maxatoms*sizeof(double *),
                                              "neighlist:firstdouble");
}

/*-----------------------------------------------------------------------
  grow elements arrays to allow for nmax elements
  triggered by more elements on a processor
  caller knows if this list stores neighs of local elements or local+ghost
-------------------------------------------------------------------------*/

void NeighList::egrow(int nmax)
{
   // skip if this list is already long enough to store nmax elements
   if (nmax <= maxeles) return;
   maxeles = nmax;

   memory->destroy(eilist);
   memory->destroy(enumneigh);
   memory->sfree(efirstneigh);

   memory->create(eilist,maxeles,"neighlist:eilist");
   memory->create(enumneigh,maxeles,"neighlist:enumneigh");
   efirstneigh = (int **) memory->smalloc(maxeles*sizeof(int *),
		                          "neighlist:firstneigh");

}

/* ----------------------------------------------------------------------
   insure stencils are large enough for smax bins
   style = BIN or MULTI
------------------------------------------------------------------------- */

void NeighList::stencil_allocate(int smax, int style)
{

  if (style == BIN) {
    if (smax > maxstencil) {
      maxstencil = smax;
      memory->destroy(stencil);
      memory->create(stencil,maxstencil,"neighlist:stencil");
    }

  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   if growflag = 0, maxatoms & maxpage will also be 0
   if stencilflag = 0, maxstencil * maxstencil_multi will also be 0
------------------------------------------------------------------------- */

bigint NeighList::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(ilist,maxatoms);
  bytes += memory->usage(numneigh,maxatoms);
  bytes += memory->usage(eilist, maxeles);
  bytes += memory->usage(enumneigh,maxeles);
  bytes += maxatoms * sizeof(int *);
  bytes += maxeles*sizeof(int *);

  int nmypage = comm->nthreads;

  if (ipage) {
    for (int i = 0; i < nmypage; i++) {
      bytes += ipage[i].size();
      bytes += epage[i].size();
    }
  }

  if (dnum && dpage) {
    for (int i = 0; i < nmypage; i++) {
      bytes += maxatoms * sizeof(double *);
      bytes += dpage[i].size();
    }
  }

  if (maxstencil) bytes += memory->usage(stencil,maxstencil);

  if (maxstencil_multi) {
    bytes += memory->usage(stencil_multi,atom->ntypes,maxstencil_multi);
    bytes += memory->usage(distsq_multi,atom->ntypes,maxstencil_multi);
  }

  return bytes;
}
