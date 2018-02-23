#include "node_neigh_list.h"
#include "neigh_list.h"
#include "element.h"
#include "comm.h"
#include "neighbor.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define PGDELTA 1

NodeNeighList::NodeNeighList(CAC *cac):
   Pointers(cac)
{
   maxnodes = 0;
   maxeles = 0;

   ninum = 0;
   neinum = 0;
   elist = NULL;
   enlist = NULL;
   numneigh = NULL;
   firstneigh = NULL;
   nfirstneigh = NULL;

   dnum = 0;

   npage = NULL;
   nnpage = NULL;
}

NodeNeighList::~NodeNeighList()
{
   //error->all(FLERR, "test from ~NodeNeighList");
   memory->destroy(elist);
   memory->destroy(enlist);
   memory->destroy(numneigh);
   memory->sfree(firstneigh);
   memory->sfree(nfirstneigh);
 
   delete [] npage;
   delete [] nnpage;
}

/*------------------------------------------------------------------*/

void NodeNeighList::grow(int nmaxele, int nmaxinp)
{
   if (nmaxele<=maxeles && nmaxinp<=maxnodes) return;
   maxeles = nmaxele;
   maxnodes = nmaxinp;

  memory->destroy(elist);
  memory->destroy(enlist);
  memory->destroy(numneigh);
  memory->sfree(firstneigh);
  memory->sfree(nfirstneigh);

  memory->create(elist, maxeles, "node_neighlist:elist");
  memory->create(enlist,maxeles,"node_neigh_list:enlist");
  memory->create(numneigh, maxnodes, "node_neighlist:numneigh");
  firstneigh = (int **) memory->smalloc(maxnodes*sizeof(int *),
		        "node_neighlist:firstneigh");
  nfirstneigh= (int **) memory->smalloc(maxnodes*sizeof(int *),
		         "node_neighlist:nfirstneigh");
}

/*-------------------------------------------------------------------*/

void NodeNeighList::setup_pages(int pgsize_caller, int onenode_caller, 
		                int dnum_caller)
{
  pgsize = pgsize_caller;
  onenode = onenode_caller;
  dnum = dnum_caller;

  int nmypage = comm->nthreads;
  npage = new MyPage<int>[nmypage];
  nnpage = new MyPage<int>[nmypage];

  for (int i = 0; i < nmypage; i++) {
    npage[i].init(onenode,pgsize,PGDELTA);
    nnpage[i].init(onenode,pgsize,PGDELTA);
  }
}

/*-------------------------------------------------------------------------------*/

bigint NodeNeighList::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(elist, maxeles);
  bytes += memory->usage(enlist,maxeles);
  bytes += memory->usage(numneigh,maxnodes);
  bytes += maxnodes*sizeof(int *);
  bytes += maxnodes*sizeof(int *);

  int nmypage = comm->nthreads;

  if (npage) {
    for (int i = 0; i< nmypage; i++) {
      bytes += npage[i].size();
      bytes += nnpage[i].size();
    }
  }

  return bytes;
}
