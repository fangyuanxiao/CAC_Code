#include "neighbor.h"
#include "neigh_list.h"
#include "element.h"
#include "domain.h"
#include "my_page.h"
#include "comm_brick.h"
#include "error.h"

using namespace CAC_NS;

void Neighbor::element_neighbor(NeighList *list)
{
   int i,j,k,n,itype,jtype,ibin,which,imol,iele,moltemplate;
   tagint tagprev;
   double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
   int *neighptr;

   bin_elements();
  

   double **x = element->x;
   int *type = element->type;
   int *mask = element->mask;
   tagint *tag = element->tag;

   int nlocal = element->nlocal;
   int *eilist = list->eilist;
   int *enumneigh = list->enumneigh;
   int **efirstneigh = list->efirstneigh;
   int nstencil  = list->nstencil;
   int *stencil = list->stencil;
   MyPage<int> *epage = list->epage;
   int inum;
   epage->reset();
   //loop over each element, storing neighbors

   inum = 0;

   for (i=0; i<nlocal;i++) {
      n=0;
      neighptr = epage->vget();

      itype = type[i];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // loop over all elements in other bins in stencil including self
      // only store pair if 
      ibin = coord2bin(x[i]);
      
      for (k=0; k<nstencil;k++) {
       for (j = ebinhead[ibin+stencil[k]];j>=0;j=ebins[j]) {
	   if (j==i) continue;
	   
	   jtype = type[j];

	   delx = xtmp - x[j][0];
	   dely = ytmp - x[j][1];
	   delz = ztmp - x[j][2];
           rsq = delx*delx+dely*dely+delz*delz;

	   if (rsq <= cut_allsq) neighptr[n++] = j;
	}
      }

      eilist[inum++] = i;
      efirstneigh[i] = neighptr;
      enumneigh[i] = n;
      epage->vgot(n);
       if (epage->status())
	error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
   }
   list->einum = inum;
}
