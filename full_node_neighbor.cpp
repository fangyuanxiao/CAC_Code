/*-----------------------------------------------------------------------
 used for testing the full node neighbor
 -----------------------------------------------------------------------*/

#include "neighbor.h"
#include "node_neigh_list.h"
#include "element.h"
#include "neigh_list.h"
#include "comm.h"
#include "my_page.h"
#include "math.h"
#include "error.h"
#include "update.h"

using namespace CAC_NS;

#define TOL =1.0e-4 

void Neighbor::node_neighbor(NodeNeighList *nilist, NeighList *list)
{
  // define pointers to the arrays in element class
  int *type = element->type;
  int *tag = element->tag;
  int *numint = element->numint;
  double ***nodex = element->nodex;
  double **x = element->x;
  double ***Spa_int = element->Spa_int;
  int npe = element->npe;
  double ***Spa_inter_ele = element->Spa_inter_ele;
  int **list_int = element->list_int;
  int *inter_num = element->inter_num;

  double cutsq = neighbor->cutneighmaxsq;

  // define pointers to the arrays in neighlist class
  
  int einum = list->einum;
  int *eilist = list->eilist;
  int *enumneigh = list->enumneigh;
  int **efirstneigh = list->efirstneigh;

  int *ejlist;

  int i,j,ii,jj,kk,ll,mm,nn,ninter,itype,jtype,it,jt,id,isub;
  int ejnum,nint; 
  int itag, jtag;
  double ixtmp,iytmp,iztmp;
  double jxtmp,jytmp,jztmp;
  double subxtmp, subytmp, subztmp;
  double delx, dely,delz,rsq;
  int n=0;
  double npdis; // nearest projection distance between the integration point of the element
  double npdis_sub;

  // define pointers to the arrays in node_neighbor_list class

  int *nelist = nilist->elist;
  int *enlist = nilist->enlist;
  int ninum;
  int neinum;
  int *neighptr, *nneighptr;
  int *numneigh = nilist->numneigh;
  int **firstneigh = nilist->firstneigh;
  int **nfirstneigh = nilist->nfirstneigh;
  int nsubele = element->nsubele;
  double **Spa_center_subele = element->Spa_center_subele;
  int *natom_subele = element->natom_subele;
  double ***Spa_inter_subele = element->Spa_inter_subele;
  int **list_inter_subele = element->list_inter_subele;

  MyPage<int> *npage = nilist->npage;
  MyPage<int> *nnpage = nilist->nnpage;

  npage->reset();
  nnpage->reset();

// search in the list page of element neighbor
// i is the tag of the element number

  ninum = 0;
  neinum = 0;

  for (ii=0;ii<einum;ii++) {

     // got the information of the neighbor element of i 
     i = eilist[ii];
     nelist[neinum] = ninum;
     enlist[neinum] = i;

     itype = type[i];
     itag = tag[i];
     it = itype-1;
     nint = numint[itype];
     ejnum = enumneigh[i];
     ejlist = efirstneigh[i];
     
     // search node neighbors of integration points in element i
     for (jj=0; jj<nint;jj++) {

       n = 0;
       neighptr = npage->vget();
       nneighptr = nnpage->vget();
       // to have coordinate of the integration points

       ixtmp = 0;
       iytmp = 0;
       iztmp = 0;

       // interpolation process

       for (kk=0; kk<npe; kk++) {
	ixtmp +=Spa_int[it][jj][kk]*nodex[i][kk][0];
	iytmp +=Spa_int[it][jj][kk]*nodex[i][kk][1];
	iztmp +=Spa_int[it][jj][kk]*nodex[i][kk][2]; 
       }


       //if (comm->me==0) fprintf(logfile, "%d, %d, %f, %f, %f\n",itag, jj, ixtmp, iytmp, iztmp);
       // search around all the neighbor element
       // all atoms in each element is interpolated to jxtmp, jytmp, jztmp 
       // then the interpolated atoms is checked if it is the neighbor of ixtmp or not

       for (ll=0;ll<ejnum;ll++) {
	 j = ejlist[ll];
         jtype = type[j];
	 jt = jtype-1;
	 jtag = tag[j];
	 ninter = inter_num[jt];
	
	 //loop over ninter atoms inside j element
	 npdis = (ixtmp - x[j][0])*(ixtmp - x[j][0])+(iytmp - x[j][1])*(iytmp-x[j][1])+(iztmp-x[j][2])*(iztmp-x[j][2]);
	
	if (npdis <= cut_elesq) {
          for (isub = 0; isub < nsubele; isub++) {
	   // interpolate the center of the subelement
	   subxtmp = 0.0;
	   subytmp = 0.0;
	   subztmp = 0.0;

	   for (kk = 0; kk < npe; kk++) {
              subxtmp += Spa_center_subele[isub][kk]*nodex[j][kk][0];
	      subytmp += Spa_center_subele[isub][kk]*nodex[j][kk][1];
	      subztmp += Spa_center_subele[isub][kk]*nodex[j][kk][2];
	   }

	   npdis_sub = (ixtmp - subxtmp)*(ixtmp - subxtmp)+(iytmp-subytmp)*(iytmp-subytmp)+(iztmp-subztmp)*(iztmp-subztmp);
	  
	   if (npdis_sub <= cut_subelesq) {
            
           for (mm=0; mm<natom_subele[isub]; mm++) {

             jxtmp = 0;
	     jytmp = 0;
	     jztmp = 0;

	     for (kk=0; kk<npe;kk++) {
	       jxtmp +=Spa_inter_subele[isub][mm][kk]*nodex[j][kk][0];
	       jytmp +=Spa_inter_subele[isub][mm][kk]*nodex[j][kk][1];
	       jztmp +=Spa_inter_subele[isub][mm][kk]*nodex[j][kk][2]; 
             }
	     
             delx = jxtmp - ixtmp;
	     dely = jytmp - iytmp;
	     delz = jztmp - iztmp;

	     rsq = delx*delx+dely*dely+delz*delz;
	     if (rsq <= cutsq) {
	       neighptr[n] = j;
	       nneighptr[n] = list_inter_subele[isub][mm];
	       n++;
	     }
	   }
           }
         }
	}
       }
       
       // search around self element

       id = list_int[it][jj];
       ninter = inter_num[it];

      for (isub =0; isub < nsubele; isub++) {
       // initialize the subxtmp
       subxtmp = 0.0;
       subytmp = 0.0;
       subztmp = 0.0;
       
       for (kk = 0; kk < npe; kk++) {
         subxtmp += Spa_center_subele[isub][kk]*nodex[i][kk][0];
	 subytmp += Spa_center_subele[isub][kk]*nodex[i][kk][1];
	 subztmp += Spa_center_subele[isub][kk]*nodex[i][kk][2]; 
       }

       npdis_sub = (ixtmp - subxtmp)*(ixtmp - subxtmp)+(iytmp-subytmp)*(iytmp-subytmp)+(iztmp-subztmp)*(iztmp-subztmp);	 
       if (npdis_sub <= cut_subelesq) {
       for (mm=0;mm<natom_subele[isub];mm++) {

	 jxtmp = 0;
	 jytmp = 0;
	 jztmp = 0;

	 for (kk=0; kk<npe;kk++) {
           jxtmp += Spa_inter_subele[isub][mm][kk]*nodex[i][kk][0];
	   jytmp += Spa_inter_subele[isub][mm][kk]*nodex[i][kk][1];
	   jztmp += Spa_inter_subele[isub][mm][kk]*nodex[i][kk][2];
         }

	 delx = jxtmp-ixtmp;
	 dely = jytmp-iytmp;
	 delz = jztmp-iztmp;

	 rsq = delx*delx+dely*dely+delz*delz;
	// if (comm->me) fprintf(logfile, "rsq: %f\n",rsq);

	 if (rsq <= cutsq && rsq > 1.0e-4) {
	   neighptr[n] = i;
	   nneighptr[n] = list_inter_subele[isub][mm];
	   n++;
         }
       }
       }
      }

       firstneigh[ninum] = neighptr;
       nfirstneigh[ninum] = nneighptr;
       numneigh[ninum] = n;
       npage->vgot(n);
       nnpage->vgot(n);
       ninum++;
       //if (comm->me==0) fprintf(logfile, "n:%d\n",n);
     }
     neinum++;
  }

  nilist->ninum = ninum;
  nilist->neinum = neinum;

}

