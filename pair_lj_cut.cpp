#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "universe.h"
#include "memory.h"
#include "error.h"
#include "node_neigh_list.h"
#include "element.h"
#include "update.h"
#include "neighbor.h"

using namespace CAC_NS;
using namespace MathConst;

PairLJCut::PairLJCut(CAC *cac) : Pair(cac)
{
  respa_enable = 1;
  writedata = 1;
}

PairLJCut::~PairLJCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
}

/* ----------------------------------------------------------------------
  allocate all arrays
------------------------------------------------------------------------- */

void PairLJCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = universe->numeric(FLERR,arg[0]);

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCut::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();
  
  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = universe->numeric(FLERR,arg[2]);
  double sigma_one = universe->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = universe->numeric(FLERR,arg[4]);

  int count = 0;

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
  proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;

    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCut::init_style()
{
   int irequest;

   irequest = neighbor->request(this,instance_me);
   cut_respa = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJCut::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ---------------------------------------------------------------------- */

void PairLJCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  //if (comm->me==0) fprintf(screen, "force->special_lj = %f, %f, %f\n",
  //                    special_lj[1], special_lj[2],special_lj[3]);
  //if (comm->me==0) fprintf(screen, "newton_pair = %d\n", newton_pair); 

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        
      } 
    }
  }
}

/*----------------------------------------------------------------------
 Compute the force on each integration point and then reallocate the force
 to the ndoe through the shape function array of each integration point
------------------------------------------------------------------------*/

void PairLJCut::ecompute()
{
  // define pointers to the nilist
  int *nelist = nilist->elist;
  int *enlist = nilist->enlist;

  int ninum = nilist->ninum;
  int neinum = nilist->neinum;
  int *numneigh = nilist->numneigh;

  int **firstneigh = nilist->firstneigh;
  int **nfirstneigh = nilist->nfirstneigh;

  int *jlist, *njlist;
  int jnum;

  int *type = element->type;
  int *numint = element->numint;
  double ***nodex = element->nodex;
  double ***nodef = element->nodef;
  double ***Spa_int = element->Spa_int;
  int npe = element->npe;
  double ***Spa_inter_ele = element->Spa_inter_ele;
  double **weight = element->weight; 
  int nall = element->nlocal + element->nghost; 
  int nlocal = element->nlocal;

  int i,j,nj,ii,jj,kk,ll,mm,nn,itype,jtype,it,jt,nint,np;
  double ixtmp,iytmp,iztmp;
  double jxtmp,jytmp,jztmp;
  double delx, dely, delz, rsq,fpair;

  double ecutsq = cutsq[1][1];
  
  double fx,fy,fz;
  double elj1 = lj1[1][1];
  double elj2 = lj2[1][1];
  double r2inv,r6inv,forcelj;
  bigint ntimestep = update->ntimestep;
  
 // if (comm->me==14 && ntimestep == 2375) {
 //   for (ii=0; ii<neinum;ii++) {
 //     fprintf(screen, "%d,%d\n",ii,enlist[ii]);
 //   }
 //   error->all(FLERR, "test from pair");
 // }



  for (ii=0;ii<neinum;ii++) {

  // i = enlist[ii];
   np = nelist[ii];
   itype = type[ii];
   it = itype-1;
   nint = numint[itype];
  // if (i >= nlocal) error->all(FLERR,"test1");
   for (jj=0;jj<nint;jj++) {
      
     // interpolate jj integration point in i element
     ixtmp = 0;
     iytmp = 0;
     iztmp = 0;

     fx = 0;
     fy = 0;
     fz = 0;

     for (kk=0; kk<npe;kk++) {
      ixtmp +=Spa_int[it][jj][kk]*nodex[ii][kk][0];
      iytmp +=Spa_int[it][jj][kk]*nodex[ii][kk][1];
      iztmp +=Spa_int[it][jj][kk]*nodex[ii][kk][2];
     }

     //if(comm->me==0) fprintf(logfile, "i:%f,%f,%F\n",ixtmp,iytmp,iztmp);
     // get the neighbor information

     jnum = numneigh[np];
     jlist = firstneigh[np];
     njlist = nfirstneigh[np];
     
     // search for all neighbors
     for (ll=0;ll<jnum;ll++) {
      // nj atom inside j element 
       
      j = jlist[ll];
      nj = njlist[ll];
   //   if (j >= nall) error->all(FLERR,"test2");
   //   if (nj >= 2197) error->all(FLERR,"test3");
      jtype = type[j];
      jt = jtype - 1;

      // interpolation of nj atom inside j element
      jxtmp = 0;
      jytmp = 0;
      jztmp = 0;

      for (kk=0;kk<npe;kk++) {
        jxtmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][0];
	jytmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][1];
	jztmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][2];
      }

      delx = ixtmp - jxtmp;
      dely = iytmp - jytmp;
      delz = iztmp - jztmp;
      rsq = delx*delx+dely*dely+delz*delz;

     // if (comm->me==0) fprintf(logfile, "%f,%f\n", rsq, ecutsq);

      if (rsq < ecutsq) {
       
        r2inv = 1.0/rsq;
	r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv*(elj1*r6inv-elj2);

	fpair = forcelj*r2inv;
        //if (comm->me==0) fprintf(logfile,"%f\n",fpair);
	fx += delx*fpair;
	fy += dely*fpair;
	fz += delz*fpair;

      }
      //if(comm->me==0) fprintf(logfile,"%f,%f,%f\n",jxtmp,jytmp,jztmp);
     }
      // error->all(FLERR, "test from pair_lj_cut");
     // reallocate the force calculated to the nodes of the element through shape function
     //if (comm->me==0) fprintf(logfile,"%f,%f,%f\n",fx,fy,fz);

     for (kk=0;kk<npe;kk++) {
        nodef[ii][kk][0] += weight[it][jj]*Spa_int[it][jj][kk]*fx;
	nodef[ii][kk][1] += weight[it][jj]*Spa_int[it][jj][kk]*fy;
	nodef[ii][kk][2] += weight[it][jj]*Spa_int[it][jj][kk]*fz;
     }

     np++;

   }
  }

// int nlocal = element->nlocal;
//if (comm->me==0) {
//  for (ii = 0;ii<nlocal;ii++) 
//    for(jj=0;jj<npe;jj++)
//      fprintf(logfile, "%f,%f,%f\n", nodef[ii][jj][0], nodef[ii][jj][1],nodef[ii][jj][2]);
//}
//  error->all(FLERR, "test from ecompute");
}
