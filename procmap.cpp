#include "procmap.h"
#include "universe.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

#include <map>
#include <string>

using namespace CAC_NS;

#define MAXLINE 128

enum{MULTIPLE};

ProcMap::ProcMap(CAC *cac) : Pointers(cac) {}

/* ----------------------------------------------------------------------
   create a one-level 3d grid of procs
------------------------------------------------------------------------- */

void ProcMap::onelevel_grid(int nprocs, int *user_procgrid, int *procgrid,
                            int otherflag, int other_style,
                            int *other_procgrid, int *other_coregrid)
{
  int **factors;
  int npossible = factor(nprocs,NULL);
  memory->create(factors,npossible,3,"procmap:factors");
  npossible = factor(nprocs,factors);

  // constrain by 2d, user request, other partition

  if (domain->dimension == 2) npossible = cull_2d(npossible,factors,3);
  npossible = cull_user(npossible,factors,3,user_procgrid);
  if (otherflag) npossible = cull_other(npossible,factors,3,
                                        other_style,other_procgrid,
                                        other_coregrid);

  if (npossible == 0)
    error->all(FLERR,"Could not create 3d grid of processors");

  best_factors(npossible,factors,procgrid,1,1,1);

  memory->destroy(factors);
}

/* ----------------------------------------------------------------------
  generate all possible 3-integer factorizations of N
  store them in factors if non-NULL
  return # of factorizations
 ------------------------------------------------------------------------- */

int ProcMap::factor(int n, int **factors)
{
  int i,j,nyz;

  int m = 0;
  for (i = 1; i <= n; i++) {
    if (n % i) continue;
    nyz = n/i;
    for (j = 1; j <= nyz; j++) {
      if (nyz % j) continue;
      if (factors) {
        factors[m][0] = i;
        factors[m][1] = j;
        factors[m][2] = nyz/j;
      }
      m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   remove any factors where Pz != 1 for 2d
  ------------------------------------------------------------------------- */

int ProcMap::cull_2d(int n, int **factors, int m)
{
  int i = 0;
  while (i < n) {
    if (factors[i][2] != 1) {
      for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
      n--;
    } else i++;
  }
  return n;
}

/* ----------------------------------------------------------------------
  remove any factors that do not match non-zero user_factors Px,Py,Pz
 ------------------------------------------------------------------------- */

int ProcMap::cull_user(int n, int **factors, int m, int *user_factors)
{
  int i = 0;
  while (i < n) {
    int flag = 0;
    if (user_factors[0] && factors[i][0] != user_factors[0]) flag = 1;
    if (user_factors[1] && factors[i][1] != user_factors[1]) flag = 1;
    if (user_factors[2] && factors[i][2] != user_factors[2]) flag = 1;
    if (flag) {
      for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
      n--;
    } else i++;
  }
  return n;
}

/* ----------------------------------------------------------------------
   remove any factors that do not match settings from other partition
   MULTIPLE = other Nx,Ny,Nz must be multiple of my Px,Py,Pz
              where Nx,Ny,Nz = node grid = procgrid/coregrid
------------------------------------------------------------------------- */

int ProcMap::cull_other(int n, int **factors, int m,
                        int other_style, int *other_procgrid,
                        int *other_coregrid)
{
  int i = 0;
  while (i < n) {
    if (other_style == MULTIPLE) {
      int flag = 0;
      if ((other_procgrid[0]/other_coregrid[0]) % factors[i][0]) flag = 1;
      if ((other_procgrid[1]/other_coregrid[1]) % factors[i][1]) flag = 1;
      if ((other_procgrid[2]/other_coregrid[2]) % factors[i][2]) flag = 1;
      if (flag) {
        for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
        n--;
      } else i++;
    }
  }
  return n;
}

/* ----------------------------------------------------------------------
   choose best factors from list of Npossible factors
   best = minimal surface area of sub-domain
   return best = 3 factors
   return index of best factors in factors
------------------------------------------------------------------------- */

int ProcMap::best_factors(int npossible, int **factors, int *best,
                          const int sx, const int sy, const int sz)
{
  // determine cross-sectional areas for orthogonal and triclinic boxes
  // for triclinic, area = cross product of 2 edge vectors stored in h matrix
  // area[3] = surface area 3 box faces divided by sx,sy,sz
  // area[0] = xy, area[1] = xz, area[2] = yz

  double area[3];
  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd / (sx*sy);
    area[1] = domain->xprd * domain->zprd / (sx*sz);
    area[2] = domain->yprd * domain->zprd / (sy*sz);
  } else {
    double *h = domain->h;
    double a[3],b[3],c[3];
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[5]; b[1] = h[1]; b[2] = 0.0;
    MathExtra::cross3(a,b,c);
    area[0] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx*sy);
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[1] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx*sz);
    a[0] = h[5]; a[1] = h[1]; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[2] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sy*sz);
  }

  int index;
  double surf;
  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  for (int m = 0; m < npossible; m++) {
    surf = area[0]/factors[m][0]/factors[m][1] +
      area[1]/factors[m][0]/factors[m][2] +
      area[2]/factors[m][1]/factors[m][2];
    if (surf < bestsurf) {
      bestsurf = surf;
      best[0] = factors[m][0];
      best[1] = factors[m][1];
      best[2] = factors[m][2];
      index = m;
    }
  }

  return index;

}

/* ----------------------------------------------------------------------
   map processors to 3d grid via MPI_Cart routines
   MPI may do layout in machine-optimized fashion
------------------------------------------------------------------------- */

void ProcMap::cart_map(int reorder, int *procgrid,
                       int *myloc, int procneigh[3][2], int ***grid2proc)
{
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;

  MPI_Cart_create(world,3,procgrid,periods,reorder,&cartesian);
  MPI_Cart_get(cartesian,3,procgrid,periods,myloc);
  MPI_Cart_shift(cartesian,0,1,&procneigh[0][0],&procneigh[0][1]);
  MPI_Cart_shift(cartesian,1,1,&procneigh[1][0],&procneigh[1][1]);
  MPI_Cart_shift(cartesian,2,1,&procneigh[2][0],&procneigh[2][1]);

  int coords[3];
  int i,j,k;
  for (i = 0; i < procgrid[0]; i++)
    for (j = 0; j < procgrid[1]; j++)
      for (k = 0; k < procgrid[2]; k++) {
        coords[0] = i; coords[1] = j; coords[2] = k;
        MPI_Cart_rank(cartesian,coords,&grid2proc[i][j][k]);
      }

  MPI_Comm_free(&cartesian);

}

/* ----------------------------------------------------------------------
  output mapping of processors to 3d grid to file
------------------------------------------------------------------------- */

void ProcMap::output(char *file, int *procgrid, int ***grid2proc)
{
  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  FILE *fp;
  if (me == 0) {
    fp = fopen(file,"w");
    if (fp == NULL) error->one(FLERR,"Cannot open processors output file");
    fprintf(fp,"LAMMPS mapping of processors to 3d grid\n");
    fprintf(fp,"partition = %d\n",universe->iworld+1);
    fprintf(fp,"Px Py Pz = %d %d %d\n",procgrid[0],procgrid[1],procgrid[2]);
    fprintf(fp,"world-ID universe-ID original-ID: I J K: name\n\n");
  }

  // find me in the grid

  int ime,jme,kme;
  for (int i = 0; i < procgrid[0]; i++)
    for (int j = 0; j < procgrid[1]; j++)
      for (int k = 0; k < procgrid[2]; k++)
        if (grid2proc[i][j][k] == me) {
          ime = i; jme = j; kme = k;
        }

  // polled comm of grid mapping info from each proc to proc 0

  int tmp;
  int vec[6];
  char procname[MPI_MAX_PROCESSOR_NAME+1];

  vec[0] = me;
  vec[1] = universe->me;
  MPI_Comm_rank(universe->uorig,&vec[2]);
  vec[3] = ime + 1;
  vec[4] = jme + 1;
  vec[5] = kme + 1;

  int len;
  MPI_Get_processor_name(procname,&len);
  procname[len] = '\0';

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Recv(vec,6,MPI_INT,iproc,0,world,MPI_STATUS_IGNORE);
        MPI_Recv(procname,MPI_MAX_PROCESSOR_NAME+1,MPI_CHAR,
                 iproc,0,world,MPI_STATUS_IGNORE);
      }
      fprintf(fp,"%d %d %d: %d %d %d: %s\n",
              vec[0],vec[1],vec[2],vec[3],vec[4],vec[5],procname);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Send(vec,6,MPI_INT,0,0,world);
    MPI_Send(procname,strlen(procname)+1,MPI_CHAR,0,0,world);
  }

  if (me == 0) fclose(fp);
  //error->all(FLERR,"test from procmap");
}
