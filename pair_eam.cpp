#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_eam.h"
#include "atom.h"
#include "force.h"
#include "element.h"
#include "comm.h"
#include "neighbor.h"
#include "universe.h"
#include "neigh_list.h"
#include "node_neigh_list.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace CAC_NS;

#define MAXLINE 1024

/*--------------------------------------------------------------------------------------------*/

PairEAM::PairEAM(CAC *cac) : Pair(cac)
{
  restartinfo = 0;
  manybody_flag = 1;

  enmax = 0;
  nmax = 0;
  rho = NULL;
  fp = NULL;
  erho = NULL;

  nfuncfl = 0;
  funcfl = NULL;

  setfl = NULL;
  fs = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;

  // set comm size needed by this pair

  comm_forward = 1;
  comm_reverse = 1;
  ecomm_forward = element->npe;
}

/*----------------------------------------------------------------------------------------
  check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------------------------*/

PairEAM::~PairEAM()
{
  if (copymode) return;
 
  memory->destroy(rho);
  memory->destroy(fp);
  memory->destroy(erho);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] type2frho;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl;i++) {
      delete [] funcfl[i].file;
      memory->destroy(funcfl[i].frho);
      memory->destroy(funcfl[i].rhor);
      memory->destroy(funcfl[i].zr);
    }
    memory->sfree(funcfl);
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}


/*----------------------------------------------------------------------------------------------*/

void PairEAM::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR, "Illegal pair_style command");
}


/*------------------------------------------------------------------------------------------------
  set coeff for one or more type pairs
  read DYNAMO funcfl file
-------------------------------------------------------------------------------------------------*/

void PairEAM::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);
  // for now the default results for ilo ihi jlo jhi are 1 1 1 1

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2],funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *)
      memory->srealloc(funcfl,nfuncfl*sizeof(Funcfl),"pair:funcfl");
    read_file(arg[2]);
    int n = strlen(arg[2]) + 1;
    funcfl[ifuncfl].file = new char[n];
    strcpy(funcfl[ifuncfl].file,arg[2]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i=j

  int count = 0;
  for (int i = ilo; i <= ihi;i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (i==j) {
        setflag[i][i] = 1;
	map[i] = ifuncfl;
	atom->set_mass(i,funcfl[ifuncfl].mass);
	element->set_mass(funcfl[ifuncfl].mass);
	count++;
      }
    }
  }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
  //error->all(FLERR, "test from pair_eam");
}

/*------------------------------------------------------------------------------
  allocate all arrays
-------------------------------------------------------------------------------*/

void PairEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <=n; i++)
    for (int j = i; j<=n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i<=n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
}

/*---------------------------------------------------------------------------------
  read potential values from a DYNAMO single element funcfl file
-----------------------------------------------------------------------------------*/

void PairEAM::read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl-1];

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  int tmp;
  if (me==0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg", &tmp, &file->mass);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %d %lg %lg",
	   &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->mass,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  memory->create(file->frho,(file->nrho+1),"pair:frho");
  memory->create(file->rhor,(file->nr+1),"pair:rhor");
  memory->create(file->zr,(file->nr+1),"pair:zr");

  if (me == 0) grab(fptr,file->nrho,&file->frho[1]);
  MPI_Bcast(&file->frho[1],file->nrho,MPI_DOUBLE,0,world);

 // if (me==0) grab(fptr,file->nr,&file->zr[1]);
 // MPI_Bcast(&file->zr[1],file->nr,MPI_DOUBLE,0,world);

  if(me==0) grab(fptr,file->nr,&file->rhor[1]);
  MPI_Bcast(&file->rhor[1],file->nr,MPI_DOUBLE,0,world);

  if(me==0) grab(fptr,file->nr,&file->zr[1]);
  MPI_Bcast(&file->zr[1],file->nr,MPI_DOUBLE,0,world);

  if (me== 0) fclose(fptr); 
}

/*------------------------------------------------------------------------------
  grab n values from file fp and put them in list
  values can be several to a line
  only called by proc 0
--------------------------------------------------------------------------------*/

void PairEAM::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fptr);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/*----------------------------------------------------------------------------------
  init specific to this pair style
------------------------------------------------------------------------------------*/

void PairEAM::init_style()
{
  // convert read-in file(s) to arrays and spline them
  file2array();
  array2spline();
  neighbor->request(this,instance_me);
  //error->all(FLERR, "test from init_style");
  
}

/*----------------------------------------------------------------------------------
  init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------------------*/

double PairEAM::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax,funcfl[m].cut);
  } else if (setfl) cutmax = setfl->cut;
  else if (fs) cutmax = fs->cut;

  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/*-----------------------------------------------------------------------------------
  convert read_in funcfl potential(s) to standard array format
  interpolate all file values to a single grid and cutoff
-------------------------------------------------------------------------------------*/

void PairEAM::file2array()
{
  int i,j,k,m,n;
  int ntypes = atom->ntypes;
  double sixth = 1.0/6.0;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr,file->dr);
    drho = MAX(drho,file->drho);
    rmax = MAX(rmax,(file->nr-1)*file->dr);
    rhomax = MAX(rhomax,(file->nrho-1)*file->drho);
  }

  // set nr, nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);

  //if (comm->me==0) fprintf(screen, "nr = %d, nrho = %d\n", nr,nrho);

  //--------------------------------------------------------------------
  // setup frho arrays
  //--------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r, p, cof1,cof2,cof3,cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m-1)*drho;
      p = r/file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      frho[n][m] = cof1*file->frho[k-1] + cof2*file->frho[k] +
	           cof3*file->frho[k+1] + cof4*file->frho[k+2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  //--------------------------------------------------------------------------
  // setup rhor arrays
  //--------------------------------------------------------------------------
 
  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m-1)*dr;
      p = r/file->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      rhor[n][m] = cof1*file->rhor[k-1] + cof2*file->rhor[k] +
	      cof3*file->rhor[k+1] + cof4*file->rhor[k+2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used
  
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  //--------------------------------------------------------------------
  // setup z2r arrays
  //--------------------------------------------------------------------
 
  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl*(nfuncfl+1)/2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r,nr+1,"pair:z2r");

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff

  double zri,zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
	r = (m-1)*dr;

//	p = r/ifile->dr + 1.0;
//	k = static_cast<int> (p);
//	k = MIN(k,ifile->nr-2);
//	k = MAX(k,2);
//	p -= k;
//	p = MIN(p,2.0);
//	cof1 = -sixth*p*(p-1.0)*(p-2.0);
//	cof2 = 0.5*(p*p-1.0)*(p-2.0);
//	cof3 = -0.5*p*(p+1.0)*(p-2.0);
//	cof4 = sixth*p*(p*p-1.0);
//	zri = cof1*ifile->zr[k-1] + cof2*ifile->zr[k] +
//            cof3*ifile->zr[k+1] + cof4*ifile->zr[k+2];


	p = r/jfile->dr + 1.0;
	k = static_cast<int> (p);
        k = MIN(k,jfile->nr-2);
	k = MAX(k,2);
	p -= k;
	p = MIN(p,2.0);
	cof1 = -sixth*p*(p-1.0)*(p-2.0);
	cof2 = 0.5*(p*p-1.0)*(p-2.0);
	cof3 = -0.5*p*(p+1.0)*(p-2.0);
        cof4 = sixth*p*(p*p-1.0);
	zrj = cof1*jfile->zr[k-1] + cof2*jfile->zr[k] +
	      cof3*jfile->zr[k+1] + cof4*jfile->zr[k+2];
	//z2r[n][m] = 27.2*0.529*zri*zrj;
	z2r[n][m] = zrj;

      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid)
  //  type2z2r is not used by non-opt
  //  but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol ==-1) {
	type2z2r[i][j] = 0;
	continue;
      }
      if (irow < icol) {
	irow = map[j];
	icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n+=m+1;
      n+=icol;
      type2z2r[i][j] = n;
    }
  }
  //if (comm->me==0) fprintf(screen,"type2z2r[1][1] = %d\n",type2z2r[1][1]);
  //error->all(FLERR, "test from file2array");
}

/*----------------------------------------------------------------------------------------*/

void PairEAM::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline,nfrho,nrho+1,7,"pair:frho");
  memory->create(rhor_spline,nrhor,nr+1,7,"pair:rhor");
  memory->create(z2r_spline,nz2r,nr+1,7,"pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr, dr,rhor[i],rhor_spline[i]);

  for (int i=0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/*-----------------------------------------------------------------------------------------*/

void PairEAM::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5*(spline[3][6] - spline[1][6]);
  spline[n-1][5] = 0.5*(spline[n][6] - spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6] - spline[m+2][6]) +
		    8.0*(spline[m+1][6]-spline[m-1][6]))/12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6] - spline[m][6])-
	    2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] - 
	    2.0*(spline[m+1][6] - spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m<=n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/*-----------------------------------------------------------------------------------------*/

void PairEAM::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
  tagint *tag = atom->tag;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
      for (i = 0; i < nall; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;
  
  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for ( jj = 0; jj < jnum;jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx +dely*dely+delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	p = sqrt(rsq)*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rho[i] += ((coeff[3]*p+coeff[4])*p+coeff[5])*p+coeff[6];
      }
    }
  }
  
  //if (comm->me==0) {
  //  for (ii = 0; ii < nlocal; ii++) fprintf(logfile, "%d,%f\n", tag[ii],rho[ii]);
  //}

  //error->all(FLERR, "test from compute in pair eam");
  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //  will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p+coeff[1])*p + coeff[2];

  }
 
  // communicate derivative of embedding function 
  comm->forward_comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j&=NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	r = sqrt(rsq);
        p = r*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
        
	// rhoip = derivative of (density at atom j due to atom i)
	// rhojp = derivative of (density at atom j due to atom j)
	// phi = pair potential energy
	// phip = phi'
	// z2 = phi*r
	// z2p = (phi * r)' = (phi' r) + phi
	// psip needs both fp[i] and fp[j] terms since r_ij appears in two 
	//   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	//   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
	
	coeff = rhor_spline[type2rhor[itype][jtype]][m];
	rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rhojp = (coeff[0]*p + coeff[1])*p+coeff[2];
	coeff = z2r_spline[type2z2r[itype][jtype]][m];
	z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
	z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

	recip = 1.0/r;
	phi = z2*recip;
	phip = z2p*recip - phi*recip;
	psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	//psip = phip;
	fpair = -psip*recip;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;
      }
    }
  }
}

/*----------------------------------------------------------------------------------------
  compute the force on the node force on each element:
  the calculation process is as follows:
  1. calculate the electron density on the nodes
  2. interpolate the electron at each integration point
  3. calculate the derivative of fp at the integration point
  4. interpolate the electron density at the neighbor of the integration point
  5. calculate the derivative of fp at the neighbor of the integartion point
  6. calculate the derivative of the electron density at the neighbor of the integration point
  7. calculate the pair interaction
  8. Got the force
------------------------------------------------------------------------------------------*/

void PairEAM::ecompute()
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

  // define pointers related to element class
  int *type = element->type;
  int *numint = element->numint;
  double ***nodex = element->nodex;
  double ***nodef = element->nodef;
  double ***Spa_int = element->Spa_int;
  int npe = element->npe;
  double ***Spa_inter_ele = element->Spa_inter_ele;
  double **weight = element->weight;
  int **node_list = element->node_list;
  int nlocal = element->nlocal;
  tagint *tag = element->tag;

  // variable used in calculate the force

  int i,j,nj,ii,jj,kk,ll,mm,nn,itype,jtype,it,jt,nint,np,nt,nl;
  double ixtmp,iytmp,iztmp;
  double jxtmp,jytmp,jztmp;
  double delx,dely,delz,rsq,fpair,r,recip,evdwl;

  double fx,fy,fz;

  double p;
  int m;
  double *coeff;

  // rhoi: the electron density at integration point
  // rhoj: the electron density at the neighbor site of the integration point
  double rhoi, rhoj;

  //fpi: derivative of embedding energy at integration point
  //fpj: derivative of embedding energy at the neighbor site of the integration point
  double fpi, fpj;

  // rhoip: derivative of electron density with distance
  // since there is only one atom type in our simulation 
  // rhoip = rhojp
  double rhoip;

  // variable related with pair potential in eam potential

  double phi, phip, z2, z2p, psip,rhotmp;
  
  evdwl = 0.0;
//  ev_setup();

  // allocate arrays to store the electron density at nodes
  // the size need to be element->nmax*npe 
  bigint ntimestep = update->ntimestep;
MPI_Barrier(world);
if (comm->me==0 && ntimestep > 334) fprintf(screen, "test1 from pair eam\n");
  if (element->nmax > enmax) {
    memory->destroy(erho);
    enmax = element->nmax;
    memory->create(erho,enmax,npe,"pair:erho");
    //error->all(FLERR, "test from pair eam");
  }

  for (i = 0; i < enmax; i++) {
    for (j = 0; j < npe; j++) {
      erho[i][j] = 0.0;
    }
  }
  //if (ntimestep >333) fprintf(screen, "%d, %d,%d,%d\n", enmax, element->nlocal,element->nghost,element->nmax);
MPI_Barrier(world);
if (comm->me==0 && ntimestep > 334) fprintf(screen, "test2 from pair eam,%f\n",cutforcesq);
//if (comm->me==0 && ntimestep > 334) fprintf(screen, "%f, %d\n",rdr,nr); 
  // calculate the electron density at element node


if (ntimestep > 334) {
  //fprintf(screen, "neinum = %d, nlocal = %d\n",neinum, nlocal);
  for (ii=0; ii < neinum; ii++) {
    i = enlist[ii]; 
    np = nelist[ii];
    itype = type[i];
    it = itype -1;
    nint = numint[itype];
   // if (nint < npe) error->all(FLERR,"nint < npe");
   // if (i <0 || i >= enmax) error->all(FLERR,"test from pair eam");

    for(jj=0; jj<nint; jj++) {

     // ixtmp = 0.0;
     // iytmp = 0.0;
     // iztmp = 0.0;
     // for (kk = 0; kk < npe; kk++) {
     //	ixtmp += Spa_int[it][jj][kk]*nodex[i][kk][0];
     //	iytmp += Spa_int[it][jj][kk]*nodex[i][kk][1];
     //	iztmp += Spa_int[it][jj][kk]*nodex[i][kk][2];
     // }

      ixtmp = nodex[i][jj][0];
      iytmp = nodex[i][jj][1];
      iztmp = nodex[i][jj][2];

      jnum = numneigh[np];
      jlist = firstneigh[np];
      njlist = nfirstneigh[np];
      fprintf(screen, "tag = %d,%d jnum = %d\n",tag[ii],jj,jnum);
      for (ll = 0; ll < jnum; ll++) {
        j = jlist[ll];
	nj = njlist[ll];
	jtype = type[j];
	jt = jtype-1;

	// interpolation of nj atom inside j element
	jxtmp = 0;
	jytmp = 0;
	jztmp = 0;

	for (kk = 0; kk < npe; kk++) {
          jxtmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][0];
	  jytmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][1];
	  jztmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][2];
        }

	delx = ixtmp - jxtmp;
	dely = iytmp - jytmp;
	delz = iztmp - jztmp;
        rsq = delx*delx+dely*dely+delz*delz;
        //fprintf(logfile, "rsq = %f\n",rsq);
	//if (rsq < 0.0001) error->all(FLERR, "test from pair eam");
       //if (rsq < cutforcesq) {
	 // p = sqrt(rsq)*rdr + 1.0;
	//  m = static_cast<int> (p);
	//  m = MIN(m,nr-1);
	//  p -= m;
	// p = MIN(p,1.0);
	//  coeff = rhor_spline[type2rhor[1][1]][m];
	//  rhotmp = ((coeff[3]*p+coeff[4])*p + coeff[5])*p+coeff[6];
        //}


      }
      np++;
    }
    
  }
}

//if (ntimestep > 334) fprintf(screen, "ii = %d, neinum = %d, me = %d\n", ii, neinum, comm->me);


MPI_Barrier(world);
if (comm->me==0 && ntimestep > 334) {
  fprintf(screen, "test d from pair eam\n");
  error->all(FLERR, "test from pair eam");
}

  for (ii=0; ii <neinum;ii++) {

    i = enlist[ii];
    np = nelist[ii];
    itype = type[i];
    it = itype - 1;
    nint = numint[itype];

    for (jj=0; jj <npe; jj++) {
      
      nl = node_list[itype][jj];
      ixtmp = nodex[i][jj][0];
      iytmp = nodex[i][jj][1];
      iztmp = nodex[i][jj][2];
      
      nt = np + nl;

      jnum = numneigh[nt];
      jlist = firstneigh[nt];
      njlist = nfirstneigh[nt];

      // search for all neighbors 

      for (ll = 0; ll < jnum; ll++) {
	// nj atom inside j element
	
        j = jlist[ll];
	nj = njlist[ll];
	jtype = type[j];
	jt = jtype-1;

	// interpolation of nj atom inside j element
	jxtmp = 0;
	jytmp = 0;
	jztmp = 0;

	for (kk = 0; kk < npe; kk++) {
          jxtmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][0];
	  jytmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][1];
	  jztmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][2];
        }

	delx = ixtmp - jxtmp;
	dely = iytmp - jytmp;
	delz = iztmp - jztmp;
        rsq = delx*delx+dely*dely+delz*delz;

	if (rsq < cutforcesq) {
	  p = sqrt(rsq)*rdr + 1.0;
	  m = static_cast<int> (p);
	  m = MIN(m,nr-1);
	  p -= m;
	  p = MIN(p,1.0);
	  coeff = rhor_spline[type2rhor[1][1]][m];
	  erho[i][jj] += ((coeff[3]*p+coeff[4])*p + coeff[5])*p+coeff[6];
        }

      }
    }
  }

MPI_Barrier(world);
if (comm->me==0 && ntimestep > 333) fprintf(screen, "test3 from pair eam\n");
 
  comm->eforward_comm_pair(this);
MPI_Barrier(world);
if (comm->me==0 && ntimestep > 333) fprintf(screen, "test4 from pair eam\n");
 

  // for finest mesh the electron density calculated here is the same as those calculated by md

 // if (comm->me==0) {
   // for (ii = 0; ii <nlocal; ii++) {
     // for (jj = 0; jj<npe; jj++) fprintf(logfile, "%d, %d,%f\n", tag[ii],jj,erho[ii][jj]);
   // }
  //}
  //error->all(FLERR, "test from pair ecompute");
 
  // For each integration point, first the electron density at that integration site is got by 
  // interpolation. Then the derivative of embedding energy at the integration site is calculated.
  // Then search around the neighbor of the integration site. First the electron density at the neighbor would also be
  // interpolated. Then the derivative of embedding energy is calculated at the neighbor site.
  // Then the derivative of density is calculated
  // pair potential calculated
  // the force on the integration point is calculated
  
   for (ii=0; ii < neinum; ii++) {
     i = enlist[ii];
     np = nelist[ii];
     itype = type[i];
     it = itype - 1;
     nint = numint[itype];

     for (jj=0; jj<nint;jj++) {

       // interpolate jj integration point in i element
       ixtmp = 0;
       iytmp = 0;
       iztmp = 0;

       fx = 0;
       fy = 0;
       fz = 0;

       rhoi = 0;

       for (kk=0; kk<npe; kk++) {
	 ixtmp += Spa_int[it][jj][kk]*nodex[i][kk][0];
	 iytmp += Spa_int[it][jj][kk]*nodex[i][kk][1];
	 iztmp += Spa_int[it][jj][kk]*nodex[i][kk][2];
	 rhoi += Spa_int[it][jj][kk]*erho[i][kk];
       }

       // fpi = derivative of embedding energy at this integration point using the rhoi calculated above

        p = rhoi*rdrho + 1.0;
        m = static_cast<int> (p);
        m = MAX(1, MIN(m, nrho-1));
        p -= m;
        p = MIN(p,1.0);
        coeff = frho_spline[type2frho[1]][m];
        fpi = (coeff[0]*p + coeff[1])*p + coeff[2];

        // get the neighbor information
	
        jnum = numneigh[np];
        jlist = firstneigh[np];
        njlist = nfirstneigh[np];

        // search for all neighbors
	for (ll=0; ll<jnum; ll++) {
          
          // nj atom inside j element
	  j = jlist[ll];
	  nj = njlist[ll];
	  jtype = type[j];
	  jt = jtype-1;

	  // interpolation of nj atom inside j element
	  jxtmp = 0;
	  jytmp = 0;
	  jztmp = 0;
	  rhoj = 0;

	  for (kk=0; kk<npe; kk++) {
	    jxtmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][0];
	    jytmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][1];
	    jztmp += Spa_inter_ele[jt][nj][kk]*nodex[j][kk][2];
	    rhoj += Spa_inter_ele[jt][nj][kk]*erho[j][kk];
	  }

	  delx = ixtmp - jxtmp;
	  dely = iytmp - jytmp;
	  delz = iztmp - jztmp;
	  rsq = delx*delx + dely*dely + delz*delz;

	  if (rsq < cutforcesq) {
	    
	    // first calculate fpj = derivative of embedding energy at this neighbor of integration point using rhoj calculated
	    p = rhoj*rdrho + 1.0;
	    m = static_cast<int> (p);
	    m = MAX(1,MIN(m,nrho-1));
	    p -= m;
	    p = MIN(p,1.0);
	    coeff = frho_spline[type2frho[1]][m];
	    fpj = (coeff[0]*p + coeff[1])*p + coeff[2];

	    // then using the distance between the two atoms calculate:
	    // 1. rhoip = derivative of (density at atom j due to atom i)
	    // 2. pair interaction
	    // phi = pair potential energy
	    // phip = phi'
	    // z2 = phi*r
	    // z2p = (phi * r)' = (phi' r) + phi

	    r = sqrt(rsq);
	    p = r*rdr + 1.0;
	    m = static_cast<int> (p);
	    m = MIN(m,nr-1);
	    p -= m;
	    p = MIN(p,1.0);

	    coeff = rhor_spline[type2rhor[1][1]][m];
	    rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
	    coeff = z2r_spline[type2z2r[1][1]][m];
	    z2p =  (coeff[0]*p + coeff[1])*p + coeff[2];
	    z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

	    recip = 1.0/r;
	    phi = z2*recip;
	    phip = z2p*recip - phi*recip;
	    psip = (fpi + fpj)*rhoip + phip;
	   // psip = phip;
            fpair = -psip*recip;

	    fx += delx*fpair;
	    fy += dely*fpair;
	    fz += delz*fpair;


          }
        }

        for (kk=0; kk <npe; kk++) {
          nodef[i][kk][0] += weight[it][jj]*Spa_int[it][jj][kk]*fx;
	  nodef[i][kk][1] += weight[it][jj]*Spa_int[it][jj][kk]*fy;
	  nodef[i][kk][2] += weight[it][jj]*Spa_int[it][jj][kk]*fz;
	}	
	np++;
     }
   }

//if (comm->me==0 && ntimestep > 333) fprintf(screen, "test5 from pair eam");
 
}	

/*-----------------------------------------------------------------------------------------------*/

int PairEAM::pack_forward_comm(int n, int *list, double *buf,
		               int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
  }
  return m;
}

/*---------------------------------------------------------------------------------------------*/

int PairEAM::epack_forward_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int npe = element->npe;
  m = 0;

  for (i=0; i<n; i++) {
    j = list[i];
    for (int k=0; k <npe; k++) {
      buf[m++] = erho[j][k];
    }
  }
  return m;
}

/*------------------------------------------------------------------------------------------*/

void PairEAM::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp[i] = buf[m++];
}

/*----------------------------------------------------------------------------------------------*/

void PairEAM::eunpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  int npe = element->npe;
  m = 0;
  last = first + n;
  for (i=first; i <last; i++) {
    for (int k = 0; k < npe; k++) {
      erho[i][k] = buf[m++];
    }
  }
}
