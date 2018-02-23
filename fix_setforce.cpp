#include "string.h"
#include "stdlib.h"
#include "fix_setforce.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
//#include "comm.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/*-------------------------------------------------------*/

FixSetForce::FixSetForce(CAC *cac, int narg, char **arg) :
  Fix(cac,narg,arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix setforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  sforce = NULL;
  
  xstr = ystr = zstr = NULL;

  if (strcmp(arg[3],"NULL")==0) {
    xstyle = NONE;
  } else {
    xvalue = universe->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }

  if (strcmp(arg[4],"NULL")==0) {
    ystyle = NONE;
  } else {
    yvalue = universe->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }

  if (strcmp(arg[5],"NULL")==0) {
    zstyle = NONE;
  } else {
    zvalue = universe->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  
  varflag = CONSTANT;

  maxatom = atom->nmax;
  if (maxatom) memory->create(sforce,maxatom,3,"setforce:sforce");
}

/*----------------------------------------------------------------------------*/

FixSetForce::~FixSetForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/*--------------------------------------------------------------------------*/

int FixSetForce::setmask()
{
  //error->all(FLERR, "test from fix_setforce");
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/*---------------------------------------------------------------------------*/

void FixSetForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    if (atom->nlocal) post_force(vflag);
    if (element->nlocal) epost_force(vflag);
  }
}	

/*----------------------------------------------------------------------------*/

void FixSetForce::post_force(int vflag)
{
  //if (comm->me==0) fprintf(screen, "vflag = %d\n", vflag);
  //error->all(FLERR, "test from post_force in FixSetForce");
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (int i=0; i <  nlocal; i++)
      if (mask[i] & groupbit) {
       foriginal[0] += f[i][0];
       foriginal[1] += f[i][1];
       foriginal[2] += f[i][2];
       if (xstyle) f[i][0] = xvalue;
       if (ystyle) f[i][1] = yvalue;
       if (zstyle) f[i][2] = zvalue;
      }
  }

}

/*---------------------------------------------------------------------------*/

void FixSetForce::epost_force(int vflag)
{
  double **x = element->x;
  double ***nodef = element->nodef;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int npe = element->npe;
  int i,j;

  eforiginal[0] = eforiginal[1] = eforiginal[2] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (i=0;i<nlocal;i++) 
      if (mask[i] & groupbit) {
       for (j=0;j<npe;j++) {
	eforiginal[0] += nodef[i][j][0];
	eforiginal[1] += nodef[i][j][1];
	eforiginal[2] += nodef[i][j][2];

	if (xstyle) nodef[i][j][0] = xvalue;
	if (ystyle) nodef[i][j][1] = yvalue;
	if (zstyle) nodef[i][j][2] = zvalue;
       }
      }
  }
}
