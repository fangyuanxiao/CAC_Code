#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "velocity.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "domain.h"
#include "input.h"
#include "universe.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

enum{CREATE,SET,SCALE,RAMP,ZERO};
enum{ALL,LOCAL,GEOM};
enum{NONE,CONSTANT,EQUAL,ATOM};

#define WARMUP 100
#define SMALL 0.001

/*----------------------------------------------------------------------------*/

Velocity::Velocity(CAC *cac) : Pointers(cac) {}

/*----------------------------------------------------------------------------*/

void Velocity::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal velocity command");

  if (domain->box_exist==0)
    error->all(FLERR,"Velocity command before simulation box is defined");
  if (atom->natoms == 0 && cac->element_flag==0)
    error->all(FLERR,"Velocity command with no atoms and elements existing");
  
  // atoms masses must all be set

  atom->check_mass();

  //identify group

  igroup = group->find(arg[0]);
  if (igroup==-1) error->all(FLERR,"Cound not find velocity group ID");
  groupbit = group->bitmask[igroup];

  // identify style

  if (strcmp(arg[1],"set")==0) style = SET;
  else if (strcmp(arg[1],"ramp")==0) style = RAMP;
  else error->all(FLERR, "Illegal velocity command");

  // set defaults

  temperature = NULL;
  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 0;
  rotation_flag = 0;
  bias_flag = 0;
  loop_flag = ALL;
  scale_flag = 1;
  rfix = -1;

  // read options from end of input line
  // change defaults as options specify

  if (style == SET) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);

  // initialize velocities based on style
  // create() invoked differently, so can be called externally

  if (style == SET) set(narg-2,&arg[2]);
  else if (style == RAMP) ramp(narg-2,&arg[2]);



  //error->all(FLERR, "test from velocity command");
}

/*------------------------------------------------------------------
  parse optional parameters at end of velocity input line
--------------------------------------------------------------------*/

void Velocity::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal velocity command");

  int iarg = 0;
  while (iarg < narg) {
   if (strcmp(arg[iarg],"sum")==0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
     if (strcmp(arg[iarg+1],"no")==0) sum_flag = 0;
     else if (strcmp(arg[iarg+1],"yes")==0) sum_flag = 1;
     else error->all(FLERR, "Illegal velocity command");
     iarg+=2;
   }else if (strcmp(arg[iarg],"units")==0) {
     if (iarg+2 > narg) error->all(FLERR, "Illegal velocity command");
     if (strcmp(arg[iarg+1],"box")==0) scale_flag = 0;
     else if (strcmp(arg[iarg+1],"lattice")==0) scale_flag = 1;
     else error->all(FLERR, "Illegal velocity command");
     iarg += 2;
   } else error->all(FLERR, "Illegal velocity command");
  }


}

/*-------------------------------------------------------------------*/

void Velocity::set(int narg, char **arg)
{
  int xstyle, ystyle, zstyle, varflag;
  double vx, vy, vz;
  char *xstr, *ystr, *zstr;
  int xvar,yvar,zvar;

  xstyle = ystyle = zstyle = CONSTANT;
  xstr = ystr = zstr = NULL;

  if (strcmp(arg[0],"NULL") == 0) xstyle = NONE;
  else vx = universe->numeric(FLERR,arg[0]);

  if (strcmp(arg[1],"NULL")==0) ystyle = NONE;
  else vy = universe->numeric(FLERR, arg[1]);

  if (strcmp(arg[2],"NULL")==0) zstyle = NONE;
  else vz = universe->numeric(FLERR,arg[2]);

  // set and apply scale factors

  xscale = yscale = zscale = 1.0;

  if (xstyle && !xstr) {
    vx *= xscale;
  }

  if (ystyle && !ystr) {
    vy *= yscale;
  }

  if (zstyle && !zstr) {
    vz *= zscale;
  }

  varflag = CONSTANT;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (varflag == CONSTANT) {
    for (int i = 0; i <nlocal; i++) {
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
	  if (xstyle) v[i][0] = vx;
	  if (ystyle) v[i][1] = vy;
	  if (zstyle) v[i][2] = vz;
        } else {
	  if (xstyle) v[i][0] += vx;
	  if (ystyle) v[i][1] += vy;
	  if (zstyle) v[i][2] += vz;
	}
      }
    }
  }

  double ***nodev = element->nodev;
  int *emask = element->mask;
  int enlocal = element->nlocal;
  int npe = element->npe;

  if (varflag == CONSTANT) {
    for (int i=0; i<enlocal;i++) {
      if (emask[i] & groupbit) {
        if (sum_flag ==0) {
          for (int j=0;j<npe;j++) {
             if (xstyle) nodev[i][j][0] = vx;
	     if (ystyle) nodev[i][j][1] = vy;
	     if (zstyle) nodev[i][j][2] = vz;
	  }
        } else {
	  for (int j=0; j<npe;j++) {
	    if (xstyle) nodev[i][j][0] +=vx;
	    if (ystyle) nodev[i][j][1] +=vy;
	    if (zstyle) nodev[i][j][2] +=vz;
	  }
	}
      }
    }
  }

  // clean up

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
}

/*---------------------------------------------------------------------------
   apply a ramped set of velocities
-----------------------------------------------------------------------------*/

void Velocity::ramp(int narg, char **arg)
{
  // set scale factors
  xscale = yscale = zscale = 1.0;

  int v_dim;
  if (strcmp(arg[0],"vx")==0) v_dim = 0;
  else if (strcmp(arg[0],"vy")==0) v_dim=1;
  else if (strcmp(arg[0],"vz")==0) v_dim=2;
  else error->all(FLERR, "Illegal velocity command");

  if (v_dim == 2 && domain->dimension == 2)
    error->all(FLERR, "Velocity ramp in z for a 2d problem");

  double v_lo, v_hi;
  if (v_dim==0) {
    v_lo = xscale*universe->numeric(FLERR,arg[1]);
    v_hi = xscale*universe->numeric(FLERR,arg[2]);
  } else if (v_dim==1) {
    v_lo = yscale*universe->numeric(FLERR,arg[1]);
    v_hi = yscale*universe->numeric(FLERR,arg[2]);
  } else if (v_dim==2) {
    v_lo = zscale*universe->numeric(FLERR, arg[1]);
    v_hi = zscale*universe->numeric(FLERR, arg[2]);
  }

  int coord_dim;
  if (strcmp(arg[3],"x")==0) coord_dim=0;
  else if (strcmp(arg[3],"y")==0) coord_dim=1;
  else if (strcmp(arg[3],"z")==0) coord_dim=2;
  else error->all(FLERR, "Illegal velocity command");

  double coord_lo,coord_hi;
  if (coord_dim==0) {
    coord_lo = xscale*universe->numeric(FLERR,arg[4]);
    coord_hi = xscale*universe->numeric(FLERR,arg[5]);
  } else if (coord_dim==1) {
    coord_lo = yscale*universe->numeric(FLERR,arg[4]);
    coord_hi = yscale*universe->numeric(FLERR,arg[5]);
  } else if (coord_dim==2) {
    coord_lo = zscale*universe->numeric(FLERR,arg[4]);
    coord_hi = zscale*universe->numeric(FLERR,arg[5]);
  }

  // vramp = ramped velocity component for v_dim
  // add or set based on sum_flag

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fraction, vramp;

  for (int i=0; i<nlocal;i++) {
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim]-coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi-v_lo);
      if (sum_flag) v[i][v_dim] += vramp;
      else v[i][v_dim] = vramp;
    }
  }

  // give velocity to element
  double ***nodev = element->nodev;
  double ***nodex = element->nodex;
  int *emask = element->mask;
  int enlocal = element->nlocal;
  int npe = element->npe;

  for (int i=0; i<enlocal; i++) {
    if (emask[i] & groupbit) {
      for (int j=0; j<npe;j++) {
	fraction = (nodex[i][j][coord_dim]-coord_lo) / (coord_hi - coord_lo);
	fraction = MAX(fraction,0.0);
	fraction = MIN(fraction,1.0);
	vramp = v_lo + fraction*(v_hi-v_lo);
	if (sum_flag) nodev[i][j][v_dim] += vramp;
	else nodev[i][j][v_dim] = vramp;
      }
    }
  }
  //error->all(FLERR, "test from ramp in velocity");
}

