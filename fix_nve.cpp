#include "stdio.h"
#include "string.h"
#include "fix_nve.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

/*-----------------------------------------------------------------*/

FixNVE::FixNVE(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVE::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

//  if (strstr(update->integrate_style,"respa"))
//    step_respa = ((Respa *) update->integrate)->step;

}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVE::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  //if (igroup == atom->firstgroup) 
    //if (comm->me==0) fprintf(screen, "testing\n");
    if (rmass) {
    //if (comm->me==0) fprintf(screen, "testing\n");
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  } else {
    //if (comm->me==0) fprintf(screen, "testing\n");
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        //if (comm->me==0) fprintf(screen, "testing\n");
        dtfm = dtf / mass[type[i]];
	if (comm->me==0) fprintf(screen, "dtf = %g,mass = %g\n",dtf, mass[type[i]]);
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

void FixNVE::einitial_integrate(int vflag) {
   double dtfm;
   
   // update v and x of elements in group
   
   double **x = element->x;
   double ***nodex = element->nodex;
   double ***nodev = element->nodev;
   double ***nodef = element->nodef;
   double *mass = element->mass;
   int *type = element->type;
   int *mask = element->mask;
   //int *tag = element->tag;
   //int itag;
   int nlocal = element->nlocal;
   // double *rmass = atom->rmass;
   int npe = element->npe;
   double inv_npe = 1/static_cast<double>(npe);
   int i,j;
   //if (comm->me==0) fprintf(logfile, "inv_npe:%f, %e\n",inv_npe, dtf);
  
   for (i = 0; i<nlocal; i++) {
    if (mask[i]&groupbit) {
     // calculate the distributed mass on each node
     dtfm = dtf / (inv_npe*mass[type[i]]);
    // if (comm->me==0) fprintf(screen, "dtf = %g, mass = \n",dtf, inv_npe);
     for (j = 0; j<npe; j++) {
      //if (comm->me==0) fprintf(logfile, "f: %d, %d, %e, %e, %e\n", itag,j,nodef[i][j][0],nodef[i][j][1],nodef[i][j][2]);
      nodev[i][j][0] += dtfm * nodef[i][j][0];
      nodev[i][j][1] += dtfm * nodef[i][j][1];
      nodev[i][j][2] += dtfm * nodef[i][j][2];
      //if (comm->me==0) fprintf(logfile, "v: %e, %e, %e\n", nodev[i][j][0],nodev[i][j][1],nodev[i][j][2]);
      nodex[i][j][0] += dtv * nodev[i][j][0];
      nodex[i][j][1] += dtv * nodev[i][j][1];
      nodex[i][j][2] += dtv * nodev[i][j][2];
      //if (comm->me==0) fprintf(logfile, "x: %e, %e, %e\n", nodex[i][j][0],nodex[i][j][1],nodex[i][j][2]);
     }
    }
   }
   // update center point of elements
   element_center_clear();
   //error->all(FLERR, "test from enitial_integrate");
   for ( i = 0; i<nlocal;i++) {
   // itag = tag[i];
    for ( j = 0; j<npe;j++) {
     x[i][0] += nodex[i][j][0];
     x[i][1] += nodex[i][j][1];
     x[i][2] += nodex[i][j][2];
    }
    x[i][0] = x[i][0]/npe;
    x[i][1] = x[i][1]/npe;
    x[i][2] = x[i][2]/npe;
   // if (itag==5) fprintf(logfile,"x: %d, %e, %e, %e\n", itag, x[i][0],x[i][1],x[i][2]);
   }

   //error->all(FLERR, "test from einitial_integrate");

}

/*-----------------------------------------------------------------------*/
void FixNVE::element_center_clear()
{
   int i;
   size_t nbytes;

   int nlocal = element->nlocal;

   nbytes = sizeof(double) *nlocal;
   if (nbytes) {
     memset(&element->x[0][0], 0,3*nbytes);
   }
//   double **x = element->x;
//   if (comm->me==0) 
//      for (int ii=0; ii<nlocal; ii++)
//	fprintf(logfile, "%f,%f,%f\n",x[ii][0],x[ii][1],x[ii][2]);

//   error->all(FLERR, "test from element_center_clear");
}

/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }
}

/*------------------------------------------------------------------*/

void FixNVE::efinal_integrate()
{
  double dtfm;

  // update v of nodes in group

  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  double *mass = element->mass;
  int *type = element->type;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int npe = element->npe;
  double inv_npe = 1/static_cast<double>(npe);

  for (int i = 0; i<nlocal; i++) {
   if (mask[i]&groupbit) {
     dtfm = dtf / (inv_npe*mass[type[i]]);
     for (int j = 0; j<npe; j++) {
       nodev[i][j][0] += dtfm * nodef[i][j][0];
       nodev[i][j][1] += dtfm * nodef[i][j][1];
       nodev[i][j][2] += dtfm * nodef[i][j][2];
     }
   }	    
  }

}
