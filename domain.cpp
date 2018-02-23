#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "math.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "output.h"
#include "modify.h"
#include "fix.h"
#include "region.h"
#include "region_block.h"
#include "update.h"
#include "thermo.h"
#include "force.h"
#include "universe.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "element.h"

using namespace CAC_NS;
using namespace MathConst;

enum{NO_REMAP,X_REMAP,V_REMAP};    // same as fix_deform.cpp
enum{IGNORE,WARN,ERROR};           // same as thermo.cpp
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

#define BIG   1.0e20
#define SMALL 1.0e-4
#define DELTAREGION 4
#define BONDSTRETCH 1.1

/* ----------------------------------------------------------------------
   default is periodic
------------------------------------------------------------------------- */

Domain::Domain(CAC *cac) : Pointers(cac)
{
  box_exist = 0;

  dimension = 3;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;

  minxlo = minxhi = 0.0;
  minylo = minyhi = 0.0;
  minzlo = minzhi = 0.0;

  triclinic = 0;
  tiltsmall = 1;

  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  xy = xz = yz = 0.0;

  h[3] = h[4] = h[5] = 0.0;
  h_inv[3] = h_inv[4] = h_inv[5] = 0.0;
  h_rate[0] = h_rate[1] = h_rate[2] =
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;

  prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
  prd_half_lamda[0] = prd_half_lamda[1] = prd_half_lamda[2] = 0.5;
  boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
  boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;

//  lattice = NULL;
  char **args = new char*[2];
  args[0] = (char *) "none";
  args[1] = (char *) "1.0";
//  set_lattice(2,args);
  delete [] args;

  nregion = maxregion = 0;
  regions = NULL;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
//  delete lattice;
 for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  // set box_change flags if box size/shape/sub-domains ever change
  // due to shrink-wrapping or fixes that change box size/shape/sub-domains

  box_change_size = box_change_shape = box_change_domain = 0;

  if (nonperiodic == 2) box_change_size = 1;
  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->box_change_size) box_change_size = 1;
    if (modify->fix[i]->box_change_shape) box_change_shape = 1;
    if (modify->fix[i]->box_change_domain) box_change_domain = 1;
  }

  box_change = 0;
  if (box_change_size || box_change_shape || box_change_domain) box_change = 1;

  // check for fix deform

  deform_flag = deform_vremap = deform_groupbit = 0;
  
}

/*-------------------------------------------------------------------------
 create a new region
 -------------------------------------------------------------------------*/

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal region command");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTAREGION;
    regions = (Region **)
	    memory->srealloc(regions, maxregion*sizeof(Region *),"domain:regions");
  }

  // create the region

//  if (strcmp(arg[1],"none")==0) error->all(FLERR, "Unknown region style");

//  #define REGION_CLASS
//  #define RegionStyle(key,Class) \ 
//    else if (strcmp(arg[1],#key)==0) \
//	regions[nregion] = new Class(cac, narg,arg);
//  #include "style_region.h"
//  #undef REGION_CLASS

//    else error->all(FLERR,"Unknown region style");
  regions[nregion] = new RegBlock(cac,narg,arg);

  // initialize any region variables via init()
  // in case region is used between runs, e.g. to print a variable
  
  regions[nregion]->init();
  nregion++;
}

/* ----------------------------------------------------------------------
  (re)set boundary settings
  flag = 0, called from the input script
  flag = 1, called from change box command
 ------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg, int flag)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  char c;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];

      if (c == 'p') boundary[idim][iside] = 0;
      else if (c == 'f') boundary[idim][iside] = 1;
      else if (c == 's') boundary[idim][iside] = 2;
      else if (c == 'm') boundary[idim][iside] = 3;
      else {
        if (flag == 0) error->all(FLERR,"Illegal boundary command");
        if (flag == 1) error->all(FLERR,"Illegal change_box command");
      }
    }

  for (int idim = 0; idim < 3; idim++)
    if ((boundary[idim][0] == 0 && boundary[idim][1]) ||
        (boundary[idim][0] && boundary[idim][1] == 0))
      error->all(FLERR,"Both sides of boundary must be periodic");

  if (boundary[0][0] == 0) xperiodic = 1;
  else xperiodic = 0;
  if (boundary[1][0] == 0) yperiodic = 1;
  else yperiodic = 0;
  if (boundary[2][0] == 0) zperiodic = 1;
  else zperiodic = 0;

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
    nonperiodic = 1;
    if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
        boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
        boundary[2][0] >= 2 || boundary[2][1] >= 2) nonperiodic = 2;
  }
}

/* ----------------------------------------------------------------------
   print box info, orthogonal or triclinic
------------------------------------------------------------------------- */

void Domain::print_box(const char *str)
{
  if (comm->me == 0) {
    if (screen) {
      if (triclinic == 0)
        fprintf(screen,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
      else {
        char *format = (char *)
          "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
        fprintf(screen,format,
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2],
                xy,xz,yz);
      }
    }
    if (logfile) {
      if (triclinic == 0)
        fprintf(logfile,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
      else {
        char *format = (char *)
          "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
        fprintf(logfile,format,
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2],
                xy,xz,yz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi and triclinic tilts are already set
   expandflag = 1 if need to expand box in shrink-wrapped dims
   not invoked by read_restart since box is already expanded
   if don't prevent further expansion, restarted triclinic box
     with unchanged tilt factors can become a box with atoms outside the box
------------------------------------------------------------------------- */

void Domain::set_initial_box(int expandflag)
{
   // error checks for orthogonal and triclinic domains

   if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one(FLERR,"Box bounds are invalid");

  if (domain->dimension == 2 && (xz != 0.0 || yz != 0.0))
    error->all(FLERR,"Cannot skew triclinic box in z for 2d simulation");

  // error check or warning on triclinic tilt factors

  if (triclinic) {
    if ((fabs(xy/(boxhi[0]-boxlo[0])) > 0.5 && xperiodic) || 
        (fabs(xz/(boxhi[0]-boxlo[0])) > 0.5 && xperiodic) ||
        (fabs(yz/(boxhi[1]-boxlo[1])) > 0.5 && yperiodic)) {
      if (tiltsmall)
        error->all(FLERR,"Triclinic box skew is too large");
      else if (comm->me == 0)
        error->warning(FLERR,"Triclinic box skew is large");
    }
  }

  // set small based on box size and SMALL
  // this works for any unit system

  small[0] = SMALL * (boxhi[0] - boxlo[0]);
  small[1] = SMALL * (boxhi[1] - boxlo[1]);
  small[2] = SMALL * (boxhi[2] - boxlo[2]);

  // if expanding, adjust box lo/hi for shrink-wrapped dims

  if (!expandflag) return;
  
  if (boundary[0][0] == 2) boxlo[0] -= small[0];
  else if (boundary[0][0] == 3) minxlo = boxlo[0];
  if (boundary[0][1] == 2) boxhi[0] += small[0];
  else if (boundary[0][1] == 3) minxhi = boxhi[0];

  if (boundary[1][0] == 2) boxlo[1] -= small[1];
  else if (boundary[1][0] == 3) minylo = boxlo[1];
  if (boundary[1][1] == 2) boxhi[1] += small[1];
  else if (boundary[1][1] == 3) minyhi = boxhi[1];

  if (boundary[2][0] == 2) boxlo[2] -= small[2];
  else if (boundary[2][0] == 3) minzlo = boxlo[2];
  if (boundary[2][1] == 2) boxhi[2] += small[2];
  else if (boundary[2][1] == 3) minzhi = boxhi[2];
  
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];

  prd_half[0] = xprd_half = 0.5*xprd;
  prd_half[1] = yprd_half = 0.5*yprd;
  prd_half[2] = zprd_half = 0.5*zprd;

  if (triclinic) {
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);

    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
    boxlo_bound[2] = boxlo[2];

    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
    boxhi_bound[2] = boxhi[2];
  }
}

/* ----------------------------------------------------------------------
    set local subbox params for orthogonal boxes
   assumes global box is defined and proc assignment has been made
   uses comm->xyz_split or comm->mysplit
     to define subbox boundaries in consistent manner
   insure subhi[max] = boxhi
------------------------------------------------------------------------- */

void Domain::set_local_box()
{
  if (triclinic) return;

  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;

  sublo[0] = boxlo[0] + xprd*xsplit[myloc[0]];
  if (myloc[0] < procgrid[0]-1) subhi[0] = boxlo[0] + xprd*xsplit[myloc[0]+1];
  else subhi[0] = boxhi[0];

  sublo[1] = boxlo[1] + yprd*ysplit[myloc[1]];
  if (myloc[1] < procgrid[1]-1) subhi[1] = boxlo[1] + yprd*ysplit[myloc[1]+1];
  else subhi[1] = boxhi[1];

  sublo[2] = boxlo[2] + zprd*zsplit[myloc[2]];
  if (myloc[2] < procgrid[2]-1) subhi[2] = boxlo[2] + zprd*zsplit[myloc[2]+1];
  else subhi[2] = boxhi[2];

}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   adjust 3 image flags encoded in image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before doing remap
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::remap(double *x, imageint &image)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];
  imageint idim,otherdims;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }

  if (xperiodic) {
    while (coord[0] < lo[0]) {
      coord[0] += period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim--;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    while (coord[0] >= hi[0]) {
      coord[0] -= period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim++;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) {
      coord[1] += period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    while (coord[1] >= hi[1]) {
      coord[1] -= period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) {
      coord[2] += period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    while (coord[2] >= hi[2]) {
      coord[2] -= period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);

}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for one atom
   x = H lamda + x0;
   lamda and x can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::lamda2x(double *lamda, double *x)
{
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::x2lamda(double *x, double *lamda)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];

  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
}

/* ----------------------------------------------------------------------
   check warn if any proc's subbox is smaller than thresh
    since may lead to lost atoms in comm->exchange()
  current callers set thresh = neighbor skin
------------------------------------------------------------------------- */

void Domain::subbox_too_small_check(double thresh)
{
  int flag = 0;
  if (!triclinic) {
    if (subhi[0]-sublo[0] < thresh || subhi[1]-sublo[1] < thresh) flag = 1;
    if (dimension == 3 && subhi[2]-sublo[2] < thresh) flag = 1;
  }
  
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0) 
    error->warning(FLERR,"Proc sub-domain size < neighbor skin, "
                   "could lead to lost atoms");
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 or 20 bits for each dimension depending on sizeof(imageint)
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::pbc()
{
  int i;
  imageint idim,otherdims;
  double *lo,*hi,*period;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;

  lo = boxlo;
  hi = boxhi;
  period = prd;

  for (i = 0; i < nlocal; i++) {
    //if (comm->me==0) fprintf(screen, "xperiodic = %d, yperiodic = %d, zperiodic = %d\n", 
    //                                                   xperiodic, yperiodic, zperiodic);
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
        x[i][0] += period[0];
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
        x[i][0] = MAX(x[i][0],lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] -= h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }
    
    if (yperiodic) {
      if (x[i][1] < lo[1]) {
        x[i][1] += period[1];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[5];
          v[i][1] += h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i][1] >= hi[1]) {
        x[i][1] -= period[1];
        x[i][1] = MAX(x[i][1],lo[1]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[5];
          v[i][1] -= h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }
    
    if (zperiodic) {
      if (x[i][2] < lo[2]) {
        x[i][2] += period[2];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[4];
          v[i][1] += h_rate[3];
          v[i][2] += h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x[i][2] >= hi[2]) {
        x[i][2] -= period[2];
        x[i][2] = MAX(x[i][2],lo[2]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[4];
          v[i][1] -= h_rate[3];
          v[i][2] -= h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
    }
  }
}


/*--------------------------------------------------------------------------
  enforce PBC for element and modify box image flags for each element
  called every reneighboring and by other commands that change elements
  resulting coord must satisfy lo <= coord < hi
  MAX is important since coord - prd < lo can happen when coord = hi
  if fix deform, remap velocity of fix group atoms by box edge velocities 
  for triclinic, elements must be in lamda coords (0-1) before pbc is called
  image = 10 or 20 bits for each dimension depending of sizeof(imageint)
  increment/decrement in wrap- around fashion
---------------------------------------------------------------------------*/

void Domain::elepbc()
{
  int i,m;
  imageint idim, otherdims;
  double *lo, *hi, *period;
  int nlocal = element->nlocal;
  double **x = element->x;
  double ***nodex = element->nodex;
  int *mask = element->mask;
  imageint *image = element->image;
  int npe = element->npe;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  }

  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
	 x[i][0] +=period[0];
         for (m = 0; m<npe; m++) { 
           nodex[i][m][0] +=period[0];
	 }

	 idim = image[i] & IMGMASK;
	 otherdims = image[i] ^idim;
	 idim--;
	 idim &= IMGMASK;
	 image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
	for (m=0; m<npe;m++) {
	  nodex[i][m][0] +=period[0];
	}
	x[i][0] = MAX(x[i][0],lo[0]);
	idim = image[i] & IMGMASK;
	otherdims = image[i] ^ idim;
	idim++;
	idim &= IMGMASK;
	image[i] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[i][1] < lo[1]) {
	x[i][1] += period[1];
        for (m=0; m<npe;m++) {
	  nodex[i][m][1] +=period[1];
	}
	idim = (image[i] >> IMGBITS) & IMGMASK;
	otherdims = image[i] ^ (idim << IMGBITS);
	idim--;
	idim &= IMGMASK;
	image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i][1] >= hi[1]) {
	x[i][1] -= period[1];
	for (m=0; m<npe; m++) {
	  nodex[i][m][1] -= period[1];
	}
	x[i][1] = MAX(x[i][1],lo[1]);
        idim = (image[i] >> IMGBITS) & IMGMASK;
	otherdims = image[i] ^ (idim << IMGBITS);
	idim++;
	idim &= IMGMASK;
	image[i] = otherdims | (idim << IMGBITS);
      }	
    }

    if (zperiodic) {
       if (x[i][2] < lo[2]) {
	 x[i][2] +=period[2];
	 for (m=0; m<npe;m++) {
	   nodex[i][m][2] += period[2];
         }
	 idim = image[i] >> IMG2BITS;
	 otherdims = image[i] ^ (idim << IMG2BITS);
	 idim--;
	 idim &= IMGMASK;
	 image[i] = otherdims | (idim << IMG2BITS);
       }
       if (x[i][2] >= hi[2]) {
	 x[i][2] -= period[2];
	 for (m=0; m<npe; m++) {
	    nodex[i][m][2] -= period[2];
	 }
	 x[i][2] = MAX(x[i][2],lo[2]);
	 idim = image[i] >> IMG2BITS;
	 otherdims = image[i] ^ (idim << IMG2BITS);
	 idim++;
	 idim &= IMGMASK;
	 image[i] = otherdims | (idim << IMG2BITS);
       }
    }
  }
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   for triclinic, atoms must be in lamda coords (0-1) before reset_box is called
------------------------------------------------------------------------- */

void Domain::reset_box()
{
  // perform shrink-wrapping
  // compute extent of atoms on this proc
  // for triclinic, this is done in lamda space

  //if (comm->me==0) fprintf(screen, "nonperiodic = %d\n", nonperiodic);
  if (nonperiodic == 2) {
    double extent[3][2],all[3][2];

    extent[2][0] = extent[1][0] = extent[0][0] = BIG;
    extent[2][1] = extent[1][1] = extent[0][1] = -BIG;

    double **x = atom->x;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      extent[0][0] = MIN(extent[0][0],x[i][0]);
      extent[0][1] = MAX(extent[0][1],x[i][0]);
      extent[1][0] = MIN(extent[1][0],x[i][1]);
      extent[1][1] = MAX(extent[1][1],x[i][1]);
      extent[2][0] = MIN(extent[2][0],x[i][2]);
      extent[2][1] = MAX(extent[2][1],x[i][2]);
    }

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // in shrink-wrapped dims, set box by atom extent
    // if minimum set, enforce min box size settings
    // for triclinic, convert lamda extent to box coords, then set box lo/hi
    // decided NOT to do the next comment - don't want to sneakily change tilt
    // for triclinic, adjust tilt factors if 2nd dim is shrink-wrapped,
    //   so that displacement in 1st dim stays the same

    if (triclinic == 0) {
      if (xperiodic == 0) {
        if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = MIN(-all[0][0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = all[0][1] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = MIN(-all[1][0]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = all[1][1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
      }
      if (zperiodic == 0) {
        if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = MIN(-all[2][0]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = all[2][1] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
      }
    } 
  }
  
  // reset box whether shrink-wrapping or not

  set_global_box();
  set_local_box();

}

/* ----------------------------------------------------------------------
   warn if image flags of any bonded atoms are inconsistent
   could be a problem when using replicate or fix rigid
------------------------------------------------------------------------- */

void Domain::image_check()
{
  if (!atom->molecular) return;
}

void Domain::box_too_small_check()
{
  if (!atom->molecular) return;
}

/* ----------------------------------------------------------------------
   format boundary string for output
   assume str is 9 chars or more in length
------------------------------------------------------------------------- */

void Domain::boundary_string(char *str)
{
  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (boundary[idim][iside] == 0) str[m++] = 'p';
      else if (boundary[idim][iside] == 1) str[m++] = 'f';
      else if (boundary[idim][iside] == 2) str[m++] = 's';
      else if (boundary[idim][iside] == 3) str[m++] = 'm';
    }
    str[m++] = ' ';
  }
  str[8] = '\0';
}

/*-------------------------------------------------------------------------
  return region index if name matches existing region ID
  return -1 if no such region
---------------------------------------------------------------------------*/

int Domain::find_region(char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
}
