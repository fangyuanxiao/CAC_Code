#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_viscous.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

/*---------------------------------------------------------------------------*/

FixViscous::FixViscous(CAC *cac,int narg, char **arg) :
  Fix(cac,narg,arg)
{
  //error->all(FLERR, "test from FixViscous");
  if (narg < 4) error->all(FLERR,"Illegal fix viscous command");

  double gamma_one = universe->numeric(FLERR,arg[3]);

  egamma = new double[element->ntypes+1];
  

  gamma = new double[atom->ntypes+1];
  for (int i = 1; i<=atom->ntypes;i++) gamma[i] = gamma_one;
  //error->all(FLERR,"test from FixViscous");
  // for element the coefficient is multiplied by the propotion mass on each node
  for (int i=0; i<element->ntypes;i++) egamma[i+1] = gamma_one;
 // error->all(FLERR, "test from FixViscous");
}

/*---------------------------------------------------------------------------*/

FixViscous::~FixViscous()
{
  delete [] gamma;
  delete [] egamma;
}

/*---------------------------------------------------------------------------*/

int FixViscous::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask; 
}

/*----------------------------------------------------------------------------*/

void FixViscous::init()
{

}

/*----------------------------------------------------------------------------*/

void FixViscous::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    if (atom->nlocal) post_force(vflag);
    if (element->nlocal) epost_force(vflag);
  }
}

/*---------------------------------------------------------------------------*/

void FixViscous::post_force(int vflag)
{
  // apply drag force to atoms in group
  // diretion is opposed to velocity vector
  // magnitude depends on atom type

  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double drag;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      drag = gamma[type[i]];
      f[i][0] -= drag*v[i][0];
      f[i][1] -= drag*v[i][1];
      f[i][2] -= drag*v[i][2];
    }
}

/*----------------------------------------------------------------------------*/
void FixViscous::epost_force(int vflag)
{
   // apply drag force to elements in group
   // direction is opposed to velocity vector
   // magnitude depends on element type

   double ***nodev = element->nodev;
   double ***nodef = element->nodef;
   int *mask = element->mask;
   int *type = element->type;
   int nlocal = element->nlocal;
   int npe = element->npe;
   int i,j;
   double drag;

   for (i = 0; i<nlocal;i++) {
     if (mask[i] & groupbit) {
       drag = egamma[type[i]];
       for (j=0;j<npe;j++) {
	nodef[i][j][0] -= drag*nodev[i][j][0];
	nodef[i][j][1] -= drag*nodev[i][j][1];
	nodef[i][j][2] -= drag*nodev[i][j][2];
       }
     }
   }

}
