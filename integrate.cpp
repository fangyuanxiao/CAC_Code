#include "stdlib.h"
#include "integrate.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "compute.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

Integrate::Integrate(CAC *cac, int narg, char **arg) : Pointers(cac)
{
  elist_global = elist_atom = NULL;
  vlist_global = vlist_atom = NULL;
  external_force_clear = 0;
}

/* ---------------------------------------------------------------------- */

Integrate::~Integrate()
{
  delete [] elist_global;
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
}

/* ---------------------------------------------------------------------- */

void Integrate::init()
{
  update->atimestep = update->ntimestep;

  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0; 

}

/* ----------------------------------------------------------------------
   setup lists of computes for global and per-atom PE and pressure
------------------------------------------------------------------------- */

void Integrate::ev_setup()
{
  delete [] elist_global;
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  elist_global = elist_atom = NULL;
  vlist_global = vlist_atom = NULL;

  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peflag) nelist_global++;
    if (modify->compute[i]->peatomflag) nelist_atom++;
    if (modify->compute[i]->pressflag) nvlist_global++;
    if (modify->compute[i]->pressatomflag) nvlist_atom++;
  }  

  if (nelist_global) elist_global = new Compute*[nelist_global];
  if (nelist_atom) elist_atom = new Compute*[nelist_atom];
  if (nvlist_global) vlist_global = new Compute*[nvlist_global];
  if (nvlist_atom) vlist_atom = new Compute*[nvlist_atom];

  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peflag)
      elist_global[nelist_global++] = modify->compute[i];
    if (modify->compute[i]->peatomflag)
      elist_atom[nelist_atom++] = modify->compute[i];
    if (modify->compute[i]->pressflag)
      vlist_global[nvlist_global++] = modify->compute[i];
    if (modify->compute[i]->pressatomflag)
      vlist_atom[nvlist_atom++] = modify->compute[i];
  }
}

/* ----------------------------------------------------------------------
   set eflag,vflag for current iteration
   invoke matchstep() on all timestep-dependent computes to clear their arrays
   eflag/vflag based on computes that need info on this ntimestep
   eflag = 0 = no energy computation
   eflag = 1 = global energy only
   eflag = 2 = per-atom energy only
   eflag = 3 = both global and per-atom energy
   vflag = 0 = no virial computation (pressure)
   vflag = 1 = global virial with pair portion via sum of pairwise interactions
   vflag = 2 = global virial with pair portion via F dot r including ghosts
   vflag = 4 = per-atom virial only
   vflag = 5 or 6 = both global and per-atom virial
------------------------------------------------------------------------- */

void Integrate::ev_set(bigint ntimestep)
{
  int i,flag;
  flag = 0;
  int eflag_global = 0;
  
  flag = 0;
  int eflag_atom = 0;
  
  eflag = eflag_global + eflag_atom;
  flag = 0;
  int vflag_global = 0;

  flag = 0;
  int vflag_atom = 0;
 
  vflag = vflag_global + vflag_atom;
}
