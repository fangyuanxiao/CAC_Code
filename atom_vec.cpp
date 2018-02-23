#include "string.h"
#include "stdlib.h"
#include "atom_vec.h"
#include "atom.h"
#include "domain.h"
#include "error.h"

using namespace CAC_NS;

#define DELTA 16384
#define DELTA_BONUS 8192

/*----------------------------------------------------------------------*/

AtomVec::AtomVec(CAC *cac) : Pointers(cac)
{
  nmax = 0;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 0;
  mass_type = dipole_type = 0;
  forceclearflag = 0;
  size_data_bonus = 0;
  cudable = kokkosable = 0;
  nargcopy = 0;
  argcopy = NULL;
}

/*-----------------------------------------------------------------------*/

AtomVec::~AtomVec()
{
  for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
  delete [] argcopy;
}

/* ----------------------------------------------------------------------
     make copy of args for use by restart & replicate
    -------------------------------------------------------------------------*/

void AtomVec::store_args(int narg, char **arg)
{
  nargcopy = narg;
  argcopy = new char*[nargcopy];
  for (int i = 0; i < nargcopy; i++) {
    int n = strlen(arg[i]) + 1;
    argcopy[i] = new char[n];
    strcpy(argcopy[i],arg[i]);
  }
}

/* ----------------------------------------------------------------------
    no additional args by default
  ------------------------------------------------------------------------- */

void AtomVec::process_args(int narg, char **arg)
{
  if (narg) error->all(FLERR,"Invalid atom_style command");
}

/* ----------------------------------------------------------------------
  grow nmax so it is a multiple of DELTA
 ------------------------------------------------------------------------- */

void AtomVec::grow_nmax()
{
  nmax = nmax/DELTA * DELTA;
  nmax += DELTA;
}

/* ----------------------------------------------------------------------
   copy of velocity remap settings from Domain
------------------------------------------------------------------------- */

void AtomVec::init()
{
  deform_vremap = domain->deform_vremap;
  deform_groupbit = domain->deform_groupbit;
  h_rate = domain->h_rate;
}


