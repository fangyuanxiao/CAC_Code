#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "compute.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace CAC_NS;

#define DELTA 4
#define BIG MAXTAGINT

// allocate space for static class instance variable and initialize it

int Compute::instance_total = 0;

/* ---------------------------------------------------------------------- */

Compute::Compute(CAC *cac, int narg, char **arg) : Pointers(cac)
{
  instance_me = instance_total++;

  if (narg < 3) error->all(FLERR,"Illegal compute command");

  // compute ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
                 "Compute ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find compute group ID");
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;
  size_vector_variable = size_array_rows_variable = 0;

  tempflag = pressflag = peflag = 0;
  pressatomflag = peatomflag = 0;
  create_attribute = 0;
  tempbias = 0;

  timeflag = 0;
  comm_forward = comm_reverse = 0;
  dynamic = 0;
  dynamic_group_allow = 1;
  cudable = 0;

  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_peratom = invoked_local = -1;
  invoked_flag = 0;

  // set modify defaults

  extra_dof = domain->dimension;
  dynamic_user = 0;
  fix_dof = 0;

  // setup list of timesteps

  ntime = maxtime = 0;
  tlist = NULL;

  // data masks

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;

  // force init to zero in case these are used as logicals

  vector = vector_atom = vector_local = NULL;
  array = array_atom = array_local = NULL;

}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  delete [] id;
  delete [] style;
  memory->destroy(tlist);
}

/* ----------------------------------------------------------------------
  reset extra_dof to its default value
 ------------------------------------------------------------------------- */

void Compute::reset_extra_dof()
{
  extra_dof = domain->dimension;
}

/* ----------------------------------------------------------------------
   add ntimestep to list of timesteps the compute will be called on
   do not add if already in list
   search from top downward, since list of times is in decreasing order
------------------------------------------------------------------------- */

void Compute::addstep(bigint ntimestep)
{
  // i = location in list to insert ntimestep

  int i;
  for (i = ntime-1; i >= 0; i--) {
    if (ntimestep == tlist[i]) return;
    if (ntimestep < tlist[i]) break;
  }
  i++;

  // extend list as needed

  if (ntime == maxtime) {
    maxtime += DELTA;
    memory->grow(tlist,maxtime,"compute:tlist");
  }

  // move remainder of list upward and insert ntimestep

  for (int j = ntime-1; j >= i; j--) tlist[j+1] = tlist[j];
  tlist[i] = ntimestep;
  ntime++;

}
