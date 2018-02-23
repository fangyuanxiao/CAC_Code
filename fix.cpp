#include "string.h"
#include "ctype.h"
#include "fix.h"
#include "atom.h"
#include "group.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

// allocate space for static class instance variable and initialize it

  int Fix::instance_total = 0;

/*----------------------------------------------------------------------*/

Fix::Fix(CAC *cac, int narg, char **arg) : Pointers(cac)
{
  instance_me = instance_total++;
  
  // fix ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Fix ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  restart_global = restart_peratom = restart_file = 0;
  force_reneighbor = 0;
  box_change_size = box_change_shape = box_change_domain = 0;
  thermo_energy = 0;
  rigid_flag = 0;
  virial_flag = 0;
  no_change_box = 0;
  time_integrate = 0;
  time_depend = 0;
  create_attribute = 0;
  restart_pbc = 0;
  wd_header = wd_section = 0;
  dynamic_group_allow = 0;
  dof_flag = 0;
  special_alter_flag = 0;
  cudable_comm = 0;

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;
  size_vector_variable = size_array_rows_variable = 0;

  comm_forward = comm_reverse = comm_border = 0;
  restart_reset = 0;

  // reasonable defaults
  // however, each fix that uses these values should explicitly set them

  nevery = 1;
  global_freq = 1;

  // per-atom virial
  // set vflag_atom = 0 b/c some fixes grow vatom in grow_arrays()
  //   which may occur outside of timestepping

  maxvatom = 0;
  vatom = NULL;
  vflag_atom = 0;

  copymode = 0;

}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
  if (copymode) return;

  delete [] id;
  delete [] style;
  memory->destroy(vatom);
}
