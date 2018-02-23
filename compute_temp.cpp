#include "mpi.h"
#include "string.h"
#include "compute_temp.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeTemp::ComputeTemp(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp command");
  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTemp::~ComputeTemp()
{
  delete [] vector;
}
