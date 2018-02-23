#include "mpi.h"
#include "string.h"
#include "compute_pe.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "domain.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputePE::ComputePE(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute pe command");
  if (igroup) error->all(FLERR,"Compute pe must use group all");

  scalar_flag = 1;
  extscalar = 1;
  peflag = 1;
  timeflag = 1;

  if (narg == 3) {
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = 1;
    thermoflag = 1;
  } else {
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = 0;
    thermoflag = 0;
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else error->all(FLERR,"Illegal compute pe command");
      iarg++;
    }
  }
}


