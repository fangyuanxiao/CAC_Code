#include "string.h"
#include "stdlib.h"
#include "update.h"
#include "integrate.h"
#include "style_integrate.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
//#include "region.h"
#include "compute.h"
#include "output.h"
#include "comm.h"
#include "universe.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

Update::Update(CAC *cac) : Pointers(cac)
{

  char *str;

  ntimestep = 0;
  atime = 0.0;
  atimestep = 0;
  first_update = 0;

  whichflag = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  setupflag = 0;
  multireplica = 0;

  restrict_output = 0;

  eflag_global = vflag_global = -1;

  unit_style = NULL;
  set_units("metal");

  integrate_style = NULL;
  integrate = NULL;

  str = (char *) "verlet";
  create_integrate(1,&str,1);
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  delete [] unit_style;

  delete [] integrate_style;
  delete integrate;
}

/* ---------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg, int trysuffix)
{
  if (narg < 1) error->all(FLERR,"Illegal run_style command");

  delete [] integrate_style;
  delete integrate;

  int sflag;
  new_integrate(arg[0],narg-1,&arg[1],trysuffix,sflag);

  int n = strlen(arg[0]) + 1;
  integrate_style = new char[n];
  strcpy(integrate_style,arg[0]);
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{

   if (strcmp(style,"lj") == 0) {

    force->boltz = 1.0;
    force->hplanck = 0.18292026;  // using LJ parameters for argon
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 0.005;
    neighbor->skin = 0.3;
   } else if (strcmp(style,"metal") == 0) {
    force->boltz = 8.617343e-5;
    force->hplanck = 4.135667403e-3;
    force->mvv2e = 1.0364269e-4;
    force->ftm2v = 1.0 / 1.0364269e-4;
    force->mv2d = 1.0 / 0.602214129;
    force->nktv2p = 1.6021765e6;
    force->qqr2e = 14.399645;
    force->qe2f = 1.0;
    force->vxmu2f = 0.6241509647;
    force->xxt2kmu = 1.0e-4;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0e-3;
    force->qelectron = 1.0;
    dt = 0.001;
    neighbor->skin = 2.0;
  }else error->all(FLERR,"Illegal units command");
 
  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ----------------------------------------------------------------------
   create the Integrate style, first with suffix appended
 ------------------------------------------------------------------------- */

void Update::new_integrate(char *style, int narg, char **arg,
                           int trysuffix, int &sflag)
{
  sflag = 0;
  if (0) return;

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
  else if (strcmp(style,#key) == 0) integrate = new Class(cac,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS

  else error->all(FLERR,"Illegal integrate style");

}

/* ----------------------------------------------------------------------
   reset timestep as called from input script
------------------------------------------------------------------------- */

//void Update::reset_timestep(int narg, char **arg)
//{
//  if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
//  bigint newstep = universe->bnumeric(FLERR,arg[0]);
//  if (comm->me==0) fprintf(screen, "%d\n",newstep);
//  reset_timestep(newstep);
//}

/* ----------------------------------------------------------------------
   reset timestep
   called from rerun command and input script (indirectly)
------------------------------------------------------------------------- */

//void Update::reset_timestep(bigint newstep)
//{
//  ntimestep = newstep;
//  if (ntimestep < 0) error->all(FLERR,"Timestep must be >= 0");
//  if (ntimestep > MAXBIGINT) error->all(FLERR,"Too big a timestep");

  // set atimestep to new timestep
  // so future update_time() calls will be correct

//  atimestep = ntimestep;

  // trigger reset of timestep for output
  // do not allow any timestep-dependent fixes to be already defined

//  output->reset_timestep(ntimestep);

//}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // if USER-CUDA mode is enabled:
  // integrate/minimize style must be CUDA variant
  //if (comm->me==0) fprintf(screen, "whichflag = %d\n", whichflag);
  if (whichflag == 0) return;
  if (whichflag == 1) integrate->init();

}

bigint Update::memory_usage()
{
  bigint bytes = 0;
  bytes += integrate->memory_usage();
  return bytes;
}

/* ----------------------------------------------------------------------
   update elapsed simulation time
   called at end of runs or when timestep size changes
------------------------------------------------------------------------- */

void Update::update_time()
{
  atime += (ntimestep-atimestep) * dt;
  atimestep = ntimestep;
}
