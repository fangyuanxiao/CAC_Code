#include "stdio.h"
#include "string.h"
#include "modify.h"
#include "style_compute.h"
#include "style_fix.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "fix.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

#define DELTA 4
#define BIG 1.0e20
#define NEXCEPT 5       // change when add to exceptions in add_fix()

/*------------------------------------------------------------------------------*/

Modify::Modify(CAC *cac) : Pointers(cac)
{

  nfix = maxfix = 0;
  n_initial_integrate = n_post_integrate = 0;
  n_pre_exchange = n_pre_neighbor = 0;
  n_pre_force = n_post_force = 0;
  n_final_integrate = n_end_of_step = n_thermo_energy = 0;
  n_initial_integrate_respa = n_post_integrate_respa = 0;
  n_pre_force_respa = n_post_force_respa = n_final_integrate_respa = 0;
  n_min_pre_exchange = n_min_pre_force = n_min_post_force = n_min_energy = 0;

  fix = NULL;
  fmask = NULL;
  list_initial_integrate = list_post_integrate = NULL;
  list_pre_exchange = list_pre_neighbor = NULL;
  list_pre_force = list_post_force = NULL;
  list_final_integrate = list_end_of_step = NULL;
  list_thermo_energy = NULL;
  list_initial_integrate_respa = list_post_integrate_respa = NULL;
  list_pre_force_respa = list_post_force_respa = NULL;
  list_final_integrate_respa = NULL;
  list_min_pre_exchange = list_min_pre_neighbor = NULL;
  list_min_pre_force = list_min_post_force = NULL;
  list_min_energy = NULL;

  end_of_step_every = NULL;

  list_timeflag = NULL;

  nfix_restart_global = 0;
  id_restart_global = style_restart_global = state_restart_global = NULL;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = NULL;
  index_restart_peratom = NULL;

  ncompute = maxcompute = 0;
  compute = NULL;

  // fill map with fixes listed in style_fix.h

  fix_map = new std::map<std::string,FixCreator>();

#define FIX_CLASS
#define FixStyle(key,Class) \
  (*fix_map)[#key] = &fix_creator<Class>;
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS

  // fill map with computes listed in style_compute.h

  compute_map = new std::map<std::string,ComputeCreator>();

#define COMPUTE_CLASS
#define ComputeStyle(key,Class) \
  (*compute_map)[#key] = &compute_creator<Class>;
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS

}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
   // delete all fixes
   // do it via delete_fix() so callbacks in Atom are also updated correctly

  while (nfix) delete_fix(fix[0]->id);
  memory->sfree(fix);
  memory->destroy(fmask);

  // delete all computes

  for (int i = 0; i < ncompute; i++) delete compute[i];
  memory->sfree(compute);

  delete [] list_initial_integrate;
  delete [] list_post_integrate;
  delete [] list_pre_exchange;
  delete [] list_pre_neighbor;
  delete [] list_pre_force;
  delete [] list_post_force;
  delete [] list_final_integrate;
  delete [] list_end_of_step;
  delete [] list_thermo_energy;
  delete [] list_initial_integrate_respa;
  delete [] list_post_integrate_respa;
  delete [] list_pre_force_respa;
  delete [] list_post_force_respa;
  delete [] list_final_integrate_respa;
  delete [] list_min_pre_exchange;
  delete [] list_min_pre_neighbor;
  delete [] list_min_pre_force;
  delete [] list_min_post_force;
  delete [] list_min_energy;

  delete [] end_of_step_every;
  delete [] list_timeflag;

  restart_deallocate();

  delete compute_map;
  delete fix_map;
}

/* ----------------------------------------------------------------------
one instance per fix in style_fix.h
------------------------------------------------------------------------- */

template <typename T>
Fix *Modify::fix_creator(CAC *cac, int narg, char **arg)
{
  return new T(cac,narg,arg);
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
   Atom class must update indices in its list of callbacks to fixes
 ------------------------------------------------------------------------- */

void Modify::delete_fix(const char *id)
{
  int ifix = find_fix(id);
  if (ifix < 0) error->all(FLERR,"Could not find fix ID to delete");
  delete fix[ifix];
  atom->update_callback(ifix);

  // move other Fixes and fmask down in list one slot
  
  for (int i = ifix+1; i < nfix; i++) fix[i-1] = fix[i];
  for (int i = ifix+1; i < nfix; i++) fmask[i-1] = fmask[i];
  nfix--;
}

/* ----------------------------------------------------------------------
  find a fix by ID
  return index of fix or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_fix(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) return -1;
  return ifix;
}

/* ----------------------------------------------------------------------
   one instance per compute in style_compute.h
------------------------------------------------------------------------- */

template <typename T>
Compute *Modify::compute_creator(CAC *cac, int narg, char **arg)
{
  return new T(cac,narg,arg);
}

/* ----------------------------------------------------------------------
   delete all lists of restart file Fix info
------------------------------------------------------------------------- */

void Modify::restart_deallocate()
{
  if (nfix_restart_global) {
    for (int i = 0; i < nfix_restart_global; i++) {
      delete [] id_restart_global[i];
      delete [] style_restart_global[i];
      delete [] state_restart_global[i];
    }
    delete [] id_restart_global;
    delete [] style_restart_global;
    delete [] state_restart_global;
  }

  if (nfix_restart_peratom) {
    for (int i = 0; i < nfix_restart_peratom; i++) {
      delete [] id_restart_peratom[i];
      delete [] style_restart_peratom[i];
    }
    delete [] id_restart_peratom;
    delete [] style_restart_peratom;
    delete [] index_restart_peratom;
  }

  nfix_restart_global = nfix_restart_peratom = 0;
}

/* ----------------------------------------------------------------------
  add a new compute
 ------------------------------------------------------------------------- */
void Modify::add_compute(int narg, char **arg, int trysuffix)
{
  if (narg < 3) error->all(FLERR,"Illegal compute command");

  // error check

    for (int icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(arg[0],compute[icompute]->id) == 0)
      error->all(FLERR,"Reuse of compute ID");

  // extend Compute list if necessary

  if (ncompute == maxcompute) {
    maxcompute += DELTA;
    compute = (Compute **)
      memory->srealloc(compute,maxcompute*sizeof(Compute *),"modify:compute");
  }

  // create the Compute
  // try first with suffix appended

  compute[ncompute] = NULL;

  if (compute[ncompute] == NULL && 
      compute_map->find(arg[2]) != compute_map->end()) {
    ComputeCreator compute_creator = (*compute_map)[arg[2]];
    compute[ncompute] = compute_creator(cac,narg,arg);
  }

  if (compute[ncompute] == NULL) error->all(FLERR,"Unknown compute style");

  ncompute++;
}

/* ----------------------------------------------------------------------
  add a new fix or replace one with same ID
------------------------------------------------------------------------- */

void Modify::add_fix(int narg, char **arg, int trysuffix)
{
  if (narg < 3) error->all(FLERR,"Illegal fix command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID");

 int ifix,newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break; 

  if (ifix < nfix) {
  
    newflag = 0;
    int match = 0;
    if (strcmp(arg[2],fix[ifix]->style) == 0) match = 1;
    if (!match) error->all(FLERR,"Replacing a fix, but new style != old style");

    if (fix[ifix]->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Replacing a fix, but new group != old group");

    delete fix[ifix];
    fix[ifix] = NULL;

  }else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **) memory->srealloc(fix,maxfix*sizeof(Fix *),"modify:fix");
      memory->grow(fmask,maxfix,"modify:fmask");
    }
  }

  // create the Fix
  // try first with suffix appended

  fix[ifix] = NULL;

  if (fix[ifix] == NULL && fix_map->find(arg[2]) != fix_map->end()) {
    FixCreator fix_creator = (*fix_map)[arg[2]];
    fix[ifix] = fix_creator(cac,narg,arg);
  }
  
  if (fix[ifix] == NULL) error->all(FLERR,"Unknown fix style");

  if (newflag) nfix++;
  fmask[ifix] = fix[ifix]->setmask();
  fix[ifix]->post_constructor();
}

/* ----------------------------------------------------------------------
   find a compute by ID
   return index of compute or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_compute(const char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,compute[icompute]->id) == 0) break;
  if (icompute == ncompute) return -1;
  return icompute;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
  int i,j;

  // delete storage of restart info since it is not valid after 1st run

  restart_deallocate();

  // create lists of fixes to call at each stage of run

  list_init(INITIAL_INTEGRATE,n_initial_integrate,list_initial_integrate);
  list_init(POST_INTEGRATE,n_post_integrate,list_post_integrate);
  list_init(PRE_EXCHANGE,n_pre_exchange,list_pre_exchange);
  list_init(PRE_NEIGHBOR,n_pre_neighbor,list_pre_neighbor);
  list_init(PRE_FORCE,n_pre_force,list_pre_force);
  list_init(POST_FORCE,n_post_force,list_post_force);
  list_init(FINAL_INTEGRATE,n_final_integrate,list_final_integrate);
  list_init_end_of_step(END_OF_STEP,n_end_of_step,list_end_of_step);
  list_init_thermo_energy(THERMO_ENERGY,n_thermo_energy,list_thermo_energy);

  list_init(INITIAL_INTEGRATE_RESPA,
            n_initial_integrate_respa,list_initial_integrate_respa);
  list_init(POST_INTEGRATE_RESPA,
            n_post_integrate_respa,list_post_integrate_respa);
  list_init(POST_FORCE_RESPA,
            n_post_force_respa,list_post_force_respa);
  list_init(PRE_FORCE_RESPA,
            n_pre_force_respa,list_pre_force_respa);
  list_init(FINAL_INTEGRATE_RESPA,
            n_final_integrate_respa,list_final_integrate_respa);

  list_init(MIN_PRE_EXCHANGE,n_min_pre_exchange,list_min_pre_exchange);
  list_init(MIN_PRE_NEIGHBOR,n_min_pre_neighbor,list_min_pre_neighbor);
  list_init(MIN_PRE_FORCE,n_min_pre_force,list_min_pre_force);
  list_init(MIN_POST_FORCE,n_min_post_force,list_min_post_force);
  list_init(MIN_ENERGY,n_min_energy,list_min_energy);

  // init each fix
  // not sure if now needs to come before compute init
  // used to b/c temperature computes called fix->dof() in their init,
  // and fix rigid required its own init before its dof() could be called,
  // but computes now do their DOF in setup()
  
 // if (comm->me==0) fprintf(screen, "fix number = %d\n", nfix);
  for (i = 0; i < nfix; i++) fix[i]->init();

  // set global flag if any fix has its restart_pbc flag set

  restart_pbc_any = 0;
  for (i = 0; i < nfix; i++)
    if (fix[i]->restart_pbc) restart_pbc_any = 1;

  // create list of computes that store invocation times

  list_init_compute();

  // init each compute
  // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
  // add initial timestep to all computes that store invocation times
  //   since any of them may be invoked by initial thermo
  // do not clear out invocation times stored within a compute,
  //   b/c some may be holdovers from previous run, like for ave fixes

  for (i = 0; i < ncompute; i++) {
    compute[i]->init();
    compute[i]->invoked_scalar = -1;
    compute[i]->invoked_vector = -1;
    compute[i]->invoked_array = -1;
    compute[i]->invoked_peratom = -1;
    compute[i]->invoked_local = -1;
  }

  addstep_compute_all(update->ntimestep);

  // error if any fix or compute is using a dynamic group when not allowed

  for (i = 0; i < nfix; i++)
    if (!fix[i]->dynamic_group_allow && group->dynamic[fix[i]->igroup]) {
      char str[128];
      sprintf(str,"Fix %s does not allow use of dynamic group",fix[i]->id);
      error->all(FLERR,str);
    }

  for (i = 0; i < ncompute; i++)
    if (!compute[i]->dynamic_group_allow && 
        group->dynamic[compute[i]->igroup]) {
      char str[128];
      sprintf(str,"Compute %s does not allow use of dynamic group",fix[i]->id);
      error->all(FLERR,str);
    }

 // warn if any particle is time integrated more than once
 
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int *flag = new int[nlocal];
  for (i = 0; i < nlocal; i++) flag[i] = 0;
  
  int groupbit;
  for (i = 0; i < nfix; i++) {
    if (fix[i]->time_integrate == 0) continue;
    groupbit = fix[i]->groupbit;
    for (j = 0; j < nlocal; j++)
      if (mask[j] & groupbit) flag[j]++;
  }

  int check = 0;
  for (i = 0; i < nlocal; i++)
    if (flag[i] > 1) check = 1;

  delete [] flag;

  int checkall;
  MPI_Allreduce(&check,&checkall,1,MPI_INT,MPI_SUM,world);

  if (comm->me == 0 && checkall)
    error->warning(FLERR,
                   "One or more atoms are time integrated more than once"); 
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) list[n++] = i;

}

/* ----------------------------------------------------------------------
   create list of fix indices for end_of_step fixes
   also create end_of_step_every[]
------------------------------------------------------------------------- */

void Modify::list_init_end_of_step(int mask, int &n, int *&list)
{
  delete [] list;
  delete [] end_of_step_every;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];
  end_of_step_every = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every[n++] = fix[i]->nevery;
    } 
}

/* ----------------------------------------------------------------------
   create list of fix indices for thermo energy fixes
   only added to list if fix has THERMO_ENERGY mask
   and its thermo_energy flag was set via fix_modify
------------------------------------------------------------------------- */

void Modify::list_init_thermo_energy(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask && fix[i]->thermo_energy) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask && fix[i]->thermo_energy) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of compute indices for computes which store invocation times
------------------------------------------------------------------------- */

void Modify::list_init_compute()
{
  delete [] list_timeflag;

  n_timeflag = 0;
//  if (comm->me==0) fprintf(screen, "compute number = %d\n",ncompute);

  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) n_timeflag++;
  list_timeflag = new int[n_timeflag];

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) list_timeflag[n_timeflag++] = i;

}

/* ----------------------------------------------------------------------
   loop over all computes
   schedule next invocation for those that store invocation times
   called when not sure what computes will be needed on newstep
   do not loop only over n_timeflag, since may not be set yet
------------------------------------------------------------------------- */

void Modify::addstep_compute_all(bigint newstep)
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    if (compute[icompute]->timeflag) compute[icompute]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void Modify::setup_pre_exchange()
{
//    if (comm->me==0) fprintf(screen, "n_pre_exchange = %d\n", n_pre_exchange);
//  if (update->whichflag <= 1)
//    for (int i = 0; i < n_pre_exchange; i++)
//      fix[list_pre_exchange[i]]->setup_pre_exchange();
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void Modify::setup_pre_neighbor()
{
//  if (update->whichflag == 1)
//    for (int i = 0; i < n_pre_neighbor; i++)
//      fix[list_pre_neighbor[i]]->setup_pre_neighbor();
//  else if (update->whichflag == 2)
//    for (int i = 0; i < n_min_pre_neighbor; i++)
//      fix[list_min_pre_neighbor[i]]->min_setup_pre_neighbor();
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup_pre_force(int vflag)
{
  //if (update->whichflag == 1)
    //if (comm->me==0) fprintf(screen, "n_pre_force = %d\n", n_pre_force);
    
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup(int vflag)
{
  // compute setup needs to come before fix setup
  //   b/c NH fixes need DOF of temperature computes
  // fix group setup() is special case since populates a dynamic group
  //   needs to be done before temperature compute setup

  //for (int i = 0; i < nfix; i++)
   // if (strcmp(fix[i]->style,"GROUP") == 0)
     // if (comm->me==0) fprintf(screen, "testing\n");
 
  if (update->whichflag == 1)
    for (int i = 0; i < nfix; i++) fix[i]->setup(vflag);
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate(int vflag)
{
  //if (comm->me==0) fprintf(screen, "vflag = %d, n_initial_integrate = %d\n", vflag, n_initial_integrate);
  for (int i = 0; i < n_initial_integrate; i++) {
    if (atom->nlocal) fix[list_initial_integrate[i]]->initial_integrate(vflag);
    if (cac->element_flag) fix[list_initial_integrate[i]]->einitial_integrate(vflag);
  }
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate()
{
  for (int i = 0; i < n_final_integrate; i++) {
   if (atom->nlocal) fix[list_final_integrate[i]]->final_integrate();
   if (cac->element_flag) fix[list_final_integrate[i]]->efinal_integrate(); 
  }
}

bigint Modify::memory_usage()
{
  bigint bytes = 0;
  for (int i = 0; i < nfix; i++)
    bytes += static_cast<bigint> (fix[i]->memory_usage());
  for (int i = 0; i < ncompute; i++)
    bytes += static_cast<bigint> (compute[i]->memory_usage());
  return bytes;
}

/* ----------------------------------------------------------------------
   loop over computes that store invocation times
   if its invoked flag set on this timestep, schedule next invocation
   called everywhere that computes are used, after computes are invoked
------------------------------------------------------------------------- */

void Modify::addstep_compute(bigint newstep)
{
  for (int icompute = 0; icompute < n_timeflag; icompute++)
    if (compute[list_timeflag[icompute]]->invoked_flag)
      compute[list_timeflag[icompute]]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void Modify::post_run()
{
  for (int i = 0; i < nfix; i++) fix[i]->post_run();
}

/*------------------------------------------------------------------------
  post_force call, only for relevant fixes
--------------------------------------------------------------------------*/

void Modify::post_force(int vflag)
{
  for (int i = 0; i < n_post_force;i++) {
   if (atom->nlocal) fix[list_post_force[i]]->post_force(vflag);
   if (element->nlocal) fix[list_post_force[i]]->epost_force(vflag);
  }
}
