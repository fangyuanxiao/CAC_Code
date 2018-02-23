#include "string.h"
#include "verlet.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
//#include "bond.h"
//#include "angle.h"
//#include "dihedral.h"
//#include "improper.h"
//#include "kspace.h"
#include "output.h"
#include "dump.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

Verlet::Verlet(CAC *cac, int narg, char **arg) :
  Integrate(cac, narg, arg) {}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void Verlet::init()
{
  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR,"No fixes defined, atoms won't move");

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_fdotr_compute via sum over ghosts

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // detect if fix omp is present for clearing force arrays


  // set flags for arrays to clear in force_clear()
  

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void Verlet::setup()
{
  if (comm->me == 0 && screen) {
    fprintf(screen,"Setting up Verlet run ...\n");
    fprintf(screen,"  Unit style  : %s\n", update->unit_style);
    fprintf(screen,"  Current step: " BIGINT_FORMAT "\n", update->ntimestep);
    fprintf(screen,"  Time step   : %g\n", update->dt);
  }

  update->setupflag = 1;

  comm->setup();

  if (neighbor->style) {
    neighbor->setup_bins();
  }
 
  comm->exchange();
  comm->borders();

  // interpolation of the atoms
  element->Interpolate_to_atom_vec();
  cac->element_flag = 0; 
  neighbor->build();
  output->centro_atom(neighbor->lists[0]);
  output->dump[0]->write(); 
  MPI_Barrier(world);
  error->all(FLERR,"using error to stop here");
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void Verlet::force_clear()
{
  int i;
  size_t nbytes;

  if (external_force_clear) return;
  
  // clear force on all particles
  // if either newton flag is set, also include ghosts
  // when using threads always clear all forces.

  int nlocal = atom->nlocal;

  if (neighbor->includegroup == 0) {
    nbytes = sizeof(double) * nlocal;
    if (nbytes) {
      memset(&atom->f[0][0],0,3*nbytes);
    }
  }
}

/*-----------------------------------------------------------------------
  clear force on own & ghost nodes
  clear other arrays as needed
-------------------------------------------------------------------------*/

void Verlet::node_force_clear()
{
  int i;
  size_t nbytes;

  int nlocal = element->nlocal;
  int npe = element->npe;

  nbytes = sizeof(double)*nlocal*npe;
  if (nbytes) {
    memset(&element->nodef[0][0][0],0,3*nbytes);
  }

}
/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void Verlet::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;
  for (int i = 0; i < n; i++) {
    
    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);
    
    // initial time integration
    modify->initial_integrate(vflag);    

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      if (atom->nlocal) comm->forward_comm();
      if (cac->element_flag) comm->eforward_comm();
      timer->stamp(TIME_COMM);
    } else {
      if (atom->nlocal) domain->pbc();
      if (cac->element_flag) domain->elepbc();
      timer->stamp();
      if (cac->element_flag) element->element_size();
      neighbor->setup_cut();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
      comm->exchange();
      comm->borders();
      if (atom->nlocal) neighbor->build();
      if (cac->element_flag) neighbor->ebuild(); 
      if (cac->element_flag) neighbor->nbuild();
      timer->stamp(TIME_NEIGHBOR);
    }

    if (cac->element_flag) node_force_clear();

    timer->stamp();
   
    if (pair_compute_flag) {
      if (cac->element_flag) force->pair->ecompute();
      timer->stamp(TIME_PAIR);
    }
 
    if (n_post_force) modify->post_force(vflag);
    modify->final_integrate();
    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }

    if (ntimestep%100==0) {
      if (comm->me==0) fprintf(screen, "%d\n", ntimestep);
      if (comm->me==0) fprintf(logfile, "%d\n", ntimestep);
    }

  }
}

/* ---------------------------------------------------------------------- */

void Verlet::cleanup()
{
  modify->post_run();
  domain->box_too_small_check();
  update->update_time();
}
