#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "style_dump.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "input.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "thermo.h"
#include "modify.h"
#include "force.h"
#include "dump.h"
#include "update.h"
//#include "write_restart.h"
#include "universe.h"
#include "memory.h"
#include "error.h"
#include "element.h"
using namespace CAC_NS;

#define DELTA 1

/* ----------------------------------------------------------------------
   initialize all output
------------------------------------------------------------------------- */

Output::Output(CAC *cac) : Pointers(cac)
{
   // create default computes for temp,pressure,pe

  char **newarg = new char*[4];
//  newarg[0] = (char *) "thermo_temp";
//  newarg[1] = (char *) "all";
//  newarg[2] = (char *) "temp";
//  modify->add_compute(3,newarg,1);

//  newarg[0] = (char *) "thermo_press";
//  newarg[1] = (char *) "all";
//  newarg[2] = (char *) "pressure";
//  newarg[3] = (char *) "thermo_temp";
//  modify->add_compute(4,newarg,1);

//  newarg[0] = (char *) "thermo_pe";
//  newarg[1] = (char *) "all";
//  newarg[2] = (char *) "pe";
//  modify->add_compute(3,newarg,1);

  delete [] newarg;

  // create default Thermo class

//  newarg = new char*[1];
//  newarg[0] = (char *) "one";
//  thermo = new Thermo(cac,1,newarg);
 // delete [] newarg;
  thermo = NULL;
  thermo_every = 0;
  var_thermo = NULL;

  ndump = 0;
  max_dump = 0;
  every_dump = NULL;
  next_dump = NULL;
  last_dump = NULL;
  var_dump = NULL;
  ivar_dump = NULL;
  dump = NULL;

  restart_flag = restart_flag_single = restart_flag_double = 0;
//  restart_every_single = restart_every_double = 0;
//  last_restart = -1;
//  restart1 = restart2a = restart2b = NULL;
//  var_restart_single = var_restart_double = NULL;
//  restart = NULL;
  nmax = 0;
  centro = NULL;
  distsq = NULL;
  nearest = NULL;
  maxneigh = 0;
  nnn = 12;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

Output::~Output()
{
  if (thermo) delete thermo;
  delete [] var_thermo;

  memory->destroy(every_dump);
  memory->destroy(next_dump);
  memory->destroy(last_dump);
  memory->destroy(centro);
  memory->destroy(distsq);
  memory->destroy(nearest);
  for (int i = 0; i < ndump; i++) delete [] var_dump[i];
  memory->sfree(var_dump);
  memory->destroy(ivar_dump);
  for (int i = 0; i < ndump; i++) delete dump[i];
  memory->sfree(dump);

//  delete [] restart1;
//  delete [] restart2a;
//  delete [] restart2b;
//  delete [] var_restart_single;
//  delete [] var_restart_double;
//  delete restart;
}

/* ----------------------------------------------------------------------
   set thermo output frequency from input script
------------------------------------------------------------------------- */

void Output::set_thermo(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal thermo command");

  thermo_every = universe->inumeric(FLERR,arg[0]);
  if (thermo_every < 0) error->all(FLERR,"Illegal thermo command");
}

/* ----------------------------------------------------------------------
   new Thermo style
------------------------------------------------------------------------- */

void Output::create_thermo(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal thermo_style command");
  
  if (domain->box_exist == 0)
    error->all(FLERR,"Thermo_style command before simulation box is defined");

  if (thermo->modified && comm->me == 0)
    error->warning(FLERR,"New thermo_style command, "
                   "previous thermo_modify settings will be lost");

  delete thermo;
  thermo = NULL;
  thermo = new Thermo(cac,narg,arg);
}

/* ----------------------------------------------------------------------
   add a Dump to list of Dumps
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump command");

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0)
      error->all(FLERR,"Reuse of dump ID");
  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find dump group ID");
  if (universe->inumeric(FLERR,arg[3]) <= 0) 
    error->all(FLERR,"Invalid dump frequency");

  if (ndump == max_dump) {
    max_dump += DELTA;
    dump = (Dump **)
      memory->srealloc(dump,max_dump*sizeof(Dump *),"output:dump");
    memory->grow(every_dump,max_dump,"output:every_dump");
    memory->grow(next_dump,max_dump,"output:next_dump");
    memory->grow(last_dump,max_dump,"output:last_dump");
    var_dump = (char **)
      memory->srealloc(var_dump,max_dump*sizeof(char *),"output:var_dump");
    memory->grow(ivar_dump,max_dump,"output:ivar_dump");
  }

  if (0) return;

  #define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) dump[ndump] = new Class(cac,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Unknown dump style");

  every_dump[ndump] = universe->inumeric(FLERR,arg[3]);
  if (every_dump[ndump] <= 0) error->all(FLERR,"Illegal dump command");
  last_dump[ndump] = -1;
  var_dump[ndump] = NULL;
  ndump++;
}

/* ----------------------------------------------------------------------
   timestep is being changed, called by update->reset_timestep()
   reset next timestep values for dumps, restart, thermo output
   reset to smallest value >= new timestep
   if next timestep set by variable evaluation,
     eval for ntimestep-1, so current ntimestep can be returned if needed
     no guarantee that variable can be evaluated for ntimestep-1
       if it depends on computes, but live with that rare case for now
------------------------------------------------------------------------- */

//void Output::reset_timestep(bigint ntimestep)
//{
// next_dump_any = MAXBIGINT;
// for (int idump = 0; idump < ndump; idump++) {
//  if (every_dump[idump]) {
//    next_dump[idump] = (ntimestep/every_dump[idump])*every_dump[idump];
//  }
// } 

//}

/* ---------------------------------------------------------------------- */

void Output::init()
{
  //thermo->init();
  for (int i = 0; i < ndump; i++) dump[i]->init();  
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do thermo last, so will print after memory_usage
   memflag = 0/1 for printing out memory usage
------------------------------------------------------------------------- */

void Output::setup(int memflag)
{
  bigint ntimestep = update->ntimestep;
  //if (comm->me==0) fprintf(screen, "testing ntimestep = %d\n", ntimestep);
  // perform dump at start of run only if:
  //   current timestep is multiple of every and last dump not >= this step
  //   this is first run after dump created and firstflag is set
  //   note that variable freq will not write unless triggered by firstflag
  // set next_dump to multiple of every or variable value
  // set next_dump_any to smallest next_dump
  // wrap dumps that invoke computes and variable eval with clear/add
  // if dump not written now, use addstep_compute_all() since don't know
  //   what computes the dump write would invoke
  // if no dumps, set next_dump_any to last+1 so will not influence next

  int writeflag;

   

  if (ndump && update->restrict_output == 0) {
    //if (comm->me==0) fprintf(screen, "testing in output set up\n"); 
    for (int idump = 0; idump < ndump; idump++) {
      writeflag = 0;
      if (every_dump[idump] && ntimestep % every_dump[idump] == 0 &&
          last_dump[idump] != ntimestep) writeflag = 1;
      if (last_dump[idump] < 0 && dump[idump]->first_flag == 1) writeflag = 1;
      //if (comm->me==0) fprintf(screen, "output_writeflag = %d\n", writeflag);  
      if (writeflag) {
        dump[idump]->write();
//	error->all(FLERR, "test from output setup");
        last_dump[idump] = ntimestep;
      }
      if (every_dump[idump])
        next_dump[idump] =
          (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      //if (comm->me==0) fprintf(screen, "nex_dump = %d\n", next_dump[idump]);
        //if (dump[idump]->clearstep || every_dump[idump] == 0) {
        //   if (comm->me==0) fprintf(screen, "testing in output setup\n");
        //}
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  } else next_dump_any = update->laststep + 1;

  //if (comm->me==0) fprintf(screen, "next_dump_any = %d\n", next_dump_any); 
  
  if (memflag) memory_usage();

  if (thermo!=NULL) { 
     thermo->header(); 
     thermo->compute(0);
     last_thermo = ntimestep;
  } 

  //if (var_thermo) {
  //  if (comm->me==0) fprintf(screen, "testing\n");
  //}

  //if (thermo_every) {
  //  next_thermo = (ntimestep/thermo_every)*thermo_every + thermo_every;
  //  next_thermo = MIN(next_thermo,update->laststep);
  //} else next_thermo = update->laststep;

  //modify->addstep_compute(next_thermo);

  //next = MIN(next_dump_any,next_thermo);
  next = next_dump_any;
   
 // error->all(FLERR, "test from output setup");
  
}

/* ----------------------------------------------------------------------
   sum and print memory usage
   result is only memory on proc 0, not averaged across procs
------------------------------------------------------------------------- */

void Output::memory_usage()
{
  bigint bytes = 0;
  bytes += atom->memory_usage();
  bytes += neighbor->memory_usage();
  bytes += comm->memory_usage();
  bytes += update->memory_usage();
  bytes += force->memory_usage();
  bytes += modify->memory_usage(); 
  bytes += element->memory_usage();
  double mbytes = bytes/1024.0/1024.0;
  for (int i = 0; i < ndump; i++) bytes += dump[i]->memory_usage();

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Memory usage per processor = %g Mbytes\n",mbytes);
    if (logfile)
      fprintf(logfile,"Memory usage per processor = %g Mbytes\n",mbytes);
  } 
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last output doesn't
   do dump/restart before thermo so thermo CPU time will include them
------------------------------------------------------------------------- */

void Output::write(bigint ntimestep)
{
  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep) {
        if (last_dump[idump] != ntimestep) {
          dump[idump]->write();
          last_dump[idump] = ntimestep;
        }
        if (every_dump[idump]) next_dump[idump] += every_dump[idump];
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }
  
  next = next_dump_any;
}

void Output::centro_atom(NeighList *list)
{
  int i,j,k,ii,jj,kk,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,value;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (atom->nlocal > nmax) {
    memory->destroy(centro);
    nmax = atom->nmax;
    memory->create(centro,nmax,"centro/atom:centro");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int nhalf = nnn/2;
  int npairs = nnn*(nnn-1)/2;
  double *pairs = new double[npairs];

  double **x = atom->x;
  double cutsq = neighbor->cutneighmaxsq;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
   
    if (jnum > maxneigh) {
      memory->destroy(distsq);
      memory->destroy(nearest);
      maxneigh = jnum;
      memory->create(distsq,maxneigh,"centro/atom:distsq");
      memory->create(nearest,maxneigh,"centro/atom:nearest");
    } 

    // loop over list of all neighbors within force cutoff
    // distsq[] = distance sq to each
    // nearest[] = atom indices of neighbors

    n = 0;
    for (jj = 0; jj < jnum; jj++) {
       j = jlist[jj];

       delx = xtmp - x[j][0];
       dely = ytmp - x[j][1];
       delz = ztmp - x[j][2];
       rsq = delx*delx + dely*dely + delz*delz;
       if (rsq < cutsq) {
	 distsq[n] = rsq;
	 nearest[n++] = j;
       } 
    }

    // if not nnn neighbors, centro = 0.0

    if (n < nnn) {
      centro[i] = 0.0;
      continue;
    }

    // store nnn nearest neighs in 1st nnn locations of distsq and nearest

    select2(nnn,n,distsq,nearest);

    // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
    // pairs = squared length of each R

    n = 0;
    for (j = 0; j < nnn; j++) {
      jj = nearest[j];
      for (k = j+1; k < nnn; k++) {
	kk = nearest[k];
	delx = x[jj][0] + x[kk][0] - 2.0*xtmp;
	dely = x[jj][1] + x[kk][1] - 2.0*ytmp;
	delz = x[jj][2] + x[kk][2] - 2.0*ztmp;
	pairs[n++] = delx*delx + dely*dely + delz*delz;
      }
    }

    // store nhalf smallest pair distances in 1st nhalf locations of pairs

    select(nhalf,npairs,pairs);

    // centrosymmetry = sum of nhalf smallest squared values

    value = 0.0;
    for (j=0; j < nhalf; j++) value += pairs[j];
    centro[i] = value;
  }

}

/*----------------------------------------------------------------------------------------
  2 select routines from Numerical Recipes (slightly modified)
  find k smallest values in array of length n
  2nd routine sorts auxillary array at same time
-------------------------------------------------------------------------------------------*/

#define SWAP(a,b) tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

void Output::select(int k, int n, double *arr)
{
  int i,ir,j,l,mid;
  double a,tmp;

  arr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
      }
      return;
     } else {
       mid=(l+ir) >> 1;
       SWAP(arr[mid],arr[l+1])
       if (arr[l] > arr[ir]) {
	 SWAP(arr[l],arr[ir])
       }
       if (arr[l+1] > arr[ir]) {
         SWAP(arr[l+1],arr[ir])
       }
       if (arr[l] > arr[l+1]) {
	 SWAP(arr[l],arr[l+1])
       }
       i = l+1;
       j = ir;
       a = arr[l+1];
       for (;;) {
         do i++; while (arr[i] < a);
	 do j--; while (arr[j] > a);
	 if (j < i) break;
	 SWAP(arr[i],arr[j])
       }
       arr[l+1] = arr[j];
       arr[j] = a;
       if (j >= k) ir = j-1;
       if (j <= k) l = i;
    }
  }
}

/*----------------------------------------------------------------------------------*/

void Output::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      ISWAP(iarr[mid],iarr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	ISWAP(iarr[l],iarr[ir])
       }
       if (arr[l+1] > arr[ir]) {
	 SWAP(arr[l+1],arr[ir])
         ISWAP(iarr[l+1],iarr[ir])
       }
       if (arr[l] > arr[l+1]) {
	 SWAP(arr[l],arr[l+1])
         ISWAP(iarr[l],iarr[l+1])
       }
       i = l+1;
       j = ir;
       a = arr[l+1];
       ia = iarr[l+1];
       for (;;) {
	 do i++; while (arr[i] < a);
	 do j--; while (arr[j] > a);
	 if (j < i) break;
	 SWAP(arr[i],arr[j])
	 ISWAP(iarr[i],iarr[j])
       }
       arr[l+1] = arr[j];
       arr[j] = a;
       iarr[l+1] = iarr[j];
       iarr[j] = ia;
       if (j >= k) ir = j-1;
       if (j <= k) l = i;
    }
  }
}
