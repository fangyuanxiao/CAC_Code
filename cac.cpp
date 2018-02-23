#include "mpi.h"
#include "string.h"
#include "ctype.h"
#include "cac.h"
#include "style_atom.h"
#include "style_command.h"
#include "style_pair.h"
#include "style_fix.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "input.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "neighbor.h"
#include "comm.h"
#include "comm_brick.h"
#include "group.h"
#include "output.h"
#include "domain.h"
#include "force.h"
#include "timer.h"
#include "stdlib.h"

using namespace CAC_NS;

CAC::CAC(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this, communicator);

  screen = NULL;
  logfile = NULL;
  infile = NULL;

  initclock = MPI_Wtime();
  
  // parse input switches 

  int inflag = 0;
  element_flag = 0;
  suffix = suffix2 = NULL;
  suffix_enable = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-in")==0 ||
        strcmp(arg[iarg],"-i")==0) {
     if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
        inflag = iarg + 1;
        iarg +=2;
     }
  }
   
  universe->uscreen = stdout;
  universe->ulogfile = fopen("log.cac","w");
  if (universe->ulogfile == NULL)
       error->universe_warn(FLERR,"Cannot open log.cac for writing");

  screen = universe->uscreen;
  logfile = universe->ulogfile;
  world = universe->uworld;
  
  if (universe->me ==0){
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
         char str[128];
         sprintf(str, "Cannot open input script %s", arg[inflag]);
         error->one(FLERR,str);
      }
   }
  
   // check consistency of datatype settings in cactype.h
   //
   if (sizeof(smallint) != sizeof(int))
     error->all(FLERR,"Smallint setting in cactype.h is invalid");
   if (sizeof(imageint) < sizeof(smallint))
     error->all(FLERR,"Imageint setting in cactype.h is invalid");
   if (sizeof(tagint) < sizeof(smallint))
     error->all(FLERR,"Tagint setting in cactype.h is invalid");
   if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
     error->all(FLERR,"Bigint setting in cactype.h is invalid");

   int mpisize;
   MPI_Type_size(MPI_LMP_TAGINT,&mpisize);
   if (mpisize != sizeof(tagint))
     error->all(FLERR,"MPI_LMP_TAGINT and tagint in "
                "cactype.h are not compatible");
   MPI_Type_size(MPI_LMP_BIGINT,&mpisize);
   if (mpisize != sizeof(bigint))
      error->all(FLERR,"MPI_LMP_BIGINT and bigint in "
                 "cactype.h are not compatible");

   #ifdef LAMMPS_SMALLBIG
     if (sizeof(smallint) != 4 || sizeof(imageint) != 4 || 
          sizeof(tagint) != 4 || sizeof(bigint) != 8)
          error->all(FLERR,"Small to big integers are not sized correctly");
   #endif
   #ifdef LAMMPS_BIGBIG
     if (sizeof(smallint) != 4 || sizeof(imageint) != 8 || 
         sizeof(tagint) != 8 || sizeof(bigint) != 8)
        error->all(FLERR,"Small to big integers are not sized correctly");
   #endif
   #ifdef LAMMPS_SMALLSMALL
     if (sizeof(smallint) != 4 || sizeof(imageint) != 4 || 
         sizeof(tagint) != 4 || sizeof(bigint) != 4)
         error->all(FLERR,"Small to big integers are not sized correctly");
   #endif

   input = new Input(this, narg,arg);
   create();
   //int nr = static_cast<int> (499.3+0.5);
   //if (comm->me==0) fprintf(screen, "nr = %d\n",nr);
   //error->all(FLERR, "test from cac");
    
}

CAC::~CAC()
{
   //error->all(FLERR, "test from ~CAC");
   destroy();
  // error->all(FLERR, "test from ~CAC");

   if (logfile) fclose(logfile);
   if (screen && screen != stdout) fclose(screen);
   logfile = NULL;
   if (screen != stdout) screen = NULL;
   //error->all(FLERR, "test from ~CAC");
   if (infile && infile != stdin) fclose(infile);
   
   if (world != universe->uworld) MPI_Comm_free(&world);
   //error->all(FLERR,"test from ~CAC");
   delete input;
   //error->all(FLERR, "test from ~CAC");
   delete universe;
   //error->all(FLERR, "test from ~CAC");
  // class Error *ee;
  // ee->all(FLERR,"test from ~CAC");
   //error->done();
   delete error;
   //MPI_Barrier(world);
   //MPI_Finalize();
   //exit(1);
   //ee->all(FLERR,"test from ~CAC");
   delete memory;
   //MPI_Barrier(world);
   //MPI_Finalize();
   //exit(1);
   //ee->all(FLERR, "test from ~CAC");
   
}

void CAC::create()
{
   //error->all(FLERR, "test from create in cac");
   comm = new CommBrick(this);
   neighbor = new Neighbor(this);
   domain = new Domain(this);
   atom = new Atom(this);
   atom->create_avec("atomic",0,NULL,1);
   
   element = new Element(this); 
   group = new Group(this);
   force = new Force(this);

  modify = new Modify(this);
  // error->all(FLERR, "test from create in cac"); 
  output = new Output(this);
  
  update = new Update(this);
  timer = new Timer(this); 
  //error->all(FLERR, "test from create in cac");
}

void CAC::destroy()
{
  delete update;
  delete neighbor;
  delete comm;
  delete force;
  delete group;
  delete output;
  delete modify;

  delete domain;
  delete atom;
  delete element;
  delete timer;

  modify = NULL;
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void CAC::init()
{
   update->init();
   force->init();
   domain->init();
   element->init(); // must before neighbor for element size
   atom->init();
   modify->init();
   neighbor->init();
   comm->init();
   output->init();
}
