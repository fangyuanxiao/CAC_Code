#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "domain.h"
#include "neighbor.h"
#include "force.h"
#include "output.h"
#include "modify.h"
#include "pair.h"
#include "dump.h"
#include "compute.h"
#include "fix.h"
#include "universe.h"
#include "procmap.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define BUFMIN 1000

enum{SINGLE,MULTI};
enum{MULTIPLE};
enum{ONELEVEL,TWOLEVEL,NUMA,CUSTOM};
enum{CART,CARTREORDER,XYZ};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};

Comm::Comm(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  mode = 0;
  bordergroup = 0;
  cutghostuser = 0.0;
  ghost_velocity = 0;

  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  coregrid[0] = coregrid[1] = coregrid[2] = 1;
  gridflag = ONELEVEL;
  mapflag = CART;
  customfile = NULL;
  outfile = NULL;
  recv_from_partition = send_to_partition = -1;
  otherflag = 0;
  maxexchange_atom = maxexchange_fix = 0;

  grid2proc = NULL;
  xsplit = ysplit = zsplit = NULL;
  rcbnew = 0;
  
  nthreads = 1;
}

Comm::~Comm()
{ 
  memory->destroy(grid2proc);
  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);
  delete [] customfile;
  delete [] outfile;
}

/* ----------------------------------------------------------------------
   create a 3d grid of procs based on Nprocs and box size & shape
   map processors to grid, setup xyz split for a uniform grid
------------------------------------------------------------------------- */

void Comm::set_proc_grid(int outflag)
{
  ProcMap *pmap = new ProcMap(cac);

  pmap->onelevel_grid(nprocs,user_procgrid,procgrid,
                        otherflag,other_style,other_procgrid,other_coregrid);

  if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs)
    error->all(FLERR,"Bad grid of processors");
  if (domain->dimension == 2 && procgrid[2] != 1)
    error->all(FLERR,"Processor count in z must be 1 for 2d simulation");

  // grid2proc[i][j][k] = proc that owns i,j,k location in 3d grid
  
  if (grid2proc) memory->destroy(grid2proc);
  memory->create(grid2proc,procgrid[0],procgrid[1],procgrid[2],
                 "comm:grid2proc");

  // map processor IDs to 3d processor grid
  // produces myloc, procneigh, grid2proc

  pmap->cart_map(0,procgrid,myloc,procneigh,grid2proc);

  if (outflag && me == 0) {
    if (screen) {
      fprintf(screen,"  %d by %d by %d MPI processor grid\n",
              procgrid[0],procgrid[1],procgrid[2]);
      if (gridflag == NUMA || gridflag == TWOLEVEL)
        fprintf(screen,"  %d by %d by %d core grid within node\n",
                coregrid[0],coregrid[1],coregrid[2]);
    }
    if (logfile) {
      fprintf(logfile,"  %d by %d by %d MPI processor grid\n",
              procgrid[0],procgrid[1],procgrid[2]);
      if (gridflag == NUMA || gridflag == TWOLEVEL)
        fprintf(logfile,"  %d by %d by %d core grid within node\n",
                coregrid[0],coregrid[1],coregrid[2]);
    }
  }

  if(outfile) pmap->output(outfile,procgrid,grid2proc);
  
  delete pmap;

  // set xsplit,ysplit,zsplit for uniform spacings

  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);

  memory->create(xsplit,procgrid[0]+1,"comm:xsplit");
  memory->create(ysplit,procgrid[1]+1,"comm:ysplit");
  memory->create(zsplit,procgrid[2]+1,"comm:zsplit");

  for (int i = 0; i < procgrid[0]; i++) xsplit[i] = i * 1.0/procgrid[0];
  for (int i = 0; i < procgrid[1]; i++) ysplit[i] = i * 1.0/procgrid[1];
  for (int i = 0; i < procgrid[2]; i++) zsplit[i] = i * 1.0/procgrid[2];

  xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;  
}

/* ----------------------------------------------------------------------
   proc 0 reads Nlines from file into buf and bcasts buf to all procs
   caller allocates buf to max size needed
   each line is terminated by newline, even if last line in file is not
   return 0 if successful, 1 if get EOF error before read is complete
------------------------------------------------------------------------- */

int Comm::read_lines_from_file(FILE *fp, int nlines, int maxline, char *buf)
{
  int m;
  
    if (me == 0) {
    m = 0;
    for (int i = 0; i < nlines; i++) {
      if (!fgets(&buf[m],maxline,fp)) {
	m = 0;
	break;
      }
      m += strlen(&buf[m]);
    }
    if (m) {
      if (buf[m-1] != '\n') strcpy(&buf[m++],"\n");
      m++;
    }
  }

  MPI_Bcast(&m,1,MPI_INT,0,world);
  if (m == 0) return 1;
  MPI_Bcast(buf,m,MPI_CHAR,0,world);
  return 0;

}

/* ----------------------------------------------------------------------
   common to all Comm styles
------------------------------------------------------------------------- */

void Comm::init()
{
  triclinic = domain->triclinic;
  map_style = atom->map_style;

  // check warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in exchange()
  // really should check every exchange() in case box size is shrinking
  //   but seems overkill to do that (fix balance does perform this check)
  // if (me==0) fprintf(screen, "skin = %f\n", neighbor->skin);  

  domain->subbox_too_small_check(neighbor->skin);

  // comm_only = 1 if only x,f are exchanged in forward/reverse comm
  // comm_x_only = 0 if ghost_velocity since velocities are added

  comm_x_only = atom->avec->comm_x_only;
  comm_f_only = atom->avec->comm_f_only;
  if (ghost_velocity) comm_x_only = 0;

  // set per-atom sizes for forward/reverse/border comm
  // augment by velocity and fix quantities if needed

  size_forward = atom->avec->size_forward;
  size_reverse = atom->avec->size_reverse;
  size_border = atom->avec->size_border;

  esize_forward = element->size_forward;
  esize_reverse = element->size_reverse;
  esize_border = element->size_border;  
 
  if (ghost_velocity) size_forward += atom->avec->size_velocity;
  if (ghost_velocity) size_border += atom->avec->size_velocity;

  for (int i = 0; i < modify->nfix; i++)
    size_border += modify->fix[i]->comm_border; 

  // per-atom limits for communication
  // maxexchange = max # of datums in exchange comm, set in exchange()
  // maxforward = # of datums in largest forward comm
  // maxreverse = # of datums in largest reverse comm
  // query pair,fix,compute,dump for their requirements
  // pair style can force reverse comm even if newton off

  maxforward = MAX(size_forward,size_border);
  maxreverse = size_reverse;
  emaxforward = MAX(esize_forward,esize_border);
  emaxreverse = esize_reverse;

  if (force->pair) maxforward = MAX(maxforward,force->pair->comm_forward);
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse);

  for (int i = 0; i < modify->nfix; i++) {
    maxforward = MAX(maxforward,modify->fix[i]->comm_forward);
    maxreverse = MAX(maxreverse,modify->fix[i]->comm_reverse);
  }

  for (int i = 0; i < modify->ncompute; i++) {
    maxforward = MAX(maxforward,modify->compute[i]->comm_forward);
    maxreverse = MAX(maxreverse,modify->compute[i]->comm_reverse);
  }

  for (int i = 0; i < output->ndump; i++) {
    maxforward = MAX(maxforward,output->dump[i]->comm_forward);
    maxreverse = MAX(maxreverse,output->dump[i]->comm_reverse);
  }

  if (force->newton == 0) maxreverse = 0;
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse_off);

  //if (me==0) fprintf(screen, "maxforward = %d, maxreverse = %d\n", maxforward, maxreverse);
  //if (me==0) fprintf(screen, "size_forward = %d, size_reverse = %d\n", size_forward, size_reverse);
}

/*----------------------------------------------------------------------------------------------
   set dimensions for 3d grid of processors, and associated flags
   invoked from input script by processors command
-----------------------------------------------------------------------------------------------*/

void Comm::set_processors(int narg, char **arg)
{
   int iarg=0;
   while (iarg < narg) {
     if (strcmp(arg[iarg],"file")==0) {
       if (iarg+2 > narg) error->all(FLERR, "Illegal processors command");
       delete [] outfile;
       int n = strlen(arg[iarg+1]) + 1;
       outfile = new char[n];
       strcpy(outfile, arg[iarg+1]);
       iarg +=2;
     } else error->all(FLERR,"Illegal processors command");
   }
}
