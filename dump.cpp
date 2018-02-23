#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "dump.h"
#include "atom.h"
#include "element.h"
//#include "irregular.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "universe.h"

using namespace CAC_NS;

Dump *Dump::dumpptr;

#define BIG 1.0e20
#define EPSILON 1.0e-6

enum{ASCEND,DESCEND};

/* ---------------------------------------------------------------------- */

Dump::Dump(CAC *cac, int narg, char **arg) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  n = strlen(arg[4]) + 1;
  filename = new char[n];
  strcpy(filename,arg[4]);

  comm_forward = comm_reverse = 0;

  first_flag = 0;
  flush_flag = 1;
  format = NULL;
  eformat = NULL;
  format_user = NULL;
  format_default = NULL;
  clearstep = 0;
  sort_flag = 0;
  append_flag = 0;
  buffer_allow = 0;
  buffer_flag = 0;
  padflag = 0;
  nstep = 0;

  emaxbuf = maxbuf = maxids = maxsort = maxproc = 0;
  buf = bufsort = NULL;
  ids = idsort = NULL;
  index = proclist = NULL;
//  irregular = NULL;

  maxsbuf = 0;
  emaxsbuf = 0;
  sbuf = NULL;

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //  if ends in .bin = binary file
  //  else if ends in .gz = gzipped text file
  //  else ASCII text file

  fp = NULL;
  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;

  multiproc = 0;
  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;
  multiname = NULL;

  char *ptr;
  if ((ptr = strchr(filename,'%'))) {
    if (strstr(style,"mpiio"))
      error->all(FLERR,
                 "Dump file MPI-IO output not allowed with % in filename");
    multiproc = 1;
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    MPI_Comm_split(world,me,0,&clustercomm);
    multiname = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",filename,me,ptr+1);
    *ptr = '%';
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;

}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;
  delete [] multiname;

  delete [] format;
  delete [] eformat;
  delete [] format_default;
  delete [] format_user;

  memory->destroy(buf);
  memory->destroy(bufsort);
  memory->destroy(ids);
  memory->destroy(idsort);
  memory->destroy(index);
  memory->destroy(proclist);
//  delete irregular;
  
  memory->destroy(sbuf);

  if (multiproc) MPI_Comm_free(&clustercomm);

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  init_style();
  
  if (!sort_flag) {
    //if (me==0) fprintf(screen, "dump init sort_flag testing\n");
    memory->destroy(bufsort);
    memory->destroy(ids);
    memory->destroy(idsort);
    memory->destroy(index);
    memory->destroy(proclist);
   
    maxids = maxsort = maxproc = 0;
    bufsort = NULL;
    ids = idsort = NULL;
    index = proclist = NULL; 
  } 
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return
  
  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  char *filecurrent = filename;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    sprintf(filecurrent, "%s" BIGINT_FORMAT "%s",
		    filestar,update->ntimestep, ptr+1);
    *ptr = '*';

  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    fp = fopen(filecurrent,"w");
    
    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  } else fp = NULL;
  
  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void Dump::write()
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();
  
  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
  
  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_LMP_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader);
 
  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

 pack(NULL);

 // if buffering, convert doubles into strings
 // insure sbuf is sized for communicating
 // cannot buffer if output is to binary file
 
 if (buffer_flag && !binary) {
   nsme = convert_string(nme,buf);
   int nsmin, nsmax;
   MPI_Allreduce(&nsme,&nsmin,1,MPI_INT,MPI_MIN,world);
   if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
   if (multiproc != nprocs)
     MPI_Allreduce(&nsme,&nsmax,1,MPI_INT,MPI_MAX,world);
   else nsmax = nsme;

   if (nsmax>maxsbuf) {
     maxsbuf = nsmax;
     memory->grow(sbuf,maxsbuf,"dump:sbuf");
   }
 }


  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines,nchars;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles
  
  if (buffer_flag == 0 || binary) {
    //if (me==0) fprintf(screen, "testing\n");
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
        } else nlines = nme;
        
        write_data(nlines,buf);
      } 
    } else {
       MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
       MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
    }
  } else {
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs;iproc++) {
	if (iproc) {
          MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
	  MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
	  MPI_Wait(&request,&status);
	  MPI_Get_count(&status, MPI_CHAR,&nchars);
        } else nchars = nsme;

	write_data(nchars,(double *) sbuf);
      } 
      if (flush_flag) fflush(fp);
    } else {
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(sbuf,nsme,MPI_CHAR,fileproc,0,world);
    }
  }
  // since the tag changes and even the tag don't reprensent the real connectivity
  // so here I just change to output the connectivity one by one

  // if file per timestep,close file if I am filewriter
  nstep++;

  if(filewriter) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int Dump::count()
{
  int nn = 0;
  int i;
  int nlocal = atom->nlocal;
  double *centro = output->centro;

  if (atom->nlocal) { 
    for (i=0;i<nlocal;i++) {
      if (centro[i]>1 && centro[i]<10) nn++;
    }
  }
  return nn;
}

bigint Dump::memory_usage()
{
  bigint bytes = memory->usage(buf,size_one*maxbuf);
  bytes += memory->usage(sbuf,maxsbuf);
  return bytes;
}
