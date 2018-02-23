#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "finish.h"
#include "cac.h"
#include "universe.h"
//#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_vec.h"
//#include "molecule.h"
#include "comm.h"
#include "force.h"
//#include "kspace.h"
#include "update.h"
//#include "min.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "timer.h"
#include "output.h"
#include "memory.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

Finish::Finish(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void Finish::end(int flag)
{
  int i,m,nneigh,nneighfull;
  int histo[10];
  int loopflag,minflag,prdflag,tadflag,timeflag,fftflag,histoflag,neighflag;
  double time,tmp,ave,max,min;
  double time_loop,time_other;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs); 

  bigint nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&element->nelements,1,MPI_LMP_BIGINT,MPI_SUM,world);

  loopflag = 1;
  minflag = prdflag = tadflag = timeflag = fftflag = histoflag = neighflag = 0;

  if (flag == 1) {
    timeflag = histoflag = 1;
    neighflag = 1;
    if (update->whichflag == 1 &&
        strcmp(update->integrate_style,"verlet/split") == 0 &&
        universe->iworld == 1) neighflag = 0;
    //if (me==0) fprintf(screen, "neighflag = %d\n",neighflag);
  }

  if (loopflag) {
    time_other = timer->array[TIME_LOOP] -
      (timer->array[TIME_PAIR] + timer->array[TIME_BOND] +
       timer->array[TIME_KSPACE] + timer->array[TIME_NEIGHBOR] +
       timer->array[TIME_COMM] + timer->array[TIME_OUTPUT]);

    time_loop = timer->array[TIME_LOOP];
    MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time_loop = tmp/nprocs;

    if (me == 0) {
      if (screen) fprintf(screen,
                          "Loop time of %g on %d procs for %d steps with "
                          BIGINT_FORMAT " elements\n",
                          time_loop,nprocs,update->nsteps,element->nelements);
      if (logfile) fprintf(logfile,
                           "Loop time of %g on %d procs for %d steps with "
                           BIGINT_FORMAT " elements\n",
                           time_loop,nprocs,update->nsteps,element->nelements);
    }

   if (time_loop == 0.0) time_loop = 1.0;
  } 

  if (timeflag) {
    if (me==0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    time = timer->array[TIME_PAIR];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
	fprintf(screen,"Pair time (%%) = %g (%g)\n",
		time, time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Pair time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }

    time = timer->array[TIME_NEIGHBOR];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me==0) {
      if (screen)
	fprintf(screen,"Neigh time (%%) = %g (%g)\n",
		 time, time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Neigh time (%%) = %g (%g)\n",
	        time, time/time_loop*100.0);
    }

    time = timer->array[TIME_COMM];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me==0) {
      if (screen)
	 fprintf(screen,"Comm time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Comm time (%%) = %g (%g)\n",
		 time,time/time_loop*100.0);
    }

    time = timer->array[TIME_OUTPUT];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me==0) {
      if (screen)
	fprintf(screen,"Output time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Output time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me==0) {
      if (screen)
	fprintf(screen,"Other time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Other time (%%) = %g (%g)\n",
	         time, time/time_loop*100.0);
    }
  }

  if (histoflag) {
    if (me==0) {
       if (screen) fprintf(screen,"\n");
       if (logfile) fprintf(logfile,"\n");
    }

    tmp = element->nlocal;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me==0) {
      if (me==0) {
	if (screen) {
          fprintf(screen,"Nlocal:   %g ave %g max %g min\n",ave,max,min);
	  fprintf(screen, "Histogram:");
	  for (i = 0; i < 10; i++) fprintf(screen, " %d",histo[i]);
	  fprintf(screen,"\n");
        }
	if (logfile) {
	  fprintf(logfile,"Nlocal:    %g ave %g max %g min\n",ave, max, min);
	  fprintf(logfile, "Histogram:");
	  for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
	  fprintf(logfile,"\n");
	}
      }
    }

    tmp = element->nghost;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
	fprintf(screen,"Nghost:    %g ave %g max %g min\n",ave,max,min);
	fprintf(screen,"Histogram:");
	for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
	fprintf(screen,"\n");
      }
      if (logfile) {
	fprintf(logfile,"Nghost:    %g ave %g max %g min\n",ave,max,min);
	fprintf(logfile,"Histogram:");
	for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
	fprintf(logfile,"\n");
      }
    }
  }

  if (logfile) fflush(logfile);

}

/*---------------------------------------------------------------------------------*/

void Finish::stats(int n, double *data,
		   double *pave, double *pmax, double *pmin,
		   int nhisto,int *histo)
{
  int i,m;
  int *histotmp;
  
  double min = 1.0e20;
  double max = -1.0e20;
  double ave = 0.0;
  for (i = 0; i < n; i++) {
    ave += data[i];
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  int ntotal;
  MPI_Allreduce(&n,&ntotal,1,MPI_INT,MPI_SUM,world);
  double tmp;
  MPI_Allreduce(&ave,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  ave = tmp/ntotal;
  MPI_Allreduce(&min,&tmp,1,MPI_DOUBLE,MPI_MIN,world);
  min = tmp;
  MPI_Allreduce(&max,&tmp,1,MPI_DOUBLE,MPI_MAX,world);
  max = tmp;

  for (i = 0; i < nhisto;i++) histo[i] = 0;

  double del = max - min;
  for (i = 0; i < n; i++) {
    if (del == 0.0) m = 0;
    else m = static_cast<int> ((data[i]-min)/del * nhisto);
    if (m > nhisto-1) m = nhisto-1;
    histo[m]++;
  }

  memory->create(histotmp,nhisto,"finish:histotmp");
  MPI_Allreduce(histo,histotmp,nhisto,MPI_INT,MPI_SUM,world);
  for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
  memory->destroy(histotmp);

  *pave = ave;
  *pmax = max;
  *pmin = min;
}
