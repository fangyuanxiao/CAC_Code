/*--------------------------------------------------------------------------
 *  Simple code for testing the CAC code from lammps. 
 *  Want to move the function from lammps into CAC
 *  to make CAC more effcient     Hao Chen 01-12-2016
---------------------------------------------------------------------------- */  

#include "mpi.h"
#include "cac.h"
#include "input.h"
#include "string.h"
#include "error.h"

using namespace CAC_NS;

int main(int argc, char **argv)
{
   MPI_Init(&argc,&argv);

   CAC *my_cac = new CAC(argc,argv,MPI_COMM_WORLD);
   //class Error *ee;
   my_cac->input->file();
   //my_cac->error->all(FLERR, "test from main");
   delete my_cac;

   //ee->all(FLERR, "test from main");
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize(); 
}
