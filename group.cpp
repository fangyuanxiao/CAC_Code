#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "input.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace CAC_NS;

#define MAX_GROUP 32

enum{TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};

#define BIG 1.0e20

Group *Group::cptr;

Group::Group(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);

  names = new char*[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];

  for (int i = 0; i < MAX_GROUP; i++) names[i] = NULL;
  for (int i = 0; i < MAX_GROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++) inversemask[i] = bitmask[i] ^ ~0;
  for (int i = 0; i < MAX_GROUP; i++) dynamic[i] = 0;

  // create "all" group

  char *str = (char *) "all";
  int n = strlen(str) + 1;
  names[0] = new char[n];
  strcpy(names[0],str);
  ngroup = 1;
}

/* ----------------------------------------------------------------------
 free all memory
------------------------------------------------------------------------- */

Group::~Group()
{
  for (int i = 0; i < MAX_GROUP; i++) delete [] names[i];
  delete [] names;
  delete [] bitmask;
  delete [] inversemask;
  delete [] dynamic;
}

/* ----------------------------------------------------------------------
  return group index if name matches existing group, -1 if no such group
------------------------------------------------------------------------- */

int Group::find(const char *name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && strcmp(name,names[igroup]) == 0) return igroup;
  return -1;
}

/*------------------------------------------------------------------------
  assign atoms to a new or existing group
--------------------------------------------------------------------------*/

void Group::assign(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR, "Group command before simulation box is defined");
  if (narg < 2) error->all(FLERR, "Illegal group command");

  int igroup = find(arg[0]);
  if (igroup == -1) {
    if (ngroup == MAX_GROUP) error->all(FLERR, "Too many groups");
    igroup = find_unused(); 
    int n = strlen(arg[0]) + 1;
    names[igroup] = new char[n];
    strcpy(names[igroup],arg[0]);
    ngroup++;
  }
 // if (comm->me==0) fprintf(screen, "igroup = %d\n",igroup);
 // error->all(FLERR, "test from assign in group");
 
 double **x = atom->x;
 double **ex = element->x;
 int *mask = atom->mask;
 int *emask = element->mask;
 int nlocal = atom->nlocal;
 int enlocal = element->nlocal;
 int bit = bitmask[igroup];

 // style = region
 // add to group if atom and element center is in region
 
 if (strcmp(arg[1],"region") == 0) {
   
   if (narg !=3) error->all(FLERR, "Illegal group command");

   int iregion = domain->find_region(arg[2]);
   //if (comm->me==0) fprintf(screen, "%d\n", iregion);
   //error->all(FLERR, "test from assign in group");

   if (iregion == -1) error->all(FLERR, "Group region ID does not exist");
   domain->regions[iregion]->init();
   domain->regions[iregion]->prematch();
   
   for(i = 0; i <nlocal; i++) {
     if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
       mask[i] |= bit;
   }

   for(i=0;i<enlocal;i++) {
      if (domain->regions[iregion]->match(ex[i][0],ex[i][1],ex[i][2]))
	emask[i] |= bit;
   }
 } else if (strcmp(arg[1],"subtract")==0) {
  
   if (narg < 4) error->all(FLERR,"Illegal group command");

   int length = narg-2;
   int *list = new int[length];

   int jgroup;
   for (int iarg = 2; iarg < narg ;iarg++) {
     jgroup = find(arg[iarg]);
     if (jgroup==-1) error->all(FLERR,"Group ID does not exist");
     list[iarg-2] = jgroup;
   }

   // add to group if in 1st group in list

   int otherbit = bitmask[list[0]];

   for (i = 0; i < nlocal; i++)
     if (mask[i] & otherbit) mask[i] |= bit;

   for (i=0; i<enlocal;i++)
     if (emask[i] & otherbit) emask[i] |= bit;

   //remove atoms if they are in any of the other groups
   //AND with inverse mask removes the atom from group

   int inverse = inversemask[igroup];

   for (int ilist = 1; ilist<length;ilist++) {
     otherbit = bitmask[list[ilist]];
     for (i=0;i<nlocal;i++)
       if (mask[i]&otherbit) mask[i] &= inverse;

     for (i=0;i<enlocal;i++)
       if (emask[i]&otherbit) emask[i] &= inverse;
   }

   delete [] list;
   // style = union
 } else if (strcmp(arg[1],"union")==0) {
   
   if (narg < 3) error->all(FLERR,"Illegal group command");

   int length = narg-2;
   int *list = new int[length];

   int jgroup;
   for (int iarg = 2; iarg < narg; iarg++) {
     jgroup = find(arg[iarg]);
     if (jgroup==-1) error->all(FLERR,"Group ID does not exist");
     list[iarg-2] = jgroup;
   }

   // add to group if in any other group in list
   
   int otherbit;

   for (int ilist = 0; ilist < length;ilist++) {
     otherbit = bitmask[list[ilist]];
     for (i=0;i<nlocal;i++)
       if (mask[i] & otherbit) mask[i] |= bit;

     for (i=0; i<enlocal;i++)
       if (emask[i] & otherbit) emask[i] |= bit;
   }

   delete [] list;
 } else error->all(FLERR, "Illegal group command");

 int n;
 n = 0;
 for (i = 0; i<nlocal; i++) if (mask[i] &bit) n++;

 int m;
 m=0;
 for (i=0; i<enlocal; i++) if (emask[i] &bit) m++;

 double rlocal = n;
 double all;

 MPI_Allreduce(&rlocal,&all,1,MPI_DOUBLE,MPI_SUM,world);

 double erlocal = m;
 double eall;

 MPI_Allreduce(&erlocal,&eall,1,MPI_DOUBLE,MPI_SUM,world);

 if (me==0) {
   if (screen && atom->nlocal)
       fprintf(screen,"%.15g atoms in group %s\n",all,names[igroup]);
   if (screen && cac->element_flag)
       fprintf(screen,"%.15g elements in group %s\n",eall,names[igroup]);
   if (logfile && atom->nlocal)
       fprintf(logfile,"%.15g atoms in group %s\n",all,names[igroup]);
   if (logfile && cac->element_flag)
       fprintf(logfile,"%.15g elements in group %s\n",eall,names[igroup]);
 }
}

/*-------------------------------------------------------------------------
   return index of first available group
   should never be called when group limit has been reached
---------------------------------------------------------------------------*/

int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == NULL) return igroup;
  return -1;
}
