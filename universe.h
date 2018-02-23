#ifndef CAC_UNIVERSE_H
#define CAC_UNIVERSE_H

#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class Universe : protected Pointers {
 public:

  MPI_Comm uworld;        // communicator for entire universe
  int me,nprocs;          // my place in universe

  FILE *uscreen;          // universe screen output
  FILE *ulogfile;         // universe logfile

  int existflag;          // 1 if universe exists due to -partition flag
  int nworlds;            // # of worlds in universe
  int iworld;             // which world I am in
  int *procs_per_world;   // # of procs in each world
  int *root_proc;         // root proc in each world

  MPI_Comm uorig;         // original communicator passed to LAMMPS instance
  int *uni2orig;          // proc I in universe uworld is
                          // proc uni2orig[I] in original communicator

  double numeric(const char *, int, char *);
  int inumeric(const char *, int, char *);
  bigint bnumeric(const char *, int, char *);
  tagint tnumeric(const char *, int, char *);

  Universe(class CAC *, MPI_Comm);
  ~Universe();
  void reorder(char *, char *);
  void add_world(char *);
  int consistent();
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid -reorder N value

Self-explanatory.

E: Nprocs not a multiple of N for -reorder

Self-explanatory.

E: Cannot open -reorder file

Self-explanatory.

E: Unexpected end of -reorder file

Self-explanatory.

E: Invalid entry in -reorder file

Self-explanatory.

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch LAMMPS.

*/
