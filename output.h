#ifndef CAC_OUTPUT_H
#define CAC_OUTPUT_H

#include "pointers.h"

namespace CAC_NS {

class Output : protected Pointers {
 public:
  bigint next;                 // next timestep for any kind of output

  bigint next_thermo;          // next timestep for thermo output
  int thermo_every;            // output freq for thermo, 0 if first/last only
  bigint last_thermo;          // last timestep thermo was output
  char *var_thermo;            // variable name for thermo freq, NULL if every
  int ivar_thermo;             // variable index for thermo frequency
  class Thermo *thermo;        // Thermodynamic computations

  int ndump;                   // # of Dumps defined
  int max_dump;                // max size of Dump list
  bigint next_dump_any;        // next timestep for any Dump
  int *every_dump;             // write freq for each Dump, 0 if var
  bigint *next_dump;           // next timestep to do each Dump
  bigint *last_dump;           // last timestep each snapshot was output
  char **var_dump;             // variable name for dump frequency
  int *ivar_dump;              // variable index for dump frequency
  class Dump **dump;           // list of defined Dumps

  int restart_flag;            // 1 if any restart files are written
  int restart_flag_single;     // 1 if single restart files are written
  int restart_flag_double;     // 1 if double restart files are written
//  bigint next_restart;         // next timestep to write any restart file
//  bigint next_restart_single;  // next timestep to write a single restart file
//  bigint next_restart_double;  // next timestep to write a double restart file
//  int restart_every_single;    // single restart file write freq, 0 if var
//  int restart_every_double;    // double restart file write freq, 0 if var
//  bigint last_restart;         // last timestep any restart file was output
//  int restart_toggle;          // 0 if use restart2a as prefix, 1 if restart2b
//  char *var_restart_single;    // variable name for single restart freq
//  char *var_restart_double;    // variable name for double restart freq
//  int ivar_restart_single;     // index of var_restart_single
//  int ivar_restart_double;     // index of var_restart_double
//  char *restart1;              // name single restart file
//  char *restart2a,*restart2b;  // names of double restart files
//  class WriteRestart *restart; // class for writing restart files
  int nmax,maxneigh,nnn;
  double *distsq;
  int *nearest;
  double *centro;

  Output(class CAC *);
  ~Output();
  void init();
  void setup(int memflag = 1);       // initial output before run/min
  void set_thermo(int, char **);     // set thermo output freqquency
  void create_thermo(int, char **);  // create a thermo style
  void write(bigint);                // output for current timestep
  void add_dump(int, char **);       // add a Dump to Dump list
  void memory_usage();               // print out memory usage
  void centro_atom(class NeighList *); // centro symmetric parameter calculated
  void select(int, int, double *);
  void select2(int, int, double *, int *);
//  void reset_timestep(bigint);       // reset next timestep for all output

};

}

#endif
