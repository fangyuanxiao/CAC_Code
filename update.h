#ifndef CAC_UPDATE_H
#define CAC_UPDATE_H

#include "pointers.h"

namespace CAC_NS {

class Update : protected Pointers {
 public:
  double dt;                      // timestep
  double etol,ftol;               // minimizer tolerances on energy/force
  bigint ntimestep;               // current step (dynamics or min iterations)
  int nsteps;                     // # of steps to run (dynamics or min iter)
  int whichflag;                  // 0 for unset, 1 for dynamics, 2 for min
  double atime;                   // simulation time at atime_step
  bigint atimestep;               // last timestep atime was updated
  bigint firststep,laststep;      // 1st & last step of this run
  bigint beginstep,endstep;       // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after
  int max_eval;                   // max force evaluations for minimizer
  int restrict_output;            // 1 if output should not write dump/restart
  int setupflag;                  // set when setup() is computing forces
  int multireplica;               // 1 if min across replicas, else 0

  bigint eflag_global,eflag_atom;  // timestep global/peratom eng is tallied on
  bigint vflag_global,vflag_atom;  // ditto for virial

  char *unit_style;

  class Integrate *integrate;
  char *integrate_style;

  Update(class CAC *);
  ~Update();
  void init();
  void set_units(const char *);
  void create_integrate(int, char **, int);
//  void reset_timestep(int, char **);
//  void reset_timestep(bigint);
  bigint memory_usage();
  void update_time();

private:
  void new_integrate(char *, int, char **, int, int &);

};

}

#endif
