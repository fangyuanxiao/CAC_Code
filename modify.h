#ifndef CAC_MODIFY_H
#define CAC_MODIFY_H

#include "stdio.h"
#include "pointers.h"
#include <map>
#include <string>

namespace CAC_NS {

class Modify : protected Pointers {
 public:
  int nfix,maxfix;
  int n_initial_integrate,n_post_integrate,n_pre_exchange,n_pre_neighbor;
  int n_pre_force,n_post_force;
  int n_final_integrate,n_end_of_step,n_thermo_energy;
  int n_initial_integrate_respa,n_post_integrate_respa;
  int n_pre_force_respa,n_post_force_respa,n_final_integrate_respa;
  int n_min_pre_exchange,n_min_pre_neighbor;
  int n_min_pre_force,n_min_post_force,n_min_energy;

  int restart_pbc_any;       // 1 if any fix sets restart_pbc
  int nfix_restart_global;   // stored fix global info from restart file
  int nfix_restart_peratom;  // stored fix peratom info from restart file

  class Fix **fix;           // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  int ncompute,maxcompute;   // list of computes
  class Compute **compute;

  Modify(class CAC *);
  virtual ~Modify();
  virtual void init();
  virtual void setup(int);
  virtual void setup_pre_exchange();
  virtual void setup_pre_neighbor();
  virtual void setup_pre_force(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void post_force(int);
  void add_fix(int, char **, int trysuffix=0); 
 
  void add_compute(int, char **, int trysuffix=0);
  void delete_fix(const char *);
  int  find_fix(const char *);
  void restart_deallocate();
  int find_compute(const char *);

  void addstep_compute(bigint);
  void addstep_compute_all(bigint);
  bigint memory_usage();
  virtual void post_run();

protected:

  // lists of fixes to apply at different stages of timestep

  int *list_initial_integrate,*list_post_integrate;
  int *list_pre_exchange,*list_pre_neighbor;
  int *list_pre_force,*list_post_force;
  int *list_final_integrate,*list_end_of_step,*list_thermo_energy;
  int *list_initial_integrate_respa,*list_post_integrate_respa;
  int *list_pre_force_respa,*list_post_force_respa;
  int *list_final_integrate_respa;
  int *list_min_pre_exchange,*list_min_pre_neighbor;
  int *list_min_pre_force,*list_min_post_force;
  int *list_min_energy;

  int *end_of_step_every;

  int n_timeflag;            // list of computes that store time invocation
  int *list_timeflag;

  char **id_restart_global;           // stored fix global info
  char **style_restart_global;        // from read-in restart file
  char **state_restart_global;

  char **id_restart_peratom;          // stored fix peratom info
  char **style_restart_peratom;       // from read-in restart file
  int *index_restart_peratom;

  int index_permanent;        // fix/compute index returned to library call

  void list_init(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
  void list_init_thermo_energy(int, int &, int *&);
  void list_init_compute();

private:
  typedef Compute *(*ComputeCreator)(CAC *, int, char **);
  std::map<std::string,ComputeCreator> *compute_map;

  typedef Fix *(*FixCreator)(CAC *, int, char **);
  std::map<std::string,FixCreator> *fix_map;

  template <typename T> static Compute *compute_creator(CAC *, int, char **);
  template <typename T> static Fix *fix_creator(CAC *, int, char **);

};

}

#endif
