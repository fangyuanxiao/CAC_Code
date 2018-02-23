#ifndef CAC_FIX_H
#define CAC_FIX_H

#include "pointers.h"

namespace CAC_NS {

class Fix : protected Pointers {
public:
  static int instance_total;     // # of Fix classes ever instantiated

  char *id,*style;
  int igroup,groupbit;

  int restart_global;            // 1 if Fix saves global state, 0 if not
  int restart_peratom;           // 1 if Fix saves peratom state, 0 if not
  int restart_file;              // 1 if Fix writes own restart file, 0 if not
  int force_reneighbor;          // 1 if Fix forces reneighboring, 0 if not

  int box_change_size;           // 1 if Fix changes box size, 0 if not
  int box_change_shape;          // 1 if Fix changes box shape, 0 if not
  int box_change_domain;         // 1 if Fix changes proc sub-domains, 0 if not

  bigint next_reneighbor;        // next timestep to force a reneighboring
  int thermo_energy;             // 1 if fix_modify enabled ThEng, 0 if not
  int nevery;                    // how often to call an end_of_step fix
  int rigid_flag;                // 1 if Fix integrates rigid bodies, 0 if not
  int virial_flag;               // 1 if Fix contributes to virial, 0 if not
  int no_change_box;             // 1 if cannot swap ortho <-> triclinic
  int time_integrate;            // 1 if fix performs time integration, 0 if no
  int time_depend;               // 1 if requires continuous timestepping
  int create_attribute;          // 1 if fix stores attributes that need
                                 //      setting when a new atom is created
  int restart_pbc;               // 1 if fix moves atoms (except integrate)
                                 //       so write_restart must remap to PBC
  int wd_header;                 // # of header values fix writes to data file
  int wd_section;                // # of sections fix writes to data file
  int dynamic_group_allow;       // 1 if can be used with dynamic group, else 0
  int dof_flag;                  // 1 if has dof() method (not min_dof())
  int special_alter_flag;        // 1 if has special_alter() meth for spec lists
  int cudable_comm;              // 1 if fix has CUDA-enabled communication

  int scalar_flag;               // 0/1 if compute_scalar() function exists
  int vector_flag;               // 0/1 if compute_vector() function exists
  int array_flag;                // 0/1 if compute_array() function exists
  int size_vector;               // length of global vector
  int size_array_rows;           // rows in global array
  int size_array_cols;           // columns in global array
  int size_vector_variable;      // 1 if vec length is unknown in advance
  int size_array_rows_variable;  // 1 if array rows is unknown in advance
  int global_freq;               // frequency s/v data is available at

  int peratom_flag;              // 0/1 if per-atom data is stored
  int size_peratom_cols;         // 0 = vector, N = columns in peratom array
  int peratom_freq;              // frequency per-atom data is available at

  int local_flag;                // 0/1 if local data is stored
  int size_local_rows;           // rows in local vector or array
  int size_local_cols;           // 0 = vector, N = columns in local array
  int local_freq;                // frequency local data is available at

  int extscalar;            // 0/1 if global scalar is intensive/extensive
  int extvector;            // 0/1/-1 if global vector is all int/ext/extlist
  int *extlist;             // list of 0/1 int/ext for each vec component
  int extarray;             // 0/1 if global array is intensive/extensive

  double *vector_atom;           // computed per-atom vector
  double **array_atom;           // computed per-atom array
  double *vector_local;          // computed local vector
  double **array_local;          // computed local array

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)
  int comm_border;               // size of border communication (0 if none)

  double virial[6];              // accumlated virial
  double **vatom;                // accumulated per-atom virial

  int restart_reset;             // 1 if restart just re-initialized fix

  Fix(class CAC *, int, char **);
  virtual ~Fix();
  virtual void init() {}
  virtual void setup(int) {}
  virtual int setmask() = 0;
  virtual void initial_integrate(int) {}
  virtual void einitial_integrate(int) {}
  virtual void final_integrate() {}
  virtual void efinal_integrate() {}
  virtual void post_constructor() {}
  virtual void post_run() {}
  virtual void post_force(int) {}
  virtual void epost_force(int) {}
  virtual double memory_usage() {return 0.0;}
 protected:
  int instance_me;        // which Fix class instantiation I am

  int evflag;
  int vflag_global,vflag_atom;
  int maxvatom;

  int copymode;   // if set, do not deallocate during destruction
                  // required when classes are used as functors by Kokkos


  // union data struct for packing 32-bit and 64-bit ints into double bufs
  // see atom_vec.h for documentation

  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };
};

namespace FixConst {
  static const int INITIAL_INTEGRATE =       1<<0;
  static const int POST_INTEGRATE =          1<<1;
  static const int PRE_EXCHANGE =            1<<2;
  static const int PRE_NEIGHBOR =            1<<3;
  static const int PRE_FORCE =               1<<4;
  static const int POST_FORCE =              1<<5;
  static const int FINAL_INTEGRATE =         1<<6;
  static const int END_OF_STEP =             1<<7;
  static const int THERMO_ENERGY =           1<<8;
  static const int INITIAL_INTEGRATE_RESPA = 1<<9;
  static const int POST_INTEGRATE_RESPA =    1<<10;
  static const int PRE_FORCE_RESPA =         1<<11;
  static const int POST_FORCE_RESPA =        1<<12;
  static const int FINAL_INTEGRATE_RESPA =   1<<13;
  static const int MIN_PRE_EXCHANGE =        1<<14;
  static const int MIN_PRE_NEIGHBOR =        1<<15;
  static const int MIN_PRE_FORCE =           1<<16;
  static const int MIN_POST_FORCE =          1<<17;
  static const int MIN_ENERGY =              1<<18;
  static const int POST_RUN =                1<<19;
  static const int FIX_CONST_LAST =          1<<20;
}

}

#endif
