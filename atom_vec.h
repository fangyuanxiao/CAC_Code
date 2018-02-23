#ifndef CAC_ATOM_VEC_H
#define CAC_ATOM_VEC_H

#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class AtomVec : protected Pointers {
public:
  int molecular;                       // 0 = atomic, 1 = molecular system
  int bonds_allow,angles_allow;        // 1 if bonds, angles are used
  int dihedrals_allow,impropers_allow; // 1 if dihedrals, impropers used
  int mass_type;                       // 1 if per-type masses
  int dipole_type;                     // 1 if per-type dipole moments
  int forceclearflag;                  // 1 if has forceclear() method
  int comm_x_only;                     // 1 if only exchange x in forward comm
  int comm_f_only;                     // 1 if only exchange f in reverse comm
  int size_forward;                    // # of values per atom in comm
  int size_reverse;                    // # in reverse comm
  int size_border;                     // # in border comm
  int size_velocity;                   // # of velocity based quantities
  int size_data_atom;                  // number of values in Atom line
  int size_data_vel;                   // number of values in Velocity line
  int size_data_bonus;                 // number of values in Bonus line
  int xcol_data;                       // column (1-N) where x is in Atom line
  class Molecule **onemols;            // list of molecules for style template
  int nset;                            // # of molecules in list
  int cudable;                         // 1 if atom style is CUDA-enabled
  int kokkosable;                      // 1 if atom style is KOKKOS-enabled
  int *maxsend;                        // CUDA-specific variable
  int nargcopy;          // copy of command-line args for atom_style command
  char **argcopy;        // used when AtomVec is realloced (restart,replicate)

  AtomVec(class CAC *);
  virtual ~AtomVec();
  void store_args(int, char **);
  virtual void process_args(int, char**);

  virtual void grow(int) = 0;
  virtual void init();
  virtual void copy(int, int, int) = 0;
  virtual void clear_bonus() {}
  virtual int pack_exchange(int, double *) = 0;
  virtual int unpack_exchange(double *) = 0;
  virtual int pack_border_vel(int, int *, double *, int, int *) = 0;
  virtual void unpack_border(int, int, double *) = 0;
  virtual void unpack_border_vel(int, int, double *) = 0;
  virtual int pack_border(int, int *, double *, int, int *) = 0;
  virtual int pack_comm(int, int *, double *, int, int *) = 0;
  virtual int pack_comm_vel(int, int *, double *, int, int *) = 0;
  virtual void unpack_comm(int, int, double *) = 0;
  virtual void unpack_comm_vel(int, int, double *) = 0;
  virtual void data_atom(double *, imageint, char **) = 0;

  virtual bigint memory_usage() = 0;
 protected:
  int nmax;                             // local copy of atom->nmax
  int deform_vremap;                    // local copy of domain properties
  int deform_groupbit;
  double *h_rate;

  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };
  
  void grow_nmax();
};

}

#endif
