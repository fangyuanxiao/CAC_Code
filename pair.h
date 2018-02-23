#ifndef CAC_PAIR_H
#define CAC_PAIR_H

#include "pointers.h"

namespace CAC_NS {

class Pair : protected Pointers {

 public:
  static int instance_total;     // # of Pair classes ever instantiated

  double eng_vdwl,eng_coul;      // accumulated energies
  double virial[6];              // accumulated virial
  double virialtmp[6];           // template virial stress
  double *eatom,**vatom;         // accumulated per-atom energy/virial
  double **enode,***vnode;       // accumulated per-node energy/virial

  double cutforce;               // max cutoff for all atom pairs
  double **cutsq;                // cutoff sq for each atom pair
  int **setflag;                 // 0/1 = whether each i,j has been set

  int comm_forward;              // size of forward communication (0 if none)
  int ecomm_forward;
  int comm_reverse;              // size of reverse communication (0 if none)
  int comm_reverse_off;          // size of reverse comm even if newton off

  int single_enable;             // 1 if single() routine exists
  int restartinfo;               // 1 if pair style writes restart info
  int respa_enable;              // 1 if inner/middle/outer rRESPA routines
  int one_coeff;                 // 1 if allows only one coeff * * call
  int manybody_flag;             // 1 if a manybody potential
  int no_virial_fdotr_compute;   // 1 if does not invoke virial_fdotr_compute()
  int writedata;                 // 1 if writes coeffs to data file
  int ghostneigh;                // 1 if pair style needs neighbors of ghosts
  double **cutghost;             // cutoff for each ghost pair

  int ewaldflag;                 // 1 if compatible with Ewald solver
  int pppmflag;                  // 1 if compatible with PPPM solver
  int msmflag;                   // 1 if compatible with MSM solver
  int dispersionflag;            // 1 if compatible with LJ/dispersion solver
  int tip4pflag;                 // 1 if compatible with TIP4P solver
  int dipoleflag;                // 1 if compatible with dipole solver
  int reinitflag;                // 1 if compatible with fix adapt and alike

  int tail_flag;                 // pair_modify flag for LJ tail correction
  double etail,ptail;            // energy/pressure tail corrections
  double etail_ij,ptail_ij;

  int evflag;                    // energy,virial settings
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;

  int ncoultablebits;            // size of Coulomb table, accessed by KSpace
  int ndisptablebits;            // size of dispersion table
  double tabinnersq;
  double tabinnerdispsq;
  double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
  double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
  double *rdisptable, *drdisptable, *fdisptable, *dfdisptable;
  double *edisptable, *dedisptable;
  int ncoulshiftbits,ncoulmask;
  int ndispshiftbits, ndispmask;

  int nextra;                    // # of extra quantities pair style calculates
  double *pvector;               // vector of extra pair quantities

  int single_extra;              // number of extra single values calculated
  double *svector;               // vector of extra single quantities

  class NeighList *list;         // standard neighbor list used by most pairs
  class NeighList *listinner;    // rRESPA lists used by some pairs
  class NeighList *listmiddle;
  class NeighList *listouter;
  class NodeNeighList *nilist;

  unsigned int datamask;
  unsigned int datamask_ext;

  int compute_flag;              // 0 if skip compute()

  Pair(class CAC *);
  virtual ~Pair();

  void init();

  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;
  virtual void write_data_all(FILE *) {}
  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual void init_nilist(class NodeNeighList *);
  virtual double init_one(int, int) {return 0.0;}
  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
  virtual void compute(int, int) = 0;  
  virtual void ecompute() = 0;
  void compute_dummy(int, int);
  virtual double memory_usage();

  virtual int pack_forward_comm(int, int *, double *, int, int *) {return 0;}
  virtual int epack_forward_comm(int, int *, double *) {return 0;}
  virtual void unpack_forward_comm(int,int,double *) {}
  virtual void eunpack_forward_comm(int, int, double *) {}
  virtual void ev_setup();

 protected:
  int instance_me;        // which Pair class instantiation I am

  enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};   // mixing options

  int special_lj[4];           // copied from force->special_lj for Kokkos

  int allocated;               // 0/1 = whether arrays are allocated
  int suffix_flag;             // suffix compatibility flag

  int offset_flag,mix_flag;            // flags for offset and mixing
  double tabinner;                     // inner cutoff for Coulomb table
  double tabinner_disp;                 // inner cutoff for dispersion table

  // custom data type for accessing Coulomb tables

  typedef union {int i; float f;} union_int_float_t;

  double THIRD;

  int vflag_fdotr;
  int maxeatom,maxvatom;
  int maxelement;

  int copymode;   // if set, do not deallocate during destruction
                  // required when classes are used as functors by Kokkos

  inline int sbmask(int j) {
    return j >> SBBITS & 3;
  }
                                  

};

}

#endif
