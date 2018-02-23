#ifndef CAC_THERMO_H
#define CAC_THERMO_H

#include "pointers.h"

namespace CAC_NS {

class Thermo : protected Pointers {
 public:
  char *style;
  int normflag;          // 0 if do not normalize by atoms, 1 if normalize
  int modified;          // 1 if thermo_modify has been used, else 0
  int cudable;           // 1 if all computes used are cudable
  int lostflag;          // IGNORE,WARN,ERROR
  int lostbond;          // ditto for atoms in bonds

  Thermo(class CAC *, int, char **);
  ~Thermo();
  void init();
  void header();
  void compute(int);

 private:
  char *line;
  char **keyword;
  int *vtype;

  int nfield,nfield_initial;
  int me;

  char **format,**format_user;
  char *format_float_one_def,*format_float_multi_def;
  char *format_int_one_def,*format_int_multi_def;
  char *format_float_user,*format_int_user,*format_bigint_user;
  char format_multi[128];
  char format_bigint_one_def[8],format_bigint_multi_def[8];

  int normvalue;         // use this for normflag unless natoms = 0
  int normuserflag;      // 0 if user has not set, 1 if has
  int normuser;

  int firststep;
  int lostbefore;
  int flushflag,lineflag;

  double last_tpcpu,last_spcpu;
  double last_time;
  bigint last_step;

  bigint natoms;

  int ivalue;            // integer value to print
  double dvalue;         // double value to print
  bigint bivalue;        // big integer value to print
  int ifield;            // which field in thermo output is being computed
  int *field2index;      // which compute,fix,variable calcs this field
  int *argindex1;        // indices into compute,fix scalar,vector
  int *argindex2;

  int index_temp,index_press_scalar,index_press_vector,index_pe;
  char *id_temp,*id_press,*id_pe;
  class Compute *temperature,*pressure,*pe;

  int ncompute;                // # of Compute objects called by thermo
  char **id_compute;           // their IDs
  int *compute_which;          // 0/1/2 if should call scalar,vector,array
  class Compute **computes;    // list of ptrs to the Compute objects

  int nfix;                    // # of Fix objects called by thermo
  char **id_fix;               // their IDs
  class Fix **fixes;           // list of ptrs to the Fix objects

  int nvariable;               // # of variables evaulated by thermo
  char **id_variable;          // list of variable names
  int *variables;              // list of Variable indices
  
  void allocate();
  void deallocate();

  void parse_fields(char *);

  typedef void (Thermo::*FnPtr)();
  void addfield(const char *, FnPtr, int);
  FnPtr *vfunc;  

  void compute_step();
};
}

#endif
