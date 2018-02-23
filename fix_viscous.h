#ifdef FIX_CLASS

FixStyle(viscous, FixViscous)

#else

#ifndef CAC_FIX_VISCOUS_H
#define CAC_FIX_VISCOUS_H

#include "fix.h"

namespace CAC_NS {

class FixViscous : public Fix {
public:
  FixViscous(class CAC *, int, char **);
  virtual ~FixViscous();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void epost_force(int);

protected:
  double *gamma;
  double *egamma;
  int nlevels_respa;

};

}

#endif
#endif
