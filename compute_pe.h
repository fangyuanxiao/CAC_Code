#ifdef COMPUTE_CLASS

ComputeStyle(pe,ComputePE)

#else

#ifndef CAC_COMPUTE_PE_H
#define CAC_COMPUTE_PE_H

#include "compute.h"

namespace CAC_NS {

class ComputePE : public Compute {
 public:
  ComputePE(class CAC *, int, char **);
  ~ComputePE() {}
  void init() {}
  double compute_scalar();

 private:
  int pairflag,bondflag,angleflag,dihedralflag,improperflag,kspaceflag;
};

}

#endif
#endif
