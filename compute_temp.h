#ifdef COMPUTE_CLASS

ComputeStyle(temp,ComputeTemp)

#else

#ifndef CAC_COMPUTE_TEMP_H
#define CAC_COMPUTE_TEMP_H

#include "compute.h"

namespace CAC_NS {

class ComputeTemp : public Compute {
 public:
  ComputeTemp(class CAC *, int, char **);
  virtual ~ComputeTemp();
  void init() {}

 protected:
  double tfactor;

};

}

#endif

#endif
