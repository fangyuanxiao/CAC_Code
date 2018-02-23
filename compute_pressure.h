#ifdef COMPUTE_CLASS

ComputeStyle(pressure,ComputePressure)

#else

#ifndef CAC_COMPUTE_PRESSURE_H
#define CAC_COMPUTE_PRESSURE_H

#include "compute.h"

namespace CAC_NS {

class ComputePressure : public Compute {
 public:
  ComputePressure(class CAC *, int, char **);
  virtual ~ComputePressure();
  void init();

 protected:
  double boltz,nktv2p,inv_volume;
  int nvirial,dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int fixflag,kspaceflag;

};
}

#endif
#endif
