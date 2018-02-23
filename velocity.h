#ifdef COMMAND_CLASS

CommandStyle(velocity,Velocity)

#else

#ifndef CAC_VELOCITY_H
#define CAC_VELOCITY_H

#include "pointers.h"

namespace CAC_NS {

class Velocity : protected Pointers {
  public: 
   Velocity(class CAC *);
   void command(int,char **);
   void options(int,char **);

  private:
   int igroup,groupbit;
   int style;
   int dist_flag,sum_flag,momentum_flag,rotation_flag;
   int bias_flag, loop_flag, scale_flag,rfix;
   double xscale,yscale,zscale;
   class Compute *temperature;

   void set(int, char **);
   void ramp(int, char **);

};
}

#endif
#endif
