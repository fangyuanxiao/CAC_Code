#ifndef CAC_TIMER_H
#define CAC_TIMER_H

#include "pointers.h"

enum{TIME_LOOP,TIME_PAIR,TIME_BOND,TIME_KSPACE,TIME_NEIGHBOR,
     TIME_COMM,TIME_OUTPUT,TIME_N};

namespace CAC_NS {

class Timer : protected Pointers {
 public:
  double *array;

  Timer(class CAC *);
  ~Timer();
  void init();
  void stamp();
  void stamp(int);
  void barrier_start(int);
  void barrier_stop(int);
  double elapsed(int);

 private:
  double previous_time;
};

}

#endif
