#ifndef CAC_GROUP_H
#define CAC_GROUP_H

#include "stdio.h"
#include "pointers.h"
#include <map>

namespace CAC_NS {

class Group : protected Pointers {

public:
  int ngroup;                  // # of defined groups
  char **names;                // name of each group
  int *bitmask;                // one-bit mask for each group
  int *inversemask;            // inverse mask for each group
  int *dynamic;                // 1 if dynamic, 0 if not
  
  Group(class CAC *);
  ~Group();
  int find(const char *);       // lookup name in list of groups
  void assign(int,char **); 

private:
  int me;
  std::map<tagint,int> *hash;
  
  int find_unused();

  static Group *cptr;
  int molbit;
  

};

}

#endif
