#ifndef CAC_NODE_NEIGH_LIST_H
#define CAC_NODE_NEIGH_LIST_H

#include "pointers.h"
#include "my_page.h"

namespace CAC_NS {

class NodeNeighList :  protected Pointers {

public: 
  int index;                                  // index of which neigh list it is 
                                              // needed when a class invokes it directly
					      // also indexes the request it came from
 int buildflag;                               // 1 if pair_build invoked every reneigh
 int growflag;                                // 1 if stores atom-based arrays & pages
 int stencilflag;                             // 1 if it stores neighbors of ghosts

 // data structs to store neighbor pairs I, J and associated values
 
 int ninum;                                   // # of I integration points are stored for
 int neinum;                                   // # of I elements are stored for
 int *elist;                                  // local indices of I elements
 int *enlist;                                 // correspondance between elements and nodes
 int *numneigh;                               // # of J neighbors for each I integration points
 int **firstneigh;                            // ptr to 1st J int value of each I integration points
 int **nfirstneigh;                           // ptr to 1st J int value in each element of each I integration points
 
 int pgsize;                                  // size of each page
 int onenode;                                 // max size for one integration points
 int dnum;                                    // # of doubles per neighbor, 0 if none
 MyPage<int> *npage;                          // pages of integration point neighbor indices
 MyPage<int> *nnpage;                         // pages of integration point neighbor indices along x direction

 NodeNeighList(class CAC *);
 virtual ~NodeNeighList();
 void grow(int,int);                          // first int is number of elements and the second one is number of integration points
 void setup_pages(int, int, int);             // setup page data structures 
 bigint memory_usage();
 protected:
   int maxnodes;                              // size of allocated node arrays
   int maxeles;                               // size of allocated element arrays 

};
}

#endif
