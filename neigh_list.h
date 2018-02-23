#ifndef CAC_NEIGH_LIST_H
#define CAC_NEIGH_LIST_H

#include "pointers.h"
#include "my_page.h"

namespace CAC_NS {

class NeighList : protected Pointers {

public:
  int index;                       // index of which neigh list it is
                                   // needed when a class invokes it directly
                                   // also indexes the request it came from

  int buildflag;                   // 1 if pair_build invoked every reneigh
  int growflag;                    // 1 if stores atom-based arrays & pages
  int stencilflag;                 // 1 if stores stencil arrays
  int ghostflag;                   // 1 if it stores neighbors of ghosts

  // data structs to store neighbor pairs I,J and associated values
  
  int inum;                        // # of I atoms neighbors are stored for
  int einum;                       // # of I element neighbors are stored for
  int gnum;                        // # of ghost atoms neighbors are stored for
  int *ilist;                      // local indices of I atoms
  int *numneigh;                   // # of J neighbors for each I atom
  int **firstneigh;                // ptr to 1st J int value of each I atom
  double **firstdouble;            // ptr to 1st J double value of each I atom

  int *eilist;                     // local indices of I elements
  int *enumneigh;                  // # of J neighbors for each I element
  int **efirstneigh;                // ptr to 1st J int value of each I atom 

  int pgsize;                      // size of each page
  int oneatom;                     // max size for one atom
  int dnum;                        // # of doubles per neighbor, 0 if none
  MyPage<int> *ipage;              // pages of atom neighbor indices
  MyPage<int> *epage;              // pages of element neighbor indices
  MyPage<double> *dpage;           // pages of neighbor doubles, if dnum > 0

  bigint last_build;           // timestep of last build for occasional lists

  // atom types to skip when building list
  // iskip,ijskip are just ptrs to corresponding request
  
  int *iskip;         // iskip[i] = 1 if atoms of type I are not in list
  int **ijskip;       // ijskip[i][j] = 1 if pairs of type I,J are not in list

  int respamiddle;              // 1 if this respaouter has middle list
  NeighList *listinner;         // me = respaouter, point to respainner
  NeighList *listmiddle;        // me = respaouter, point to respamiddle
  NeighList *listfull;          // me = half list, point to full I derive from
  NeighList *listcopy;          // me = copy list, point to list I copy from
  NeighList *listskip;          // me = skip list, point to list I skip from

  int maxstencil;                  // max size of stencil
  int nstencil;                    // # of bins in stencil
  int *stencil;                    // list of bin offsets
  int **stencilxyz;                // bin offsets in xyz dims

  int maxstencil_multi;            // max sizes of stencils
  int *nstencil_multi;             // # bins in each type-based multi stencil
  int **stencil_multi;             // list of bin offsets in each stencil
  double **distsq_multi;           // sq distances to bins in each stencil

  NeighList(class CAC *);
  virtual ~NeighList();
  void setup_pages(int, int, int);      // setup page data structures
  void grow(int);                       // grow maxlocal for atoms
  void egrow(int);                      // grow maxlocal for elements
  void stencil_allocate(int, int);      // allocate stencil arrays
  void copy_skip_info(int *, int **);   // copy skip info from a neigh request
  bigint memory_usage();

  protected:
    int maxatoms;                  // size of allocated atom arrays
    int maxeles;                    // size of allocated element arrays

};
   
}

#endif
