#ifndef CAC_COMM_BRICK_H
#define CAC_COMM_BRICK_H

#include "comm.h"

namespace CAC_NS{

class CommBrick : public Comm {
public:
  CommBrick(class CAC *);

  virtual ~CommBrick();

  virtual void init();
  virtual void setup();                        // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0);    // forward comm of atom coords
  virtual void eforward_comm(int dummy = 0);   // forward comm of element and node coords 
  virtual void exchange();                     // move atoms to new procs
  virtual void borders();                      // setup list of atoms to comm

  virtual void forward_comm_pair(class Pair *);  // forward comm from a pair
  virtual void eforward_comm_pair(class Pair *);  // forward comm of elements from a pair
  virtual bigint memory_usage();
protected:  
  int nswap;                       // # of swaps to perform = sum of maxneed
  int recvneed[3][2];              // # of procs away I recv atoms from
  int sendneed[3][2];              // # of procs away I send atoms to
  int maxneed[3];                  // max procs away any proc needs, per dim
  int maxswap;                     // max # of swaps memory is allocated for
  int *sendnum, *recvnum;          // # of atoms to send/recv in each swap
  int *sendproc, *recvproc;        // proc to send/recv to/from at each swap
  int *size_forward_recv;          // # of values to recv in each forward comm
  int *size_reverse_send;          // # to send in each reverse comm
  int *size_reverse_recv;          // # to recv in each reverse comm
  double *slablo,*slabhi;          // bounds of slab to send at each swap
  double **multilo, **multihi;     // bounds of slabs for multi-type swap
  double **cutghostmulti;          // cutghost on a per-type basis
  int *pbc_flag;                   // general flag for sending atoms thru PBC
  int **pbc;                       // dimension flags for PBC adjustments

  int *firstrecv;                  // where to put 1st recv atom in each swap
  int **sendlist;                  // list of atoms to send in each swap
  int *maxsendlist;                // max size of send list for each swap
   
  int *esendnum, *erecvnum;        // # of elements to send/recv in each swap 
  int *esize_forward_recv;         // # of values to recv in each forward comm for elements
  int *esize_reverse_send;         // # to send in each reverse comm for elements
  int *esize_reverse_recv;         // # to recv in each reverse comm for elements
  
  int *efirstrecv;                 // where to put 1st recv elements in each swap
  int **esendlist;                 // list of elements to send in each swap 
  int *emaxsendlist;               // max size of send list of elments for each swap 

  double *buf_send;                // send buffer for all comm
  double *buf_recv;                // recv buffer for all comm
  int maxsend, maxrecv;            // current size of send/recv buffer
  int bufextra;                    // extra space beyond maxsend in send buffer
  int smax, rmax;             // max size in atoms of single borders send/recv
  int esmax, ermax;           // max size in elements of single borders send/recv

  void init_buffers();
  virtual void grow_send(int, int);         // reallocate send buffer
  virtual void grow_recv(int);              // free/allocate recv buffer
  virtual void grow_swap(int);              // grow swap and multi arrays
  virtual void grow_list(int, int);         // reallocate one sendlist
  virtual void egrow_list(int,int);         // reallocate one esendlist
  virtual void allocate_swap(int);
  virtual void allocate_multi(int);         // allocate multi arrays
  virtual void free_swap();
  virtual void free_multi();                // free multi arrays
};

}

#endif
