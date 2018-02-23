#ifndef CAC_DUMP_H
#define CAC_DUMP_H

#include "mpi.h"
#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump
  char *filename;            // user-specified file
  int igroup,groupbit;       // group that Dump is performed on

  int first_flag;            // 0 if no initial dump, 1 if yes initial dump
  int clearstep;             // 1 if dump invokes computes, 0 if not

  int comm_forward;          // size of forward communication (0 if none)
  int comm_reverse;          // size of reverse communication (0 if none)
  int nstep;

  static Dump *dumpptr;         // holds a ptr to Dump currently being used

  Dump(class CAC *, int, char **);
  virtual ~Dump();
  void init();
  virtual void write();
  virtual bigint memory_usage();

 protected:
  int me,nprocs;             // proc info

  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all,

  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  char *multiname;           // filename with % converted to cluster ID
  MPI_Comm clustercomm;      // MPI communicator within my cluster of procs

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int sort_flag;             // 1 if sorted output
  int append_flag;           // 1 if open file in append mode, 0 if not
  int buffer_allow;          // 1 if style allows for buffer_flag, 0 if not
  int buffer_flag;           // 1 if buffer output as one big string, 0 if not
  int padflag;               // timestep padding in filename
  int singlefile_opened;     // 1 = one big file, already opened, else 0
  int sortcol;               // 0 to sort on ID, 1-N on columns
  int sortcolm1;             // sortcol - 1
  int sortorder;             // ASCEND or DESCEND

  char boundstr[9];          // encoding of boundary flags
  char *format_default;      // default format string
  char *format_user;         // format string set by user
  char *format;              // format string for the file write
  char *eformat;       
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one atom
  int esize_one;
  int nme;                   // # of atoms in this dump from me
  int enme;                  // # of elements in this dump for me
  int ensme;                 // # of chars in string output from me
  int nsme;                  // # of chars in string output from me

  double boxxlo,boxxhi;      // local copies of domain values
  double boxylo,boxyhi;      // lo/hi are bounding box for triclinic
  double boxzlo,boxzhi;
  double boxxy,boxxz,boxyz;

  bigint ntotal;             // total # of per-atom lines in snapshot
  bigint entotal;            // total # of element lines in snapshot
  int reorderflag;           // 1 if OK to reorder instead of sort
  int ntotal_reorder;        // # of atoms that must be in snapshot
  int nme_reorder;           // # of atoms I must own in snapshot
  tagint idlo;               // lowest ID I own when reordering

  int maxbuf;                // size of buf
  int emaxbuf;
  double *buf;               // memory for atom quantities
  int maxsbuf;               // size of sbuf
  int emaxsbuf;
  char *sbuf;                // memory for atom quantities in string format

  int maxids;                // size of ids
  int maxsort;               // size of bufsort, idsort, index
  int maxproc;               // size of proclist
  tagint *ids;               // list of atom IDs, if sorting on IDs
  double *bufsort;
  tagint *idsort;
  int *index,*proclist;
 
  virtual void init_style() = 0;
  virtual void openfile();
  virtual void write_header(bigint) = 0;
  virtual int count();
  virtual void pack(tagint *) = 0;
  virtual void epack() = 0;
  virtual int convert_string(int, double *) {return 0;}
  virtual int econvert_string(int, double *) {return 0;}
  virtual void write_data(int, double *) = 0;
  virtual void ewrite_data(int, double *) = 0;
};
}

#endif
