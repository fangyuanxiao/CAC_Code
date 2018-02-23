#ifdef DUMP_CLASS

DumpStyle(atom,DumpAtom)

#else

#ifndef CAC_DUMP_ATOM_H
#define CAC_DUMP_ATOM_H

#include "dump.h"

namespace CAC_NS {

class DumpAtom : public Dump {
 public:
  DumpAtom(CAC *, int, char**);

 protected:
  int scale_flag;            // 1 if atom coords are scaled, 0 if no
  int image_flag;            // 1 if append box count to atom coords, 0 if no

  char *columns;             // column labels
  void init_style();
  void write_header(bigint);
  void pack(tagint *);
  void epack () {}
  int convert_string(int, double *);
  void write_data(int, double *);
  void ewrite_data(int,double *) {}
  
  typedef void (DumpAtom::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_item(bigint);
  
  typedef void (DumpAtom::*FnPtrPack)(tagint *);
  FnPtrPack pack_choice;               // ptr to pack functions
  void pack_noscale_noimage(tagint *);
  void epack_noscale_noimage();

  typedef int (DumpAtom::*FnPtrConvert)(int, double *);
  FnPtrConvert convert_choice;          // ptr to convert data functions
  int convert_noimage(int, double *);

  typedef void (DumpAtom::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;              // ptr to write data functions
  void write_string(int, double *);
  void write_lines_noimage(int, double *);
};
}
#endif
#endif
