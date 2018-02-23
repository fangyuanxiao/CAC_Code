#ifdef DUMP_CLASS

DumpStyle(xyz, DumpXyz)

#else

#ifndef CAC_DUMP_XYZ_H
#define CAC_DUMP_XYZ_H

#include "dump.h"

namespace CAC_NS {

class DumpXyz : public Dump {
  public:
   DumpXyz(CAC *, int, char **);

  protected:
   int scale_flag;         // 1 if atom coords are scaled, 0 if no
   int image_flag;         // 1 if append box count to atom coords, 0 if no

   char *columns;          // column labels
   void init_style();
   void write_header(bigint);
   void pack(tagint *);
   void epack() {}
   void write_data(int, double *);
   void ewrite_data(int, double *) {}

   void header_item(bigint);
   void pack_noscale_noimage(tagint *);
   void epack_noscale_noimage();
   void write_lines_noimage(int, double *);

};
}
#endif
#endif
