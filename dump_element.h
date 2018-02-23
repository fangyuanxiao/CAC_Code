#ifdef DUMP_CLASS

DumpStyle(element, DumpElement)

#else

#ifndef CAC_DUMP_ELEMENT_H
#define CAC_DUMP_ELEMENT_H

#include "dump.h"

namespace CAC_NS {

class DumpElement : public Dump {
  public:
   DumpElement(CAC *, int,char**);

 protected:
   int scale_flag;                // 1 if element coords are scaled, 0 if no
   int image_flag;                // 1 if append box count to element coords, 0 if no
   //char eformat; 
   char *columns;                 // column labels
   void init_style();
   void write_header(bigint);
   void pack(tagint *);
   void epack();
   int convert_string(int, double *);
   int econvert_string(int, double *);

   void write_data(int, double *);
   void ewrite_data(int, double *);

   typedef void (DumpElement:: *FnPtrHeader) (bigint);
   FnPtrHeader header_choice;                       // ptr to write header functions
   void header_item(bigint);

   typedef void (DumpElement::*FnPtrWrite)(int,double *);
   FnPtrWrite write_choice;
   FnPtrWrite ewrite_choice;
   void write_string(int, double *);
   void write_lines_noimage(int, double *);
   void ewrite_string(int,double *);
   void ewrite_lines_noimage(int, double *);

   void epack_noscale_noimage();
   int convert_noimage(int, double *);
   int econvert_noimage(int, double *);

};
}

#endif
#endif
