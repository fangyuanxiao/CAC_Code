#include "string.h"
#include "dump_xyz.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "output.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

/*----------------------------------------------------------------------------*/

DumpXyz::DumpXyz(CAC *cac, int narg, char **arg) : Dump(cac, narg,arg)
{
  if (narg != 5) error->all(FLERR, "Illegal dump xyz command");

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 0;
  buffer_flag = 0;
  format_default = NULL;
}

/*------------------------------------------------------------------------------------*/

void DumpXyz::init_style()
{
  if (image_flag == 0) size_one=4;
  else size_one = 6;

  delete [] format;
  char *str;
  str = (char *) "Cu %g %g %g %g";
  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format, str);
  strcat(format,"\n");

  if (multifile == 0) openfile();
}

/*----------------------------------------------------------------------------------------*/

void DumpXyz::header_item(bigint ndump)
{
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"XYZ\n");
}

/*------------------------------------------------------------------------------------------*/

void DumpXyz::pack_noscale_noimage(tagint *ids)
{
   int m,n;

   double **x = atom->x;
   double *centro = output->centro; 
   int nlocal = atom->nlocal;

   m = n = 0;
   for (int i = 0; i<nlocal;i++) {
     if (centro[i]>1 && centro[i]<10) {
       buf[m++] = x[i][0];
       buf[m++] = x[i][1];
       buf[m++] = x[i][2];
       buf[m++] = centro[i];
     }
  }
}

/*----------------------------------------------------------------------------------------*/

void DumpXyz::epack_noscale_noimage()
{
  int m,n;

  int *mask = element->mask;
  double ***nodex = element->nodex;

  int nlocal = element->nlocal;
  int npe = element->npe;

  m = n = 0;

  for (int i = 0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j=0; j<npe; j++) {
	buf[m++] = nodex[i][j][0];
	buf[m++] = nodex[i][j][1];
	buf[m++] = nodex[i][j][2];
	//error->all(FLERR, "test from dump xyz");
      }
    }
  }
  //error->all(FLERR, "test from dump xyz");
}

/*------------------------------------------------------------------------------------------*/

void DumpXyz::pack(tagint *ids)
{
   if (atom->nlocal) pack_noscale_noimage(ids);
}
/*------------------------------------------------------------------------------------------*/

void DumpXyz::write_header(bigint ndump)
{
   header_item(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpXyz::write_data(int n, double *mybuf)
{
   write_lines_noimage(n, mybuf);
}

/*------------------------------------------------------------------------------------------*/

void DumpXyz::write_lines_noimage(int n, double *mybuf)
{

  double ixtmp, iytmp, iztmp, icentro; 
  int m = 0;
  for (int i = 0; i<n;i++) {
   
    ixtmp = buf[m];
    iytmp = buf[m+1];
    iztmp = buf[m+2];
    icentro = buf[m+3];
    m += size_one;
    fprintf(fp,format,ixtmp,iytmp,iztmp,icentro);
  }

}


