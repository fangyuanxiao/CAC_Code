#include "string.h"
#include "dump_element.h"
#include "domain.h"
#include "element.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

/*--------------------------------------------------------------------------------*/

DumpElement::DumpElement(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
   if (narg != 5) error->all(FLERR,"Illegal dump element command");

   scale_flag = 0;
   image_flag = 0;
   buffer_allow = 1;
   buffer_flag = 1;
   format_default = NULL;
  // error->all(FLERR, "test from dump element");
}

/*--------------------------------------------------------------------------------*/

void DumpElement::init_style()
{
  if (image_flag == 0) size_one = 9;
  else size_one = 12;

  esize_one = element->npe;

  delete [] format;
  char *str;
  str = (char *) "%g %g %g %g %g %g %g %g %g";
  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  delete [] eformat;
  str = (char *) BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT " " BIGINT_FORMAT;
  n = strlen(str) + 2;
  eformat = new char[n];
  strcpy(eformat, str);
  strcat(eformat,"\n");

  columns = (char *) "\"x\", \"y\", \"z\", \"vx\", \"vy\", \"vz\", \"fx\", \"fy\", \"fz\"";

  header_choice = &DumpElement::header_item;
  //error->all(FLERR,"test from init_style in dump element");

  //pack_choice = &DumpElement::pack_noscale_noimage;
  
  if (buffer_flag ==1) write_choice = &DumpElement::write_string;
  else write_choice = &DumpElement::write_lines_noimage;

  if (buffer_flag ==1) ewrite_choice = &DumpElement::ewrite_string;
  else ewrite_choice = &DumpElement::ewrite_lines_noimage;


}

/*------------------------------------------------------------------------------------------*/

void DumpElement::header_item(bigint ndump)
{
  int npe = element->npe;
  bigint nele = ndump/npe;
  fprintf(fp, "title=\"results for cac\"\n");
  fprintf(fp, "variables=%s\n", columns);
  fprintf(fp,"zone t=\"load step "BIGINT_FORMAT "\",", update->ntimestep);
  fprintf(fp," n=  " BIGINT_FORMAT" e =  " BIGINT_FORMAT, ndump, nele);
  fprintf(fp," datapacking=point, zonetype=febrick\n");
  //error->all(FLERR, "test from header_item in dump element");
}

/*-------------------------------------------------------------------------------------------*/

void DumpElement::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice) (ndump);
}

/*-----------------------------------------------------------------------------------------*/

void DumpElement::pack(tagint *ids)
{
   epack_noscale_noimage();
}

/*------------------------------------------------------------------------------------------*/

void DumpElement::epack_noscale_noimage()
{
  int m,n;

  int *mask = element->mask;
  double ***nodex = element->nodex;
  double ***nodef = element->nodef;
  double ***nodev = element->nodev;
  int nlocal = element->nlocal;
  int npe = element->npe;

  m=n=0;
  for (int i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j=0; j<npe; j++) {
	buf[m++] = nodex[i][j][0];
	buf[m++] = nodex[i][j][1];
	buf[m++] = nodex[i][j][2];
	buf[m++] = nodev[i][j][0];
	buf[m++] = nodev[i][j][1];
	buf[m++] = nodev[i][j][2];
	buf[m++] = nodef[i][j][0];
	buf[m++] = nodef[i][j][1];
	buf[m++] = nodef[i][j][2];
      }
    }
  }
}

void DumpElement::epack()
{
  int m;

  tagint **nodetag = element->nodetag;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int npe = element->npe;

  m = 0;
  for (int i = 0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j<npe; j++) {
	buf[m++] = nodetag[i][j];
      }
    }
  }
}

/*--------------------------------------------------------------------------------------------*/

int DumpElement::convert_string(int n, double *mybuf)
{
  return convert_noimage(n,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpElement::convert_noimage(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i<n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset], format,
		    mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],
		    mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],
		    mybuf[m+8]);
    m += size_one;
  }
  return offset;
}

/*----------------------------------------------------------------------------------------------------*/
int DumpElement::econvert_string(int n, double *mybuf)
{
  return econvert_noimage(n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

int DumpElement::econvert_noimage(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i <n;i++) {
    offset += sprintf(&sbuf[offset],eformat,
		    static_cast<tagint>(mybuf[m]),
		    static_cast<tagint>(mybuf[m+1]),
		    static_cast<tagint>(mybuf[m+2]),
		    static_cast<tagint>(mybuf[m+3]),
		    static_cast<tagint>(mybuf[m+4]),
		    static_cast<tagint>(mybuf[m+5]),
		    static_cast<tagint>(mybuf[m+6]),
		    static_cast<tagint>(mybuf[m+7]));
    m += esize_one;
  }
  return offset;
}	

/*----------------------------------------------------------------------------------------------------*/

void DumpElement::write_data(int n, double *mybuf)
{
   (this->*write_choice) (n,mybuf);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpElement::write_string(int n, double *mybuf)
{
  int m;
  fwrite(mybuf, sizeof(char),n,fp);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpElement::ewrite_string(int n, double *mybuf)
{
   int m;
   fwrite(mybuf,sizeof(char),n,fp);
}

/*-------------------------------------------------------------------------------------------------*/

void DumpElement::write_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format, mybuf[m], mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],
		    mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8]);
    m += size_one;
  }
}

/*----------------------------------------------------------------------------------------------------*/

void DumpElement::ewrite_lines_noimage(int n, double *mybuf)
{
   int m = 0;
   for (int i = 0; i<n; i++) {
     fprintf(fp, eformat, static_cast<tagint>(mybuf[m]),
		          static_cast<tagint>(mybuf[m+1]),
			  static_cast<tagint>(mybuf[m+2]),
			  static_cast<tagint>(mybuf[m+3]),
			  static_cast<tagint>(mybuf[m+4]),
		          static_cast<tagint>(mybuf[m+5]),
			  static_cast<tagint>(mybuf[m+6]),
			  static_cast<tagint>(mybuf[m+7]));
     m += esize_one;
   }
}

/*----------------------------------------------------------------------------------------------------*/

void DumpElement::ewrite_data(int n, double *mybuf)
{
   (this->*ewrite_choice) (n, mybuf);
}
