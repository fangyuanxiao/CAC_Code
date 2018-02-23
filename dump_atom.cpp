#include "string.h"
#include "dump_atom.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpAtom::DumpAtom(CAC *cac, int narg, char **arg) : Dump(cac, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump atom command");

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  format_default = NULL;

}

/* ---------------------------------------------------------------------- */

void DumpAtom::init_style()
{
  if (image_flag == 0) size_one = 11;
  else size_one = 8;

  //if (comm->me==0) fprintf(screen, "size_one = %d\n", size_one);  
  // default format depends on image flags

  delete [] format;
  char *str;
  str = (char *) TAGINT_FORMAT " %d %g %g %g %g %g %g %g %g %g";
  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

 // setup boundary string
 
  domain->boundary_string(boundstr); 
  
  // setup column string

    columns = (char *) "id type x y z fx fy fz vx vy vz";
 

  // setup function ptrs

  //if (comm->me==0) fprintf(screen, "binary = %d\n", binary);
  header_choice = &DumpAtom::header_item;

  //if (comm->me==0) fprintf(screen, "scale_flag = %d, image_flag = %d\n", scale_flag, image_flag);
  pack_choice = &DumpAtom::pack_noscale_noimage;

  convert_choice = &DumpAtom::convert_noimage;
  
  if (buffer_flag == 1) write_choice = &DumpAtom::write_string;
  else write_choice = &DumpAtom::write_lines_noimage;

  // open single file, one time only

  if (multifile == 0) openfile();   
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_noscale_noimage(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = f[i][0];
      buf[m++] = f[i][1];
      buf[m++] = f[i][2];
      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];
      if (ids) ids[n++] = tag[i];
    }
}

/*-----------------------------------------------------------------------*/

void DumpAtom::epack_noscale_noimage()
{
   int m,n;

   tagint **nodetag = element->nodetag;
   int *mask = element->mask;
   double ***nodex = element->nodex;
   double ***nodef = element->nodef;
   double ***nodev = element->nodev;
   int nlocal = element->nlocal;
   int npe = element->npe;

   m = n =0;
   for (int i =0; i <nlocal; i++) {
     if (mask[i] & groupbit) {
       for (int j=0;j<npe;j++) {
	 buf[m++] = nodetag[i][j];
	 buf[m++] = 1;
	 buf[m++] = nodex[i][j][0];
	 buf[m++] = nodex[i][j][1];
         buf[m++] = nodex[i][j][2];
	 buf[m++] = nodef[i][j][0];
	 buf[m++] = nodef[i][j][1];
	 buf[m++] = nodef[i][j][2];
	 buf[m++] = nodev[i][j][0];
	 buf[m++] = nodev[i][j][1];
	 buf[m++] = nodev[i][j][2];
       }
     }
   }

 // error->all(FLERR, "test from epack_noscale_noimage in pack_atom");
}

/* ---------------------------------------------------------------------- */

int DumpAtom::convert_noimage(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
                      static_cast<tagint> (mybuf[m]),
                      static_cast<int> (mybuf[m+1]),
                      mybuf[m+2],mybuf[m+3],mybuf[m+4],
		      mybuf[m+5],mybuf[m+6],mybuf[m+7],
		      mybuf[m+8],mybuf[m+9],mybuf[m+10]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_string(int n, double *mybuf)
{
  int m;
  fwrite(mybuf,sizeof(char),n,fp);
  //m = 0;
  //if (me==0)
  //  fprintf(screen, "atom->nlocal = %d\n", atom->nlocal);
  // for (int i=0; i<atom->nlocal; i++)
  //   fprintf(screen, "%d, %d, %f, %f, %f\n", buf[m++], buf[m++], buf[m++], buf[m++], buf[m++]);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpAtom::convert_string(int n, double *mybuf)
{
  return (this->*convert_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            static_cast<tagint> (mybuf[m]), static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6],
	    mybuf[m+7],mybuf[m+8],mybuf[m+9],mybuf[m+10]);
    m += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack(tagint *ids)
{
  if (atom->nlocal) (this->*pack_choice)(ids);
  if (cac->element_flag) epack_noscale_noimage();
}
