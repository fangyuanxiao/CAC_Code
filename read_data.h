#ifdef COMMAND_CLASS

CommandStyle(read_data,ReadData)

#else

#ifndef CAC_READ_DATA_H
#define CAC_READ_DATA_H

#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class ReadData : protected Pointers {
 public:
  ReadData(class CAC *);
  ~ReadData();
  void command(int, char **);

  private:
  int me,compressed;
  char *line,*keyword,*buffer,*style;
  FILE *fp;
  char **arg;
  int narg,maxarg;
  char argoffset1[8],argoffset2[8];

  bigint id_offset;

  bigint natoms;
  bigint nelements;
  bigint nnodes;
  bigint nbonds,nangles,ndihedrals,nimpropers;
  int ntypes;
  int entypes;
  int npe;
  int nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;

  double boxlo[3],boxhi[3];
  double xy,xz,yz;
  int triclinic;

  int addflag,offsetflag,shiftflag;
  tagint addvalue;
  int toffset,boffset,aoffset,doffset,ioffset;
  double shift[3];
  int extra_atom_types,extra_bond_types,extra_angle_types;
  int extra_dihedral_types,extra_improper_types;
  int groupbit;

  int totalnumint;
  int nfix;         
  int *fix_index;
  char **fix_header;
  char **fix_section;

  void open(char *);
  void header();
  void parse_keyword(int);
  void skip_lines(bigint);
  int style_match(const char *, const char *);
  
  void atoms();
  void mass();
  void elements();
  void nodes();
  void naae();
  void numint();
  void intp();
};
}

#endif
#endif

