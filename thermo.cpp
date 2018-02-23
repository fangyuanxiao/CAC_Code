#include "cactype.h" 
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "thermo.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "universe.h"
//#include "lattice.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "input.h"
//#include "variable.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
//#include "bond.h"
//#include "angle.h"
//#include "dihedral.h"
//#include "improper.h"
//#include "kspace.h"
#include "output.h"
#include "timer.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;
using namespace MathConst;

#define ONE "step"
#define MULTI "etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press"

enum{IGNORE,WARN,ERROR};           // same as several files
enum{ONELINE,MULTILINE};
enum{INT,FLOAT,BIGINT};
enum{SCALAR,VECTOR,ARRAY};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

#define DELTA 8

/* ---------------------------------------------------------------------- */

Thermo::Thermo(CAC *cac, int narg, char **arg) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  modified = 0;
  normuserflag = 0;
  lineflag = ONELINE;
  lostflag = lostbond = ERROR;
  lostbefore = 0;
  flushflag = 0;

  if (strcmp(style,"one") == 0) {
    line = new char[256+6*64];
    strcpy(line,ONE);
  } else if (strcmp(style,"custom") == 0) {
    if (narg == 1) error->all(FLERR,"Illegal thermo style custom command");
    line = new char[256+narg*64];
    line[0] = '\0';
    for (int iarg = 1; iarg < narg; iarg++) {
      strcat(line,arg[iarg]);
      strcat(line," ");
    }
    line[strlen(line)-1] = '\0';
  } else error->all(FLERR,"Illegal thermo style command");

  temperature = NULL;
  pressure = NULL;
  pe = NULL;

  index_temp = index_press_scalar = index_press_vector = index_pe = -1;

  id_temp = (char *) "thermo_temp";
  id_press = (char *) "thermo_press";
  id_pe = (char *) "thermo_pe";

  nfield_initial = atom->count_words(line);
  allocate();
  parse_fields(line);

  char *bigint_format = (char *) BIGINT_FORMAT;
  char *fformat_multi = (char *) "---------------- Step %%8%s ----- "
    "CPU = %%11.4f (sec) ----------------";

  sprintf(format_multi,fformat_multi,&bigint_format[1]);
  format_float_one_def = (char *) "%12.8g";
  format_float_multi_def = (char *) "%14.4f";
  format_int_one_def = (char *) "%8d";
  format_int_multi_def = (char *) "%14d";
  sprintf(format_bigint_one_def,"%%8%s",&bigint_format[1]);
  sprintf(format_bigint_multi_def,"%%14%s",&bigint_format[1]);

  format_float_user = NULL;
  format_int_user = NULL;
  format_bigint_user = NULL;
}

/* ---------------------------------------------------------------------- */

Thermo::~Thermo()
{
  delete [] style;
  delete [] line;

  deallocate();

  delete [] format_float_user;
  delete [] format_int_user;
  delete [] format_bigint_user;
}

/* ----------------------------------------------------------------------
    allocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::allocate()
{
  int n = nfield_initial + 1;

  keyword = new char*[n];
  for (int i = 0; i < n; i++) keyword[i] = new char[32];
  vfunc = new FnPtr[n];
  vtype = new int[n]; 

  format = new char*[n];
  for (int i = 0; i < n; i++) format[i] = new char[32];
  format_user = new char*[n];
  for (int i = 0; i < n; i++) format_user[i] = NULL;

  field2index = new int[n];
  argindex1 = new int[n];
  argindex2 = new int[n];

  ncompute = 0;
  id_compute = new char*[3*n];
  compute_which = new int[3*n];
  computes = new Compute*[3*n];

  nfix = 0;
  id_fix = new char*[n];
  fixes = new Fix*[n];  
}

/* ----------------------------------------------------------------------
   deallocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::deallocate()
{
  int n = nfield_initial + 1;

  for (int i = 0; i < n; i++) delete [] keyword[i];
  delete [] keyword;
  delete [] vfunc;
  delete [] vtype;

  for (int i = 0; i < n; i++) delete [] format[i];
  delete [] format;
  for (int i = 0; i < n; i++) delete [] format_user[i];
  delete [] format_user;

  delete [] field2index;
  delete [] argindex1;
  delete [] argindex2;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  delete [] id_compute;
  delete [] compute_which;
  delete [] computes;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  delete [] id_fix;
  delete [] fixes;
}

/* ----------------------------------------------------------------------
   parse list of thermo keywords from str
   set compute flags (temp, press, pe, etc)
------------------------------------------------------------------------- */

void Thermo::parse_fields(char *str)
{
  nfield = 0;

  char *word = strtok(str," \0");

  while (word) {
   
    if (strcmp(word,"step") == 0) {
      addfield("Step",&Thermo::compute_step,BIGINT);
    } else error->all(FLERR,"Unknown keyword in thermo_style custom command");
   
    word = strtok(NULL," \0"); 
  } 
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Thermo::addfield(const char *key, FnPtr func, int typeflag)
{
  strcpy(keyword[nfield],key);
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

void Thermo::compute_step()
{
  bivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Thermo::init()
{
  int i,n;

  // set normvalue to default setting unless user has specified it

  if (normuserflag) normvalue = normuser;
  else if (strcmp(update->unit_style,"lj") == 0) normvalue = 1;
  else normvalue = 0;

  // add Volume field if volume changes and not style = custom
  // this check must come after domain init, so box_change is set

  nfield = nfield_initial;

  // set format string for each field
  // include keyword if lineflag = MULTILINE
  // add '/n' every 3 values if lineflag = MULTILINE
  // add trailing '/n' to last value

  char *ptr;
  for (i = 0; i < nfield; i++) {
   format[i][0] = '\0';
   if (vtype[i] == FLOAT) {
      ptr = format_float_one_def;
    } else if (vtype[i] == INT) {
       ptr = format_int_one_def;
    } else if (vtype[i] == BIGINT) {
      ptr = format_bigint_one_def;
    }
   n = strlen(format[i]);
   sprintf(&format[i][n],"%s ",ptr);
   if (i == nfield-1) strcat(format[i],"\n");    
  } 
  
  // find current ptr for each Compute ID
  // cudable = 0 if any compute used by Thermo is non-CUDA

  cudable = 1;
  //if (me==0) fprintf(screen, "ncompute = %d\n", ncompute);
  int icompute;
  for (i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find thermo compute ID");
    computes[i] = modify->compute[icompute];
    cudable = cudable && computes[i]->cudable;
  }

  // find current ptr for each Fix ID
  // check that fix frequency is acceptable with thermo output frequency
  //if (me==0) fprintf(screen, "nfix = %d\n", nfix);
  int ifix;
  for (i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find thermo fix ID");
    fixes[i] = modify->fix[ifix];
    if (output->thermo_every % fixes[i]->global_freq)
      error->all(FLERR,"Thermo and fix not computed at compatible times");
  }

  // set ptrs to keyword-specific Compute objects

  if (index_temp >= 0) temperature = computes[index_temp];
  if (index_press_scalar >= 0) pressure = computes[index_press_scalar];
  if (index_press_vector >= 0) pressure = computes[index_press_vector];
  if (index_pe >= 0) pe = computes[index_pe];
}

void Thermo::header()
{
  if (lineflag == MULTILINE) return;
  int loc = 0;
  for (int i = 0; i < nfield; i++)
    loc += sprintf(&line[loc],"%s ",keyword[i]);
  sprintf(&line[loc],"\n");
  if (me == 0) {
    if (screen) fprintf(screen,"%s",line);
    if (logfile) fprintf(logfile,"%s",line);
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute(int flag)
{
  int i;
  firststep = flag;
  bigint ntimestep = update->ntimestep;

  int loc = 0;  

  for (ifield = 0; ifield < nfield; ifield++) {
    (this->*vfunc[ifield])();
    if (vtype[ifield] == FLOAT)
      loc += sprintf(&line[loc],format[ifield],dvalue);
    else if (vtype[ifield] == INT)
      loc += sprintf(&line[loc],format[ifield],ivalue);
    else if (vtype[ifield] == BIGINT) {
      loc += sprintf(&line[loc],format[ifield],bivalue);
    }
  }

  if (me == 0) {
    if (screen) fprintf(screen,"%s",line);
    if (logfile) {
      fprintf(logfile,"%s",line);
      if (flushflag) fflush(logfile);
    }
  }
}
