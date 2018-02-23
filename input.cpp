#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "sys/stat.h"
#include "input.h"
#include "style_command.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "comm_brick.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "compute.h"
#include "update.h"
#include "neighbor.h"
#include "universe.h"
#include "error.h"
#include "memory.h"
#include "element.h"
using namespace CAC_NS;

#define DELTALINE 256
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(CAC *cac, int argc, char **argv) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);

  maxline = maxcopy = maxwork = 0;
  line = copy = work = NULL;
  narg = maxarg = 0;
  arg = NULL;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = NULL;
  jump_skip = 0;
  ifthenelse_flag = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *),"input:infiles");
    infiles[0] = infile;
  } else infiles = NULL;

   command_map = new std::map<std::string,CommandCreator>();

   #define COMMAND_CLASS
#define CommandStyle(key,Class) \
  (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS

}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  memory->sfree(line);
  memory->sfree(copy);
  memory->sfree(work);
  if (labelstr) delete [] labelstr;
  memory->sfree(arg);
  memory->sfree(infiles);

}

/* ---------------------------------------------------------------------- */

template <typename T>
void Input::command_creator(CAC *cac, int narg, char **arg)
{
  T cmd(cac);
  cmd.command(narg,arg);
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int m,n;
   
  while (1) {

  // read a line from input script
  // n = length of line including str terminator, 0 if end of file
  // if line ends in continuation char '&', concatenate next line

  if (me == 0) {
     m = 0;
     while (1) {
       if (maxline-m < 2) reallocate(line,maxline,0);

  // end of file reached, so break
  //  n == 0 if nothing read, else n = line with str terminator

      if (fgets(&line[m],maxline-m,infile) == NULL) {
          if (m) n = strlen(line) + 1;
          else n = 0;
          break;
        }

   // continue if last char read was not a newline
   // could happen if line is very long

       m = strlen(line);
       if (line[m-1] != '\n') continue;
    
   // continue reading if final printable char is & char
   // or if odd number of triple quotes
   // else break with n = line with str terminator

      m--;
      while (m >= 0 && isspace(line[m])) m--;
      if (m < 0 || line[m] != '&') {
     if (numtriple(line) % 2) {
        m += 2;
	continue;
     }
         line[m+1] = '\0';
          n = m+2;
          break;
     }
    } 
   }
    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
       if (me == 0) {
	//fprintf(screen,"nfile = %d\n",nfile);
        if (infile != stdin) {
          fclose(infile);
          infile = NULL;
        }
        nfile--;
	//fprintf(screen, "nfile = %d\n",nfile);
      }
      MPI_Bcast(&nfile,1,MPI_INT,0,world);
      //fprintf(screen, "me = %d,nfile = %d\n",me,nfile);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }
    
    if (n > maxline) reallocate(line,maxline,n);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s\n",line);
      if (echo_log && logfile) fprintf(logfile,"%s\n",line);
    }
   
   // parse the line
   // if no command, skip to next line in input script
   
   parse();
   if (command == NULL) continue;
//   if (me == 0) fprintf(logfile,"%s\n",command);
   
   // execute the command

   if (execute_command()) {
      char *str = new char[maxline+32];
      sprintf(str,"Unknown command: %s",line);
      error->all(FLERR,str);
    }
    
  }

  //error->all(FLERR, "test from file in input");
}

/* ----------------------------------------------------------------------
  parse copy of command line by inserting string terminators
  strip comment = all chars from # on
  replace all $ via variable substitution except within quotes
  command = first word
  narg = # of args
  arg[] = individual args
  treat text between single/double/triple quotes as one arg via nextword()
 ------------------------------------------------------------------------- */

void Input::parse()
{
   // duplicate line into copy string to break into words

  int n = strlen(line) + 1;
  if (n > maxcopy) reallocate(copy,maxcopy,n);
  strcpy(copy,line);

  // strip any # comment by replacing it with 0
  // do not strip from a # inside single/double/triple quotes
  // quoteflag = 1,2,3 when encounter first single/double,triple quote
  // quoteflag = 0 when encounter matching single/double,triple quote

  int quoteflag = 0;
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && !quoteflag) {
      *ptr = '\0';
      break;
    }
   if (quoteflag == 0) {
      if (strstr(ptr,"\"\"\"") == ptr) {
	quoteflag = 3;
	ptr += 2;
      }
      else if (*ptr == '"') quoteflag = 2;
      else if (*ptr == '\'') quoteflag = 1;
    } else {
      if (quoteflag == 3 && strstr(ptr,"\"\"\"") == ptr) {
	quoteflag = 0;
	ptr += 2;
      }
      else if (quoteflag == 2 && *ptr == '"') quoteflag = 0;
      else if (quoteflag == 1 && *ptr == '\'') quoteflag = 0;
    }
   ptr++;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  // command = 1st arg in copy string
  
  char *next;
  command = nextword(copy,&next);
  if (command == NULL) return;

  // point arg[] at each subsequent arg in copy string
  // nextword() inserts string terminators into copy string to delimit args
  // nextword() treats text between single/double/triple quotes as one arg
  
  narg = 0;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = nextword(ptr,&next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
}

/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   treat text between single/double/triple quotes as one arg
   matching quote must be followed by whitespace char if not end of string
   strip quotes from returned word
   return ptr to start of word or NULL if no word in string
   also return next = ptr after word
------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start,*stop;

  // start = first non-whitespace char
  
  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return NULL;

  // if start is single/double/triple quote:
  // start = first char beyond quote
  // stop = first char of matching quote
  // next = first char beyond matching quote
  // next must be NULL or whitespace
  // if start is not single/double/triple quote:
  // stop = first whitespace char after start
  // next = char after stop, or stop itself if stop is NULL
  
  if (strstr(start,"\"\"\"") == start) {
    stop = strstr(&start[3],"\"\"\"");
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start += 3;
    *next = stop+3;
    if (**next && !isspace(**next))
      error->all(FLERR,"Input line quote not followed by whitespace");
  } else if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1],*start);
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start++;
    *next = stop+1;
    if (**next && !isspace(**next))
      error->all(FLERR,"Input line quote not followed by whitespace");
  } else {
    stop = &start[strcspn(start," \t\n\v\f\r")];
    if (*stop == '\0') *next = stop;
    else *next = stop+1;
  }
  
  // set stop to NULL to terminate word

  *stop = '\0';
  return start;
}

/* ----------------------------------------------------------------------
  return number of triple quotes in line
------------------------------------------------------------------------- */

int Input::numtriple(char *line)
{
  int count = 0;
  char *ptr = line;
  while ((ptr = strstr(ptr,"\"\"\""))) {
    ptr += 3;
    count++;
  }
  return count;
}


/* ----------------------------------------------------------------------
  rellocate a string
  if n > 0: set max >= n in increments of DELTALINE
  if n = 0: just increment max by DELTALINE
------------------------------------------------------------------------- */

void Input::reallocate(char *&str, int &max, int n)
{
  if (n) {
    while (n > max) max += DELTALINE;
  } else max += DELTALINE;
  
  str = (char *) memory->srealloc(str,max*sizeof(char),"input:str");
}

/* ----------------------------------------------------------------------
  process a single parsed command
  return 0 if successful, -1 if did not recognize command
 ------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  if (!strcmp(command,"clear")) clear();
  else if (!strcmp(command,"atom_style")) atom_style();
  else if (!strcmp(command,"units")) units();
  else if (!strcmp(command,"dimension")) dimension();
  else if (!strcmp(command,"boundary")) boundary();
  else if (!strcmp(command,"newton")) newton();
  else if (!strcmp(command,"neighbor")) neighbor_command();
  else if (!strcmp(command,"pair_coeff")) pair_coeff();
  else if (!strcmp(command,"pair_style")) pair_style();
  else if (!strcmp(command,"atom_style")) atom_style(); 
  else if (!strcmp(command,"neigh_modify")) neigh_modify();
  else if (!strcmp(command,"fix")) fix();
  else if (!strcmp(command,"region")) region();
  else if (!strcmp(command,"group")) group_command();
  else if (!strcmp(command,"thermo")) thermo();
  else if (!strcmp(command,"thermo_style")) thermo_style();
  else if (!strcmp(command,"dump")) dump();
  else if (!strcmp(command,"timestep")) timestep();
  else if (!strcmp(command,"processors")) processors();
  else if (!strcmp(command,"elementsplit")) element_split();

  else flag = 0;

  // return if command was listed above
  
  if (flag) return 0;

  // invoke commands added via style_command.h

  if (command_map->find(command) != command_map->end()) {
    CommandCreator command_creator = (*command_map)[command];
    command_creator(cac,narg,arg);
    return 0;
  }

  // unrecognized command

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */



void Input::region()
{
   domain->add_region(narg,arg);
   //error->all(FLERR, "test from region of input");
}

/*------------------------------------------------------------------------*/

void Input::group_command()
{
  group->assign(narg,arg);
}	

/*------------------------------------------------------------------------*/

void Input::atom_style()
{
  if (narg < 1) error->all(FLERR,"Illegal atom_style command");
  if (domain->box_exist)
    error->all(FLERR,"Atom_style command after simulation box is defined");
  atom->create_avec(arg[0],narg-1,&arg[1],1);
}

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all(FLERR,"Illegal clear command");
// if (me==0) fprintf(logfile,"clear \n");
  cac->destroy();
  cac->create();
//  lmp->post_create(0,NULL,NULL,NULL);
}

/* ---------------------------------------------------------------------- */

void Input::units()
{
  if (narg != 1) error->all(FLERR,"Illegal units command");
  if (domain->box_exist)
    error->all(FLERR,"Units command after simulation box is defined");
  update->set_units(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (narg != 1) error->all(FLERR,"Illegal dimension command");
  if (domain->box_exist)
    error->all(FLERR,"Dimension command after simulation box is defined");
  domain->dimension = universe->inumeric(FLERR,arg[0]);
  if (domain->dimension != 2 && domain->dimension != 3)
    error->all(FLERR,"Illegal dimension command");

  // must reset default extra_dof of all computes
  // since some were created before dimension command is encounteredi

  for (int i = 0; i < modify->ncompute; i++)
    modify->compute[i]->reset_extra_dof();
  
  if (domain->dimension == 2) {
     element->dimension = domain->dimension;
     element->npe = 4; 
     element->size_border = 18;
     element->size_forward = 12;
  }

  
  // if (me==0) fprintf(screen, "npe = %d\n", element->npe);
//   if (me==0) fprintf(logfile,"dimension \n");
 // error->all(FLERR, "testing on setting the dimension");
}

/*-----------------------------------------------------------------------*/

void Input::element_split()
{
  if (narg != 1) error->all(FLERR, "Illegal element_split command");
  element->nsplit_ele = universe->inumeric(FLERR,arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  if (domain->box_exist)
    error->all(FLERR,"Boundary command after simulation box is defined");
  domain->set_boundary(narg,arg,0);

   if (me==0) fprintf(logfile,"boundary \n");
}

/* ---------------------------------------------------------------------- */

void Input::newton()
{
   int newton_pair=1,newton_bond=1;

  if (narg == 1) {
    if (strcmp(arg[0],"off") == 0) newton_pair = newton_bond = 0;
    else if (strcmp(arg[0],"on") == 0) newton_pair = newton_bond = 1;
    else error->all(FLERR,"Illegal newton command");
  } else if (narg == 2) {
    if (strcmp(arg[0],"off") == 0) newton_pair = 0;
    else if (strcmp(arg[0],"on") == 0) newton_pair= 1;
    else error->all(FLERR,"Illegal newton command");
    if (strcmp(arg[1],"off") == 0) newton_bond = 0;
    else if (strcmp(arg[1],"on") == 0) newton_bond = 1;
    else error->all(FLERR,"Illegal newton command");
  } else error->all(FLERR,"Illegal newton command");

  force->newton_pair = newton_pair;

  if (domain->box_exist && (newton_bond != force->newton_bond))
    error->all(FLERR,"Newton bond change after simulation box is defined");
  force->newton_bond = newton_bond;

  if (newton_pair || newton_bond) force->newton = 1;
  else force->newton = 0;

}

/* ---------------------------------------------------------------------- */

void Input::neighbor_command()
{
  neighbor->set(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neigh_modify()
{
  neighbor->modify_params(narg,arg);
}

/* ----------------------------------------------------------------------
   if old pair style exists and new style is same, just change settings
   else create new pair class
------------------------------------------------------------------------- */

void Input::pair_style()
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");
  if (force->pair) {
    int match = 0;
    if (strcmp(arg[0],force->pair_style) == 0) match = 1;
    if (!match && cac->suffix_enable) {
      char estyle[256];
      if (cac->suffix) {
        sprintf(estyle,"%s/%s",arg[0],cac->suffix);
        if (strcmp(estyle,force->pair_style) == 0) match = 1;
      }
      if (cac->suffix2) {
        sprintf(estyle,"%s/%s",arg[0],cac->suffix2);
        if (strcmp(estyle,force->pair_style) == 0) match = 1;
      }
    }
    if (match) {
      force->pair->settings(narg-1,&arg[1]);
      return;
    }
  }

  force->create_pair(arg[0],1);
  if (force->pair) force->pair->settings(narg-1,&arg[1]);  
}

/* ---------------------------------------------------------------------- */

void Input::pair_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Pair_coeff command before simulation box is defined");
  if (force->pair == NULL)
    error->all(FLERR,"Pair_coeff command before pair_style is defined");
  force->pair->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{	
  modify->add_fix(narg,arg,1);
}

/* ---------------------------------------------------------------------- */

void Input::thermo()
{
  output->set_thermo(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_style()
{
  output->create_thermo(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  output->add_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (narg != 1) error->all(FLERR,"Illegal timestep command");
  update->dt = universe->numeric(FLERR,arg[0]);
}

/*-----------------------------------------------------------------------*/

void Input::processors()
{
  if (domain->box_exist)
    error->all(FLERR,"Processors command after simulation box is defined");
  comm->set_processors(narg,arg);
}
