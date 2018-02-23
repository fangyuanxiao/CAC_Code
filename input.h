#ifndef CAC_INPUT_H
#define CAC_INPUT_H

#include "stdio.h"
#include "pointers.h"
#include <map>
#include <string>

namespace CAC_NS {

class Input : protected Pointers {
public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  class Variable *variable;    // defined variables

  Input(class CAC *, int, char **);
  ~Input();
  void file();

private:
  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line,*copy,*work;      // input line & copy and work string
  int maxline,maxcopy,maxwork; // max lengths of char strings
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile,maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise
  int ifthenelse_flag;         // 1 if executing commands inside an if-then-else

  FILE **infiles;              // list of open input files

  typedef void (*CommandCreator)(CAC *, int, char **);
  std::map<std::string,CommandCreator> *command_map;

  template <typename T> static void command_creator(CAC *, int, char **);
 
  void parse();                          // parse an input text line
  char *nextword(char *, char **);       // find next word in string with quotes
  int numtriple(char *);                 // count number of triple quotes
  void reallocate(char *&, int &, int);  // reallocate a char string
  int execute_command();                 // execute a single command

  void clear();                // input script commands
  void dimension();
  void units();
  void boundary();
  void newton(); 
  void atom_style();
  void neighbor_command();
  void neigh_modify();
  void pair_style();
  void pair_coeff();
  void fix();
  void thermo();
  void thermo_style();
  void dump();
  void timestep();
  void processors();
  void region();
  void group_command();
  void element_split();

};

}

#endif
