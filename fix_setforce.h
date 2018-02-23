# ifdef FIX_CLASS

 FixStyle(setforce, FixSetForce)

#else

#ifndef CAC_FIX_SET_FORCE_H
#define CAC_FIX_SET_FORCE_H

#include "fix.h"

namespace CAC_NS {

class FixSetForce : public Fix {
public: 
 FixSetForce(class CAC *, int, char **);
 ~FixSetForce();
 int setmask();
 void setup(int);
 void post_force(int);
 void epost_force(int);

private:
 double xvalue,yvalue,zvalue;
 int varflag,iregion;
 char *xstr, *ystr, *zstr;
 char * idregion;
 int xvar,yvar,zvar,xstyle,ystyle,zstyle;
 double foriginal[3],foriginal_all[3],eforiginal[3];
 int force_flag;
 int nlevels_respa;

 int maxatom;
 double **sforce;

};
}
#endif
#endif
