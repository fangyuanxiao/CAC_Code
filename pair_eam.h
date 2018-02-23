#ifdef PAIR_CLASS

PairStyle(eam,PairEAM)

#else

#ifndef CAC_PAIR_EAM_H
#define CAC_PAIR_EAM_H

#include "stdio.h"
#include "pair.h"

namespace CAC_NS {

class PairEAM : public Pair {
  public:

   double cutmax;

   // potentials as array data

   int nrho, nr;
   int nfrho, nrhor,nz2r;
   double **frho,**rhor,**z2r;
   int *type2frho,**type2rhor,**type2z2r;

   // potentials in spline form used for force computation
   double dr, rdr, drho, rdrho, rhomax;
   double ***rhor_spline, ***frho_spline, ***z2r_spline;

   PairEAM(class CAC*);
   virtual ~PairEAM();
   void settings(int, char **);
   virtual void coeff(int, char **);
   virtual void compute(int,int);
   virtual void ecompute();
   void init_style();
   double init_one(int,int);

   virtual int pack_forward_comm(int, int *, double *, int, int *);
   virtual int epack_forward_comm(int,int *, double *);
   virtual void unpack_forward_comm(int, int, double *);
   virtual void eunpack_forward_comm(int, int, double *);

 protected:
  int nmax;                        // allocated size of per-atom arrays
  int enmax;
  double cutforcesq;

  // per-atom arrays

  double *rho,*fp;
  double **erho;

  // potentials as file data

  int *map;                        // which element each atom type maps to

  struct Funcfl {
    char *file;
    int nrho,nr;
    double drho,dr,cut,mass;
    double *frho,*rhor,*zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;	
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,**rhor,***z2r;
  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,***rhor,***z2r;
  };
  Fs *fs;

  virtual void allocate();
  virtual void array2spline();
  void interpolate(int,double,double *, double **);
  void grab(FILE *, int, double *);

  virtual void read_file(char *);
  virtual void file2array();

};
}

#endif
#endif
