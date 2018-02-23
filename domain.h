#ifndef CAC_DOMAIN_H
#define CAC_DOMAIN_H

#include "math.h"
#include "pointers.h"

namespace CAC_NS {

class Domain : protected Pointers {
 public:
  int box_exist;                         // 0 = not yet created, 1 = exists
  int dimension;                         // 2 = 2d, 3 = 3d
  int nonperiodic;                       // 0 = periodic in all 3 dims
                                         // 1 = periodic or fixed in all 6
                                         // 2 = shrink-wrap in any of 6
  int xperiodic,yperiodic,zperiodic;     // 0 = non-periodic, 1 = periodic
  int periodicity[3];                    // xyz periodicity as array

  int boundary[3][2];                    // settings for 6 boundaries
                                         // 0 = periodic
                                         // 1 = fixed non-periodic
                                         // 2 = shrink-wrap non-periodic
                                         // 3 = shrink-wrap non-per w/ min

  int triclinic;                         // 0 = orthog box, 1 = triclinic
  int tiltsmall;                         // 1 if limit tilt, else 0

                                         // orthogonal box
  double xprd,yprd,zprd;                 // global box dimensions
  double xprd_half,yprd_half,zprd_half;  // half dimensions
  double prd[3];                         // array form of dimensions
  double prd_half[3];                    // array form of half dimensions

                                         // triclinic box
                                         // xprd,xprd_half,prd,prd_half =
                                         // same as if untilted
  double prd_lamda[3];                   // lamda box = (1,1,1)
  double prd_half_lamda[3];              // lamda half box = (0.5,0.5,0.5)

  double boxlo[3],boxhi[3];              // orthogonal box global bounds

                                         // triclinic box
                                         // boxlo/hi = same as if untilted
  double boxlo_lamda[3],boxhi_lamda[3];  // lamda box = (0,1)
  double boxlo_bound[3],boxhi_bound[3];  // bounding box of tilted domain
  double corners[8][3];                  // 8 corner points

                                         // orthogonal box & triclinic box
  double minxlo,minxhi;                  // minimum size of global box
  double minylo,minyhi;                  // when shrink-wrapping
  double minzlo,minzhi;                  // tri only possible for non-skew dims

                                         // orthogonal box
  double sublo[3],subhi[3];              // sub-box bounds on this proc

                                         // triclinic box
                                         // sublo/hi = undefined
  double sublo_lamda[3],subhi_lamda[3];  // bounds of subbox in lamda

                                         // triclinic box
  double xy,xz,yz;                       // 3 tilt factors
  double h[6],h_inv[6];                  // shape matrix in Voigt notation
  double h_rate[6],h_ratelo[3];          // rate of box size/shape change

  int box_change;                // 1 if any of next 3 flags are set, else 0
  int box_change_size;           // 1 if box size changes, 0 if not
  int box_change_shape;          // 1 if box shape changes, 0 if not
  int box_change_domain;         // 1 if proc sub-domains change, 0 if not

  int deform_flag;                // 1 if fix deform exist, else 0
  int deform_vremap;              // 1 if fix deform remaps v, else 0
  int deform_groupbit;            // atom group to perform v remap for

//  class Lattice *lattice;                  // user-defined lattice

  int nregion;                             // # of defined Regions
  int maxregion;                           // max # list can hold
  class Region **regions;                  // list of defined Regions

  Domain(class CAC *);
  virtual ~Domain();
  virtual void init();
  void set_initial_box(int expandflag=1);
  virtual void set_global_box();
  virtual void set_local_box();
  virtual void reset_box();
  virtual void pbc();
  virtual void elepbc();
  void set_boundary(int, char **, int);
  void print_box(const char *);
  void remap(double *, imageint &);
  void box_too_small_check();
  void subbox_too_small_check(double);
  void boundary_string(char *);
  void image_check();
  virtual void lamda2x(double *, double *);
  virtual void x2lamda(double *, double *);
  void add_region(int,char **);
  int find_region(char *); 

 protected:
  double small[3];                  // fractions of box lengths
};

}

#endif

