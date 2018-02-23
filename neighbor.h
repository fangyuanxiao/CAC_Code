#ifndef CAC_NEIGHBOR_H
#define CAC_NEIGHBOR_H

#include "pointers.h"

namespace CAC_NS {

class Neighbor : protected Pointers {

public:
   int style;                      // 0,1,2 = nsq, bin, multi
   int every;                      // build every this many steps
   int delay;                      // delay build for this many steps
   int dist_check;                 // 0 = always build, 1 = only if 1/2 dist
   int ago;                        // how many steps ago neighboring occurred
   int pgsize;                     // size of neighbor page
   int oneatom;                    // max # of neighbors for one atom
   int includegroup;               // only build pairwise lists for this group
   int build_once;                 // 1 if only build lists once per run
   int cudable;                    // GPU <-> CPU communication flag for CUDA

   double skin;                    // skin distance
   double cutneighmin;             // min neighbor cutoff for all type pairs
   double cutneighmax;             // max neighbor cutoff for all type pairs
   double *cuttype;                // for each type, max neigh cut w/ others
   double cut_all;

   int binsizeflag;                // user-chosen bin size
   double binsize_user;            // set externally by some accalerator pkgs

   bigint ncalls;                  // # of times build has been called
   bigint ndanger;                 // # of dangerous builds
   bigint lastcall;                // timestep of last neighbor::build() call
  
   int nrequest;                    // requests for pairwise neighbor lists
   class NeighRequest **requests;   // from Pair, Fix, Compute, Command classes
   int maxrequest;

   int old_style, old_nrequest;    // previous run info to avoid
   int old_triclinic,old_pgsize;   // re-creation of pairwise neighbor lists
   int old_oneatom, old_every;
   int old_delay,old_check;
   double old_cutoff;

   class NeighRequest **old_requests;   

   int nlist;                       // pairwise neighbor lists
   class NeighList **lists;

   class NodeNeighList *nilist;      // node pairwise neighbor based on element lists

   int nbondlist;                   // list of bonds to compute
   int **bondlist;                
   int nanglelist;                  // list of angles to compute
   int **anglelist;  
   int ndihedrallist;               // list of dihedral to compute
   int **dihedrallist;
   int nimproperlist;               // list of improper to compute
   int **improperlist;
   double cutneighmaxsq;

   Neighbor(class CAC*);
   virtual ~Neighbor();
   virtual void init();
   int decide();                     // decide whether to build or not
   virtual int check_distance();     // check max distance moved since last build
   int request(void *, int instance=0);  // another class requests a neigh list
   void set(int, char **);           // set neighbor style and skin distance
   void modify_params(int, char**);  // modify parameters that control builds
   bigint memory_usage();
   void setup_bins();                // setup bins based on box and cutoff
   void setup_cut();                 // update cut off for changed element size
   int cluster_check;               // 1 if check bond/angle/etc satisfies minimg
   virtual void build(int topoflag=1);  // create all neighbor lists (pair,bond)
   virtual void ebuild(int topoflag=1);  // create element neighbor lists
   virtual void nbuild(int topoflag=1);  // create node neighbor lists

protected:
   int me,nprocs;

   int maxatom;                    // size of atom-based NeighList arrays
   int maxele;
   int maxinp;
   int maxbond,maxangle,maxdihedral,maximproper;   // size of bond lists
   int maxwt;                      // max weighting factor applied + 1

   int must_check;                 // 1 if must check other classes to reneigh
   int restart_check;              // 1 if restart enabled, 0 if no
   int fix_check;                  // # of fixes that induce reneigh
   int *fixchecklist;              // which fixes to check

   double **cutneighsq;            // neighbor cutneigh sq for each type pair
   double **cutneighghostsq;       // neighbor cutnsq for each ghost type pair
   double *cuttypesq;              // cuttype squared
   double cut_allsq;
   double cut_ele;
   double cut_elesq;
   double cut_subele;
   double cut_subelesq;
   double triggersq;               // trigger = build when atom moves this dist             
   int ngrow_flag;

   double **xhold;                 // atoms coords at last neighbor build
   int maxhold;                    // size of xhold array
   double ***nodexhold;            // node coords at last neighbor build
   int emaxhold;                   // size of element array
   int boxcheck;                   // 1 if need to store box size
   double boxlo_hold[3],boxhi_hold[3]; //box size at last neighbor build
   double corners_hold[8][3];          //box corners at last neighbor build

   int binatomflag;                 // bin atoms or not when build neigh list 
                                   // turned off by build_one()
   int nbinx,nbiny,nbinz;           // # of global bins
   int *bins;                       // ptr to next atom in each bin
   int *ebins;                      // ptr to next element in each bin
   int maxbin;                      // size of bins array
   int emaxbin;                     // size of bins array

   int *binhead;                    // ptr to 1st atom in each bin
   int *ebinhead;
   int maxhead;                     // size of binhead array
   
   int mbins;                       // # of local bins and offset
   int mbinx, mbiny, mbinz;
   int mbinxlo,mbinylo,mbinzlo;
   
   double binsizex,binsizey,binsizez;
   double bininvx,bininvy,bininvz;

   int sx,sy,sz,smax;                // bin stencil extents

   int dimension;                    // 2/3 for 2d/3d
   int triclinic;                    // 0 if domain is orthog, 1 if triclinic
   int newton_pair;                  // 0 if newton off, 1 if on for pairwise

   double *bboxlo,*bboxhi;           // ptrs to full domain bounding box
   double (*corners)[3];             // ptr to 8 corners of triclinic box

   double inner[2],middle[2];        // rRESPA cutoffs for extra lists
   double cut_inner_sq;              // outer cutoff for inner neighbor list
   double cut_middle_sq;             // outer cutoff for middle neighbor list
   double cut_middle_inside_sq;      //inner cutoff for middle neighbor list

   int special_flag[4];              // flags for 1-2, 1-3, 1-4 neighbors

   int anyghostlist;                 // 1 if any non-occational list
                                     // stores neighbors of ghosts
                         
   int exclude;                      // 0 if no type/group exclusions, 1 if yes

   int nex_type;                     // # of entries in type exclusion list
   int maxex_type;                   // max # in type list
   int *ex1_type,*ex2_type;          // pairs of types to exclude
   int **ex_type;                    // 2d array of excluded type pairs

   int nex_group;                    // # of entries in group exclusion list
   int maxex_group;                  // max # in molecule list
   int *ex1_group,*ex2_group;        // pairs of group #'s to exclude
   int *ex1_bit,*ex2_bit;            // pairs of group bits to exclude

   int nex_mol;                      // # of entries in molecule exclusion list
   int maxex_mol;                    // max # in molecule list
   int *ex_mol_group;                // molecule group #'s to exclude
   int *ex_mol_bit;                  // molecule group bits to exclude

   int nblist, nglist,nslist;        // # of entries in molecule exclustion list
   int *blist;                       // lists to build every reneighboring
   int *glist;                       // lists to grow atom arrays every reneigh
   int *slist;                       // lists to grow stencil arrays every reneigh

   void bin_atoms();                     // bin all atoms
   void bin_elements();                  // bin all elements
   double bin_distance(int, int, int);   // distance between binx
   int coord2bin(double *);              // mapping atom coord to a bin

   int copymode;

   virtual void choose_build(int, class NeighRequest *);
   void choose_stencil(int, class NeighRequest *); 
   // pairwise build functions   

   typedef void (Neighbor::*PairPtr)(class NeighList *);
   PairPtr *pair_build;
   void copy_from(class NeighList *);
   void half_bin_no_newton(class NeighList *);

   // pairwise stencil creation functions
   void element_neighbor(class NeighList *);
   void node_neighbor(class NodeNeighList *, class NeighList *);
   void check_node_neighbor(class NodeNeighList *);
   void output_ele_neigh(class NeighList*);
  typedef void (Neighbor::*StencilPtr)(class NeighList *, int, int, int);
  StencilPtr *stencil_create;

  void stencil_half_bin_2d_no_newton(class NeighList *, int, int, int);
  void stencil_half_bin_3d_no_newton(class NeighList *, int, int, int);

};

}

#endif
