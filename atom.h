#ifndef CAC_ATOM_H
#define CAC_ATOM_H

#include "pointers.h"

namespace CAC_NS {

class Atom : protected Pointers {
public:
   char *atom_style;
   class AtomVec *avec;
   
   //atom counts

   bigint natoms;                     //total # of atoms in systems, could be 0
                                      //natoms may not be current if atoms lost

   int nlocal,nghost;                 //# of owned and ghost atoms on this proc
   int nmax;                          // max # of owned+ghost in arrays on this proc
   int tag_enable;                    // 0/1 if atom ID tags are defined
   int molecular;                     // 0 = atomic, 1 = standard molecular system,
                                      // 2 = molecule template system

   bigint nbonds,nangles,ndihedrals,nimpropers;
   int ntypes,nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;
   int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;
   int extra_bond_per_atom,extra_angle_per_atom;
   int extra_dihedral_per_atom,extra_improper_per_atom;

   int firstgroup;               // store atoms in this group first, -1 if unset
   int nfirst;                   // # of atoms in first group on this proc
   char *firstgroupname;         // group-ID to store first, NULL if unset

   // per-atom arrays
   // customize by adding new array
   
   tagint *tag;
   int *type,*mask;
   imageint *image;
   double **x,**v,**f;

   tagint *molecule;
   int *molindex,*molatom;

   double *q,**mu;
   double **omega,**angmom,**torque;
   double *radius,*rmass;
   int *ellipsoid,*line,*tri,*body;

   // PERI package
  
   double *vfrac,*s0;
   double **x0;

   // USER-EFF and USER-AWPMD packages
  
   int *spin;
   double *eradius,*ervel,*erforce,*ervelforce;
   double *cs,*csforce,*vforce;
   int *etag;

   // USER-SPH package

   double *rho,*drho,*e,*de,*cv;
   double **vest;

   // USER-SMD package
   
   double *contact_radius;
   double **smd_data_9;
   double **smd_stress;
   double *eff_plastic_strain;
   double *eff_plastic_strain_rate;
   double *damage;

  // molecular info
  
  int **nspecial;                 // 0,1,2 = cummulative # of 1-2, 1-3, 1-4 neighs
  tagint **special;               // IDs of 1-2,1-3,1-4 neighs of each atom
  int maxspecial;                 // special[nlocal][maxspecial]

  int *num_bond;
  int **bond_type;
  tagint **bond_atom;

  int *num_angle;
  int **angle_type;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;

  int *num_improper;
  int **improper_type;
  tagint **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;

  // custom arrays used by fix property/atom

  int **ivector;
  double **dvector;
  char **iname,**dname;
  int nivector,ndvector;

  // used by USER-CUDA to flag used per-atom arrays
  
  unsigned int datamask;
  unsigned int datamask_ext;

  // atom style and per-atom array existence flags
  // customize by adding new flag

  int sphere_flag,ellipsoid_flag,line_flag,tri_flag,body_flag;
  int peri_flag,electron_flag;
  int ecp_flag;
  int wavepacket_flag,sph_flag;

  int molecule_flag,molindex_flag,molatom_flag;
  int q_flag,mu_flag;
  int rmass_flag,radius_flag,omega_flag,torque_flag,angmom_flag;
  int vfrac_flag,spin_flag,eradius_flag,ervel_flag,erforce_flag;
  int cs_flag,csforce_flag,vforce_flag,ervelforce_flag,etag_flag;
  int rho_flag,e_flag,cv_flag,vest_flag;

  // USER-SMD package

  int smd_flag;
  int contact_radius_flag;
  int smd_data_9_flag;
  int smd_stress_flag;
  int x0_flag;
  int eff_plastic_strain_flag;
  int eff_plastic_strain_rate_flag;
  int damage_flag;

  // Peridynamics scale factor, used by dump cfg

  double pdscale;

  // molecule templates
  // each template can be a set of consecutive molecules
  // each with same ID (stored in molecules)
  // 1st molecule in template stores nset = # in set

  int nmolecule;
//  class Molecule **molecules;

  // extra peratom info in restart file destined for fix & diag

  double **extra;

  // per-type arrays

  double *mass;
  int *mass_setflag;

  // callback ptrs for atom arrays managed by fix classes
  
  int nextra_grow,nextra_restart,nextra_border;  // # of callbacks of each type
  int *extra_grow,*extra_restart,*extra_border;  // index of fix to callback to
  int nextra_grow_max,nextra_restart_max;        // size of callback lists
  int nextra_border_max;
  int nextra_store;

  int map_style;                  // style of atom map: 0=none, 1=array, 2=hash
  int map_user;                   // user selected style = same 0,1,2
  tagint map_tag_max;             // max atom ID that map() is setup for

  // spatial sorting of atoms

  int sortfreq;             // sort atoms every this many steps, 0 = off
  bigint nextsort;          // next timestep to sort on
  double userbinsize;       // requested sort bin size

  // indices of atoms with same ID

  int *sametag;              // sametag[I] = next atom with same ID, -1 if no more

  // functions

   Atom(class CAC *);
  ~Atom();

  void update_callback(int);

  void create_avec(const char *, int, char **, int);
  class AtomVec *new_avec(const char *, int, int &);
  void init();
  void setup();
  int count_words(const char *);

  void deallocate_topology();
  virtual void allocate_type_arrays();
  void data_atoms(int, char *, tagint, int, int, double *);
  void set_mass(const char *, int);
  void set_mass(int,double);
  void check_mass();
  virtual void sort();

  bigint memory_usage();
  int memcheck(const char *);
protected:

  // global to local ID mapping
 
  int *map_array;       // direct map via array that holds map_tag_max
  int map_maxarray;     // allocated size of map_array (1 larger than this)

  struct HashElem {     // hashed map
    tagint global;      // key to search on = global ID
    int local;          // value associated with key = local index
    int next;           // next entry in this bucket, -1 if last
  };
  int map_nhash;        // # of entries hash table can hold
  int map_nused;        // # of actual entries in hash table
  int map_free;         // ptr to 1st unused entry in hash table
  int map_nbucket;      // # of hash buckets
  int *map_bucket;      // ptr to 1st entry in each bucket
  HashElem *map_hash;   // hash table

  int max_same;         // allocated size of sametag

   // spatial sorting of atoms

  int nbins;                      // # of sorting bins
  int nbinx,nbiny,nbinz;          // bins in each dimension
  int maxbin;                     // max # of bins
  int maxnext;                    // max size of next,permute
  int *binhead;                   // 1st atom in each bin
  int *next;                      // next atom in bin
  int *permute;                   // permutation vector
  double bininvx,bininvy,bininvz; // inverse actual bin sizes
  double bboxlo[3],bboxhi[3];     // bounding box of my sub-domain

  int memlength;                  // allocated size of memstr
  char *memstr;                   // string of array names already counted

  void setup_sort_bins();
};

}

#endif
