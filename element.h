#ifndef CAC_ELEMENT_H
#define CAC_ELEMENT_H

#include "pointers.h"

namespace CAC_NS {

class Element : protected Pointers {
  public: 
    bigint nnodes;
    bigint nelements;
    bigint n_interpolated;
 
    int nlocal,nghost;
    int ntypes;

    int size_forward;
    int size_reverse;
    int size_border;
    int size_data_element;
    int dimension;
    int map_style;
    tagint map_tag_max;

    double ***nodex, ***nodev, ***nodef;
    tagint **nodetag;
    int ***intp;                       // integration points of each type of element
    double ***intpd;                   // parent coordinate of integration points
    double **weight;                   // weight of each integration point
    double ***sub_ele_pc;              // parent coordinate of sub element
    double **center_subele_pc;         // parent coordinate of the center of the sub element
    int ** list_inter_subele;          // list of inter atoms inside sub element
    int *natom_subele;                 // # of atoms inside the sub element
    int nsplit_ele;
    int nsubele;
    double ***Spa_inter_subele;       // shape function of the sub element
    double **Spa_center_subele;      // shape function of the center of the sub element

    double ***inter_in_ele;            // store the parent coordinate of all 
                                       // the atoms in each type of element
    double ***Spa_inter_ele;           // whole shape function inside each type of elements
    int *inter_num;                     // # of interpolation points in each type of
                                       // element
    int **list_int;                    // integration point # in inter_in_ele          
    int **node_list;                   // node point # in Spa_inter_ele 

    int maxnum_int;                   // store the max # of integration points for all types of element
    int nint_local;                   // # of local integration points 
    double **x;
    double **inpx;
    tagint *tag;
    int    *type, *mask;
    imageint *image;
    double *mass;
    int **plane_index;

    int     xcol_data; 
    int *sametag;
    int     inp_flag;       // flag which determine weather the element need to be 
                            // inerpolated. 1 to be interpolated. 0 not.
    bigint   num_inp;
    int      num_inp_local;
    int      num_inp_ghost;
    int   **nae;            // number of atoms along each direction (important for interpolation of the atoms inside
                            // the element)
    int *nae_setflag;      
    int *numint;            // number of integration points
    double ***Spa_int;                // shape function array
    int totalnumint;        // count total # of int points of all kinds of element
    int npe;                // number of nodes per element dependent on the dimension of the simulations
    int nmax;               // maximum element that the allocated array can store
    int nmaxinp;            // maximum # of integration point 
    double *size;
    double max_size;
    Element(class CAC *);
    ~Element();
    void grow(int);
    void grow_nae(int);   
    void data_elements(int, char *);
    void data_element(double *, imageint, char **);
    void data_nodes(int, char *);
    void data_node(int, char **);
    void data_node_vel(int,char **);
    void set_numint(const char *, int);
    void set_intp(char *, int);
    void set_mass(double);
    inline int map(tagint global) {
      if (map_style == 1) return map_array[global];
      else if (map_style == 2) return map_find_hash(global);
      else return -1;
    };

    void map_init(int check = 1);
    int map_style_set();
    void map_delete();
    void map_set();
    int map_find_hash(tagint);
    void init();
    int pack_exchange(int, double *);
    int unpack_exchange(double *);
    void copy(int, int, int);
    int pack_comm(int, int *, double *, int, int *);
    void unpack_comm(int,int,double *);
    int pack_border(int, int *, double *, int, int *);
    void unpack_border(int, int, double *);
    void allocate_type_arrays();
    void allocate_int_arrays();
    void set_nae(const char *,int);
    void num_atoms_from_element();
    void Interpolate_to_atom_vec();
    void interpolate_atom_vec();
    void count_intp();
    void element_size();
    bigint memory_usage();
  protected:
    int *map_array;
    int map_maxarray;

    struct HashElem {     // hashed map
      tagint global;      // key to search on = global ID
      int local;          // value associated with key = local index
      int next;           // next entry in this bucket, -1 if last
    };    
    int map_nhash;
    int map_nused;
    int map_free;
    int map_nbucket;
    int *map_bucket;
    HashElem *map_hash;
    
    int max_same;
    //int nmax;
    //int nmaxinp;
    void grow_nmax();
    void grow_nmaxinp();
    int next_prime(int);
   // void element_size();
    void intp_pcoord();
    void shape_array_int();
    void shape_array_inter_ele();
    void shape_array_inter_subele();
    void inter_ele();
    void find_node_list();
    void specify_plane();
    void specify_sub_ele_pc(int);
    void count_sub_element();       // count # of atoms inside the sub element 

    //void count_intp();

  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };

};
}

#endif
