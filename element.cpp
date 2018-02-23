#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "element.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec_atomic.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define NDELTA 4096
#define DELTA 1024
#define EPSILON 1.0e-6
#define EXTRA 1000
#define TOL 0.01
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};

Element::Element(CAC *cac) : Pointers(cac)
{
 dimension = 3;
 inp_flag = 0;
 xcol_data = 3;
 nmax = 0;
 nmaxinp = 0;
 nelements = 0;
 nnodes = 0;
 ntypes = 0; 
 nlocal = 0;
 nghost = 0;
 nsplit_ele = 2;

 nodex = NULL;
 nodev = NULL;
 nodef = NULL;
 nodetag = NULL;
 // 3d 8 node element default 
 npe = 8;
 //shape function of integration points
 Spa_int = NULL; 
 //shape function of all atoms in each type of elements
 Spa_inter_ele = NULL; 
 // shape function of atoms in each sub element
 Spa_inter_subele = NULL;
 center_subele_pc = NULL;
 Spa_center_subele = NULL;
 list_inter_subele = NULL;

 x = NULL;
 tag = NULL;
 type = NULL;
 mask = NULL;
 num_inp = 0;
 inpx = NULL;
 
 //ntype-length arrays
 
 nae = NULL;
 nae_setflag = NULL;
 mass = NULL;
 numint = NULL;
 intp = NULL;
 intpd = NULL;
 sub_ele_pc = NULL;
 natom_subele = NULL;
 weight = NULL;
 totalnumint = 0;
 nint_local = 0;
 

 size = NULL;
 max_size = 0.0;
 max_same = 0;
 map_array = NULL;
 map_bucket = NULL;
 map_hash = NULL;
 sametag = NULL;
 map_style = 0;
 map_tag_max = -1;
 map_maxarray = map_nhash = -1; 
 
 size_data_element = 5;
 size_forward = 24;
 size_reverse = 3;
 size_border = 30;
 
 inter_in_ele = NULL;
 inter_num = NULL;
 list_int = NULL; 
 node_list = NULL;
 plane_index = NULL;
}

Element::~Element()
{
  memory->destroy(nodex);
  memory->destroy(nodev);
  memory->destroy(nodef);
  memory->destroy(nodetag);
  memory->destroy(x);
  memory->destroy(inpx);
  memory->destroy(tag);
  memory->destroy(mask);
  memory->destroy(type);
  memory->destroy(size);
  memory->destroy(nae);
  memory->destroy(plane_index);
  memory->destroy(sub_ele_pc);
  memory->destroy(natom_subele);
  memory->destroy(center_subele_pc);
  memory->destroy(Spa_center_subele);
  delete [] nae_setflag;
  delete [] numint;
  delete [] mass;
  memory->destroy(intp);
  memory->destroy(intpd);
  memory->destroy(weight);
  map_delete();

  for (int ii=0;ii<ntypes;ii++) 
    memory->destroy(Spa_int[ii]);
  delete [] Spa_int;

  for (int ii=0;ii<ntypes;ii++)
    memory->destroy(list_int[ii]);
  delete [] list_int;

  memory->destroy(node_list);

  for (int ii=0;ii<ntypes;ii++)
    memory->destroy(inter_in_ele[ii]);
  delete [] inter_in_ele;
 
  for (int ii=0; ii<ntypes;ii++)
    memory->destroy(Spa_inter_ele[ii]);
  delete [] Spa_inter_ele;

  for (int ii=0; ii < nsubele;ii++)
    memory->destroy(Spa_inter_subele[ii]);
  delete [] Spa_inter_subele;

  for (int ii=0; ii<nsubele;ii++)
    memory->destroy(list_inter_subele[ii]);
  delete [] list_inter_subele;

  delete [] inter_num;
}

/* ----------------------------------------------------------------------
   grow element arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n (learn from the atom_vec_atomic.cpp)
-------------------------------------------------------------------------*/

void Element::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");
//  if (comm->me==0) fprintf(screen,"npe in grow is %d\n",npe);
//  error->all(FLERR, "error from grow in element\n"); 
  memory->grow(tag,nmax,"element:tag");
  memory->grow(type,nmax,"element:type");
  memory->grow(mask,nmax,"element:mask");
  memory->grow(image,nmax,"element:image");
  memory->grow(x,nmax,3,"element:x");
  memory->grow(nodex,nmax,npe,3,"node:x");
  memory->grow(nodev,nmax,npe,3,"node:v");
  memory->grow(nodef,nmax,npe,3,"node:f");
  memory->grow(nodetag, nmax, npe, "node:tag");
}

/* ----------------------------------------------------------------------
   grow nmax so it is a multiple of DELTA
------------------------------------------------------------------------- */

void Element::grow_nmax()
{
  nmax = nmax/DELTA * DELTA;
  nmax += DELTA;
}

/* ----------------------------------------------------------------------
   unpack 9*n lines from Element section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Element::data_elements(int n, char *buf)
{
  int m,xptr,iptr;
  imageint imagedata;
  double xdata[3],lamda[3];
  double *coord;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = atom->count_words(buf);
  *next = '\n';

  if (nwords != size_data_element && nwords != size_data_element+3)
    error->all(FLERR, "Incorrect element format in data file");
  
  char **values = new char*[nwords];

  double epsilon[3];

  epsilon[0] = domain->prd[0] * EPSILON;
  epsilon[1] = domain->prd[1] * EPSILON;
  epsilon[2] = domain->prd[2] * EPSILON;

  double sublo[3],subhi[3];

  sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
  sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
  sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];

  if (comm->layout != LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }
  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // xptr = which word in line starts xyz coords of element
  // iptr = which word in line starts ix, iy, iz image flags

  xptr = xcol_data - 1;
  int imageflag = 0;
  if (nwords > size_data_element) imageflag = 1;
  if (imageflag) iptr = nwords - 3; 

  
  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    values[0] = strtok(buf," \t\n\r\f");

    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL," \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR,"Incorrect atom format in data file");
    }
    
    if (imageflag)
      imagedata = ((imageint) (atoi(values[iptr]) + IMGMAX) & IMGMASK) |
        (((imageint) (atoi(values[iptr+1]) + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (atoi(values[iptr+2]) + IMGMAX) & IMGMASK) << IMG2BITS);
    else imagedata = ((imageint) IMGMAX << IMG2BITS) |
           ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    xdata[0] = atof(values[xptr]);
    xdata[1] = atof(values[xptr+1]);
    xdata[2] = atof(values[xptr+2]);
    
    domain->remap(xdata,imagedata);
    coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
       data_element(xdata,imagedata,values);
    }
    buf = next+1; 
  }
  delete [] values;
}

void Element::data_element(double *coord, imageint imagetmp, char **values)
{
  if (nlocal == nmax) grow(0);  
 
  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  //fprintf(screen, "me = %d, tag = %d\n",comm->me, tag[nlocal]);
  if (type[nlocal] <= 0 || type[nlocal] > ntypes)
    error->one(FLERR,"Invalid atom type in Elements section of data file"); 

  image[nlocal] = imagetmp;
  // 1 equal to the all group, see data_atom
  mask[nlocal] = 1;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  
  nlocal++;
}

void Element::map_init(int check)
{
  int recreate = 0;
  if (check) recreate = map_style_set();
  
  if (map_style == 1 && map_tag_max > map_maxarray) recreate = 1;
  else if (map_style == 2 && nlocal+nghost > map_nhash) recreate = 1;
  map_delete();
  //if (comm->me ==0) fprintf(screen, "recreate = %d\n",recreate);

  if (map_style == 1)
  {
    map_maxarray = map_tag_max;
    memory->create(map_array,map_maxarray+1,"atom:map_array");
    for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;
    //if (comm->me==0) fprintf(screen, "testing in creating map_array\n");
  } else {
    int nper = static_cast<int> (nelements/comm->nprocs);
    map_nhash = MAX(nper,nmax);
    map_nhash *= 2;
    map_nhash = MAX(map_nhash,1000);

    map_nbucket = next_prime(map_nhash);

    map_bucket = new int[map_nbucket];
    for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;

    map_hash = new HashElem[map_nhash];
    map_nused = 0;
    map_free = 0;
    for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
    map_hash[map_nhash-1].next = -1; 
  }
}

int Element::map_style_set()
{
  tagint max = -1;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&map_tag_max,1,MPI_LMP_TAGINT,MPI_MAX,world);
  //if (comm->me==0) fprintf(screen, "map_tag_max = %d\n", map_tag_max);
  int map_style_old = map_style; 
  if (map_tag_max > 1000000) map_style = 2;
  else map_style = 1;
  //if (comm->me==0) fprintf(screen, "map_style = %d\n", map_style);
  int recreate = 0;
  if (map_style != map_style_old) recreate = 1;
  return recreate;
}

/* ----------------------------------------------------------------------
   return next prime larger than n
------------------------------------------------------------------------- */

int Element::next_prime(int n)
{
  int factor;
 
  int nprime = n+1;
  if(nprime % 2 == 0) nprime++;
  int root = static_cast<int> (sqrt(1.0*n)) + 2;

  while (nprime <= MAXSMALLINT) {
    for (factor = 3; factor <root; factor++)
      if (nprime % factor ==0) break;
    if (factor == root) return nprime;
    nprime +=2;
  }

  return MAXSMALLINT;
}

void Element::map_delete()
{
  memory->destroy(sametag);
  sametag = NULL;
  max_same = 0;

  if (map_style == 1) {
    memory->destroy(map_array);
    map_array = NULL;
  } else {
    if (map_nhash) {
      delete [] map_bucket;
      delete [] map_hash;
      map_bucket = NULL;
      map_hash = NULL;
    }
    map_nhash = 0;
  }
}

void Element::map_set()
{
  int nall = nlocal + nghost;
   //fprintf(screen, "me = %d,nall = %d\n",comm->me,nall);
  if (map_style == 1) {

    if (nall > max_same) {
      max_same = nall + EXTRA;
      //fprintf(screen, "me=%d, max_same = %d\n", comm->me,max_same);
      memory->destroy(sametag);
      memory->create(sametag,max_same,"atom:sametag");
    }
    
    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_array[tag[i]];
      map_array[tag[i]] = i;
     // fprintf(screen, "me = %d, tag[i] = %d, map_array = %d\n",comm->me,tag[i],i);
    }
  } else {
    if (nall > max_same) {
      max_same = nall + EXTRA;
      memory->destroy(sametag);
      memory->create(sametag,max_same,"atom:sametag");
    }

   int previous,ibucket,index;
   tagint global;
   
    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_find_hash(tag[i]);

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index > -1) {
        map_hash[index].local = i;
        continue;
      }
      index = map_free;
      map_free = map_hash[map_free].next;
      if (previous == -1) map_bucket[ibucket] = index;
      else map_hash[previous].next = index;
      map_hash[index].global = global;
      map_hash[index].local = i;
      map_hash[index].next = -1;
      map_nused++;
    }
  }
}

int Element::map_find_hash(tagint global)
{
  int local = -1;
  int index = map_bucket[global % map_nbucket];
  while (index > -1) {
    if (map_hash[index].global == global) {
      local = map_hash[index].local;
      break;
    }
    index = map_hash[index].next;
  }
  return local;
}

void Element::data_nodes(int n, char *buf)
{
  int j,m=0;
  tagint tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = atom->count_words(buf);
  *next = '\n';

  char **values = new char*[nwords];

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    
    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = ATOTAGINT(values[0]);
    if ( (m = map(tagdata))>=0) {
       if (nwords==6) data_node(m,&values[1]);
       if (nwords==9) data_node_vel(m,&values[1]);
    }
    buf = next + 1; 
  }

   delete [] values;
}

void Element::data_node(int m, char **values)
{
  int nn = ATOTAGINT(values[0])-1;

  nodex[m][nn][0] = atof(values[1]);
  nodex[m][nn][1] = atof(values[2]);
  nodex[m][nn][2] = atof(values[3]); 

  nodev[m][nn][0] = 0.0;
  nodev[m][nn][1] = 0.0;
  nodev[m][nn][2] = 0.0;

  tagint mm = ATOTAGINT(values[4]);
  nodetag[m][nn] = mm;

}

void Element::data_node_vel(int m, char **values)
{
  int nn = ATOTAGINT(values[0])-1;

  nodex[m][nn][0] = atof(values[1]);
  nodex[m][nn][1] = atof(values[2]);
  nodex[m][nn][2] = atof(values[3]);

  nodev[m][nn][0] = atof(values[4]);
  nodev[m][nn][1] = atof(values[5]);
  nodev[m][nn][2] = atof(values[6]);

  tagint mm = ATOTAGINT(values[7]);
  nodetag[m][nn] = mm;
}

void Element::init()
{
  int n,nn,ii;
  n = ntypes;
  memory->create(size, n+1, "element:size");
  element_size();

  inter_in_ele = new double**[ntypes];
  inter_num = new int[ntypes];
  Spa_inter_ele = new double**[ntypes];

  for (ii=0;ii<ntypes;ii++) {
    nn = nae[ii+1][0]*nae[ii+1][1]*nae[ii+1][2];
    inter_num[ii] = nn;
    memory->create(inter_in_ele[ii],nn,3,"element:inter_in_ele");
    memory->create(Spa_inter_ele[ii],nn,npe,"element:Spa_inter_ele");
  }

  inter_ele();
  shape_array_inter_ele();
}


/*-------------------------------------------------------------*/

void Element::element_size()
{
  double *max = new double[ntypes+1];
  double length = 0.0;
  double lx,ly,lz;
  for (int i = 1; i<=ntypes; i++) max[i] = 0.0;
  // if (comm->me==0) fprintf(screen, "testing element size\n");

  for (int i=0; i<nlocal; i++) {
    //length = 0.0;
    
    int itype = type[i];
    for (int j=0; j<4; j++) {
      if (j<2){
        lx = nodex[i][j][0]-nodex[i][j+6][0];
        ly = nodex[i][j][1]-nodex[i][j+6][1];
        lz = nodex[i][j][2]-nodex[i][j+6][2];
      } else {
        lx = nodex[i][j][0]-nodex[i][j+2][0];
        ly = nodex[i][j][1]-nodex[i][j+2][1];
        lz = nodex[i][j][2]-nodex[i][j+2][2];
      }
      length = MAX(length, sqrt(lx*lx+ly*ly+lz*lz));
    }
    
    max[itype] = length;
    //if (comm->me==1) fprintf(screen, "max = %f\n", max[1]);
  }
   
  
  MPI_Allreduce(max,size,ntypes+1,MPI_DOUBLE,MPI_MAX,world);
   
  //{
  //  if (comm->me==0) fprintf(screen, "size = %f\n", size[j]);
  //}
  
  for (int j=1; j<=ntypes; j++) max_size = MAX(max_size, size[j]);
  if (comm->me==0) fprintf(screen, "max_size = %f\n", max_size); 
  delete [] max;
  double *xx = new double[3];
  double distance;
  for (int i=0; i<nlocal; i++) {
    xx[0] =xx[1] =xx[2] =0.0;
    for (int j=0; j<npe; j++) {
      xx[0] += nodex[i][j][0];
      xx[1] += nodex[i][j][1];
      xx[2] += nodex[i][j][2];
    }
    xx[0] = xx[0]/npe;
    xx[1] = xx[1]/npe;
    xx[2] = xx[2]/npe;
    distance = (xx[0]-x[i][0])*(xx[0]-x[i][0])
               +(xx[1]-x[i][1])*(xx[1]-x[i][1])
               +(xx[2]-x[i][2])*(xx[2]-x[i][2]);
    if (distance > TOL) error->all(FLERR, "element data has error");
   }
  delete [] xx;
  //error->all(FLERR, "test from element_size");
}

/* ----------------------------------------------------------------------
   pack data for element I for sending to another proc
   first 8*3 xyz
------------------------------------------------------------------------- */

int Element::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  for (int j=0; j<npe; j++)
  {
    buf[m++] = nodex[i][j][0];
    buf[m++] = nodex[i][j][1];
    buf[m++] = nodex[i][j][2];
    buf[m++] = nodev[i][j][0];
    buf[m++] = nodev[i][j][1];
    buf[m++] = nodev[i][j][2];
  }
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   copy element I info to element J
------------------------------------------------------------------------- */
void Element::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];

  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];

  for (int k=0; k<npe; k++)
  {
    nodex[j][k][0] = nodex[i][k][0];
    nodex[j][k][1] = nodex[i][k][1];
    nodex[j][k][2] = nodex[i][k][2];
    nodev[j][k][0] = nodev[i][k][0];
    nodev[j][k][1] = nodev[i][k][1];
    nodev[j][k][2] = nodev[i][k][2];
  }
}

/* ---------------------------------------------------------------------- */

int Element::unpack_exchange(double *buf)
{
  if (nlocal == nmax) grow(0);
  
  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  for (int j=0; j<npe; j++)
  {
    nodex[nlocal][j][0] = buf[m++];
    nodex[nlocal][j][1] = buf[m++];
    nodex[nlocal][j][2] = buf[m++];
    nodev[nlocal][j][0] = buf[m++];
    nodev[nlocal][j][1] = buf[m++];
    nodev[nlocal][j][2] = buf[m++];
  }

  nlocal++;
  return m;
}

/*------------------------------------------------------------------------*/

void Element::specify_plane()
{
  plane_index[0][0] = 0; plane_index[0][1] = 1; plane_index[0][2] = 2; plane_index[0][3] = 3;
  plane_index[1][0] = 0; plane_index[1][1] = 3; plane_index[1][2] = 7; plane_index[1][3] = 4;
  plane_index[2][0] = 0; plane_index[2][1] = 1; plane_index[2][2] = 5; plane_index[2][3] = 4;
  plane_index[3][0] = 1; plane_index[3][1] = 2; plane_index[3][2] = 6; plane_index[3][3] = 5;
  plane_index[4][0] = 2; plane_index[4][1] = 3; plane_index[4][2] = 7; plane_index[4][3] = 6;
  plane_index[5][0] = 4; plane_index[5][1] = 5; plane_index[5][2] = 6; plane_index[5][3] = 7;

}

/*-----------------------------------------------------------------------*/

void Element::specify_sub_ele_pc(int nn)
{
  // specify the first sub element parent coordinate
  double ls = 2.0/nsplit_ele;
  sub_ele_pc[0][0][0] = -1.0; sub_ele_pc[0][0][1] = -1.0; sub_ele_pc[0][0][2] = -1.0;
  sub_ele_pc[0][1][0]=ls-1.0;sub_ele_pc[0][1][1] = -1.0; sub_ele_pc[0][1][2] = -1.0;
  sub_ele_pc[0][2][0] = ls-1.0;sub_ele_pc[0][2][1]=ls-1.0;  sub_ele_pc[0][2][2] = -1.0;
  sub_ele_pc[0][3][0] = -1.0; sub_ele_pc[0][3][1] =ls-1.0;sub_ele_pc[0][3][2] = -1.0;
  sub_ele_pc[0][4][0] = -1.0; sub_ele_pc[0][4][1] = -1.0; sub_ele_pc[0][4][2]=ls-1.0;
  sub_ele_pc[0][5][0]=ls-1.0;  sub_ele_pc[0][5][1] = -1.0; sub_ele_pc[0][5][2]=ls-1.0;
  sub_ele_pc[0][6][0]=ls-1.0;sub_ele_pc[0][6][1]=ls-1.0;sub_ele_pc[0][6][2]=ls-1.0;
  sub_ele_pc[0][7][0]=-1.0;sub_ele_pc[0][7][1]=ls-1.0;sub_ele_pc[0][7][2]=ls-1.0;

  // repeat the parent coordinate along the two directions
  int in = 0;
  int ii,jj,kk,ll; 
  for (ii=0; ii<nsplit_ele; ii++)
    for (jj=0; jj<nsplit_ele; jj++)
      for (kk=0; kk<nsplit_ele; kk++) {
       if (ii==0 &&jj==0 &&kk==0) continue;
       in++;
       for (ll=0; ll <npe; ll++) {
	 sub_ele_pc[in][ll][0] = sub_ele_pc[0][ll][0] + ii*ls;
	 sub_ele_pc[in][ll][1] = sub_ele_pc[0][ll][1] + jj*ls;
	 sub_ele_pc[in][ll][2] = sub_ele_pc[0][ll][2] + kk*ls;
       }
      }

  for (ii=0; ii< nsubele; ii++) {

    center_subele_pc[ii][0] = 0.0;
    center_subele_pc[ii][1] = 0.0;
    center_subele_pc[ii][2] = 0.0;

    for (jj=0; jj<npe;jj++) {
     center_subele_pc[ii][0] += sub_ele_pc[ii][jj][0];
     center_subele_pc[ii][1] += sub_ele_pc[ii][jj][1];
     center_subele_pc[ii][2] += sub_ele_pc[ii][jj][2]; 
    }

    center_subele_pc[ii][0] = center_subele_pc[ii][0]/npe;
    center_subele_pc[ii][1] = center_subele_pc[ii][1]/npe;
    center_subele_pc[ii][2] = center_subele_pc[ii][2]/npe;

  }

  //if (comm->me==0) {
   // for (ii=0;ii<nsubele;ii++) {
  //    for (jj=0;jj<npe;jj++) {
  //	fprintf(logfile, "%f, %f, %f\n", center_subele_pc[ii][0], center_subele_pc[ii][1],center_subele_pc[ii][2]);
  //    }
  //  }
  //}
  //error->all(FLERR, "test from element");
}

/*---------------------------------------------------------------------- */

void Element::count_sub_element()
{
  double xlo, xhi, ylo, yhi, zlo, zhi;
  int ninter = inter_num[0];
  double **pinter = inter_in_ele[0];

  // search around all sub element
  for (int ii = 0; ii < nsubele; ii++) {
    
     natom_subele[ii] = 0;

     xlo = sub_ele_pc[ii][0][0];
     xhi = sub_ele_pc[ii][1][0];
     if (xlo < -0.99) xlo = -1.01;
     if (xhi > 0.99) xhi = 1.01;

     ylo = sub_ele_pc[ii][0][1];
     yhi = sub_ele_pc[ii][3][1];
     if (ylo < -0.99) ylo = -1.01;
     if (yhi > 0.99) yhi = 1.01;

     zlo = sub_ele_pc[ii][0][2];
     zhi = sub_ele_pc[ii][4][2];
     if (zlo < -0.99) zlo = -1.01;
     if (zhi > 0.99) zhi = 1.01;

     for (int jj =0; jj < ninter; jj++) {
       if (pinter[jj][0] > xlo && pinter[jj][0] <= xhi &&
           pinter[jj][1] > ylo && pinter[jj][1] <= yhi &&
	   pinter[jj][2] > zlo && pinter[jj][2] <= zhi) {
	  natom_subele[ii]++;
        }
     }
  }

  int cninter = 0;;

  for (int ii = 0; ii < nsubele; ii++) {
    cninter+=natom_subele[ii];
  }

  if (cninter!=ninter) error->all(FLERR, "sub atom total donot match total atom number");

  //error->all(FLERR, "test from count_sub_element");
    
}

/*----------------------------------------------------------------------------------------*/

void Element::shape_array_inter_subele()
{
  double xlo, xhi, ylo, yhi, zlo, zhi;
  int ninter = inter_num[0];
  double **pinter = inter_in_ele[0];
  int in = 0;
  double px, py, pz;

  for (int ii = 0; ii < nsubele;ii++) {
    
    in = 0;
    
    xlo = sub_ele_pc[ii][0][0];
    xhi = sub_ele_pc[ii][1][0];
    if (xlo < -0.99) xlo = -1.01;
    if (xhi > 0.99) xhi = 1.01;

    ylo = sub_ele_pc[ii][0][1];
    yhi = sub_ele_pc[ii][3][1];
    if (ylo < -0.99) ylo = -1.01;
    if (yhi > 0.99) yhi = 1.01;

    zlo = sub_ele_pc[ii][0][2];
    zhi = sub_ele_pc[ii][4][2];
    if (zlo < -0.99) zlo = -1.01;
    if (zhi > 0.99) zhi = 1.01;

    for (int jj = 0; jj < ninter; jj++) {
      if (pinter[jj][0] > xlo && pinter[jj][0] <= xhi &&
          pinter[jj][1] > ylo && pinter[jj][1] <= yhi &&
	  pinter[jj][2] > zlo && pinter[jj][2] <= zhi) {

	 px = pinter[jj][0];
	 py = pinter[jj][1];
	 pz = pinter[jj][2];

	 Spa_inter_subele[ii][in][0] = 0.125*(1-px)*(1-py)*(1-pz);
	 Spa_inter_subele[ii][in][1] = 0.125*(1+px)*(1-py)*(1-pz);
	 Spa_inter_subele[ii][in][2] = 0.125*(1+px)*(1+py)*(1-pz);
	 Spa_inter_subele[ii][in][3] = 0.125*(1-px)*(1+py)*(1-pz);
	 Spa_inter_subele[ii][in][4] = 0.125*(1-px)*(1-py)*(1+pz);
	 Spa_inter_subele[ii][in][5] = 0.125*(1+px)*(1-py)*(1+pz);
	 Spa_inter_subele[ii][in][6] = 0.125*(1+px)*(1+py)*(1+pz);
	 Spa_inter_subele[ii][in][7] = 0.125*(1-px)*(1+py)*(1+pz);
	 list_inter_subele[ii][in] = jj;
	 //if (comm->me==0) fprintf(logfile, "%d\n",jj);
         in++;

       }
    }
  }

  for (int ii=0; ii < nsubele; ii++)
  {
    px = center_subele_pc[ii][0];
    py = center_subele_pc[ii][1];
    pz = center_subele_pc[ii][2];
     
    Spa_center_subele[ii][0] = 0.125*(1-px)*(1-py)*(1-pz);
    Spa_center_subele[ii][1] = 0.125*(1+px)*(1-py)*(1-pz);
    Spa_center_subele[ii][2] = 0.125*(1+px)*(1+py)*(1-pz);
    Spa_center_subele[ii][3] = 0.125*(1-px)*(1+py)*(1-pz);
    Spa_center_subele[ii][4] = 0.125*(1-px)*(1-py)*(1+pz);
    Spa_center_subele[ii][5] = 0.125*(1+px)*(1-py)*(1+pz);
    Spa_center_subele[ii][6] = 0.125*(1+px)*(1+py)*(1+pz);
    Spa_center_subele[ii][7] = 0.125*(1-px)*(1+py)*(1+pz);
   // if (comm->me==0) {
     // for (int jj =0 ; jj <npe; jj++) fprintf(logfile, "%f ",Spa_center_subele[ii][jj]);
     // fprintf(logfile,"\n");
    //}  

  }

  //error->all(FLERR, "test from shape_subele");
}

/* ---------------------------------------------------------------------- */

int Element::pack_border(int n, int *list, double *buf,
                         int pbc_flag, int *pbc)
{
 int i,j,m;
 double dx,dy,dz;

 m=0;
 //fprintf(screen, "me = %d, pbc_flag = %d\n", comm->me, pbc_flag); 
 if (pbc_flag==0) {
   for (i=0; i<n; i++) {
     j = list[i];
     buf[m++] = x[j][0];
     buf[m++] = x[j][1];
     buf[m++] = x[j][2];
     buf[m++] = ubuf(tag[j]).d;
     buf[m++] = ubuf(type[j]).d;
     buf[m++] = ubuf(mask[j]).d;
     for (int k=0; k<npe; k++) {
       buf[m++] = nodex[j][k][0];
       buf[m++] = nodex[j][k][1];
       buf[m++] = nodex[j][k][2];
     }
   }
 } else {
   if (domain->triclinic==0) {
     dx = pbc[0]*domain->xprd;
     dy = pbc[1]*domain->yprd;
     dz = pbc[2]*domain->zprd;
   }
   for (i=0; i<n; i++) {
     j = list[i];
     buf[m++] = x[j][0] + dx;
     buf[m++] = x[j][1] + dy;
     buf[m++] = x[j][2] + dz;
     buf[m++] = ubuf(tag[j]).d;
     buf[m++] = ubuf(type[j]).d;
     buf[m++] = ubuf(type[j]).d;
     for (int k=0; k<npe; k++) {
       buf[m++] = nodex[j][k][0]+dx;
       buf[m++] = nodex[j][k][1]+dy;
       buf[m++] = nodex[j][k][2]+dz;
     }
   }
 }
 return m;
}

/* ---------------------------------------------------------------------- */

void Element::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i=first; i < last; i++) {
    if (i==nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    for (int k=0; k<npe; k++) {
      nodex[i][k][0] = buf[m++];
      nodex[i][k][1] = buf[m++];
      nodex[i][k][2] = buf[m++];
    }
  }
}

/*----------------------------------------------------------------------------*/

int Element::pack_comm(int n, int *list, double *buf,
		        int pbc_flag, int *pbc)
{
 int i,j,m;
 double dx,dy,dz;

 m=0;

 if (pbc_flag==0) {
  for (i=0; i<n; i++) {
   j = list[i];
   for (int k=0; k<npe; k++) {
     buf[m++] = nodex[j][k][0];
     buf[m++] = nodex[j][k][1];
     buf[m++] = nodex[j][k][2];
   }
  }
 } else {
  if (domain->triclinic==0) {
   dx = pbc[0]*domain->xprd;
   dy = pbc[1]*domain->yprd;
   dz = pbc[2]*domain->zprd;
  }
  for (i=0; i<n; i++) {
   j = list[i];
   for (int k=0; k<npe; k++) {
    buf[m++] = nodex[j][k][0]+dx;
    buf[m++] = nodex[j][k][1]+dy;
    buf[m++] = nodex[j][k][2]+dz;
   }
  }
 }
 return m;
}

/*---------------------------------------------------------------------------*/
void Element::unpack_comm(int n, int first, double *buf) 
{
 int i,m,last;

 m = 0;
 last = first + n;
 for (i=first; i < last; i++) {
  for (int k=0; k<npe; k++) {
   nodex[i][k][0] = buf[m++];
   nodex[i][k][1] = buf[m++];
   nodex[i][k][2] = buf[m++];
  }
 }
}
/*----------------------------------------------------------------------------*/

void Element::grow_nae(int n)
{
  if (n==0) grow_nmaxinp();
  else nmaxinp = n;
  if (nmaxinp < 0)
    error->one(FLERR, "Per-processor system is too big");

  memory->grow(inpx,nmaxinp,3,"element:inpx");


}

/*----------------------------------------------------------------------------*/
void Element::grow_nmaxinp()
{
  nmaxinp = nmaxinp/NDELTA *NDELTA;
  nmaxinp += NDELTA;
}
/*----------------------------------------------------------------------------*/

void Element::allocate_type_arrays()
{
  memory->create(nae,(ntypes+1),3,"element:nae");
  mass = new double[ntypes+1];
  nae_setflag = new int[ntypes+1];
  for (int itype = 1; itype<=ntypes;itype++) nae_setflag[itype] = 0;
  numint = new int[ntypes+1];
  for (int itype = 1; itype<=ntypes; itype++) numint[itype] = 0; 
}

/*----------------------------------------------------------------------------
 Allocate array to store the integration points
------------------------------------------------------------------------------*/

void Element::allocate_int_arrays()
{
  maxnum_int = 0;
  // Find out which type of element has the largest number of 
  for (int itype = 1; itype <=ntypes; itype++) {
    if (maxnum_int < numint[itype]) maxnum_int = numint[itype];
  }
  if (comm->me==0) fprintf(logfile, "Largest num of integration points: %d\n",maxnum_int);
  // allocate arrays to store integration points of each type of element

  memory->create(intp, ntypes,maxnum_int,3,"element:integration points");
  memory->create(intpd, ntypes, maxnum_int,3,"element:integration points parent coordinate");
  memory->create(weight, ntypes, maxnum_int, "element:integration weight");
}

/*----------------------------------------------------------------------------*/

void Element::set_nae(const char *str, int type_offset)
{
  if (nae==NULL) error->all(FLERR,"Cannot set nae for this element style");
  
  //error->all(FLERR,"testing from set_nae in element");

  int itype;
  int nx,ny,nz;
  int n = sscanf(str,"%d %d %d %d", &itype,&nx,&ny,&nz);
  if (n != 4) error->all(FLERR, "Invalid nae line in data file");
  itype += type_offset;

  if (itype < 1 || itype > ntypes)
    error->all(FLERR,"Invalid type for nae set");

  nae[itype][0] = nx;
  nae[itype][1] = ny;
  nae[itype][2] = nz;
  mass[itype] = nx*ny*nz*atom->mass[1];
  //if (comm->me == 0) fprintf(logfile, "mass: %f", mass[itype]);
  //if (comm->me==0) fprintf(screen, "nx = %d, ny = %d, nz = %d\n", nx, ny, nz);
  if (nx<2||ny <2||nz <2) error->all(FLERR, "Invalid naae values");
  nae_setflag[itype] = 1;
 
  if (domain->dimension ==2 && nz != 0)
    error->all(FLERR, "repeat along third direction should be zero for 2d dimension\n"); 

  if (nae[itype][0] <=0 || nae[itype][1] <=0 || nae[itype][2]<=0)
    error->all(FLERR, "Invalid nae value");
}

void Element::set_mass(double m)
{

  for (int itype = 1; itype<=ntypes; itype++) {
    mass[itype] = nae[itype][0]*nae[itype][1]*nae[itype][2]*m;
  }
}

/*------------------------------------------------------------------------*/
void Element::set_numint(const char *str, int type_offset)
{
  int itype;
  int num;
  int n = sscanf(str,"%d %d", &itype,&num);
  if (n!=2) error->all(FLERR, "Invalid numint line in data file");
  itype +=type_offset;

  if (itype <1 || itype >ntypes)
    error->all(FLERR,"Invalid type for numint set");

  numint[itype] = num;
  //if (comm->me==0) fprintf(screen, "numint = %d\n", num);
  //error->all(FLERR, "test from numint ");
}

/*----------------------------------------------------------------------*/

void Element::set_intp(char *str, int type_offset)
{
  int itype, m;
  int x, y, z;
  double w;
  char **values = new char*[6];

  values[0]= strtok(str," \t\n\r\f");
  for (int ii=1; ii<6;ii++) values[ii] = strtok(NULL," \t\n\r\f");

  itype = atoi(values[0])-1;
  m = atoi(values[1])-1;
  x = atoi(values[2]);
  y = atoi(values[3]);
  z = atoi(values[4]);
  w = atof(values[5]);
 
  intp[itype][m][0] = x;
  intp[itype][m][1] = y;
  intp[itype][m][2] = z;
  weight[itype][m] = w;
  delete [] values;
}

/*-----------------------------------------------------------------------*/
 
void Element::intp_pcoord()
{
  int ii, jj, nint, itype;
  double dx, dy, dz;

  for (ii = 0; ii<ntypes; ii++) {
      itype = ii+1;
      nint = numint[itype];
      dx = 2/(static_cast<double>(nae[itype][0])-1);
      dy = 2/(static_cast<double>(nae[itype][1])-1);
      dz = 2/(static_cast<double>(nae[itype][2])-1);
      for (jj = 0; jj<nint; jj++) {
         intpd[ii][jj][0] = dx*(static_cast<double>(intp[ii][jj][0]))-1;
	 intpd[ii][jj][1] = dy*(static_cast<double>(intp[ii][jj][1]))-1;
	 intpd[ii][jj][2] = dz*(static_cast<double>(intp[ii][jj][2]))-1;
	 //if (comm->me==0) fprintf(screen, "%f, %f, %f\n", 
	 //		                   intpd[ii][jj][0],intpd[ii][jj][1], intpd[ii][jj][2]);
      }
  }


}

/*--------------------------------------------------------------------
   To count the total atoms when the elements are all interpolated 
----------------------------------------------------------------------*/

void Element::num_atoms_from_element()
{ 
  int i,itype;
  bigint num_element, num_total;

  num_element = nlocal+nghost;
  num_inp = 0; 
  num_inp_local = 0;
  num_inp_ghost = 0;  
  num_total = 0;
  //fprintf(screen, "me = %d, num_element = %d\n", comm->me, num_element);
  //error->all(FLERR,"testing from element");
  for (i=0; i<num_element; i++) {
    itype = type[i];
    if (nae_setflag[itype]) {
      num_inp += nae[itype][0]*nae[itype][1]*nae[itype][2];
      if (i<nlocal) num_inp_local += nae[itype][0]*nae[itype][1]*nae[itype][2];
      if (i>=nlocal) num_inp_ghost += nae[itype][0]*nae[itype][1]*nae[itype][2]; 
    } else { 
      error->all(FLERR, "element type nae not setted");
    }
  }
  
}

/*--------------------------------------------------------------------------
 Put all the interpolated atoms into the vec of atoms
---------------------------------------------------------------------------*/

void Element::Interpolate_to_atom_vec()
{
  num_atoms_from_element();
  // grow the x vector in atom vec
  atom->avec->grow(num_inp);
  atom->nlocal = num_inp_local;
  atom->nghost = num_inp_ghost;
  interpolate_atom_vec(); 
}

/*-------------------------------------------------------------------------
  do the interpolation and put the interpolated into the atom vec 
--------------------------------------------------------------------------*/

void Element::interpolate_atom_vec()
{
    int i,j,k,l,itype,nin_ele;
    int num_element;
    int nn;
    double xtmp, ytmp, ztmp;
    double **x = atom->x;
    num_element = nlocal+nghost;

    nn = 0;

    for (i=0; i<num_element;i++) {
      
       itype = type[i]-1;
       nin_ele = inter_num[itype]; 
       for (j=0;j<nin_ele; j++) {

	 xtmp = 0.0;
	 ytmp = 0.0;
	 ztmp = 0.0;

	 for (k=0;k<npe;k++) {
	    xtmp += Spa_inter_ele[itype][j][k]*nodex[i][k][0];
	    ytmp += Spa_inter_ele[itype][j][k]*nodex[i][k][1];
	    ztmp += Spa_inter_ele[itype][j][k]*nodex[i][k][2];
         }
	 x[nn][0] = xtmp;
	 x[nn][1] = ytmp;
	 x[nn][2] = ztmp;
	 nn++;
       }
    }
}
 
/*--------------------------------------------------------------------------
  Count the # of local integration points
---------------------------------------------------------------------------*/

void Element::count_intp()
{
  nint_local = 0;

  int i,itype;

  for (i=0; i<nlocal; i++) {
    itype = type[i];
    nint_local += numint[itype];
  }
  nmaxinp = nint_local;
  grow_nmaxinp();
  if (comm->me==0) fprintf(screen, "local # of integration point: %d, %d\n", nint_local, nmaxinp);

  //error->all(FLERR, "test from element count intp");
}

/*--------------------------------------------------------------------------------------*/
void Element::find_node_list()
{
  int ii,jj,kk,nint,it;

  for (ii=1; ii <= ntypes; ii++) {
    nint = numint[ii];
    it = ii-1;
    for (jj = 0; jj <npe;jj++) {
      for (kk = 0; kk <nint;kk++) {
	if (fabs(Spa_int[it][kk][jj]-1.0)<TOL) { 
	  node_list[ii][jj] = kk;
	  break;
	}
      }
      if (kk == nint) error->all(FLERR, "didn't contain all node integration point");
    }
  } 

  for (ii = 1; ii <= ntypes; ii++)
    if (comm->me==0) fprintf(logfile, "%d,%d,%d,%d,%d,%d,%d,%d\n",node_list[ii][0],
		                       node_list[ii][1],node_list[ii][2],node_list[ii][3],node_list[ii][4],
				       node_list[ii][5],node_list[ii][6],node_list[ii][7]);
  //error->all(FLERR, "test from element find_node_list");
}

/*-----------------------------------------------------------------------------------
  Calculate the shape function arrays of the integration points
  Linear interpolation function is used now
------------------------------------------------------------------------------------*/

void Element::shape_array_int()
{
   int ii,nint,jj,kk;
   double px,py,pz;
   //if (comm->me==0) fprintf(logfile, "shape funciton of integration points\n");
   if (dimension==3) {
      for (ii=0; ii<ntypes; ii++) {
	nint = numint[ii+1];

	 for (jj=0; jj< nint; jj++) {
	   px = intpd[ii][jj][0];
	   py = intpd[ii][jj][1];
	   pz = intpd[ii][jj][2];

	   Spa_int[ii][jj][0] = 0.125*(1-px)*(1-py)*(1-pz);
	   Spa_int[ii][jj][1] = 0.125*(1+px)*(1-py)*(1-pz);
	   Spa_int[ii][jj][2] = 0.125*(1+px)*(1+py)*(1-pz);
	   Spa_int[ii][jj][3] = 0.125*(1-px)*(1+py)*(1-pz);
	   Spa_int[ii][jj][4] = 0.125*(1-px)*(1-py)*(1+pz);
	   Spa_int[ii][jj][5] = 0.125*(1+px)*(1-py)*(1+pz);
	   Spa_int[ii][jj][6] = 0.125*(1+px)*(1+py)*(1+pz);
	   Spa_int[ii][jj][7] = 0.125*(1-px)*(1+py)*(1+pz);
	   for (kk=0; kk<8;kk++) fprintf(logfile, "%f ", Spa_int[ii][jj][kk]);
	   fprintf(logfile,"\n");
	 }

      }
   } else error->all(FLERR, "Haven't consider 2d problem yet");   
}

/*-----------------------------------------------------------------------------------
  Calculate the shape function arrays of all the interpolated atoms inside a element
  Linear interpolation function is used now
------------------------------------------------------------------------------------*/

void Element::shape_array_inter_ele()
{
  int ii,nn,jj,kk;
  double px,py,pz;
 
  if (dimension==3) {
    for (ii=0; ii<ntypes; ii++) {
       nn = inter_num[ii];

       for (jj=0;jj<nn;jj++) {
	 px = inter_in_ele[ii][jj][0];
	 py = inter_in_ele[ii][jj][1];
	 pz = inter_in_ele[ii][jj][2];

	 Spa_inter_ele[ii][jj][0] = 0.125*(1-px)*(1-py)*(1-pz);
	 Spa_inter_ele[ii][jj][1] = 0.125*(1+px)*(1-py)*(1-pz);
	 Spa_inter_ele[ii][jj][2] = 0.125*(1+px)*(1+py)*(1-pz);
	 Spa_inter_ele[ii][jj][3] = 0.125*(1-px)*(1+py)*(1-pz);
	 Spa_inter_ele[ii][jj][4] = 0.125*(1-px)*(1-py)*(1+pz);
	 Spa_inter_ele[ii][jj][5] = 0.125*(1+px)*(1-py)*(1+pz);
	 Spa_inter_ele[ii][jj][6] = 0.125*(1+px)*(1+py)*(1+pz);
	 Spa_inter_ele[ii][jj][7] = 0.125*(1-px)*(1+py)*(1+pz);
       }
    }
  } else error->all(FLERR, "Haven't consider 2d problem yet");
}


/*--------------------------------------------------------------------------------------
  Store the parent coordinate of all the atoms in each type of elements
  used for finding neighbor of integration points and make the interpolation process faster
-------------------------------------------------------------------------------------------*/

void Element::inter_ele()
{
  int ii,jj,kk,itype,iint;
  int nx, ny, nz;
  int inum;
  double inx, iny, inz;

  for (itype=0;itype<ntypes;itype++) {
    nx = nae[itype+1][0];
    ny = nae[itype+1][1];
    nz = nae[itype+1][2];
    inx = 2/(static_cast<double>(nx)-1);
    iny = 2/(static_cast<double>(ny)-1);
    inz = 2/(static_cast<double>(nz)-1);

    inum = 0;
    for (ii=0;ii<nx; ii++) {
      for (jj=0;jj<ny;jj++) {
	for (kk=0;kk<nz;kk++) {
	 // store the parent coordinate of all intterpolated atoms in each type of element
	  inter_in_ele[itype][inum][0] = inx*static_cast<double>(ii)-1;
	  inter_in_ele[itype][inum][1] = iny*static_cast<double>(jj)-1;
	  inter_in_ele[itype][inum][2] = inz*static_cast<double>(kk)-1;
	  inum++;
        }
      }
    }
  }

}

/*--------------------------------------------------------------------------------------------------------------
  return # of bytes of allocated memory
----------------------------------------------------------------------------------------------------------------*/

bigint Element::memory_usage()
{
   bigint bytes = 0;

   bytes += memory->usage(tag,nmax);
   bytes += memory->usage(type,nmax);
   bytes += memory->usage(mask, nmax);
   bytes += memory->usage(x,nmax,3);
   bytes += memory->usage(nodex,nmax,npe,3);
   bytes += memory->usage(nodef, nmax,npe,3);
   bytes += memory->usage(nodev,nmax,npe,3);
   bytes += memory->usage(nodetag,nmax,npe);
}

