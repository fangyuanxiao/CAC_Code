#include "mpi.h"
#include "ctype.h"
#include "float.h"
#include "limits.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "node_neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"
#include "element.h"

using namespace CAC_NS;

#define EWALD_F 1.12837917

enum{NONE,RLINEAR,RSQ,BMP};

int Pair::instance_total = 0;

Pair::Pair(CAC *cac) : Pointers(cac)
{
  instance_me = instance_total++;

  THIRD = 1.0/3.0;

  eng_vdwl = eng_coul = 0.0;

  comm_forward = comm_reverse = comm_reverse_off = 0;

  single_enable = 1;
  restartinfo = 1;
  respa_enable = 0;
  one_coeff = 0;
  no_virial_fdotr_compute = 0;
  writedata = 0;
  ghostneigh = 0;

  nextra = 0;
  pvector = NULL;
  single_extra = 0;
  svector = NULL;

  ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag = dipoleflag = 0;
  reinitflag = 1;

  // pair_modify settingsx
  
  compute_flag = 1;
  manybody_flag = 0;
  offset_flag = 0;
  mix_flag = GEOMETRIC;
  tail_flag = 0;
  etail = ptail = etail_ij = ptail_ij = 0.0;
  ncoultablebits = 12;
  ndisptablebits = 12;
  tabinner = sqrt(2.0);
  tabinner_disp = sqrt(2.0);

  allocated = 0;
//  suffix_flag = Suffix::NONE;

  maxeatom = maxvatom = 0;
  maxelement = 0;
  eatom = NULL;
  vatom = NULL;
  enode = NULL;
  vnode = NULL;

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;

//  execution_space = Host;
//  datamask_read = ALL_MASK;
//  datamask_modify = ALL_MASK;

  copymode = 0;

}

Pair::~Pair()
{
  if (copymode) return;

  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(enode);
  memory->destroy(vnode);
}

/* ---------------------------------------------------------------------- */

void Pair::init()
{
  int i,j;

  if (offset_flag && tail_flag)
    error->all(FLERR,"Cannot have both pair_modify shift and tail set to yes");
  if (tail_flag && domain->dimension == 2)
    error->all(FLERR,"Cannot use pair tail corrections with 2d simulations");
  if (tail_flag && domain->nonperiodic && comm->me == 0)
    error->warning(FLERR,"Using pair tail corrections with nonperiodic system");
  if (!compute_flag && tail_flag)
    error->warning(FLERR,"Using pair tail corrections with compute set to no");
  if (!compute_flag && offset_flag)
    error->warning(FLERR,"Using pair potential shift with compute set to no");

  // for manybody potentials
  // check if bonded exclusions could invalidate the neighbor list

  // I,I coeffs must be set\
  // init_one() will check if I,J is set explicitly or inferred by mixing

  if (!allocated) error->all(FLERR,"All pair coeffs are not set");

  for (i = 1; i <= atom->ntypes; i++)
    if (setflag[i][i] == 0) error->all(FLERR,"All pair coeffs are not set");

  // style-specific initialization

  init_style();
 
  // call init_one() for each I,J
  // set cutsq for each I,J, used to neighbor
  // cutforce = max of all I,J cutoffs

  cutforce = 0.0;
  etail = ptail = 0.0;
  double cut;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
      cutforce = MAX(cutforce,cut);
      if (tail_flag) {
        etail += etail_ij;
        ptail += ptail_ij;
        if (i != j) {
          etail += etail_ij;
          ptail += ptail_ij;
        }
      }
  }

}

/* ----------------------------------------------------------------------
   init specific to a pair style
   specific pair style can override this function
     if needs its own error checks
     if needs another kind of neighbor list
   request default neighbor list = half list
------------------------------------------------------------------------- */

void Pair::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
    mixing of pair potential prefactors (epsilon)
------------------------------------------------------------------------- */

double Pair::mix_energy(double eps1, double eps2, double sig1, double sig2)
{
  if (mix_flag == GEOMETRIC)
    return sqrt(eps1*eps2);
  else if (mix_flag == ARITHMETIC)
    return sqrt(eps1*eps2);
  else if (mix_flag == SIXTHPOWER)
    return (2.0 * sqrt(eps1*eps2) *
      pow(sig1,3.0) * pow(sig2,3.0) / (pow(sig1,6.0) + pow(sig2,6.0)));
  else return 0.0;
}

/* ----------------------------------------------------------------------
   mixing of pair potential distances (sigma, cutoff)
------------------------------------------------------------------------- */

double Pair::mix_distance(double sig1, double sig2)
{
  if (mix_flag == GEOMETRIC)
    return sqrt(sig1*sig2);
  else if (mix_flag == ARITHMETIC)
    return (0.5 * (sig1+sig2));
  else if (mix_flag == SIXTHPOWER)
    return pow((0.5 * (pow(sig1,6.0) + pow(sig2,6.0))),1.0/6.0);
  else return 0.0;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   specific pair style can override this function
------------------------------------------------------------------------- */

void Pair::init_list(int which, NeighList *ptr)
{
  list = ptr;
}

/*-----------------------------------------------------------------------
  Make the nilist point to the node neighbor list needed
------------------------------------------------------------------------*/

void Pair::init_nilist(NodeNeighList *ptr)
{
   nilist = ptr;
   //error->all(FLERR, "test from init_nilist");
}

/* ---------------------------------------------------------------------- */

void Pair::compute_dummy(int eflag, int vflag)
{
  evflag = 0;
}

double Pair::memory_usage()
{
  double bytes = comm->nthreads*maxeatom * sizeof(double);
  bytes += comm->nthreads*maxvatom*6 * sizeof(double);
  bytes += maxelement*element->npe*sizeof(double);
  bytes += maxelement*element->npe*6*sizeof(double);
  return bytes;
}

void Pair::ev_setup()
{
 int i,n,j;

 if (element->nmax > maxelement) {
   maxelement = element->nmax;
   memory->destroy(enode);
   memory->destroy(vnode);
   memory->create(enode, maxelement, element->npe,"pair:enode");
   memory->create(vnode, maxelement,element->npe,6,"pair:vnode");
 }	

 eng_vdwl = 0.0;
  n = element->nlocal;
  for (i = 0; i <n; i++) {
    for (j=0; j< element->npe; j++) {
      enode[i][j] = 0.0;
      vnode[i][j][0] = 0.0;
      vnode[i][j][1] = 0.0;
      vnode[i][j][2] = 0.0;
      vnode[i][j][3] = 0.0;
      vnode[i][j][4] = 0.0;
      vnode[i][j][5] = 0.0;
    }
  }

 //error->all(FLERR, "test from ev_setup");
}	

