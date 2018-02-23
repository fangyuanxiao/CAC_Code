#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "update.h"
#include "output.h"
#include "comm.h"
#include "universe.h"
#include "memory.h"
#include "error.h"
#include "node_neigh_list.h"
using namespace CAC_NS;

#define RQDELTA 1
#define EXDELTA 1

#define LB_FACTOR 1.5
#define SMALL 1.0e-6
#define BIG 1.0e20
#define CUT2BIN_RATIO 100

enum{NSQ,BIN,MULTI};     

//#define NEIGH_LIST_DEBUG 1

Neighbor::Neighbor(CAC *cac) : Pointers(cac)
{
   MPI_Comm_rank(world,&me);
   MPI_Comm_size(world,&nprocs);

   style = BIN;
   every = 1;
   delay = 10;
   dist_check = 1;
   pgsize = 100000;
   oneatom = 2000;
   binsizeflag = 0;
   build_once = 0;
   cluster_check = 0;
   binatomflag = 1;

   cutneighsq = NULL;
   cutneighghostsq = NULL;
   cuttype = NULL;
   cuttypesq = NULL;
   fixchecklist = NULL;

   // coords at last neighboring
   
   maxhold = 0;
   xhold = NULL;
   lastcall = -1;

   emaxhold = 0;
   nodexhold = NULL;
   // binning
  
   maxhead = 0;
   binhead = NULL;
   ebinhead = NULL;
   maxbin = 0;
   emaxbin = 0;
   bins = NULL;
   ebins = NULL;
   // pair exclusion list info
  
   includegroup = 0;

   nex_type = maxex_type = 0;
   ex1_type = ex2_type = NULL;
   ex_type = NULL;

   nex_group = maxex_group = 0;
   ngrow_flag = 0;
   ex1_group = ex2_group = ex1_bit = ex2_bit = NULL;

   nex_mol = maxex_mol = 0;
   ex_mol_group = ex_mol_bit = NULL;

   // pair lists

   maxatom = 0;
   maxele = 0;
   nblist = nglist = nslist = 0;
   
   nlist = 0;
   lists = NULL;
   nilist = NULL;
   pair_build = NULL;
   stencil_create = NULL;
   blist = glist = slist = NULL;
   anyghostlist = 0;

   nrequest = maxrequest = 0;
   requests = NULL;
   
   old_nrequest = 0;
   old_requests = NULL;

   old_style = style;
   old_triclinic = 0;
   old_pgsize = 0;
   old_oneatom = oneatom;
   old_every = every;
   old_delay = delay;
   old_check = dist_check;
   old_cutoff = cutneighmax;

   // bond lists
  
   maxbond = 0;
   bondlist = NULL;
   maxangle = 0;
   anglelist = NULL;
   maxdihedral = 0;
   dihedrallist = NULL;
   maximproper = 0;
   improperlist = NULL;

   copymode = 0;
}

/*---------------------------------------------------------------------------*/

Neighbor::~Neighbor()
{

   if (copymode) return;

   memory->destroy(cutneighsq);
   memory->destroy(cutneighghostsq);
   delete [] cuttype;
   delete [] cuttypesq;
   delete [] fixchecklist;

   memory->destroy(xhold);
   memory->destroy(nodexhold);

   memory->destroy(binhead);
   memory->destroy(ebinhead);
   memory->destroy(bins);
   memory->destroy(ebins);
   memory->destroy(ex1_type);
   memory->destroy(ex2_type);
   memory->destroy(ex_type);

   memory->destroy(ex1_group);
   memory->destroy(ex2_group);
   delete [] ex1_bit;
   delete [] ex2_bit;

   memory->destroy(ex_mol_group);
   delete [] ex_mol_bit;

   for (int i = 0; i < nlist; i++) delete lists[i];
   delete [] lists;
   delete [] pair_build;
   delete [] stencil_create;
   delete [] blist;
   delete [] glist;
   delete [] slist;
   //error->all(FLERR, "test from ~neighbor");
   delete nilist;
   //error->all(FLERR, "test from ~neighbor");
   for (int i = 0; i < nrequest; i++) delete requests[i];
   memory->sfree(requests);
   for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
   memory->sfree(old_requests);

   memory->destroy(bondlist);
   memory->destroy(anglelist);
   memory->destroy(dihedrallist);
   memory->destroy(improperlist);
}

/* ----------------------------------------------------------------------
  set neighbor style and skin distance
 ------------------------------------------------------------------------- */

void Neighbor::set(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal neighbor command");

  skin = universe->numeric(FLERR,arg[0]);
  if (skin < 0.0) error->all(FLERR,"Illegal neighbor command");

  if (strcmp(arg[1],"nsq") == 0) style = NSQ;
  else if (strcmp(arg[1],"bin") == 0) style = BIN;
  else if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else error->all(FLERR,"Illegal neighbor command");

}

/* ----------------------------------------------------------------------
 *    modify parameters of the pair-wise neighbor build
 *    ------------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      every = universe->inumeric(FLERR,arg[iarg+1]);
      if (every <= 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      delay = universe->inumeric(FLERR,arg[iarg+1]);
      if (delay < 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) dist_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) dist_check = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"once") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) build_once = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) build_once = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      old_pgsize = pgsize;
      pgsize = universe->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      old_oneatom = oneatom;
      oneatom = universe->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      binsize_user = universe->numeric(FLERR,arg[iarg+1]);
      if (binsize_user <= 0.0) binsizeflag = 0;
      else binsizeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"cluster") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) cluster_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) cluster_check = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal neigh_modify command");
  }
}

/* ---------------------------------------------------------------------- */

int Neighbor::request(void *requestor, int instance)
{
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **)
       memory->srealloc(requests,maxrequest*sizeof(NeighRequest *),
                       "neighbor:requests"); 
  }

  requests[nrequest] = new NeighRequest(cac);
  requests[nrequest]->requestor = requestor;
  requests[nrequest]->requestor_instance = instance;
  nrequest++;
  return nrequest-1;
}

/*-----------------------------------------------------------------------*/

void Neighbor::setup_cut()
{
  cut_all = cutneighmax + element->max_size;
  cut_allsq = cut_all*cut_all;
  cut_ele = cutneighmax + element->max_size/2;
  cut_elesq = cut_ele*cut_ele;
  cut_subele = cutneighmax + element->max_size/(2*element->nsplit_ele);
  cut_subelesq = cut_subele*cut_subele;

}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
  int i,j,m,n;

  ncalls = ndanger = 0;
  dimension = domain->dimension;
  triclinic = domain->triclinic;
 // newton_pair = force->newton_pair;

  // error check

  if (delay > 0 && (delay % every) != 0)
    error->all(FLERR,"Neighbor delay must be 0 or multiple of every setting");

  if (pgsize < 10*oneatom)
    error->all(FLERR,"Neighbor page size must be >= 10x the one atom setting");

  // ------------------------------------------------------------------
  // settings
 
  // bbox lo/hi = bounding box of entire domain, stored by Domain

  if (triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }

  // set neighbor cutoffs (force cutoff + skin)
  // trigger determines when atoms migrate and neighbor lists are rebuilt
  //   needs to be non-zero for migration distance check
  //   even if pair = NULL and no neighbor lists are used
  // cutneigh = force cutoff + skin if cutforce > 0, else cutneigh = 0
  // cutneighghost = pair cutghost if it requests it, else same as cutneigh

  triggersq = 0.25*skin*skin;
  // if (comm->me==0) fprintf(logfile, "etriggersq = %f\n", etriggersq);
  boxcheck = 0;
  if (domain->box_change && (domain->xperiodic || domain->yperiodic ||
                             (dimension == 3 && domain->zperiodic)))
      boxcheck = 1;

 n = atom->ntypes;
 if (cutneighsq == NULL) {
   memory->create(cutneighsq,n+1,n+1,"neigh:cutneighsq");
   //memory->create(cutneighghostsq,n+1,n+1,"neigh:cutneighghostsq");
   cuttype = new double[n+1];
   cuttypesq = new double[n+1];    
 } 
   
  double cutoff,delta,cut;
  cutneighmin = BIG;
  cutneighmax = 0.0;
  
  for (i = 1; i <= n; i++) {
    cuttype[i] = cuttypesq[i] = 0.0;
    for (j = 1; j <= n; j++) {

      cutoff = 0.0;
      cut = cutoff + skin;

      cutneighsq[i][j] = cut*cut;
      cuttype[i] = MAX(cuttype[i],cut);
      cuttypesq[i] = MAX(cuttypesq[i],cut*cut);
      cutneighmin = MIN(cutneighmin,cut);
      cutneighmax = MAX(cutneighmax,cut);

    }
  }
  

  cutneighmaxsq = cutneighmax * cutneighmax;
  cut_all = cutneighmax + element->max_size;
  cut_ele = cutneighmax + element->max_size/2;
  cut_allsq = cut_all*cut_all;
  cut_elesq = cut_ele*cut_ele;
 
  // check other classes that can induce reneighboring in decide()

  restart_check = 0;
  if (output->restart_flag) restart_check = 1;

  delete [] fixchecklist;
  fix_check = 0;
  must_check = 0;

  int respa = 0;

  // ------------------------------------------------------------------
  

  if (style != NSQ && atom->nmax) {
    if (maxbin == 0) {
      maxbin = atom->nmax;
      memory->create(bins,maxbin,"bins");
    }
  }

  // exclusion lists for type, group, molecule settings from neigh_modify
  // warn if exclusions used with KSpace solver

  n = atom->ntypes;

  int same = 0;
 
  // if old and new are not the same, create new pairwise lists

  if (!same) {
     
    for (i = 0; i < nlist; i++) delete lists[i];
    delete [] lists;
    delete [] pair_build;
    delete [] stencil_create;

    nrequest = 1;
    nlist = nrequest;

    lists = new NeighList*[nrequest];

   // initialize to NULL since some may be Kokkos lists
  
   for (i = 0; i < nrequest; i++) {
      lists[i] = NULL;
    }

   // create individual lists, one per request
   // pass list ptr back to requestor (except for Command class)
   // wait to allocate initial pages until copy lists are detected

   for (i = 0; i < nrequest; i++) {
     lists[i] = new NeighList(cac);
     lists[i]->index = i;
   }

  
   // allocate initial pages for each list, except if listcopy set

   for (i = 0; i < nrequest; i++) {
       lists[i]->setup_pages(pgsize,oneatom,0);
    }

    // allocate atom arrays for neighbor lists that store them

    maxatom = atom->nmax;
    for (i = 0; i < nrequest; i++) {
      if (lists[i]) {
        lists[i]->grow(maxatom);
      }
    }    
  }    
  // output neighbor list info, only first time or when info changes

  if (!same || every != old_every || delay != old_delay || 
      old_check != dist_check || old_cutoff != cutneighmax) {
    if (me == 0) {
      const double cutghost = MAX(cutneighmax,comm->cutghostuser);

      if (logfile) {
        fprintf(logfile,"Neighbor list info ...\n");
        fprintf(logfile,"  %d neighbor list requests\n",nrequest);
        fprintf(logfile,"  update every %d steps, delay %d steps, check %s\n",
                every,delay,dist_check ? "yes" : "no");
        fprintf(logfile,"  master list distance cutoff = %g\n",cutneighmax);
        fprintf(logfile,"  ghost atom cutoff = %g\n",cutghost);
      }
      if (screen) {
        fprintf(screen,"Neighbor list info ...\n");
        fprintf(screen,"  %d neighbor list requests\n",nrequest);
        fprintf(screen,"  update every %d steps, delay %d steps, check %s\n",
                every,delay,dist_check ? "yes" : "no");
        fprintf(screen,"  master list distance cutoff = %g\n",cutneighmax);
        fprintf(screen,"  ghost atom cutoff = %g\n",cutghost);
      }
    }
  }
  
  // mark all current requests as processed
  // delete old requests
  // copy current requests and style to old for next run

  memory->sfree(old_requests);
 
}

/* ----------------------------------------------------------------------
   determine which pair_build function each neigh list needs
   based on settings of neigh request
   copy -> copy_from function
   skip -> granular function if gran with granhistory,
           respa function if respaouter,
           skip_from function for everything else
   half_from_full, half, full, gran, respaouter ->
     choose by newton and rq->newton and tri settings
     style NSQ options = newton off, newton on
     style BIN options = newton off, newton on and not tri, newton on and tri
     stlye MULTI options = same options as BIN
   if none of these, ptr = NULL since pair_build is not invoked for this list
   use "else if" b/c skip,copy can be set in addition to half,full,etc
------------------------------------------------------------------------- */

void Neighbor::choose_build(int index, NeighRequest *rq)
{
  PairPtr pb = NULL;

  if (rq->omp == 0 && rq->intel == 0) {
   
    if (rq->copy) pb = &Neighbor::copy_from; 
    
    else if (rq->half) {
      if (style == BIN) {
        if (rq->newton == 0) {
          if (newton_pair == 0) { 
            if (rq->ghost == 0) pb = &Neighbor::half_bin_no_newton;
          }
        }
      }
    }
  }

  pair_build[index] = pb;
}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
------------------------------------------------------------------------- */

void Neighbor::copy_from(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->inum = listcopy->inum;
  list->gnum = listcopy->gnum;
  list->ilist = listcopy->ilist;
  list->numneigh = listcopy->numneigh;
  list->firstneigh = listcopy->firstneigh;
  list->firstdouble = listcopy->firstdouble;
  list->ipage = listcopy->ipage;
  list->dpage = listcopy->dpage;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and other bins in stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::half_bin_no_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,imol,iatom,moltemplate;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  // bin local & ghost atoms

  if (binatomflag) bin_atoms();

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  int molecular = atom->molecular;
  if (molecular == 2) moltemplate = 1;
  else moltemplate = 0;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  // loop over each atom, storing neighbors

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    // loop over all atoms in other bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    ibin = coord2bin(x[i]);

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (i == j) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighmaxsq) {
          neighptr[n++] = j;
        }
      }
    }
    
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
 
}

/*-------------------------------------------------------------------------
  bin owned and ghost atoms
--------------------------------------------------------------------------*/

void Neighbor::bin_elements()
{
  int i, ibin;

  for (i = 0; i<mbins;i++) ebinhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost elements at end of list, which is necessary

  double **x = element->x;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int nall = nlocal + element->nghost;

  for (i=nall-1; i>=0; i--) {
    ibin = coord2bin(x[i]);
    ebins[i] = ebinhead[ibin];
    ebinhead[ibin] = i;
  }
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void Neighbor::bin_atoms()
{
  int i,ibin;

  for (i = 0; i < mbins; i++) binhead[i] = -1;
  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = nall-1; i >= 0; i--) {
    ibin = coord2bin(x[i]);
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }

}

/* ----------------------------------------------------------------------
  onvert atom coords into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int Neighbor::coord2bin(double *x)
{
  int ix,iy,iz;

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
} 

/* ----------------------------------------------------------------------
   determine which stencil_create function each neigh list needs
   based on settings of neigh request, only called if style != NSQ
   skip or copy or half_from_full -> no stencil
   half, gran, respaouter, full -> choose by newton and tri and dimension
   if none of these, ptr = NULL since this list needs no stencils
   use "else if" b/c skip,copy can be set in addition to half,full,etc
------------------------------------------------------------------------- */

void Neighbor::choose_stencil(int index, NeighRequest *rq)
{
  StencilPtr sc = NULL;

  if (rq->half || rq->gran || rq->respaouter) {
    if (style == BIN) {
      if (rq->newton == 0) {
        if (newton_pair == 0) {
          if (dimension == 2) {
             sc = &Neighbor::stencil_half_bin_2d_no_newton;
          } else if (dimension == 3) {
            sc = &Neighbor::stencil_half_bin_3d_no_newton;
          }
        }
      } 
    }
  }

  stencil_create[index] = sc;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_2d_no_newton(NeighList *list,
                                             int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;

  list->nstencil = nstencil;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_3d_no_newton(NeighList *list,
                                             int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;

  list->nstencil = nstencil;
}

/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double Neighbor::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex;

  if (j > 0) dely = (j-1)*binsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey;

  if (k > 0) delz = (k-1)*binsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez;

  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   setup neighbor binning parameters
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*binsize to -binsize, etc
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
     a dim must be filled exactly by integer # of bins
     in periodic, procs on both sides of PBC must see same bin boundary
     in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
     tilted simulation box cannot contain integer # of bins
     stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
------------------------------------------------------------------------- */

void Neighbor::setup_bins()
{
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by comm->cutghost
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3],bsubboxlo[3],bsubboxhi[3];
  double *cutghost = comm->cutghost;

  if (triclinic == 0) {
    bsubboxlo[0] = domain->sublo[0] - cutghost[0]-element->max_size/2.0;
    bsubboxlo[1] = domain->sublo[1] - cutghost[1]-element->max_size/2.0;
    bsubboxlo[2] = domain->sublo[2] - cutghost[2]-element->max_size/2.0;
    bsubboxhi[0] = domain->subhi[0] + cutghost[0]+element->max_size/2.0;
    bsubboxhi[1] = domain->subhi[1] + cutghost[1]+element->max_size/2.0;
    bsubboxhi[2] = domain->subhi[2] + cutghost[2]+element->max_size/2.0;
  }

  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // for MULTI style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size
 
  double binsize_optimal;
  if (binsizeflag) binsize_optimal = binsize_user;
  else if (style == BIN) binsize_optimal = 0.5*cutneighmax;
  if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
  double binsizeinv = 1.0/binsize_optimal;
  // test for too many global bins in any dimension due to huge global domain

  if (bbox[0]*binsizeinv > MAXSMALLINT || bbox[1]*binsizeinv > MAXSMALLINT ||
      bbox[2]*binsizeinv > MAXSMALLINT)
    error->all(FLERR,"Domain too large for neighbor bins");

  // create actual bins
  // always have one bin even if cutoff > bbox
  // for 2d, nbinz = 1
  //if (comm->me==0) fprintf(screen, "bbox %f, %f, %f\n", bbox[0],bbox[1],bbox[2]);
  nbinx = static_cast<int> (bbox[0]*binsizeinv);
  nbiny = static_cast<int> (bbox[1]*binsizeinv);
  if (dimension == 3) nbinz = static_cast<int> (bbox[2]*binsizeinv);
  else nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  // compute actual bin size for nbins to fit into box exactly
  // error if actual bin size << cutoff, since will create a zillion bins
  // this happens when nbin = 1 and box size << cutoff
  // typically due to non-periodic, flat system in a particular dim
  // in that extreme case, should use NSQ not BIN neighbor style

  binsizex = bbox[0]/nbinx;
  binsizey = bbox[1]/nbiny;
  binsizez = bbox[2]/nbinz;
 
  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;
 
  if (binsize_optimal*bininvx > CUT2BIN_RATIO ||
      binsize_optimal*bininvy > CUT2BIN_RATIO ||
      binsize_optimal*bininvz > CUT2BIN_RATIO)
    error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  int mbinxhi,mbinyhi,mbinzhi;
  double coord;

  coord = bsubboxlo[0] - SMALL*bbox[0];
  mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
  if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL*bbox[0];
  mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

  coord = bsubboxlo[1] - SMALL*bbox[1];
  mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
  if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL*bbox[1];
  mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

  if (dimension == 3) {
    coord = bsubboxlo[2] - SMALL*bbox[2];
    mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
    if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
    coord = bsubboxhi[2] + SMALL*bbox[2];
    mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);
  }

  // extend bins by 1 to insure stencil extent is included
  // if 2d, only 1 bin in z

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  if (dimension == 3) {
    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;
  } else mbinzlo = mbinzhi = 0;
  mbinz = mbinzhi - mbinzlo + 1;

  // memory for bin ptrs

  bigint bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz);
  if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
  mbins = bbin;
  if (mbins > maxhead) {
    maxhead = mbins;
    memory->destroy(binhead);
    memory->create(binhead,maxhead,"neigh:binhead"); 
  }

  // create stencil of bins to search over in neighbor list construction
  // sx,sy,sz = max range of stencil in each dim
  // smax = max possible size of entire 3d stencil
  // stencil is empty if cutneighmax = 0.0
 
  sx = static_cast<int> (cutneighmax*bininvx);
  if (sx*binsizex < cutneighmax) sx++;
  sy = static_cast<int> (cutneighmax*bininvy);
  if (sy*binsizey < cutneighmax) sy++;
  sz = static_cast<int> (cutneighmax*bininvz);
  if (sz*binsizez < cutneighmax) sz++;
  if (dimension == 2) sz = 0;
  smax = (2*sx+1) * (2*sy+1) * (2*sz+1);

  // create stencils for pairwise neighbor lists
  // only done for lists with stencilflag and buildflag set
  lists[0]->stencil_allocate(smax,style);
  stencil_half_bin_3d_no_newton(lists[0],sx,sy,sz);
}

/* ----------------------------------------------------------------------
   build perpetuals neighbor lists
   called at setup and every few timesteps during run or minimization
   topology lists also built if topoflag = 1, USER-CUDA calls with topoflag = 0
------------------------------------------------------------------------- */

void Neighbor::build(int topoflag)
{

  // store current atom positions and box size if needed


  // if any lists store neighbors of ghosts:
  //   invoke grow() if nlocal+nghost exceeds previous list size
  // else only invoke grow() if nlocal exceeds previous list size
  // only for lists with growflag set and which are perpetual (glist)
  
  if (atom->nlocal > maxatom) {
    maxatom = atom->nmax;
    lists[0]->grow(maxatom);
  }
  
  // extend atom bin list if necessary

  if (style != NSQ && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->destroy(bins);
    memory->create(bins,maxbin,"bins");
  }
 
  // check that using special bond flags will not overflow neigh lists

  if (atom->nlocal+atom->nghost > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");

  // invoke building of pair and molecular topology neighbor lists
  // only for pairwise lists with buildflag set
  // blist is for standard neigh lists, otherwise is a Kokkos list

  half_bin_no_newton(lists[0]);
}

/*-------------------------------------------------------------------------------*/

void Neighbor::ebuild(int topoflag)
{
   int i,j;

   // store current node positions 
   if (dist_check) {
     double ***nodex = element->nodex;
     int nlocal = element->nlocal;
     int npe = element->npe;
     if (nlocal > emaxhold) {
       emaxhold = element->nmax;
       memory->destroy(nodexhold);
       memory->create(nodexhold, emaxhold, npe,3,"neigh:nodexhold");
     }
     for (i=0; i<nlocal; i++) {
	for (j=0; j<npe;j++) {
	  nodexhold[i][j][0] = nodex[i][j][0];
	  nodexhold[i][j][1] = nodex[i][j][1];
	  nodexhold[i][j][2] = nodex[i][j][2];
	}

     }
   }
   
   if (element->nlocal > maxele) {
     maxele = element->nmax;
     lists[0]->egrow(maxele);
     ngrow_flag = 1;
   }

   // extend element bin list if necessay
   if (element->nmax > emaxbin) {
      emaxbin = element->nmax;
      memory->destroy(ebins);
      memory->create(ebins, emaxbin, "bins");
   }

   if (element->nlocal+element->nghost > NEIGHMASK)
     error->one(FLERR, "Too many local+ghost elements for neighbor list");
   element_neighbor(lists[0]);
}

void Neighbor::output_ele_neigh(NeighList *list)
{
  int i,j,ii,jj,einum,itag,jtag,itype;
  int *eilist,*ejlist,*enumneigh,**efirstneigh;
  int *type, *tag;

  tag = element->tag;
  einum = list->einum;
  eilist = list->eilist;
  enumneigh = list->enumneigh;
  efirstneigh = list->efirstneigh;

if (comm->me==13) {
  for (ii=0;ii<einum;ii++) {
    i = eilist[ii];
    itype = type[i];
    itag = tag[i];
    einum = enumneigh[i];
    ejlist = efirstneigh[i];
    fprintf(screen,"i = %d\n",i);
  //  for (jj=0;jj<einum;jj++) {
  //    j = ejlist[jj];
  //    jtag = tag[j];
  //    fprintf(logfile, "%d\n",jtag);
  //  }
  }
}
  
}

/*------------------------------------------------------------------------*/

void Neighbor::nbuild(int topoflag)
{

  if (ngrow_flag) {
   maxele = element->nmax;
   nilist->grow(maxele, maxele*27); 
  }

   node_neighbor(nilist,lists[0]);
   ngrow_flag = 0;
}

/*-----------------------------------------------------------------------*/

void Neighbor::check_node_neighbor(NodeNeighList *nilist)
{
  int *nelist = nilist->elist;
  int *enlist = nilist->enlist;

  int ninum = nilist->ninum;
  int neinum = nilist->neinum;
  //if (comm->me==0) fprintf(logfile, "neinum = %d,%d\n", neinum,ninum);
  int *numneigh = nilist->numneigh;

  int **firstneigh = nilist->firstneigh;
  int **nfirstneigh = nilist->nfirstneigh;

  int *jlist,*njlist;
  int jnum;

  int ii,jj,kk,i,j,np,nn,itag,jtag,itype,jtype;
  int *type = element->type;
  int *numint = element->numint;
  int *tag = element->tag;
if (comm->me==0) {
  for (ii=0; ii<neinum;ii++) {
    i = enlist[ii];
    np = nelist[ii];
    itype = type[i];
    nn = numint[itype];
    for (jj=0; jj<nn; jj++) {
      jnum = numneigh[np];
      jlist = firstneigh[np];
      njlist = nfirstneigh[np];
      for (kk=0;kk<jnum;kk++) {
        fprintf(logfile, "%d,%d\n",jlist[kk],njlist[kk]);
      }
      np++;
    }
  }
}
}
/* ---------------------------------------------------------------------- */

int Neighbor::decide()
{
  //if (comm->me==0) fprintf(screen, "must_check = %d", must_check);  
  //if (comm->me==0) fprintf(screen, "ago = %d\n", ago);
  ago++;
  //if (comm->me==0) fprintf(screen, "every = %d\n", ago % every);
  //if (comm->me==0) fprintf(screen, "delay = %d\n", delay);
  //if (comm->me==0) fprintf(screen, "dist_check = %d\n", dist_check);
  if (ago >= delay && ago % every == 0) {
    if (build_once) return 0;
    if (dist_check == 0) return 1;
    return check_distance();
  } else return 0;
}

/* ----------------------------------------------------------------------
   if any atom moved trigger distance (half of neighbor skin) return 1
   shrink trigger distance if box size has changed
   conservative shrink procedure:
     compute distance each of 8 corners of box has moved since last reneighbor
     reduce skin distance by sum of 2 largest of the 8 values
     new trigger = 1/2 of reduced skin distance
   for orthogonal box, only need 2 lo/hi corners
   for triclinic, need all 8 corners since deformations can displace all 8
------------------------------------------------------------------------- */

int Neighbor::check_distance()
{
  double delx,dely,delz,rsq;
  double delta,deltasq,delta1,delta2;

  //if (comm->me==0) fprintf(screen, "boxcheck = %d\n", boxcheck);
  deltasq = triggersq;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = x[i][0] - xhold[i][0];
    dely = x[i][1] - xhold[i][1];
    delz = x[i][2] - xhold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > deltasq) flag = 1;
  }

  double ***nodex =element->nodex;
  nlocal = element->nlocal;
  int npe = element->npe;

  for (int i = 0; i<nlocal;i++) {
   for (int j = 0; j<npe; j++) {
     delx = nodex[i][j][0] - nodexhold[i][j][0];
     dely = nodex[i][j][1] - nodexhold[i][j][1];
     delz = nodex[i][j][2] - nodexhold[i][j][2];
     rsq = delx*delx + dely*dely+delz*delz;
     if (rsq > deltasq) flag = 1;
   }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  //if (comm->me==0) fprintf(logfile, "flag: %d\n", flagall);
  //error->all(FLERR, "test from check_distance in neighbor");
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
  //return 0;
}

/* ----------------------------------------------------------------------
    return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Neighbor::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(xhold,maxhold,3);
  if (style != NSQ) {
    bytes += memory->usage(bins,maxbin);
    bytes += memory->usage(binhead,maxhead);
    bytes += memory->usage(ebins, emaxbin);
    bytes += memory->usage(ebinhead,maxhead);
  }
  for (int i = 0; i < nrequest; i++) 
    if (lists[i]) bytes += lists[i]->memory_usage();

  bytes += nilist->memory_usage();
  bytes += memory->usage(bondlist,maxbond,3);
  bytes += memory->usage(anglelist,maxangle,4);
  bytes += memory->usage(dihedrallist,maxdihedral,5);
  bytes += memory->usage(improperlist,maximproper,5);
  return bytes;
}
