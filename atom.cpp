#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom_vec.h"
#include "atom.h"
#include "style_atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "update.h"
#include "atom_masks.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;
using namespace MathConst;

#define DELTA 1
#define DELTA_MEMSTR 1024
#define EPSILON 1.0e-6
#define CUDA_CHUNK 3000
#define MAXBODY 20       // max # of lines in one body, also in ReadData class

enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/*---------------------------------------------------------------------------------*/

Atom::Atom(CAC *cac) : Pointers(cac)
{
  natoms = 0;
  nlocal = nghost = nmax = 0;
  ntypes = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;

  firstgroupname = NULL;
  sortfreq = 0;
  nextsort = 0;
  userbinsize = 0.0;
  maxbin = maxnext = 0;
  binhead = NULL;
  next = permute = NULL;

  // initialize atom arrays
  // customize by adding new array

  tag = NULL;
  type = mask = NULL;
  image = NULL;
  x = v = f = NULL;

  molecule = NULL;
  molindex = molatom = NULL;
  q = NULL;
  mu = NULL;
  omega = angmom = torque = NULL;
  radius = rmass = NULL;
  ellipsoid = line = tri = body = NULL;

  vfrac = s0 = NULL;
  x0 = NULL;

  spin = NULL;
  eradius = ervel = erforce = NULL;
  cs = csforce = vforce = ervelforce = NULL;
  etag = NULL;

  rho = drho = e = de = cv = NULL;
  vest = NULL;

  // USER-SMD

  contact_radius = NULL;
  smd_data_9 = NULL;
  smd_stress = NULL;
  eff_plastic_strain = NULL;
  eff_plastic_strain_rate = NULL;
  damage = NULL;

  // molecular info

  bond_per_atom =  extra_bond_per_atom = 0;
  num_bond = NULL;
  bond_type = NULL;
  bond_atom = NULL;
  
  angle_per_atom = extra_angle_per_atom = 0;
  num_angle = NULL;
  angle_type = NULL;
  angle_atom1 = angle_atom2 = angle_atom3 = NULL;

  dihedral_per_atom = extra_dihedral_per_atom = 0;
  num_dihedral = NULL;
  dihedral_type = NULL;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = NULL;

  improper_per_atom = extra_improper_per_atom = 0;
  num_improper = NULL;
  improper_type = NULL;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = NULL;

  maxspecial = 1;
  nspecial = NULL;
  special = NULL;

  // user-defined molecules

  nmolecule = 0;
//  molecules = NULL;

  // custom atom arrays

  nivector = ndvector = 0;
  ivector = NULL;
  dvector = NULL;
  iname = dname = NULL;
  
  // initialize atom style and array existence flags
  // customize by adding new flag\

  sphere_flag = peri_flag = electron_flag = 0;
  wavepacket_flag = sph_flag = 0;

  molecule_flag = 0;
  q_flag = mu_flag = 0;
  omega_flag = torque_flag = angmom_flag = 0;
  radius_flag = rmass_flag = 0;
  ellipsoid_flag = line_flag = tri_flag = body_flag = 0;

  vfrac_flag = 0;
  spin_flag = eradius_flag = ervel_flag = erforce_flag = ervelforce_flag = 0;
  cs_flag = csforce_flag = vforce_flag = etag_flag = 0;

  rho_flag = e_flag = cv_flag = vest_flag = 0;

  // USER-SMD

  smd_flag = 0;
  contact_radius_flag = 0;
  smd_data_9_flag = 0;
  smd_stress_flag = 0;
  x0_flag = 0;
  eff_plastic_strain_flag = 0;
  eff_plastic_strain_rate_flag = 0;
  damage_flag = 0;

  // Peridynamic scale factor

  pdscale = 1.0;

  // ntype-length arrays

  mass = NULL;
  mass_setflag = NULL;

  // callback lists & extra restart info
 
  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = NULL;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = NULL;

  //default atom ID and mapping values

  tag_enable = 1;
  map_style = map_user = 0;
  map_tag_max = -1;
  map_maxarray = map_nhash = -1;

  max_same = 0;
  sametag = NULL;
  map_array = NULL;
  map_bucket = NULL;
  map_hash = NULL;

  atom_style = NULL;
  avec = NULL;

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;
}

/*----------------------------------------------------------------------------------*/

Atom::~Atom()
{

  delete [] atom_style;
  delete avec;

  delete [] firstgroupname;
  memory->destroy(binhead);
  memory->destroy(next);
  memory->destroy(permute);

  // delete atom arrays
  // customize by adding new arrays
  
  memory->destroy(tag);
  memory->destroy(type);
  memory->destroy(mask);
  memory->destroy(image);
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);

  memory->destroy(molecule);
  memory->destroy(molindex);
  memory->destroy(molatom);

  memory->destroy(q);
  memory->destroy(mu);
  memory->destroy(omega);
  memory->destroy(angmom);
  memory->destroy(torque);
  memory->destroy(radius);
  memory->destroy(rmass);
  memory->destroy(ellipsoid);
  memory->destroy(line);
  memory->destroy(tri);
  memory->destroy(body);

  memory->destroy(vfrac);
  memory->destroy(s0);
  memory->destroy(x0);

  memory->destroy(spin);
  memory->destroy(eradius);
  memory->destroy(ervel);
  memory->destroy(erforce);
  memory->destroy(ervelforce);
  memory->destroy(cs);
  memory->destroy(csforce);
  memory->destroy(vforce);
  memory->destroy(etag);

  memory->destroy(rho);
  memory->destroy(drho);
  memory->destroy(e);
  memory->destroy(de);
  memory->destroy(cv);
  memory->destroy(vest);

  memory->destroy(contact_radius);
  memory->destroy(smd_data_9);
  memory->destroy(smd_stress);
  memory->destroy(eff_plastic_strain);
  memory->destroy(eff_plastic_strain_rate);
  memory->destroy(damage);

  memory->destroy(nspecial);
  memory->destroy(special);

  memory->destroy(num_bond);
  memory->destroy(bond_type);
  memory->destroy(bond_atom);

  memory->destroy(num_angle);
  memory->destroy(angle_type);
  memory->destroy(angle_atom1);
  memory->destroy(angle_atom2);
  memory->destroy(angle_atom3);

  memory->destroy(num_dihedral);
  memory->destroy(dihedral_type);
  memory->destroy(dihedral_atom1);
  memory->destroy(dihedral_atom2);
  memory->destroy(dihedral_atom3);
  memory->destroy(dihedral_atom4);

  memory->destroy(num_improper);
  memory->destroy(improper_type);
  memory->destroy(improper_atom1);
  memory->destroy(improper_atom2);
  memory->destroy(improper_atom3);
  memory->destroy(improper_atom4);

  // delete custom atom arrays

  for (int i = 0; i < nivector; i++) {
    delete [] iname[i];
    memory->destroy(ivector[i]);
  }
  for (int i = 0; i < ndvector; i++) {
    delete [] dname[i];
    memory->destroy(dvector[i]);
  }

  memory->sfree(iname);
  memory->sfree(dname);
  memory->sfree(ivector);
  memory->sfree(dvector);

  // delete user-defined molecules
  
//  for (int i = 0; i < nmolecule; i++) delete molecules[i];
//  memory->sfree(molecules);

  // delete per-type arrays

  delete [] mass;
  delete [] mass_setflag;

  // delete extra arrays

  memory->destroy(extra_grow);
  memory->destroy(extra_restart);
  memory->destroy(extra_border);
  memory->destroy(extra);
  
  // delete mapping data structures

//  map_delete();
}

/* ----------------------------------------------------------------------
     create an AtomVec style
     called from lammps.cpp, input script, restart file, replicate
  ------------------------------------------------------------------------- */

void Atom::create_avec(const char *style, int narg, char **arg, int trysuffix)
{

  delete [] atom_style;
  if (avec) delete avec;

  // unset atom style and array existence flags
  //may have been set by old avec
  //customize by adding new flag

  sphere_flag = peri_flag = electron_flag = 0;
  wavepacket_flag = sph_flag = 0;

  molecule_flag = 0;
  q_flag = mu_flag = 0;
  omega_flag = torque_flag = angmom_flag = 0;
  radius_flag = rmass_flag = 0;
  ellipsoid_flag = line_flag = tri_flag = body_flag = 0;
  
  vfrac_flag = 0;
  spin_flag = eradius_flag = ervel_flag = erforce_flag = ervelforce_flag = 0;
  cs_flag = csforce_flag = vforce_flag = etag_flag = 0;
  
  rho_flag = e_flag = cv_flag = vest_flag = 0;

  // create instance of AtomVec
  // use grow() to initialize atom-based arrays to length 1
  // so that x[0][0] can always be referenced even if proc has no atoms

  int sflag;
  avec = new_avec(style,trysuffix,sflag);
  avec->store_args(narg,arg);
  avec->process_args(narg,arg); 
  avec->grow(1); 
}

/* ----------------------------------------------------------------------
   generate an AtomVec class, first with suffix appended
  ------------------------------------------------------------------------- */

AtomVec *Atom::new_avec(const char *style, int trysuffix, int &sflag)
{
  sflag = 0;
  if (0) return NULL;
  
#define ATOM_CLASS
#define AtomStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(cac);
#include "style_atom.h"
#undef ATOM_CLASS

  else error->all(FLERR,"Unknown atom style");
  return NULL;
}

/* ----------------------------------------------------------------------
   decrement ptrs in callback lists to fixes beyond the deleted ifix
   happens after fix is deleted
------------------------------------------------------------------------- */

void Atom::update_callback(int ifix)
{
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix) extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix) extra_restart[i]--;
  for (int i = 0; i < nextra_border; i++)
    if (extra_border[i] > ifix) extra_border[i]--;
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int Atom::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy,n,"atom:copy");
  strcpy(copy,line);

  char *ptr;
  if ((ptr = strchr(copy,'#'))) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->destroy(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->destroy(copy);
  return n;
}

/* ----------------------------------------------------------------------
  allocate arrays of length ntypes
  only done after ntypes is set
 ------------------------------------------------------------------------- */

void Atom::allocate_type_arrays()
{
  if (avec->mass_type) {
    mass = new double[ntypes+1];
    mass_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
  }
}

/*-----------------------------------------------------------------------
 * Since I don't have bond, angle, dihedral and improper types in my simulations
 * so this deallocate function is empty
-------------------------------------------------------------------------*/

void Atom::deallocate_topology()
{
}

/* ----------------------------------------------------------------------
   unpack n lines from Atom section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_atoms(int n, char *buf, tagint id_offset, int type_offset, 
                      int shiftflag, double *shift)
{
  int m,xptr,iptr;
  imageint imagedata;
  double xdata[3],lamda[3];
  double *coord;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_atom && nwords != avec->size_data_atom + 3)
    error->all(FLERR,"Incorrect atom format in data file");
  
  char **values = new char*[nwords];

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }
  
  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

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

  // xptr = which word in line starts xyz coords
  // iptr = which word in line starts ix,iy,iz image flags

  xptr = avec->xcol_data - 1;
  int imageflag = 0;
  if (nwords > avec->size_data_atom) imageflag = 1;
  if (imageflag) iptr = nwords - 3;

  // loop over lines of atom data
  // tokenize the line into values
  // extract xyz coords and image flags
  // remap atom into simulation box

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR,"Incorrect atom format in data file");
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
    if (shiftflag) {
      xdata[0] += shift[0];
      xdata[1] += shift[1];
      xdata[2] += shift[2];
    }

    domain->remap(xdata,imagedata);
    if (triclinic) {
      domain->x2lamda(xdata,lamda);
      coord = lamda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      avec->data_atom(xdata,imagedata,values);
      if (id_offset) tag[nlocal-1] += id_offset;
      if (type_offset) {
        type[nlocal-1] += type_offset;
        if (type[nlocal-1] > ntypes)
          error->one(FLERR,"Invalid atom type in Atoms section of data file");
      }
    }
   
    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
    set a mass and flag it as set
   called from reading of data file
   type_offset may be used when reading multiple data files
------------------------------------------------------------------------- */

void Atom::set_mass(const char *str, int type_offset)
{
  if (mass == NULL) error->all(FLERR,"Cannot set mass for this atom style");

  int itype;
  double mass_one;
  int n = sscanf(str,"%d %lg",&itype,&mass_one);
  if (n != 2) error->all(FLERR,"Invalid mass line in data file");
  itype += type_offset;

  if (itype < 1 || itype > ntypes)
    error->all(FLERR,"Invalid type for mass set");

  mass[itype] = mass_one;
  mass_setflag[itype] = 1;

  if (mass[itype] <= 0.0) error->all(FLERR,"Invalid mass value");

}

/*-----------------------------------------------------------------------
  set a mass and flag it as set
  called from EAM pair routine
-------------------------------------------------------------------------*/

void Atom::set_mass(int itype, double value)
{
  if (mass == NULL) error->all(FLERR,"Cannot set mass for this atom style");
  if (itype < 1 || itype > ntypes)
    error->all(FLERR,"Invalid type for mass set");

  mass[itype] = value;
  mass_setflag[itype] = 1;

  if (mass[itype] <=0.0) error->all(FLERR,"Invalid mass value");
}

/* ---------------------------------------------------------------------- */

void Atom::init()
{
  // delete extra array since it doesn't persist past first run

  if (nextra_store) {
    memory->destroy(extra);
    extra = NULL;
    nextra_store = 0;
  }

  // check arrays that are atom type in length

  check_mass();

  // setup of firstgroup

  if (firstgroupname) {
    firstgroup = group->find(firstgroupname);
    if (firstgroup < 0)
      error->all(FLERR,"Could not find atom_modify first group ID");
  } else firstgroup = -1;

  // init AtomVec
  
  avec->init();

}

/* ----------------------------------------------------------------------
   check that all masses have been set
------------------------------------------------------------------------- */

void Atom::check_mass()
{
  if (mass == NULL) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0) error->all(FLERR,"All masses are not set");
}

/* ---------------------------------------------------------------------- */

void Atom::setup()
{
  // setup bins for sorting
  // cannot do this in init() because uses neighbor cutoff

  //if (comm->me==0) fprintf(screen, "sortfreq = %d\n", sortfreq);
  if (sortfreq > 0) setup_sort_bins();
}

/* ----------------------------------------------------------------------
   setup bins for spatial sorting of atoms
------------------------------------------------------------------------- */

void Atom::setup_sort_bins()
{
  // binsize:
  // user setting if explicitly set
  // 1/2 of neighbor cutoff for non-CUDA
  // CUDA_CHUNK atoms/proc for CUDA
  // check if neighbor cutoff = 0.0

  double binsize;
  if (userbinsize > 0.0) binsize = userbinsize;
  else binsize = 0.5 * neighbor->cutneighmax;

  if (binsize == 0.0) error->all(FLERR,"Atom sorting has bin size = 0.0");

  double bininv = 1.0/binsize;

  //if(comm->me==0) fprintf(screen, "binsize = %f\n", binsize);
  // nbin xyz = local bins
  // bbox lo/hi = bounding box of my sub-domain

  bboxlo[0] = domain->sublo[0];
  bboxlo[1] = domain->sublo[1];
  bboxlo[2] = domain->sublo[2];
  bboxhi[0] = domain->subhi[0];
  bboxhi[1] = domain->subhi[1];
  bboxhi[2] = domain->subhi[2];

  nbinx = static_cast<int> ((bboxhi[0]-bboxlo[0]) * bininv);
  nbiny = static_cast<int> ((bboxhi[1]-bboxlo[1]) * bininv);
  nbinz = static_cast<int> ((bboxhi[2]-bboxlo[2]) * bininv);
  if (domain->dimension == 2) nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  bininvx = nbinx / (bboxhi[0]-bboxlo[0]);
  bininvy = nbiny / (bboxhi[1]-bboxlo[1]);
  bininvz = nbinz / (bboxhi[2]-bboxlo[2]);

  if (1.0*nbinx*nbiny*nbinz > INT_MAX)
    error->one(FLERR,"Too many atom sorting bins");

  nbins = nbinx*nbiny*nbinz;

  // reallocate per-bin memory if needed
  
  if (nbins > maxbin) {
    memory->destroy(binhead);
    maxbin = nbins;
    memory->create(binhead,maxbin,"atom:binhead");
  }
}

/* ----------------------------------------------------------------------
   perform spatial sort of atoms within my sub-domain
   always called between comm->exchange() and comm->borders()
   don't have to worry about clearing/setting atom->map since done in comm
------------------------------------------------------------------------- */

void Atom::sort()
{
  int i,m,n,ix,iy,iz,ibin,empty;

  // set next timestep for sorting to take place
   
  nextsort = (update->ntimestep/sortfreq)*sortfreq + sortfreq;
 
  //if (comm->me==0) fprintf(screen, "nextsort = %d\n", nextsort);

  // re-setup sort bins if needed

  //if (comm->me==0) fprintf(screen, "box_change = %d\n", domain->box_change);

  if (nbins == 1) return;

  // reallocate per-atom vectors if needed

  if (nlocal > maxnext) {
    memory->destroy(next);
    memory->destroy(permute);
    maxnext = atom->nmax;
    memory->create(next,maxnext,"atom:next");
    memory->create(permute,maxnext,"atom:permute");
  } 
  //if (comm->me==0) fprintf(screen, "nmax = %d\n", atom->nmax); 

  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  // bin atoms in reverse order so linked list will be in forward order

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  for (i = nlocal-1; i >= 0; i--) {
    ix = static_cast<int> ((x[i][0]-bboxlo[0])*bininvx);
    iy = static_cast<int> ((x[i][1]-bboxlo[1])*bininvy);
    iz = static_cast<int> ((x[i][2]-bboxlo[2])*bininvz);
    ix = MAX(ix,0);
    iy = MAX(iy,0);
    iz = MAX(iz,0);
    ix = MIN(ix,nbinx-1);
    iy = MIN(iy,nbiny-1);
    iz = MIN(iz,nbinz-1);
    ibin = iz*nbiny*nbinx + iy*nbinx + ix;
    next[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  // permute = desired permutation of atoms
  // permute[I] = J means Ith new atom will be Jth old atom

  n = 0;
  for (m = 0; m < nbins; m++) {
    i = binhead[m];
    while (i >= 0) {
      permute[n++] = i;
      i = next[i];
    }
  }

  // current = current permutation, just reuse next vector
  // current[I] = J means Ith current atom is Jth old atom

  int *current = next;
  for (i = 0; i < nlocal; i++) current[i] = i;

  // reorder local atom list, when done, current = permute
  // perform "in place" using copy() to extra atom location at end of list
  // inner while loop processes one cycle of the permutation
  // copy before inner-loop moves an atom to end of atom list
  // copy after inner-loop moves atom at end of list back into list
  // empty = location in atom list that is currently empty

  for (i = 0; i < nlocal; i++) {
    if (current[i] == permute[i]) continue;
    avec->copy(i,nlocal,0);
    empty = i;
    while (permute[empty] != i) {
      avec->copy(permute[empty],empty,0);
      empty = current[empty] = permute[empty];
    }
    avec->copy(nlocal,empty,0);
    current[empty] = permute[empty];
  }  
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   call to avec tallies per-atom vectors
   add in global to local mapping storage
------------------------------------------------------------------------- */

bigint Atom::memory_usage()
{
  memlength = DELTA_MEMSTR;
  memory->create(memstr,memlength,"atom:memstr");
  memstr[0] = '\0';
  bigint bytes = avec->memory_usage();
  memory->destroy(memstr);

  bytes += max_same*sizeof(int);
  if (map_style == 1)
    bytes += memory->usage(map_array,map_maxarray);
  else if (map_style == 2) {
    bytes += map_nbucket*sizeof(int);
    bytes += map_nhash*sizeof(HashElem);
  }
  if (maxnext) {
   // if (comm->me==0) fprintf(screen, "testing from atom memory_usage");
    bytes += memory->usage(next,maxnext);
    bytes += memory->usage(permute,maxnext);
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   accumulate per-atom vec names in memstr, padded by spaces
   return 1 if padded str is not already in memlist, else 0
------------------------------------------------------------------------- */

int Atom::memcheck(const char *str)
{
  int n = strlen(str) + 3;
  char *padded = new char[n];
  strcpy(padded," ");
  strcat(padded,str);
  strcat(padded," ");
  
  if (strstr(memstr,padded)) {
    delete [] padded;
    return 0;
  }

  if (strlen(memstr) + n >= memlength) {
    memlength += DELTA_MEMSTR;
    memory->grow(memstr,memlength,"atom:memstr");
  }

  strcat(memstr,padded);
  delete [] padded;
  return 1;
}
