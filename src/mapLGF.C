/*
  Program: mapLGF.C
  Author:  D. Trinkle
  Date:    2007 June 27
  Purpose: Read in:
	   1. Cell info
           2. Set of positions (cartesian) with regional "identity"; threading direction
	   3. LGF

	   and we map the LGF to the positions, to create an easy-to-use
	   LGF for direct application in other codes, regardless of induced PBC.
	   The regional "identity" is whether an atom is in the region to be updated
	   or not; we only output LGF's where at least one atom is in the update
	   region.  If a threading direction is included, this allows the use of a 
	   2D LGF; else, it will do a 3D LGF update.  There is no assumption
	   of periodicity, so this effectively couples the system to bulk.

  Param.:  <cell> <pos-ident> <LGF>
           cell:      cell file describing our lattice
	   pos-ident: atomic positions in cartesian coordl. and atom identity
	   LGF:       lattice GF

           ==== cell ====
           a0                            # Scale factor for unit cell
           a1.x a1.y a1.z                # Cartesian coord of unit cell
           a2.x a2.y a2.z
           a3.x a3.y a3.z
           crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
           Natoms                        # Number of atoms in first unit cell
           u1.1 u1.2 u1.3                # Atom locations, in direct coord.
           ...
           uN.1 uN.2 uN.3
           ==== cell ====

           ==== pos-ident ====
           N [t1 t2 t3] # number of points, optional threading direction
           r1.x r1.y r1.z ident1  # position in cart. coord, ident
           ...                    # ident=1 == in "region 2" (to be updated)
           rN.x rN.y rN.z identN  # ident=0    else
           ==== pos-ident ====

	   ==== LGF ====
	   N [t1] [t2] [t3]           # number of points + threading direction
	   u1.1 u1.2 u1.3 Gxx .. Gzz  # position (unit) and LGF
	   ...
	   uN.1 uN.2 uN.3 Gxx .. Gzz
	   ==== LGF ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   

  Output:  We output the LGF in a "ready-to-use" format based on atom indices
           instead:

	   ==== LGF.mapped ====
	   Nupdate N  # number of atoms in update region, number of atoms total
	   i_update j  Gxx .. Gzz  # atom index in update, index in remainder, G
	   ...
	   ==== LGF.mapped ====

	   The idea is that this format, while more verbose, is irrespective
	   of any further PBC imposed on the system, and is directly applicable
	   in other codes (like VASP).
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "Dij.H"
#include "supercell.H" // indexing routines for easy lookup


//***************************** SUBROUTINES ****************************

// structure to hold atomic positions, forces, and displacements:
typedef struct 
{
  double R[3];
  int ident;
} posident_type;

int read_posident (FILE* infile, int& Np, posident_type* &p,
		   int t_unit[3], int THREAD);

inline void setmax (int &a, const int &b) 
{if (a<b) a=b;}

inline void print_mat (double mat[9]) 
{
  printf("%.6le %.6le %.6le %.6le %.6le %.6le %.6le %.6le %.6le",
	 mat[0], mat[1], mat[2],
	 mat[3], mat[4], mat[5], 
	 mat[6], mat[7], mat[8]);
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <pos-ident> <LGF>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'w', 'n'}; // flag characters.

const char* ARGEXPL = 
"  cell:      cell file describing our lattice\n\
  pos-ident: atomic positions in cartesian coord. and identity (1=LGF, 0=otherwise)\n\
  LGF:       lattice GF\n\
  -w         issue \"verbose\" warnings if LGF table entries not found\n\
  -n         no threading (i.e., 3D LGF; default is periodic 2d LGF)";

const char* FILEEXPL =
"==== cell ====\n\
a0                            # Scale factor for unit cell\n\
a1.x a1.y a1.z                # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.\n\
Natoms                        # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3                # Atom locations, in direct coord.\n\
...\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== pos-ident ====\n\
N [t1 t2 t3] # number of points, optional threading direction\n\
r1.x r1.y r1.z ident1  # position in cart. coord, ident\n\
...                    # ident=1 == in \"region 2\" (to be updated)\n\
rN.x rN.y rN.z identN  # ident=0    else\n\
==== pos-ident ====\n\
\n\
==== LGF ====\n\
N                          # number of points\n\
u1.1 u1.2 u1.3 Gxx .. Gzz  # position (unit) and LGF\n\
...\n\
uN.1 uN.2 uN.3 Gxx .. Gzz\n\
==== LGF ====\n";

int main ( int argc, char **argv ) 
{
  //  int d, i, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (int i=0; i<NFLAGS; ++i) flagon[i] = 0;
  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  int WARNING = (flagon[0]) || VERBOSE;  // turn on warning for verbose.
  int THREAD = ! flagon[1]; // by default, we thread... this may change later

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char* cell_name = args[0];
  char* posident_name = args[1];
  char* LGF_name = args[2];

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input -- we read it, but don't use it.
  int Natoms = NO_ATOMS;
  double** u_atoms = NULL;
  double atomic_mass;

  //++ ==== cell ====
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms);
  // Read in the atomic mass:
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%lf", &atomic_mass);
  myclose(infile);
  //-- ==== cell ====
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }
  if (Natoms != 1) {
    fprintf(stderr, "Sorry.  Currently we can only do single atom cells.\n");
    exit(1);
  }

  //++ ==== posident ====
  int Np, t_unit[3]; // optional threading direction
  posident_type* p;

  infile = myopenr(posident_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", posident_name);
    exit(ERROR_NOFILE);
  }
  // don't bother reading forces if we're supposed to generate forces
  ERROR = read_posident(infile, Np, p, t_unit, THREAD);
  myclose(infile);
  //-- ==== posident ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", posident_name);
    exit(ERROR);
  }

  //++ ==== LGF ====
  int Nlatt;
  point_type *LGF; // set of points

  infile = myopenr(LGF_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", LGF_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Nlatt, LGF);
  myclose(infile);
  //-- ==== LGF ====

  if (ERROR) {
    delete LGF; 
    fprintf(stderr, "Error with %s\n", LGF_name);
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  // 0. Setup indexing routines for LGF--this is only to speed up the
  // implementation; also handles threading in an semi-elegant fashion
  int*** index, maxn=0;
  // How many periodic images do we include for threading?  Seems like
  // 1 should suffice, but for dislocations... it's more like 2.
  int THREAD_RANGE = 2;
  if (!THREAD) THREAD_RANGE = 0;
  for (int n=0; n<Nlatt; ++n) {
    point_type *tp = LGF+n;
    for (int t=-THREAD_RANGE; t<=THREAD_RANGE; ++t)
      for (int d=0; d<3; ++d) 
	setmax(maxn, abs(tp->Runit[d] + t*t_unit[d]));
  }
  make_index(maxn, index);
  for (int n=0; n<Nlatt; ++n) {
    int r[3];
    point_type *tp = LGF+n;
    // if we're threading, we need to include "periodic images"
    for (int t=-THREAD_RANGE; t<=THREAD_RANGE; ++t) {
      for (int d=0; d<3; ++d) r[d] = tp->Runit[d] + t*t_unit[d];
      set_index(maxn, index, r, n);
    }
  }
  
  // 1. Now we iterate over all of our position pairs
  //  int N2 = 0;
  //  for (int np=0; np<Np; ++np) N2 += p[np].ident;
  //  printf("%d %d  # N_update  N_total\n", N2, Np);
  int nlo=Np, nhi=-1, ncount=0;
  for (int np=0; np<Np; ++np) {
    if ( p[np].ident ) {
      ++ncount;
      if (np<nlo) nlo = np;	// get the lowest point
      if (np>nhi) nhi = np;	// get the highest point
    }
  }
  printf("LGFCAR indexing: Region2(start) Region2(finish) UpdateAtom(start) UpdateAtom(finish) TotalEntries\n");
  printf("%d %d %d %d %d\n", nlo+1, nhi+1, 1, Np, Np*ncount);
  
  double cart_inv[9];
  int u[3];
  careful_inverse(cart, cart_inv);
  for (int np0=0; np0<Np; ++np0) {
    posident_type *tp0 = p + np0;
    if ( tp0->ident )  // do we have an atom in our update region?
      for (int np1=0; np1<Np; ++np1) {
	posident_type *tp1 = p + np1;
	double dR[3], du[3];
	for (int d=0; d<3; ++d) dR[d] = tp1->R[d] - tp0->R[d];
	mult_vect(cart_inv, dR, du);
	for (int d=0; d<3; ++d) u[d] = lround(du[d]);
	int n = get_index(maxn, index, u); // get the LGF element
	if (n != -1) {
	  // output our element:
	  printf("%3d %3d ", np0+1, np1+1); // increment indices by 1, each
	  print_mat(LGF[n].mat);
	  printf("\n");
	}
	else {
	  if (WARNING) {
	    fprintf(stderr, "# R0(%d)= %.8lf %.8lf %.8lf  R1(%d)= %.8lf %.8lf %.8lf\n",
		    np0+1, tp0->R[0], tp0->R[1], tp0->R[2],
		    np1+1, tp1->R[0], tp1->R[1], tp1->R[2]);
	    fprintf(stderr, "# dR= %.8lf %.8lf %.8lf  u= %d %d %d\n",
		   dR[0],dR[1],dR[2], u[0],u[1],u[2]);
	    fprintf(stderr, "# Could not find entry for u in LGF.\n");
	  }
	}
      }
  }
  // garbage collection:
  free_index(maxn, index);

  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  delete[] LGF;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


int read_posident (FILE* infile, int& Np, posident_type* &p,
		   int t_unit[3], int THREAD) 
{
  int ERROR = 0;
  char dump[512];

  nextnoncomment(infile, dump, sizeof(dump));
  if (THREAD) {
    int i = sscanf(dump, "%d %d %d %d", &Np, t_unit, t_unit+1, t_unit+2);
    if ( (i != 4) || zero(t_unit) ) {
      fprintf(stderr, "Asked for threading, but bad threading direction in read_posident\n");
      return ERROR_BADFILE;
    }
  }
  else
    sscanf(dump, "%d", &Np);

  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in read_posident\n", Np);
    return ERROR_BADFILE;
  }
  //    if (p != NULL) delete[] p;
  p = new posident_type[Np];
  if (p == NULL) {
    fprintf(stderr, "Error allocating memory...?\n");
    return ERROR_MEMORY;
  }

  int n;
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    posident_type *tp = p + n; // this point
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) break;
    // Parse it.
    int i = sscanf(dump, "%lf %lf %lf %d", 
		   &(tp->R[0]), &(tp->R[1]), &(tp->R[2]), &(tp->ident));
    ERROR = (i != 4);
    if (!ERROR) ERROR = ( (tp->ident<0) || (tp->ident>1) );
  }
  if (ERROR) {
    fprintf(stderr, "Bad read line in read_posident\n");
  }
  if (n!=Np) {
    fprintf(stderr, "Not enough points in read_posident\n");
    ERROR = ERROR_BADFILE;
  }
  return ERROR;
}
