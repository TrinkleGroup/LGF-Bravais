/*
  Program: applyLGF.C
  Author:  D. Trinkle
  Date:    October 8, 2004
  Purpose: Read in:
	   1. Cell info
           2. Set of positions (cartesian) with forces; threading direction
	   3. LGF / dynamical matrix

	   and we apply the LGF to the positions (or dynamical matrix).
	   If a threading direction is included, this allows the use of a 
	   2D LGF; else, it will do a 3D LGF update.  There is no assumption
	   of periodicity, so this effectively couples the system to bulk.

  Param.:  <cell> <pos-for> <LGF>
           cell:    cell file describing our lattice
	   pos-for: atomic positions and forces in cartesian coord.
	   LGF:     lattice GF

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

           ==== pos-for ====
           N [t1 t2 t3] # number of points, optional threading direction
           r1.x r1.y r1.z f1.x f1.y f1.z  # position, force in cart. coord
           ...
           rN.x rN.y rN.z fN.x fN.y fN.z
           ==== pos-for ====

	   ==== LGF ====
	   N                          # number of points
	   u1.1 u1.2 u1.3 Gxx .. Gzz  # position (unit) and LGF
	   ...
	   uN.1 uN.2 uN.3 Gxx .. Gzz
	   ==== LGF ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   

  Output:  We output the new positions in cartesian coord. from the update.
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
  double R[3], f[3], u[3];
} posfor_type;

int read_posfor (FILE* infile, int& Np, posfor_type* &p,
		 int t_unit[3], int THREAD, int READ_FORCES);

inline void setmax (int &a, const int &b) 
{if (a<b) a=b;}



/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <pos-for> <LGF>";

const int NFLAGS = 4;
const char USERFLAGLIST[NFLAGS] = {'w', 'n', 'd', 'e'}; // flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  pos-for: atomic positions and forces in cartesian coord.\n\
  LGF:     lattice GF / dynamical matrix (see -d switch below)\n\
  -w       issue \"verbose\" warnings if LGF table entries not found\n\
  -n       no threading (i.e., 3D LGF; default is periodic 2d LGF)\n\
  -d       calculate forces from dynamical matrix\n\
  -e       calculate the energy change (work done) by update\n\
  NOTE: -d and -e cannot be both chosen";

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
==== pos-for ====\n\
N [t1 t2 t3] # number of points, optional threading direction\n\
r1.x r1.y r1.z f1.x f1.y f1.z  # position, force in cart. coord\n\
...\n\
rN.x rN.y rN.z fN.x fN.y fN.z\n\
==== pos-for ====\n\
\n\
==== LGF ====\n\
N                          # number of points\n\
u1.1 u1.2 u1.3 Gxx .. Gzz  # position (unit) and LGF\n\
...\n\
uN.1 uN.2 uN.3 Gxx .. Gzz\n\
==== LGF ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (i=0; i<NFLAGS; ++i) flagon[i] = 0;
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
  int FORCE = flagon[2];    // default is LGF update, but we can do forces too
  int ENERGY = flagon[3];
  if (FORCE && ENERGY) {
    fprintf(stderr, "Sorry... cannot compute forces from displacement *and* work done by update.\n");
    exit(ERROR_BADFLAG);
  }

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char* cell_name = args[0];
  char* posfor_name = args[1];
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

  //++ ==== posfor ====
  int Np, t_unit[3]; // optional threading direction
  posfor_type* p;

  infile = myopenr(posfor_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", posfor_name);
    exit(ERROR_NOFILE);
  }
  // don't bother reading forces if we're supposed to generate forces
  ERROR = read_posfor(infile, Np, p, t_unit, THREAD, ! FORCE);
  myclose(infile);
  //-- ==== posfor ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", posfor_name);
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
  int THREAD_RANGE = 2, t;
  if (!THREAD) THREAD_RANGE = 0;
  for (n=0; n<Nlatt; ++n) {
    point_type *tp = LGF+n;
    for (t=-THREAD_RANGE; t<=THREAD_RANGE; ++t)
      for (d=0; d<3; ++d) 
	setmax(maxn, abs(tp->Runit[d] + t*t_unit[d]));
  }
  make_index(maxn, index);
  for (n=0; n<Nlatt; ++n) {
    int r[3];
    point_type *tp = LGF+n;
    // if we're threading, we need to include "periodic images"
    for (t=-THREAD_RANGE; t<=THREAD_RANGE; ++t) {
      for (d=0; d<3; ++d) r[d] = tp->Runit[d] + t*t_unit[d];
      set_index(maxn, index, r, n);
    }
  }
  
  // 1. Now we iterate over all of our positions and forces...
  double dR[3], cart_inv[9], du[3], f_inc[3];
  int u[3];
  careful_inverse(cart, cart_inv);
  for (np=0; np<Np; ++np) {
    posfor_type *tp0 = p + np;
    if ( (! zero_vect(tp0->f)) || FORCE )
      for (int np1=0; np1<Np; ++np1) {
	posfor_type *tp1 = p + np1;
	for (d=0; d<3; ++d) dR[d] = tp1->R[d] - tp0->R[d];
	mult_vect(cart_inv, dR, du);
	for (d=0; d<3; ++d) u[d] = lround(du[d]);
	// determine displacement...
	if (FORCE) {
	  // use f_inc as a temporary variable:
	  mult_vect(cart, u, f_inc);
	  for (d=0; d<3; ++d) du[d] = dR[d] - f_inc[d];
	}
	n = get_index(maxn, index, u); // get the LGF element
	if (n != -1) {
	  if (FORCE) {
	    // update to displacement:
	    mult_vect(LGF[n].mat, du, f_inc);
	    // I have a question about the sign... but this seems to give
	    // correct results
	    for (d=0; d<3; ++d) tp0->f[d] -= f_inc[d];
	  }
	  else {
	    // update to displacement:
	    mult_vect(LGF[n].mat, tp0->f, du);
	    for (d=0; d<3; ++d) tp1->u[d] += du[d];
	  }
	}
	else {
	  if (WARNING) {
	    fprintf(stderr, "# R0(%d)= %.8lf %.8lf %.8lf  R1(%d)= %.8lf %.8lf %.8lf\n",
		    np+1, tp0->R[0], tp0->R[1], tp0->R[2],
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

  // 2. Make sure that the total displacement doesn't change the COM:
  if (! FORCE) {
    double ucom[3] = {0,0,0};
    for (np=0; np<Np; ++np)
      for (d=0; d<3; ++d) ucom[d] += p[np].u[d];
    for (d=0; d<3; ++d) ucom[d] *= -1./(double)Np;
    for (np=0; np<Np; ++np)
      for (d=0; d<3; ++d) p[np].u[d] += ucom[d];
  }
  if (ENERGY) {
    // compute and output work done.
    double W;
    W = 0;
    for (np=0; np<Np; ++np) W += dot(p[np].f, p[np].u);
    W *= 0.5; // since we're supposedly in the quadratic regime
    fprintf(stderr, "work_done= %.8le\n", W);
  }
  
  // ****************************** OUTPUT ***************************
  if (THREAD)
    printf("%d %d %d %d\n", Np, t_unit[0], t_unit[1], t_unit[2]);
  else
    printf("%d\n", Np);
  for (np=0; np<Np; ++np) {
    posfor_type* tp = p+np;
    if (FORCE)
      printf("%18.12lf %18.12lf %18.12lf  %15.8le %15.8le %15.8le\n",
	     tp->R[0],tp->R[1],tp->R[2], tp->f[0], tp->f[1], tp->f[2]);
    else
      printf("%18.12lf %18.12lf %18.12lf\n", tp->R[0]+tp->u[0],
	     tp->R[1]+tp->u[1], tp->R[2]+tp->u[2]);
  }
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  delete[] LGF;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


int read_posfor (FILE* infile, int& Np, posfor_type* &p,
		 int t_unit[3], int THREAD, int READ_FORCES) 
{
  int n, i, icount;
  int ERROR = 0;
  char dump[512];

  nextnoncomment(infile, dump, sizeof(dump));
  if (THREAD) {
    i = sscanf(dump, "%d %d %d %d", &Np, t_unit, t_unit+1, t_unit+2);
    if ( (i != 4) || zero_vect(t_unit) ) {
      fprintf(stderr, "Asked for threading, but bad threading direction in read_posfor\n");
      return ERROR_BADFILE;
    }
  }
  else
    sscanf(dump, "%d", &Np);

  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in read_posfor\n", Np);
    return ERROR_BADFILE;
  }
  //    if (p != NULL) delete[] p;
  p = new posfor_type[Np];
  if (p == NULL) {
    fprintf(stderr, "Error allocating memory...?\n");
    return ERROR_MEMORY;
  }

  if (READ_FORCES) icount = 6;
  else             icount = 3;
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    posfor_type *tp = p + n; // this point
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) break;
    // Parse it.
    if (READ_FORCES)
      i = sscanf(dump, "%lf %lf %lf %lf %lf %lf", 
		 &(tp->R[0]), &(tp->R[1]), &(tp->R[2]),
		 &(tp->f[0]), &(tp->f[1]), &(tp->f[2]));
    else
      i = sscanf(dump, "%lf %lf %lf", 
		 &(tp->R[0]), &(tp->R[1]), &(tp->R[2]));
    if (i != icount) ERROR = ERROR_BADFILE;
    for (int d=0; d<3; ++d) tp->u[d] = 0;
    if (! READ_FORCES)
      for (int d=0; d<3; ++d) tp->f[d] = 0;
  }
  if (ERROR) {
    fprintf(stderr, "Bad read line in read_posfor\n");
  }
  if (n!=Np) {
    fprintf(stderr, "Not enough points in read_posfor\n");
    ERROR = ERROR_BADFILE;
  }
  return ERROR;
}
