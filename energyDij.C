/*
  Program: energyDij.C
  Author:  D. Trinkle
  Date:    December 15, 2004
  Purpose: Read in:
	   1. Cell info
           2. Set of positions (cartesian)
	   3. Dynamical matrix
	   4. Supercell info

	   and we apply Dij to the positions (or dynamical matrix), to compute
	   the harmonic energy.

  Param.:  <cell> <pos> <Dij> <super>
           cell:    cell file describing our lattice
	   pos:     atomic positions and forces in cartesian coord.
	   Dij:     dynamical matrix
	   super:   supercell definition

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

           ==== pos ====
           N              # number of points, optional threading direction
           r1.x r1.y r1.z # position
           ...
           rN.x rN.y rN.z
           ==== pos ====

	   ==== Dij ====
	   N                          # number of points
	   u1.1 u1.2 u1.3 Dxx .. Dzz  # position (unit) and Dij
	   ...
	   uN.1 uN.2 uN.3 Dxx .. Dzz
	   ==== Dij ====


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

const int NTHREAD = 27;

void make_thread(int super[9], int thread[NTHREAD][3]) 
{
  int t = 0;
  for (int i=-1; i<=1; ++i)
    for (int j=-1; j<=1; ++j)
      for (int k=-1; k<=1; ++k, ++t)  // we increment t here!
	for (int d=0; d<3; ++d)
	  thread[t][d] = super[3*d]*i + super[3*d+1]*j + super[3*d+2]*k;
}


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
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <pos> <Dij> <super>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  pos:     atomic positions and forces in cartesian coord.\n\
  Dij:     dynamical matrix\n\
  super:   supercell definition";

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
==== pos ====\n\
N              # number of points, optional threading direction\n\
r1.x r1.y r1.z # position\n\
...\n\
rN.x rN.y rN.z\n\
==== pos ====\n\
\n\
==== Dij ====\n\
N                          # number of points\n\
u1.1 u1.2 u1.3 Dxx .. Dzz  # position (unit) and Dij\n\
...\n\
uN.1 uN.2 uN.3 Dxx .. Dzz\n\
==== Dij ====\n";

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char* cell_name = args[0];
  char* pos_name = args[1];
  char* Dij_name = args[2];
  char* super_name = args[3];

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
  nextnoncomment(infile, dump, sizeof(dump));
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

  //++ ==== pos ====
  int Np, t_unit[3]; // optional threading direction
  posfor_type* p;

  infile = myopenr(pos_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", pos_name);
    exit(ERROR_NOFILE);
  }
  // don't bother reading forces if we're supposed to generate forces
  ERROR = read_posfor(infile, Np, p, t_unit, 0, 0);
  myclose(infile);
  //-- ==== pos ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", pos_name);
    exit(ERROR);
  }

  //++ ==== Dij ====
  int Nlatt;
  point_type *Dij; // set of points

  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Nlatt, Dij);
  myclose(infile);
  //-- ==== Dij ====

  if (ERROR) {
    delete Dij; 
    fprintf(stderr, "Error with %s\n", Dij_name);
    exit(ERROR);
  }
  
  //++ ==== super ====
  int super[9];
  
  infile = myopenr(super_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", super_name);
    exit(ERROR_NOFILE);
  }
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d %d %d %d %d %d %d", super, super+1, super+2,
	 super+3, super+4, super+5, super+6, super+7, super+8);
  myclose(infile);
  //-- ==== super ====


  // ***************************** ANALYSIS **************************
  // 0. Setup indexing routines for Dij--this is only to speed up the
  // implementation; also handles periodic BC in an semi-elegant fashion
  int*** index, maxn=0;
  int thread[NTHREAD][3], t;
  make_thread(super, thread);

  for (n=0; n<Nlatt; ++n) {
    point_type *tp = Dij+n;
    for (t=0; t<NTHREAD; ++t)
      for (d=0; d<3; ++d) 
	setmax(maxn, abs(tp->Runit[d] + thread[t][d]));
  }
  make_index(maxn, index);
  for (n=0; n<Nlatt; ++n) {
    int r[3];
    point_type *tp = Dij+n;
    // we need to include "periodic images"
    for (t=0; t<NTHREAD; ++t) {
      for (d=0; d<3; ++d) r[d] = tp->Runit[d] + thread[t][d];
      set_index(maxn, index, r, n);
    }
  }
  
  // 1. Now we iterate over all of our positions and forces...
  double dR[3], cart_inv[9], du[3], f_inc[3];
  int u[3];
  careful_inverse(cart, cart_inv);

  // determine u (displacement):
  for (np=0; np<Np; ++np) {
    posfor_type *tp0 = p + np;
    mult_vect(cart_inv, tp0->R, du);
    for (d=0; d<3; ++d) u[d] = lround(du[d]);
    // use f_inc as a temporary variable:
    mult_vect(cart, u, f_inc);
    for (d=0; d<3; ++d) tp0->u[d] = tp0->R[d] - f_inc[d];
  }

  for (np=0; np<Np; ++np) {
    posfor_type *tp0 = p + np;
    for (int np1=0; np1<Np; ++np1) {
      posfor_type *tp1 = p + np1;
      for (d=0; d<3; ++d) dR[d] = tp1->R[d] - tp0->R[d];
      mult_vect(cart_inv, dR, du);
      for (d=0; d<3; ++d) u[d] = lround(du[d]);
      // determine displacement...
      // use f_inc as a temporary variable:
      mult_vect(cart, u, f_inc);
      for (d=0; d<3; ++d) du[d] = dR[d] - f_inc[d];
      n = get_index(maxn, index, u); // get the Dij element
      if (n != -1) {
	// update to displacement:
	mult_vect(Dij[n].mat, du, f_inc);
	// I have a question about the sign... but this seems to give
	// correct results
	for (d=0; d<3; ++d) tp0->f[d] += f_inc[d];
      } else {
	fprintf(stderr, "Something terrible happened...\n  your positions, supercell, and Dij may not be compatible?\n");
      }
    }
  }
  // garbage collection:
  free_index(maxn, index);

  // compute and output work done.
  double W;
  W = 0;
  for (np=0; np<Np; ++np) W += dot(p[np].f, p[np].u);
  W *= 0.5; // since we're supposedly in the quadratic regime
  
  // ****************************** OUTPUT ***************************
  printf("%12le\n", W);
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  delete[] Dij;
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
    if ( (i != 4) || zero(t_unit) ) {
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
