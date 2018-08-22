/*
  Program: extract-Dij.C
  Author:  D. Trinkle
  Date:    August 24, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. 3 OUTCAR files describing atom positions and forces

	   and we compute
	   1. the dynamical matrix evaluated at points in OUTCAR file

  Param.:  <cell> <OUTCAR.1> <OUTCAR.2> <OUTCAR.3>
           cell:     cell file describing our lattice
	   OUTCAR.i: extracted forces and positions in OUTCAR format

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

	   ==== OUTCAR ====
	   N                      # number of atoms
	   x1 y1 z1  fx1 fy1 fz1  # position and force
	   ...
	   # NOTE: only first atom can be displaced (placed at 0),
	   # and positions must not be rotated from unit cell definition
	   ==== OUTCAR ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read in our positions and forces, and convert the position
           to unit cell coord.  The forces become dynamical matrix elements
	   after inverting the displacements.

  Output:  Usual lattice function output
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"      // cell info
#include "Dij.H"       // point_type

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************
// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <OUTCAR.1> <OUTCAR.2> <OUTCAR.3>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:     cell file describing our lattice\n\
  OUTCAR.i: extracted forces and positions in OUTCAR format";

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
==== OUTCAR ====\n\
N                      # number of atoms\n\
x1 y1 z1  fx1 fy1 fz1  # position and force\n\
...\n\
# NOTE: only first atom can be displaced (placed at 0),\n\
# and positions must not be rotated from unit cell definition\n\
==== OUTCAR ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      // Note: we don't print out the elastic const. info... mainly
      // because we just ignore all that stuff anyway.
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *outcar_name[3];
  
  // Let's pull off the args:
  cell_name = args[0];
  for (i=0; i<3; ++i)
    outcar_name[i] = args[i+1];
  
  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
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
  if (TESTING) {
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    printf("# atomic mass = %.3lf\n", atomic_mass);
  }
  //-- ==== cell ====

  //++ ==== OUTCAR.i ====
  int Np[3]; // just read in each as it is
  double** outcar[3];

  for (i=0; (i<3) && (!ERROR); ++i) {
    infile = myopenr(outcar_name[i]);
    if (infile == NULL) {
      fprintf(stderr, "Couldn't open %s for reading.\n", outcar_name[i]);
      exit(ERROR_NOFILE);
    }
    // Natoms:
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%d", &(Np[i]));
    if (Np[i] < 1) {
      ERROR = ERROR_BADFILE;
      break;
    }
    outcar[i] = new double*[Np[i]];
    for (n=0; n<Np[i]; ++n) {
      outcar[i][n] = new double[6];
      nextnoncomment(infile, dump, sizeof(dump));
      sscanf(dump, "%lf %lf %lf %lf %lf %lf",
	     outcar[i][n]  , outcar[i][n]+1,  outcar[i][n]+2, 
	     outcar[i][n]+3, outcar[i][n]+4,  outcar[i][n]+5);
    }
    myclose(infile);
  }
  //-- ==== OUTCAR.i ====

  // if error already, bail out
  if (ERROR) {
    fprintf(stderr, "Error in %s -- bad Natoms value\n", outcar_name[i]);
    exit(ERROR);
  }

  // Now do some sanity checks:
  for (i=0; i<3; ++i) 
    if (Np[i] != Np[(i+1)%3]) {
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Error: Natoms in %s (%d) does not match %s (%d)\n", 
	      outcar_name[i], Np[i], outcar_name[(i+1)%3], Np[(i+1)%3]);
    }
  if (ERROR) exit(ERROR);
  
  for (i=0; i<3; ++i) {
    j = (i+1)%3;
    for (n=1; n<Np[i]; ++n)
      if (! equal_vect(outcar[i][n], outcar[j][n]) ) {
	ERROR = ERROR_BADFILE;
	fprintf(stderr, "Error: atom pos in %s does not match %s \n", 
		outcar_name[i], outcar_name[j]);
	fprintf(stderr, "  %s: %15.8lf %15.8lf %15.8lf\n",
		outcar_name[i], 
		outcar[i][n][0], outcar[i][n][1], outcar[i][n][2]);
	fprintf(stderr, "  %s: %15.8lf %15.8lf %15.8lf\n",
		outcar_name[j], 
		outcar[j][n][0], outcar[j][n][1], outcar[j][n][2]);
      }
  }
  // if (ERROR) exit(ERROR);
  ERROR = 0;

  // ***************************** ANALYSIS **************************
  // 1. Make our points
  int N = Np[0];
  point_type Dij[N];

  double cart_inv[9];
  careful_inverse(cart, cart_inv);

  for (d=0; d<3; ++d) {
    Dij[0].Runit[d] = 0;
    Dij[0].Rcart[d] = 0;
  }
  for (n=1; n<N; ++n) {
    double temp[3];
    mult_vect(cart_inv, outcar[0][n], temp);
    for (d=0; d<3; ++d) Dij[n].Runit[d] = lround(temp[d]);
    mult_vect(cart, Dij[n].Runit, Dij[n].Rcart);
  }
  
  if (TESTING) {
    printf("## N= %d\n", N);
    for (n=0; n<N; ++n) {
      printf("## %4d: %3d%3d%3d  %15.8lf %15.8lf %15.8lf\n", n+1,
	     Dij[n].Runit[0], Dij[n].Runit[1], Dij[n].Runit[2],
	     Dij[n].Rcart[0], Dij[n].Rcart[1], Dij[n].Rcart[2]);
    }
  }

  // 2.a. Calculate our delta matrix:
  double delta[9]; // matrix of displacements:
  for (i=0; i<3; ++i)
    for (d=0; d<3; ++d)
      delta[3*d+i] = outcar[i][0][d];
  double delta_inv[9], error;
  
  careful_inverse(delta, delta_inv, error);
  if (error == 1) {
    fprintf(stderr, "ERROR!  Your matrix of displacements is...\n");
    for (i=0; i<3; ++i)
      fprintf(stderr, "delta%d = %22.8le %22.8le %22.8le\n", i,
	      delta[i], delta[i+3], delta[i+6]);
    fprintf(stderr, "...which is singular.\n");
    ERROR = ERROR_BADFILE;
  }

  if (TESTING) {
    printf("## displacement matrix:\n");
    for (i=0; i<3; ++i)
      printf("## delta%d = %22.8le %22.8le %22.8le\n", i,
	     delta[i], delta[i+3], delta[i+6]);
  }

  // 2.b Construct force matrix for each, and multiply
  for (n=0; (n<N) && (!ERROR); ++n) {
    double force[9];
    for (i=0; i<3; ++i)
      for (d=0; d<3; ++d)
	force[3*d+i] = outcar[i][n][d+3];
    mult(force, delta_inv, Dij[n].mat);
    // enforce symmetrization and negate:
    for (d=0; d<9; ++d)
      Dij[n].mat[d] = -0.5*(Dij[n].mat[d] + Dij[n].mat[index(d%3,d/3)]);
  }

  // ****************************** OUTPUT ***************************
  if (!ERROR) {
    write_Dij(stdout, N, Dij, NULL);
  }
  
  // ************************* GARBAGE COLLECTION ********************
  for (i=0; i<3; ++i) {
    for (n=0; n<N; ++n)
      delete[] outcar[i][n];
    delete[] outcar[i];
  }
  
  free_cell(Cmn_list, u_atoms, 0);

  return ERROR;
}
