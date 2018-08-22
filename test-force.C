/*
  Program: test-force.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. lattice GF evaluated at a set of points

	   And applies a test-force, and outputs the displaced set of
	   atom positions.  Since our LGF is defined relative to an
	   atom at the origin, our test force is effectively applied
	   at the origin.

  Param.:  <cell> <GL(R)> <fx> <fy> <fz>
           cell:    cell file describing our lattice
	   GL(R):   lattice GF evaluated at series of lattice sites
	   fxyz:    test force, Cartesian coord.

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

	   ==== GL(R) ====
	   N                           # number of points
	   n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG
	   ...
	   nN.1 nN.2 nN.3  Gxx .. Gzz
	   ==== GL(R) ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Currently just testing the gridding of k-points

  Output: 
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?


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
const int NUMARGS = 5;
const char* ARGLIST = "<cell> <GL(R)> <fx> <fy> <fz>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  GL(R):   lattice GF evaluated at series of lattice sites\n\
  fxyz:    Cartesian components of test force";

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
==== GL(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GF\n\
...\n\
nN.1 nN.2 nN.3  Gxx .. Gzz\n\
==== GL(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n, np; // General counting variables.

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *GL_name;
  double f0[3];
  
  // Let's pull off the args:
  cell_name = args[0];
  GL_name = args[1];
  sscanf(args[2], "%lf", &(f0[0]));
  sscanf(args[3], "%lf", &(f0[1]));
  sscanf(args[4], "%lf", &(f0[2]));
  
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

  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cart_rlv[9];
  double cell_vol;
  cell_vol = inverse(cart, cart_rlv);        // invert
  self_transpose(cart_rlv);                  // transpose in place
  mult(cart_rlv, 2*M_PI/cell_vol, cart_rlv); // scale


  //++ ==== GL(R) ====
  int Np;
  int **Runit;
  double **GL;
  int Ndim = 9;
  
  infile = myopenr(GL_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", GL_name);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Np);
  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in %s\n", Np, GL_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    Runit = new int *[Np];
    GL = new double *[Np];
  }
  
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    Runit[n] = new int[3];
    GL[n] = new double[Ndim];
    sscanf(dump, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	   &(Runit[n][0]), &(Runit[n][1]), &(Runit[n][2]),
	   &(GL[n][0]), &(GL[n][1]), &(GL[n][2]),
	   &(GL[n][3]), &(GL[n][4]), &(GL[n][5]),
	   &(GL[n][6]), &(GL[n][7]), &(GL[n][8]));
  }
  myclose(infile);
  if (n!=Np) {
    fprintf(stderr, "Didn't read enough points in %s\n", GL_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== GL(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }


  // ***************************** ANALYSIS **************************
  double Rcart[Np][3]; // Cartesian components of position
  double udisp[Np][3]; // displacement

  for (np=0; np<Np; ++np) {
    mult_vect(cart, Runit[np], Rcart[np]);
    mult_vect(GL[np], f0, udisp[np]);
  }
  

  // ****************************** OUTPUT ***************************

  printf("%d  # Number of points\n", Np);
  for (np=0; np<Np; ++np) {
    printf("%19.12le %19.12le %19.12le", 
	   Rcart[np][0]+udisp[np][0], 
	   Rcart[np][1]+udisp[np][1], 
	   Rcart[np][2]+udisp[np][2]);
    if (VERBOSE) {
      printf("  # %3d %3d %3d", Runit[np][0], Runit[np][1], Runit[np][2]);
      printf(" u= %.4le", sqrt(dot(udisp[np], udisp[np])));
    }
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  for (np=0; np<Np; ++np) {
    delete[] Runit[np];
    delete[] GL[np];
  }
  delete[] Runit;
  delete[] GL;
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
