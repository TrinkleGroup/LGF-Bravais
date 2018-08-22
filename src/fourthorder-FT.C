/*
  Program: fourthorder-FT.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a dynamical matrix, and calculate:

                       1
	   [ab,cd] = - - SUM D   x x
                       2  x   ab  c d

                          1
	   [ab,cdef] = - -- SUM D   x x x x
                         24  x   ab  c d e f

	   From these, we construct the matrices:

	   lambda (k) = SUM [ab,cd]k k
	                 cd         c d

                 (2)
	   lambda  (k) =  SUM [ab,cdef]k k k k 
	                 cdef           c d e f

	   where k is normalized.  We output the matrix:
                      -1       (2)               -1
	   [lambda(k)]  [lambda  (k)] [lambda(k)]

	   for each k point in our mesh.  This has units of A^3/eV
	   (i.e., same as our GF).

  Param.:  <cell> <Ylm> <Dij(R)> <grid>
           cell:   cell file describing our lattice
	   Ylm:    spherical harmonic expansion of (lambda(k))^-1
	   Dij(R): dynamical matrix
	   grid:   "grid" of direction cosines for computing L^-1 L(2) L^-1

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

           ==== Dij(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij
           ...
           nN.1 nN.2 nN.3  Dxx .. Dzz
           ==== Dij(R) ====

	   ==== grid ====
	   N         # Number of points
	   l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)
	   ...
	   lN mN nN
	   ==== grid ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and calculate.

  Output:  We output the points in unit coord.  (Switch for lattice coord?)
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_eigen.h>  // for determining stability
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "voigt.H"
#include "Dij.H"
#include "sphere-harm.H"
#include "lambda.H"

//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <Ylm> <Dij(R)> <grid>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'y'}; // flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  Gk(lm): spherical harmonic components of elastic GF in k-space\n\
  Dij(R): dynamical matrix\n\
  grid:   \"grid\" of direction cosines for computing L^-1 L(2) L^-1\n\
  -y: use Ylm expansion for GF instead of calcing on the fly";

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
==== G(lm) ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...) \n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
# NOTE: only EVEN l values should appear!!\n\
==== G(lm) ====\n\
\n\
==== Dij(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij\n\
...\n\
nN.1 nN.2 nN.3  Dxx .. Dzz\n\
==== Dij(R) ====\n\
\n\
==== grid ====\n\
N         # Number of points\n\
l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)\n\
...\n\
lN mN nN\n\
==== grid ====";

int main ( int argc, char **argv ) 
{
  int i, j, k, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used
  int YLM = 0;

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
  
  // flags:
  YLM = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *Dij_name, *grid_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  Dij_name = args[2];
  grid_name = args[3];

  double cart[9], cart2[9];
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
  square(cart, cart2);
  //-- ==== cell ====

  //++ ==== G(lm) ====
  // NOTE: this is different than other versions, since we're *forced*
  // to only have EVEN l values.  So [l] is really 2*l.
  int Lmax;  // Maximum l value / 2
  int Ndim;  // Dimensionality -- must be 9
  double ***RYlm, ***IYlm; // Separate our real and imag. components

  infile = myopenr(Ylm_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Ylm_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Ylm_even(infile, Lmax, Ndim, RYlm, IYlm);
  myclose(infile);

  if ( ERROR || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
            Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== G(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm, IYlm, "##");
  }

  //++ ==== Dij(R) ====
  int Np;
  point_type *p; // set of points

  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, p, READ_DIJ_MAT);
  if (ERROR) {delete p; exit(ERROR);}
  //-- ==== Dij(R) ====


  //++ ==== grid ====
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }

  int Ngrid;
  double** lmn_list;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Ngrid);
  if (Ngrid <= 0) {
    fprintf(stderr, "No points listed in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  else lmn_list = new double* [Ngrid];
  for (n=0; (!ERROR) && (n<Ngrid) && (!feof(infile)); ++n) {
    double magn;
    lmn_list[n] = new double[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", 
           &(lmn_list[n][0]), &(lmn_list[n][1]), &(lmn_list[n][2]) );
    // Normalize.
    magn = sqrt(dot(lmn_list[n],lmn_list[n]));
    if (dcomp(magn, 0.)) {
      fprintf(stderr, 
              "%s at line %d: %.5le %.5le %.5le with magn = %.5le\n",
              grid_name, i+1,
              lmn_list[0][0], lmn_list[0][1], lmn_list[0][2], magn);
      ERROR = ERROR_BADFILE;
    }
    else {
      magn = 1./magn;
      for (i=0; i<3; ++i) lmn_list[n][i] *= magn;
    }
  }
  // Make sure we read enough points...
  if (n != Ngrid) ERROR = ERROR_BADFILE;
  myclose(infile);
  //-- ==== grid ====


  // ***************************** ANALYSIS **************************
  // Use Voigt notation to generate the symmetrized elements:
  int ab, cd, ef, ii;
  int a, b, c, d, e, f;
  double lambda[9][9], lambda2[9][9][9];

  make_lambda(Np, p, lambda);
  make_lambda2(Np, p, lambda2);
  
  if (TESTING) {
    printf("## lambda: [ab,cd]: --non-zero only\n");
    print_lambda(lambda, "## ");
    printf("## lambda(2): [ab,cdef]: --non-zero only\n");
    print_lambda2(lambda2, "## ");
  }
  
  // Now, to evaluate our matrices:
  printf("%d 9 # number of points\n", Ngrid);
  double lambdak[9], lambda2k[9], inv_lambdak[9];
  double temp[9], result[9];
  double scale;

  // Spherical harmonics
  SH_workspace_type * work = allocate_SH_workspace(Lmax);
  
  for (n=0; n<Ngrid; ++n) {
    double *lmn = lmn_list[n];
    if (YLM) {
      // New method: calculate using Ylm expansion:
      eval_Ylm_expansion_R(Lmax, lmn, 9, RYlm, IYlm, inv_lambdak, work);
    }
    else {
      lambda_k(lmn, lambda, lambdak);
      // invert
      scale = inverse(lambdak, inv_lambdak);
      mult(inv_lambdak, 1./scale, inv_lambdak);
    }
    
    lambda2_k(lmn, lambda2, lambda2k);
    // multiply
    mult(inv_lambdak, lambda2k, temp);
    mult(temp, inv_lambdak, result);
    // print
    printf("%15.12lf %15.12lf %15.12lf", lmn[0], lmn[1], lmn[2]);
    for (i=0; i<9; ++i) printf(" %.15le", result[i]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  free_SH_workspace(work);
  free_Ylm(Lmax, RYlm, IYlm);
  for (n=0; n<Ngrid; ++n)
    delete[] lmn_list[n];
  delete[] lmn_list;
  
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
