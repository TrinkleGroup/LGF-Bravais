/*
  Program: phonon-eigen-Dij.C
  Author:  D. Trinkle
  Date:    2010 June 1
  Purpose: Read in a dynamical matrix, and calculate frequencies, average
           displacements, and normalized eigenvectors.  Can then be used
	   to sample displacements via quantum thermodynamics.

  Param.:  <cell> <Dij(R)> <kpt>
           cell:   cell file describing our lattice
	   Dij(R): dynamical matrix
	   kpt:    list of kpoints to calculate our phonons

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

	   ==== kpts ====
	   N k0 <type>     # Number of points, scale factor, C/L for units
	   k1.1 k1.2 k1.3 w1
	   ...
	   kN.1 kN.2 kN.3 wN
	   ==== kpts ====

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
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "kpts.H"
#include "eigen.H"
#include "cell.H"
#include "Dij.H"
#include "fourier.H"

  // Our force constants are in eV/A^2.  We want the energy in eV, and the
  // average displacement in A, for masses in amu.
  // eV_scale = hbar * sqrt(1eV/(A^2 amu)) in eV
  // Ang_scale = hbar/sqrt(eV amu) in A
const double eV_scale = 0.0646541378;
const double Ang_scale = 0.0646541378;


//***************************** SUBROUTINES ****************************

inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <Dij(R)> <kpt>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'s'}; // flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  Dij(R): dynamical matrix\n\
  kpts:   list of kpoints to calculate our phonons\n\
\n\
  -s      use 2pi as scale factor";

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
==== Dij(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij\n\
...\n\
nN.1 nN.2 nN.3  Dxx .. Dzz\n\
==== Dij(R) ====\n\
\n\
==== kpts ====\n\
N k0 <type>     # Number of points, scale factor, C/L for units\n\
k1.1 k1.2 k1.3 w1\n\
...\n\
kN.1 kN.2 kN.3 wN\n\
==== kpts ====";

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used

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
  int SCALE = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Dij_name, *kpts_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Dij_name = args[1];
  kpts_name = args[2];

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
  double a0;
  ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms, &a0);
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

  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cart_rlv[9];
  double cell_vol;
  cell_vol = inverse(cart, cart_rlv);        // invert
  self_transpose(cart_rlv);                  // transpose in place
  mult(cart_rlv, 2*M_PI/cell_vol, cart_rlv); // scale

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


  //++ ==== kpts ====

  int Nkpt;   // Total number of grid points
  
  double** kpt; // Our points -- note: has dimensionality of 4 !!
  
  infile = myopenr(kpts_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", kpts_name);
    exit(ERROR_NOFILE);
  }
  double *w = NULL;
  ERROR = read_kpts(infile, cart_rlv, Nkpt, kpt, w, 
		    READ_KPTS_MAGN | READ_KPTS_CART | READ_KPTS_WEIGHTS);
  myclose(infile);
  if (ERROR) {exit(ERROR);}
  
  //-- ==== kpts ====

  if (TESTING) {
    printf("## %d kpoints.\n", Nkpt);
    printf("##  i    k_x      k_y      k_z       |k|\n");
    for (int k=0; k<Nkpt; ++k) {
      printf("## %2d %8.5lf %8.5lf %8.5lf  %8.5lf\n",
             k+1, kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
    }
  }
  

  // *********************** ANALYSIS & OUTPUT ***********************
  double kscale = 2.*M_PI/a0;
  double kscale_1 = 1./kscale;
  if (SCALE) kscale = 2*M_PI; // for compatibility with random displacement code
  printf("%d %.12lf C\t#Nkpts #scale(2pi/a) #cartesian\n", Nkpt, kscale);
  printf("#kx(A^-1)\t#ky(A^-1)\t#kz(A^-1)\t#ene(eV)\t#x0(A)\t#ux\t#uy\t#uz\t#w\n");
  // Compute the frequencies:
  for (int k=0; k<Nkpt; ++k) {
    double Dk[9]; // our FT of Dij
    fourier(Np, p, kpt[k], Dk, KPT_CART);
    
    if (TESTING) {
      printf("## Dij~(k) =");
      print_mat(Dk);
      printf("\n");
    }
    double lambda[3];
    eigen_vect(Dk, lambda); // Get the eigenvalues and vectors
    if (TESTING) {
      printf("## lambda_%d = %.5lf %.5lf %.5lf\n", k+1, lambda[0],lambda[1],lambda[2]);
      // HORRIBLE ERROR: incorrect indexing.  The nth vector is D[n+0],D[n+3],D[n+6]
      for (int n=0; n<3; ++n)
	printf("## vect = %.5lf %.5lf %.5lf\n", Dk[n+0], Dk[n+3], Dk[n+6]);
    }
    if (zero(kpt[k][3])) {
      if (TESTING){
	printf("## --setting to 0.\n");
      }
      for (int n=0; n<3; ++n) lambda[n] = 0.;
    }

    for (int n=0; n<3; ++n) {
      printf("%15.12lf %15.12lf %15.12lf ", 
	     kpt[k][0] * kscale_1, kpt[k][1] * kscale_1, kpt[k][2] * kscale_1);
      double omega=lambda[n]/atomic_mass;
      if (omega>=0) omega = eV_scale*sqrt(omega);
      else          omega = -eV_scale*sqrt(-omega);
      double x0;
      if (zero(omega)) x0=0;
      else             x0=Ang_scale/sqrt(fabs(omega)*atomic_mass);
      printf("%15.12lf %15.12lf ", omega, x0);
      // HORRIBLE ERROR: incorrect indexing.  The nth vector is D[n+0],D[n+3],D[n+6]
      for (int d=0; d<3; ++d)
	printf("%15.12lf ", Dk[n+3*d]);
      printf("%15.12le\n", w[k]);
    }
  }
  
  // ************************* GARBAGE COLLECTION ********************
  free_kpts(Nkpt, kpt);
  delete[] w;
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
