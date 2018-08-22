/*
  Program: phononDij.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a dynamical matrix, and calculate the phonon spectra
           on a given grid.

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
	   k1.1 k1.2 k1.3
	   ...
	   kN.1 kN.2 kN.3
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

// If we are fed our GF in units of A^3/eV, then after we diagonalize
// invert and square root, we multiply by 1/(2Pi)*(ec*Na/10)^(1/2)
// to get the linear frequency nu in THz, once we divide by the
// square root of the mass in amu.  Whew.
const double THz_scale = 15.6333023;

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

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'k', 'n', 'x'}; // flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  Dij(R): dynamical matrix\n\
  kpts:    list of kpoints to calculate our phonons\n\
  ** Output options: (default: -n)\n\
  -k: output k-point in Cartesian coord\n\
  -n: output number of k-point\n\
  -x: xmgrace friendly output (overrides -k and -n)";

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
k1.1 k1.2 k1.3\n\
...\n\
kN.1 kN.2 kN.3\n\
==== kpts ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used
  int KPOINT_CART, KPOINT_NUM; // How do we want to output?
  int XMGRACE;

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

  // Flags
  KPOINT_CART = flagon[0];
  KPOINT_NUM = flagon[1];
  XMGRACE = flagon[2];
  if ( (! KPOINT_CART) && (! KPOINT_NUM) )
    KPOINT_NUM = 1;

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

  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cart_rlv[9];
  double cell_vol;
  cell_vol = inverse(cart, cart_rlv);        // invert
  self_transpose(cart_rlv);                  // transpose in place
  mult(cart_rlv, 2*M_PI/cell_vol, cart_rlv); // scale

  //++ ==== Dij(R) ====
  int Np, Ndim = 9;
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
  double k0;  // Our scale
  char type_ident;  // Character specifying the kpoint type
  int LATT;   // = 1 if we've got RLV, = 0 if cartesian
  
  double** kpt; // Our points -- note: has dimensionality of 4 !!
  
  infile = myopenr(kpts_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", kpts_name);
    exit(ERROR_NOFILE);
  }
  double *w = NULL;
  ERROR = read_kpts(infile, cart_rlv, Nkpt, kpt, w, 
		    READ_KPTS_MAGN | READ_KPTS_CART | READ_KPTS_NOWEIGHTS);
  myclose(infile);
  if (ERROR) {exit(ERROR);}
  
  //-- ==== kpts ====

  if (TESTING) {
    printf("## %d kpoints.\n", Nkpt);
    printf("##  i    k_x      k_y      k_z       |k|\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## %2d %8.5lf %8.5lf %8.5lf  %8.5lf\n",
             k+1, kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
    }
  }
  

  // ***************************** ANALYSIS **************************
  // Compute the frequencies:
  // double nu[Nkpt][3];
  double** nu = new double*[Nkpt];
  for (int n=0; n<Nkpt; ++n) nu[n] = new double[3];

  if (TESTING) {
    printf("## **** Frequency calculation ****\n");
  }
  for (k=0; k<Nkpt; ++k) {
    double Dk[9]; // our FT of Dij
    fourier(Np, p, kpt[k], Dk, KPT_CART);
    
    double lambda[3];
    eigen(Dk, lambda); // Get the eigenvalues.
    if (TESTING) {
      printf("## Dij~(k) =");
      print_mat(Dk);
      printf("\n");
      printf("## lambda_%d = %.5lf %.5lf %.5lf\n", k+1, lambda[0],lambda[1],lambda[2]);
    }
    sort3(lambda);
    if (zero(kpt[k][3])) {
      if (TESTING){
	printf("## --setting to 0.\n");
      }
      for (d=0; d<3; ++d) lambda[d] = 0.;
    }
    for (d=0; d<3; ++d)
      if (lambda[d] >= 0.)
	nu[k][d] =  THz_scale*sqrt(lambda[d]/atomic_mass);
      else
	nu[k][d] = -THz_scale*sqrt(-lambda[d]/atomic_mass);
  }
  

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.

  if (XMGRACE) {
    for (i=0; i<3; ++i) {
      for (k=0; k<Nkpt; ++k)
        printf("%3d %.12le\n", k, nu[k][i]);
      if (i!=2)
        printf("&\n");
    }
  }
  else {
    printf("%d # Nkpt\n", Nkpt);
    for (k=0; k<Nkpt; ++k) {
      if (KPOINT_NUM)
        printf("%3d ", k);
      if (KPOINT_CART)
        printf("%15.12lf %15.12lf %15.12lf ", 
               kpt[k][0], kpt[k][1], kpt[k][2]);
      // Output our "frequencies":
      printf(" %15.12lf %15.12lf %15.12lf\n", 
             nu[k][0], nu[k][1], nu[k][2]);
    }
  }


  // ************************* GARBAGE COLLECTION ********************
  for (int n=0; n<Nkpt; ++n) delete nu[n];
  delete[] nu;
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
