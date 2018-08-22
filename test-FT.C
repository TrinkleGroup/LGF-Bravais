/*
  Program: test-FT.C
  Author:  D. Trinkle
  Date:    April 26, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. dynamical matrix evaluated at a set of points
	   3. kpt mesh

	   and we compute the fourier transform using two different
	   methods (which should be equivalent).

  Param.:  <cell> <Dij(R)> <kpts>
           cell:    cell file describing our lattice
	   Dij(R):  dynamical matrix evaluated at a series of points
	   kpts:    symmetrized kpt mesh (with weights)

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
           N k0 <type>       # Number of points, scale factor, C/L for units
           k1.1 k1.2 k1.3 w1 # kpoint and weight
           ...
           kN.1 kN.2 kN.3 wN
           ==== kpts ====
	   

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Currently just testing the gridding of k-points

  Output:  The discretization correction of the lattice GF.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include <gsl/gsl_sf_bessel.h>   // Spherical bessel functions
#include <gsl/gsl_integration.h> // Needed for doing integration
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "cont-int.H" // For doing the integration over the continuous piece
#include "sphere-harm.H" // Spherical harmonic evaluation
#include "pointgroup.H"
#include "Dij.H"
#include "kpts.H"
#include "shell.H"
#include "fourier.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  //  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
  printf("xx= %15.8le yy= %15.8le zz= %15.8le  xy= %15.8le yz= %15.8le zx= %15.8le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <Dij(R)> <kpts>";

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'f', 'i', 'l'}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Dij(R):  dynamical matrix evaluated at a series of points\n\
  kpts:    symmetrized kpt mesh (with weights)\n\
  -f:      use closed shells on fourier transform\n\
  -i:      use closed shells on inverse fourier transform\n\
  -l:      use a lookup table for cos(k.R); assumes uniform mesh";

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
  int d, i, j, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter
  int FOURIER = 0;
  int INVERSEFT = 0;
  int LOOKUP = 0;

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

  FOURIER = flagon[0];
  INVERSEFT = flagon[1];
  LOOKUP = flagon[2];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Dij_name, *kpt_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Dij_name = args[1];
  kpt_name = args[2];

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


  //++ ==== Dij(R) ====
  int Np;
  point_type* Dij;
  
  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Dij(infile, cart, Np, Dij, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== Dij(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", Dij_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (n=0; n<Np; ++n) {
      printf("## %3d %3d %3d ", 
	     Dij[n].Runit[0], Dij[n].Runit[1], Dij[n].Runit[2]);
      print_mat(Dij[n].mat);
      printf("\n");
    }
  }


  //++ ==== kpts ====
  int Nkpt;
  double** kpt, *w;
  
  infile = myopenr(kpt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", kpt_name);
    exit(ERROR_NOFILE);
  }
  // the last entry specifies what we read and return:
  ERROR = read_kpts(infile, cart_rlv, Nkpt, kpt, w, 
		    READ_KPTS_MAGN | READ_KPTS_CART | READ_KPTS_WEIGHTS);
  myclose(infile);
  //-- ==== kpts ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", kpt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## %d kpoints\n", Nkpt);
    for (i=0; i<Nkpt; ++i)
      printf("## %8.5lf %8.5lf %8.5lf %.5le\n", 
	     kpt[i][0], kpt[i][1], kpt[i][2], w[i]);
  }

  // ***************************** ANALYSIS **************************
  // Let's get our point group information:
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  double g_i_ginv[MAXop][Ni][Ni];
  // 1. point group operations:
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);
  // 2. g [i] ginv, where [i] is a voigt matrix
  for (n=0; n<Nop; ++n) {
    for (i=0; i<Ni; ++i) {
      gen_g_i_ginv(i, gcart[n], gcart[inv_index[n]], g_i_ginv[n][i]);
    }
  }

  // Now, we sort Dij and shell it; we use this to reduce the number of
  // Dij's we have to calculate.
  shell_type* sh = NULL;
  int Nsh;
  
  sort_points(Dij, Np);
  ERROR = gen_shell_list(Dij, Np, gunit, Nop, sh, Nsh);

  double **Dk; // our FT
  Dk = new double *[Nkpt];
  for (k=0; k<Nkpt; ++k) Dk[k] = new double[9];

  coskR_table_type coskR;
  if (LOOKUP) {
    ERROR = init_coskR_table(Nkpt, kpt, KPT_CART, cart, coskR);
    if (ERROR) {
      fprintf(stderr, "Problem creating lookup table... perhaps kpt mesh isn't uniform?\n");
      exit(ERROR);
    }
  }
  
  for (k=0; k<Nkpt; ++k) {
    double* tkpt = kpt[k]; // this kpoint
    // Do both, and compare...
    if (FOURIER) {
      if (LOOKUP)
	fourier_cs_table(Np, Dij, Nsh, sh, Nop, inv_index, g_i_ginv, 
			 tkpt, Dk[k], KPT_CART, coskR);
      else
	fourier_cs(Np, Dij, Nsh, sh, Nop, inv_index, g_i_ginv, 
		   tkpt, Dk[k], KPT_CART);
    }
    else {
      if (LOOKUP)
	fourier_table(Np, Dij, tkpt, Dk[k], KPT_CART, coskR);
      else
	fourier(Np, Dij, tkpt, Dk[k], KPT_CART);
    }
    if (TESTING) {
      printf("## kpt: %8.5lf %8.5lf %8.5lf\n", tkpt[0], tkpt[1], tkpt[2]);
      printf("##  Dk: "); print_mat(Dk[k]); printf("\n##\n");
    }
  }


  // Now all that remains is to *invert* the FT.  This takes advantage
  // of the symmetry that must have been built into our kpoint mesh
  // and the symmetry in Gd; we run through our shells one by one

  double DR[9]; // our two attempts at inverting the transformation.
  double errR[9];
  for (d=0; d<9; ++d) errR[d] = 0.;
  for (int nsh=0; nsh<Nsh; ++nsh) {
    shell_type *tsh = sh + nsh;
    // 1. first point in the shell:
    point_type *tp = Dij + tsh->elem[0];

    if (INVERSEFT) {
      if (LOOKUP) 
	inv_fourier_cs_table(Dij, tsh, DR, Nop, inv_index, g_i_ginv, Nkpt, Dk,
			     kpt, w, KPT_CART, coskR);
      else
	inv_fourier_cs(Dij, tsh, DR, Nop, inv_index, g_i_ginv, Nkpt, Dk,
		       kpt, w, KPT_CART);
    }
    else {
      if (LOOKUP)
	inv_fourier_table(Dij, tsh, DR, Nop, inv_index, gcart, Nkpt, Dk, 
			  kpt, w, KPT_CART, coskR);
      else
	inv_fourier(Dij, tsh, DR, Nop, inv_index, gcart, Nkpt, Dk, kpt, w,
		    KPT_CART);
    }
    
    if (TESTING) {
      printf("## point: %3d %3d %3d\n", 
	     tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      printf("## Dij: "); print_mat(tp->mat); printf("\n");
      printf("## DR:  "); print_mat(DR); printf("\n##\n");
    }
    for (d=0; d<9; ++d) errR[d] += fabs(tp->mat[d]-DR[d]);
  }

  // ****************************** OUTPUT ***************************

  printf("Total deviation R:   "); print_mat(errR); printf("\n");

  // ************************* GARBAGE COLLECTION ********************
  if (LOOKUP) free_coskR_table(coskR);
  free_kpts(Nkpt, kpt);

  for (n=0; n<Nkpt; ++n)
    delete[] Dk[n];
  delete[] Dk;

  free_shell(sh, Nsh);

  delete[] Dij;
    
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
