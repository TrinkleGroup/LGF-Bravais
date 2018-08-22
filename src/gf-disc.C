/*
  Program: gf_disc.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. elastic GF in k-space spherical harmonic components
	   3. lattice GF evaluated at a set of points

	   and we compute
	   1. the cutoff kmax -- maximum k value in our sphere
	   2. the discrete correction for the GF at each of the points
	      --though we don't output points with values below our
	        relative error tolerance

	   CHANGED!! We use a new cubic envelope function.  Let's hope
	   this works "better" than the old method... *shrug*
	   We've also removed the relative error tolerance BS

  Param.:  <cell> <Gk(lm)> <GL(R)>
           cell:    cell file describing our lattice
           Gk(lm):  spherical harmonic components of elastic GF in k-space
	   GL(R):   lattice GF evaluated at series of lattice sites

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

	   ==== Gk(lm) ====
	   Lmax d       # Maximum L value, dimensionality of vector
	   l1 m1 R(...) # l,m pair -- real components
	   l1 m1 I(...) # l,m pair -- imag components
	   l2 m2 R(...) 
	   l2 m2 I(...)
	   ...
	   # Read until EOF or l=-1.
	   # NOTE: only EVEN l values should appear!!
	   ==== Gk(lm) ====

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
#include "Dij.H"
#include "pointgroup.H"
#include "shell.H"

// If we are fed our GF in units of A^3/eV, then after we diagonalize
// invert and square root, we multiply by 1/(2Pi)*(ec*Na/10)^(1/2)
// to get the linear frequency nu in THz, once we divide by the
// square root of the mass in amu.  Whew.
const double THz_scale = 15.6333;


//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
double max_k_BZ (double cart_b[9]);

// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <Gk(lm)> <GL(R)>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Gk(lm):  spherical harmonic components of elastic GF in k-space\n\
  GL(R):   lattice GF evaluated at series of lattice sites";

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
  char *cell_name, *Ylm_name, *GL_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  GL_name = args[2];

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

  ERROR = read_Ylm_even (infile, Lmax, Ndim, RYlm, IYlm);

  myclose(infile);
  //-- ==== G(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm, IYlm, "##");
  }


  //++ ==== GL(R) ====
  int Np;
  point_type *GL;
  
  infile = myopenr(GL_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", GL_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij (infile, cart, Np, GL);

  myclose(infile);
  //-- ==== GL(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }


  // ***************************** ANALYSIS **************************
  // Let's get our point group information:
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  // Now, we sort GL and shell it; we use this to reduce the number of
  // GL's we have to calculate.
  shell_type* sh;
  int Nsh;
  
  sort_points(GL, Np);
  ERROR = gen_shell_list_open(GL, Np, gunit, Nop, sh, Nsh);

  // We first need to set kmax; so that when we
  // do the integration using a sphere inscribed inside our BZ, the
  // result is exact.  We no longer use an Rcut value
  double kmax;
  
  kmax = max_k_BZ(cart_rlv);
  if (TESTING) {
    printf("## ++ Analytic spherical piece integration.\n");
    printf("## kmax = %.5lf ; volume fraction = %.5lf\n", 
	   kmax, 4*M_PI*kmax*kmax*kmax/(3.*det(cart_rlv)));
  }

  // Some workspace:
  // Spherical harmonics
  SH_workspace_type * w = allocate_SH_workspace(Lmax);
  double gint[Ndim]; // Our semicont GF evaluated at a point!  
  double ***RYlm_R, ***IYlm_R;
  double Fl_val[Lmax+1];

  ERROR = alloc_Ylm(Lmax, Ndim, RYlm_R, IYlm_R);
  if (ERROR) {
    fprintf(stderr, "Error allocating our Ylm expansion...\n");
    exit(ERROR);
  }

  for (int nsh=0; nsh<Nsh; ++nsh) {
    shell_type *tsh = sh + nsh;
    
    // 1. first point in the shell:
    point_type *tGL = GL + tsh->elem[0];
    if (TESTING) {
      printf("## Evaluation point %d : %3d %3d %3d  %.5lf %.5lf %.5lf\n", 
	     np+1, 
	     tGL->Runit[0], tGL->Runit[1], tGL->Runit[2], 
	     tGL->Rcart[0], tGL->Rcart[1], tGL->Rcart[2]);
    }

    if (tGL->Rmagn > 0) {
      // Now, we need to do the numerical integration:
      
      // Using our *new* envelope function:
      Fn_int (Lmax, kmax*tGL->Rmagn, Fl_val, TESTING, 1024);
      
      if (TESTING) {
	printf("## Fl_val:\n");
	for (l=0; l<=Lmax; ++l) {
	  printf("## %d %.12le\n", l*2, Fl_val[l]);
	}
      }
      
      // Now, we'll evaluate our spherical harmonic expansion for our given
      // Rvect; the trick will be to create a "new" set of RYlm/IYlm that
      // have been multiplied by (-1)^l F_l / (4*Pi r); then we'll have
      // the final result.
      
      double g_scale = 1./(4*M_PI*tGL->Rmagn), l_scale;
      int l_sign;
      for (l=0, l_sign = 1; l<=Lmax; ++l, l_sign = -l_sign) {
	l_scale = g_scale * l_sign * Fl_val[l];
	for (m=0; m<=(4*l); ++m) {
	  for (d=0; d<Ndim; ++d) {
	    // Initialize
	    RYlm_R[l][m][d] = RYlm[l][m][d] * l_scale;
	    IYlm_R[l][m][d] = IYlm[l][m][d] * l_scale;
	  }
	}
      }
      if (TESTING) {
	// Print out *all* the elements we just calculated
	printf("## Ylm_R components:\n");
	for (l=0; l<=Lmax; ++l) {
	  for (m=0; m<=(4*l); ++m) {
	    // Real part
	    printf("## %3d %3d R ", 2*l, m-2*l);
	    print_mat(RYlm_R[l][m]); printf("\n");
	    // Imag part
	    printf("## %3d %3d I ", 2*l, m-2*l);
	    print_mat(IYlm_R[l][m]); printf("\n");
	  }
	}
      }
      
      // Get the continuum piece:
      eval_Ylm_expansion_R (Lmax, tGL->Rcart, Ndim, RYlm_R, IYlm_R, gint, w);
    }
    else {
      // For R = 0, the continuum piece is really simple:
      double g0_scale = M_2_SQRTPI/(8*M_PI*M_PI); // strange (exact) number
      g0_scale *= kmax * kpoint_envelope_int();  // and our maximum k value...
      for (d=0; d<Ndim; ++d)
	gint[d] = g0_scale * RYlm[0][0][d];
    }

    if (TESTING) {
      printf("## semicont. GF: ");
      print_mat(gint);
      printf("\n");
    }

    // Now evaluate the discrete correction term in place:
    for (d=0; d<Ndim; ++d)
      tGL->mat[d] += -gint[d];

    // 2. rotate to the rest of the points in the shell:
    for (j=1; j<tsh->Nelem; ++j) {
      np = tsh->elem[j];
      // find the group operation that takes R0 to our R:
      for (n=0; tsh->gR0[n] != np; ++n) ;
      double temp[9];
      mult(gcart[n], tGL->mat, temp);
      mult(temp, gcart[inv_index[n]], GL[np].mat);
    }

  }

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.
  if (VERBOSE) {
    // Human readable:
    for (np=0; np<Np; ++np) {
      point_type *tp = GL+np;
      printf("# |R| = %8.5lf  Rhat = %8.5lf %8.5lf %8.5lf\n", 
	     tp->Rmagn, tp->Rcart[0], tp->Rcart[1], tp->Rcart[2]);
      printf("# "); print_mat(tp->mat); printf("\n");
    }
  }

  printf("%d %.12le  # Number of points, kmax\n", Np, kmax);
  for (np=0; np<Np; ++np) {
    point_type *tp = GL+np;
    printf("%3d %3d %3d", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    for (d=0; d<Ndim; ++d)
      printf(" %.12le", tp->mat[d]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  free_SH_workspace(w);

  delete[] GL;

  free_Ylm(Lmax, RYlm, IYlm);
  free_Ylm(Lmax, RYlm_R, IYlm_R);
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
double max_k_BZ (double cart_b[9]) 
{
  int n[3];
  const int nmax = 3; // How far out should we search for nn?
  double bb[9];
  double kmin, kmagn;
  
  square(cart_b, bb);
  kmin = bb[0]; // simple starting guess.

  for (n[0] = -nmax; n[0] <= nmax; ++(n[0]))
    for (n[1] = -nmax; n[1] <= nmax; ++(n[1]))
      for (n[2] = -nmax; n[2] <= nmax; ++(n[2])) {
	if ( (n[0] == 0) && (n[1] == 0) && (n[2] == 0))
	  continue;
	kmagn = magnsq(bb, n);
	if (kmagn < kmin) kmin = kmagn;
      }
  return 0.5*sqrt(kmin);
}
