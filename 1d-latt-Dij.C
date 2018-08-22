/*
  Program: 1d-latt-Dij.C
  Author:  D. Trinkle
  Date:    2012 June 20
  Purpose: This takes in:
           1. unit cell definition
	   2. slab coordinate system
	   3. dynamical matrix evaluated at a set of points
	   4. a 1D kpoint mesh in our BZ
	   5. a set of points to evaluate the LGF

	   and we compute
	   1. the LGF for each point

	   ==== Note for 1d ====
	   This is the ONE-DIMENSIONAL version of 2d-latt-Dij, and so most
	   of the algorithm comes from there, with changes appropriate
	   for a 1-dimensional LGF.  Note: the only thing that makes
	   this "1d" is that we make our LGF periodic along two threading
	   directions t1 and t2.  This simplifies a lot of the calculation.
	   ==== Note for 2d ====

  Param.:  cell plane Dij(R) kpts latt
           cell:    cell file describing our lattice
	   plane:   planar coordinate system
	   Dij(R):  dynamical matrix evaluated at a series of points
	   kpts:    1d kpt mesh (with weights)
	   latt:    list of lattice points to eval LGF

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

           ==== plane ====
           t1.1 t1.2 t1.3  # planar vector
           t2.1 t2.2 t2.3  # planar vector
           n1.1 n1.2 n1.3  # mutual perp. vector
           # all three vectors are given in unit cell coord.
           ==== plane ====

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

           ==== latt ====
           N               # number of points
           n1.1 n1.2 n1.3  # unit cell coord
           ...
           nN.1 nN.2 nN.3
           ==== latt ====
	   

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Stolen (mostly) from 2d-latt-Dij and changed accordingly

  Output:  The full-fledged 1d LGF at each point.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>   // Spherical bessel functions
#include <gsl/gsl_integration.h> // Needed for doing integration
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "cont-int.H" // For doing the integration over the continuous piece
#include "Dij.H"
#include "shell.H" // We only need this for the "sort_points" routine
#include "kpts.H"
#include "eigen.H"
#include "fourier.H"
#include "plane.H"
#include "lambda.H" // for evaluation (on the fly) of lambda^-1 and lambda^(2)

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************
const double BETA_ENVELOPE = 5.;
const double BETA2_ENVELOPE = BETA_ENVELOPE*BETA_ENVELOPE;
const double INV_BETA_ENVELOPE = 1./(2*BETA_ENVELOPE);
const double INV_BETA2_ENVELOPE = INV_BETA_ENVELOPE*INV_BETA_ENVELOPE;
const double BETA_PI_ENVELOPE = BETA_ENVELOPE * 0.5* M_2_SQRTPI; // beta/sqrt(pi)

inline double kpoint_gaussian_envelope (const double& x) {
  return exp(-x*x*BETA2_ENVELOPE);
}

#define SMALL_U_TOL 1e-4
#define __QUAD_2  -0.5
#define __QUAD_3   0.16666666666666666667
#define __QUAD_4  -0.41666666666666666667e-1
// ##define __QUAD_5   0.83333333333333333333e-2
// ##define __QUAD_6  -0.13888888888888888889e-2

inline double kpoint_gaussian_envelope_1 (const double& x) {
  double u = x*x*BETA2_ENVELOPE;
  if (u>SMALL_U_TOL) return 1-exp(-u);
  else
    return u*(1 + u*(__QUAD_2
		     + u*(__QUAD_3
			  + u*__QUAD_4)));
}
#undef SMALL_U_TOL
#undef __QUAD_2
#undef __QUAD_3
#undef __QUAD_4

// The inverse FT, multiplied by kmax (need to divide that back out later)
inline double F1d_int(const double& kR, const int& TESTING) {
  return -0.5*kR*erf(kR*INV_BETA_ENVELOPE) - BETA_PI_ENVELOPE*exp(-kR*kR*INV_BETA2_ENVELOPE);
}

void smallk_expand(int Np, point_type* p, double kpt[3], double k_kmax, double pk[9], 
		   double lambda[9][9], const int TYPE, const coskR_table_type *table);

inline void smallk_expand(int Np, point_type* p, double kpt[3], double k_kmax, double pk[9], 
			  double lambda[9][9], const int TYPE, const coskR_table_type &table) 
{ return smallk_expand(Np, p, kpt, k_kmax, pk, lambda, TYPE, &table);}
  
inline void smallk_expand(int Np, point_type* p, double kpt[3], double k_kmax, double pk[9], 
			  double lambda[9][9], const int TYPE) 
{ return smallk_expand(Np, p, kpt, k_kmax, pk, lambda, TYPE, NULL);}


// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
double max_k_BZ (double cart_b[9]);

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
const int NUMARGS = 5;
const char* ARGLIST = "cell plane Dij(R) kpts latt";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'l', 'r'};
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  plane:    planar coordinate system\n\
  Dij(R):  dynamical matrix evaluated at a series of points\n\
  kpts:    symmetrized kpt mesh (with weights)\n\
  latt:    list of lattice points to eval Gdisc\n\
  -l: implement a lookup table for more efficient FT--highly recommended\n\
  -r: rotate tensor to t1,t2,n coordinates";

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
==== plane ====\n\
t1.1 t1.2 t1.3  # planar vector\n\
t2.1 t2.2 t2.3  # planar vector\n\
n1.1 n1.2 n1.3  # mutual perp. vector\n\
# all three vectors are given in unit cell coord.\n\
==== plane ====\n\
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
==== kpts ====\n\
\n\
==== latt ====\n\
N               # number of points\n\
n1.1 n1.2 n1.3  # unit cell coord\n\
...\n\
nN.1 nN.2 nN.3\n\
==== latt ====";

int main ( int argc, char **argv ) 
{
  //  int i, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (int i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

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
  int LOOKUP = flagon[0];  // use a lookup table for cos(k.R)--more efficient
  int ROTATE = flagon[1];  // rotate to t1,t2,n coordinates (useful?)

  if (TESTING) {
    if (LOOKUP) printf("## Using lookup tables for FT\n");
  }
  
  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char *cell_name = args[0];
  char *plane_name = args[1];
  char *Dij_name = args[2];
  char *kpt_name = args[3];
  char *latt_name = args[4];

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

  //++ ==== plane ====
  infile = myopenr(plane_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", plane_name);
    exit(ERROR_NOFILE);
  }

  int t1_unit[3], t2_unit[3], n_unit[3];
  ERROR = read_plane(infile, t1_unit, t2_unit, n_unit);
  myclose(infile);
  //-- ==== plane ====

  if (ERROR) {
    fprintf(stderr, "Problem with dislocation: %s\n", plane_name);
    exit(ERROR);
  }
  double t1_0[3], t2_0[3], n_0[3];
  double area;
  make_planecoord(cart, t1_0, t2_0, n_0, t1_unit, t2_unit, n_unit, 
		  DISLCOORD_UNIT | DISLCOORD_NOCHANGE);
  mult_vect(cart, t1_unit, t1_0);
  mult_vect(cart, t2_unit, t2_0);
  { double t3_0[3];
    crossprod(t1_0, t2_0, t3_0);
    area = sqrt(dot(t3_0, t3_0));
  }
  if (TESTING) {
    printf("## vect    input                 cartesian\n");
    printf("## t1 = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           t1_unit[0],t1_unit[1],t1_unit[2], t1_0[0],t1_0[1],t1_0[2]);
    printf("## t2 = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           t2_unit[0],t2_unit[1],t2_unit[2], t2_0[0],t2_0[1],t2_0[2]);
    printf("## n  = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           n_unit[0],n_unit[1],n_unit[2], n_0[0],n_0[1],n_0[2]);
  }

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
    for (int n=0; n<Np; ++n) {
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
  // Check to make sure the kpoint mesh is consistent with our threading
  // directions--i.e., dot(kpt, t1) = 2*Pi*n and dot(kpt, t2) = 2*Pi*n
  { int k;
    for (k=0; k<Nkpt; ++k) {
      double dottest = fabs(dot(kpt[k], t1_0));
      dottest = dottest - 2.0*M_PI*(round(dottest*0.5*M_1_PI));
      if (! zero(dottest)) break;
      dottest = fabs(dot(kpt[k], t2_0));
      dottest = dottest - 2.0*M_PI*(round(dottest*0.5*M_1_PI));
      if (! zero(dottest)) break;
    }
    if (k != Nkpt) {
      fprintf(stderr, "kpt-set from %s does not appear to be consistent with dislocation coordinate system\n", kpt_name);
      fprintf(stderr, "%d kpoint: %18.15lf %18.15lf %18.15lf  %18.15lf\n", k+1,
	      kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
      fprintf(stderr, "%d kpoint: dot(kpt,t1)= %8.5lf\n", k+1, dot(kpt[k],t1_0));
      fprintf(stderr, "%d kpoint: dot(kpt,t2)= %8.5lf\n", k+1, dot(kpt[k],t2_0));
      exit(ERROR_BADFILE);
    }
  }

  if (TESTING) {
    printf("## %d kpoints\n", Nkpt);
    for (int i=0; i<Nkpt; ++i)
      printf("## %8.5lf %8.5lf %8.5lf %.5le | t1.k= %8.5lf t2.k= %8.5lf n.k= %8.5lf\n", 
	     kpt[i][0], kpt[i][1], kpt[i][2], w[i],
	     dot(kpt[i],t1_0), dot(kpt[i],t2_0), dot(kpt[i],n_0));
  }

  //++ ==== latt ====
  // Finally, the set of lattice points for which we want to calculate
  // our discretization correction
  point_type *Gd; // our discretization correction; held in .mat
  int Nlatt; // number of lattice points

  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Nlatt, Gd, READ_DIJ_NOMAT);
  myclose(infile);
  //-- ==== latt ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", latt_name);
    exit(ERROR);
  }

  // We "redefine" what Rmagn means by keeping it to be the distance
  // perpendicular to the t1 x t2 plane (i.e, along n_0):
  for (int n=0; n<Nlatt; ++n)
    Gd[n].Rmagn = fabs(dot(Gd[n].Rcart, n_0));

  // For ease of use later, let's sort our points:
  sort_points(Gd, Nlatt);
  
  if (TESTING) {
    printf("## lattice points to calc Gd (SORTED):\n");
    printf("## NOTE: magnitude refers to distance from t1 x t2 plane, and coord. are t1,t2,n\n");
    for (int n=0; n<Nlatt; ++n)
      printf("## %3d %3d %3d  %11.6lf %11.6lf %11.6lf  %11.6lf\n", 
	     Gd[n].Runit[0], Gd[n].Runit[1], Gd[n].Runit[2],
	     dot(Gd[n].Rcart, t1_0), dot(Gd[n].Rcart, t2_0),
	     dot(Gd[n].Rcart, n_0), Gd[n].Rmagn);
  }

  // ***************************** ANALYSIS **************************
  // 0. Initial setup.
  // 0.a. lookup table for FT:
  coskR_table_type coskR;
  if (LOOKUP) {
    ERROR = init_coskR_table(Nkpt, kpt, KPT_CART, cart, coskR);
    if (ERROR) {
      fprintf(stderr, "Problem creating lookup table... perhaps kpt mesh isn't uniform?\n"
	      );
      exit(ERROR);
    }
  }

  // 0.b. Now we set kmax; so that when we do the integration using a
  // sphere inscribed inside our BZ, the result is exact.
  double kmax;
  
  kmax = max_k_BZ(cart_rlv);
  if (TESTING) {
    printf("## ++ Analytic spherical piece integration.\n");
    printf("## kmax = %.5lf ; volume fraction = %.5lf\n", 
	   kmax, 4*M_PI*kmax*kmax*kmax/(3.*det(cart_rlv)));
  }

  // Test to see if we've got enough k-points for our set of points:
  double mink=kmax, maxR = Gd[Nlatt-1].Rmagn;
  for (int k=0; k<Nkpt; ++k)
    if ( (mink > kpt[k][3]) && (! zero(kpt[k][3])) )
      mink = kpt[k][3];
  if (TESTING) {
    printf("## maxR = %.8le  mink = %.8le\n", maxR, mink);
    printf("## pi/maxR*mink = %.5lf\n", M_PI/(maxR*mink));
  }
  // Now, our criterion is that mink*maxR must be less than Pi/6:
  if ( (mink*maxR) > (M_PI/5.99) )
    fprintf(stderr, "WARNING: maxR = %.4le  mink = %.4le  pi/maxR*mink = %.5lf ; FT might not be accurate.\n", 
	    maxR, mink, M_PI/(maxR*mink));

  double krealsmall;
  krealsmall = kmax*0.1; // the limit below which we really get "clever"
  //  if (krealsmall > (kpoint_envelope_alpha*kmax))
  //    krealsmall = (kpoint_envelope_alpha*kmax);
  //  if (krealsmall > (kpoint_envelope_alpha*kmax_0))
  //    krealsmall = (kpoint_envelope_alpha*kmax_0);

  // For the planar LGF, you only need the elastic term (i.e., lambda_inv)
  // but we'll get lambda2 in order to evaluate the k=0 term.
  double lambda[9][9], lambda2[9][9][9], GE[9], Gdisc[9];
  make_lambda(Np, Dij, lambda);
  make_lambda2(Np, Dij, lambda2);
  eval_lambda_inv(lambda, n_0, GE);
  eval_lambda2(lambda, lambda2, n_0, Gdisc);

  if (TESTING) {
    //    printf("## ksmall =     %.5lf :below which we don\'t use Fourier series.\n", ksmall);
    printf("## krealsmall = %.5lf :below which we very carefully evaluated Gdk.\n", krealsmall);
    printf("## lambda: [ab,cd]: --non-zero only\n");
    print_lambda(lambda, "## ");
    printf("## elastic GF (not multiplied by 1/k^2)\n");
    print_mat(GE); printf("\n");
    printf("## lambda(2): [ab,cdef]: --non-zero only\n");
    print_lambda2(lambda2, "## ");
    printf("## discontinuity GF\n");
    print_mat(Gdisc); printf("\n");
  }
  /*
  if (TESTING) {
    // test out the evaluation of our FT...
    double dr = 0.01;
    double F1d_val;
    for (double r=0; r<=(maxR+0.5*dr); r += dr) {
      F1d_val = F1d_int(kmax*r, TESTING) / kmax;
      printf("#### %.15lf %.15lf %.15lf\n", kmax*r, r, F1d_val);
    }
  }
  */

  // 4 ways we have of dealing with our kpoints:
  const int KREALSMALL = 0;
  const int KSMALL     = 1;
  const int KENV       = 2;
  const int KLARGE     = 3;

  int ktype[Nkpt];
  for (int k=0; k<Nkpt; ++k) {
    double* tkpt = kpt[k];
    if ( (! zero( dot(tkpt, t1_0))) || (! zero( dot(tkpt, t2_0))) ||(tkpt[3] > kmax) )
      ktype[k] = KLARGE;
    else {
      // removed code to identify ksmall:
      if (tkpt[3] > krealsmall) ktype[k] = KENV;
      else ktype[k] = KREALSMALL;
    }
  }

  double** Gdk; // the FT of our discretization correction
  Gdk = new double *[Nkpt];
  for (int n=0; n<Nkpt; ++n)
    Gdk[n] = new double[9];

  if (TESTING) {
    printf("## Discretization correction FT: Gdk:\n");
  }
  for (int k=0; k<Nkpt; ++k) {
    double* tkpt = kpt[k]; // this kpoint
    if (TESTING) {
      printf("##\n## kpt: %8.5lf %8.5lf %8.5lf ", tkpt[0], tkpt[1], tkpt[2]);
      switch (ktype[k]) {
      case KLARGE:     printf("outside envelope or plane.\n"); break;
      case KENV:       printf("in envelope, but use fourier series.\n"); break;
      case KSMALL:     printf("no fourier series.\n"); break;
      case KREALSMALL: printf("extra careful evaluation.\n"); break;
      }
    }
    
    if (dcomp(tkpt[3], 0.)) {
      // Treat gamma as a special point:
      for (int d=0; d<9; ++d) Gdk[k][d] = Gdisc[d] + BETA2_ENVELOPE*GE[d]/(kmax*kmax);
    } 
    else {
      if (ktype[k] == KREALSMALL) {
	// NOTE: this is *exact* for all k... it's just that 
	// it's really tuned to work for *small* k values...
	if (LOOKUP) smallk_expand(Np, Dij, tkpt, tkpt[3]/kmax, Gdk[k], lambda,
 				  KPT_CART, coskR);
	else        smallk_expand(Np, Dij, tkpt, tkpt[3]/kmax, Gdk[k], lambda,
				  KPT_CART);
      }
      else {
	// 1. Fourier transform Dij...
	double Dk[9];
	if (LOOKUP) fourier_table(Np, Dij, tkpt, Dk, KPT_CART, coskR);
	else        fourier(Np, Dij, tkpt, Dk, KPT_CART);
	double Dkinv[9], err;
	// Use our (more stable?) symmetric inverse:
	careful_inverse(Dk, Dkinv, err);
	if (err > 1e-13) {
	  fprintf(stderr, "Somehow managed to get an appreciable error in Dkinv: err= %.5le\n", err);
	  fprintf(stderr, "this most likely indicated that Dk:\n");
	  for (int d=0; d<9; ++d) fprintf(stderr, " %.15le", Dk[d]);
	  fprintf(stderr, "\nis singular for kpt= %.12lf %.12lf %.12lf.\n",
		  tkpt[0], tkpt[1], tkpt[2]);
	  fprintf(stderr, "det(Dk) = %.15le", det(Dk));
	}
	
	for (int d=0; d<9; ++d)
	  Gdk[k][d] = Dkinv[d];
	if (TESTING) {
	  printf("##   Dk:    "); print_mat(Dk); printf("\n");
	  printf("##   Dk^-1: "); print_mat(Dkinv); printf(" err = %.5le\n",err);
	}
	
	if (ktype[k] != KLARGE) {
	  // 2. Calculate Gsc(k)
	  double fenv = kpoint_gaussian_envelope(tkpt[3] / kmax);
	  double Gsck[9];
	  if (ktype[k] == KENV) {
	    mult(GE, fenv/(tkpt[3]*tkpt[3]), Gsck);
	    if (TESTING) {
	      printf("##   Gsck:  "); print_mat(Gsck); printf("\n");
	    }
	    // 3. Subtract:
	    for (int d=0; d<9; ++d)
	      Gdk[k][d] -= Gsck[d];
	  }
	}
      }
    }
    if (TESTING) {
      printf("##   Gdk:   "); print_mat(Gdk[k]); printf("\n");
    }
  }

  // IFT:  We do this in 2 steps for efficiency.
  // IFT.a Invert the discrete piece:

  if (TESTING) {
    printf("##\n## Inverse FT -- discrete piece:\n");
  }
  for (int np=0; np<Nlatt; ++np) {
    point_type *tGd = Gd + np;
    // 1. inverse FT the discrete piece--does not use symmetry, so 
    //    points have the (implicit) weight 1/Nkpt
    if (LOOKUP) inv_fourier_table(*tGd, Nkpt, Gdk, kpt, KPT_CART, coskR);
    else        inv_fourier(*tGd, Nkpt, Gdk, kpt, KPT_CART);
    if (TESTING) {
      printf("## %3d %3d %3d  %11.6lf %11.6lf %11.6lf  %11.6lf\n", 
	     tGd->Runit[0], tGd->Runit[1], tGd->Runit[2],
	     dot(tGd->Rcart, t1_0), dot(tGd->Rcart, t2_0),
	     dot(tGd->Rcart, n_0), tGd->Rmagn);
      printf("##   Gd: "); print_mat(tGd->mat); printf("\n");
    }
  }
  
  // IFT.b Invert the EGF piece
  double Rmagn = -1, F1d_val;
  double gint[9];
  double gscale = cell_vol/area; // = V/|t1 x t2|
  if (TESTING) {
    printf("##\n## Inverse FT -- EGF piece:\n");
  }
  for (int np=0; np<Nlatt; ++np) {
    point_type *tGd = Gd + np;
    if (! dcomp(tGd->Rmagn, Rmagn) ) {
      Rmagn = tGd->Rmagn;
      if (Rmagn > 0) 
	F1d_val = F1d_int(kmax*Rmagn, TESTING)/kmax;
      else
	F1d_val = -BETA_PI_ENVELOPE/kmax;
    }
    // now scale:
    mult(GE, F1d_val*gscale, gint);

    if (TESTING) {
      printf("## %3d %3d %3d  %11.6lf %11.6lf %11.6lf  %11.6lf\n", 
	     tGd->Runit[0], tGd->Runit[1], tGd->Runit[2],
	     dot(tGd->Rcart, t1_0), dot(tGd->Rcart, t2_0),
	     dot(tGd->Rcart, n_0), tGd->Rmagn);
      printf("## discrete FT: "); print_mat(tGd->mat); printf("\n");
      printf("## planar EGF:  "); print_mat(gint); printf("\n");
    }
    // add correction to the .mat:
    for (int d=0; d<9; ++d) tGd->mat[d] += gint[d];
    if (TESTING) {
      printf("## GF sum:      "); print_mat(tGd->mat); printf("\n");
    }
  }

  // ROTATE.  Calculate rotation matrix theta = (t1,t2,n), and then
  // rotate each matrix to be in those coordinates
  if (ROTATE) {
    double theta[9], theta_inv[9];
    double vec_magn;
    vec_magn = sqrt(dot(t1_0,t1_0));
    for (int d=0; d<3; ++d) t1_0[d] /= vec_magn;
    vec_magn = sqrt(dot(t2_0,t2_0));
    for (int d=0; d<3; ++d) t2_0[d] /= vec_magn;
    for (int d=0; d<3; ++d) {
      theta[  d*3] = t1_0[d];
      theta[1+d*3] = t2_0[d];
      theta[2+d*3] = n_0[d];
    }
    careful_inverse(theta, theta_inv);
    // NOTE: this theta^-1 G theta -> G works because the multiplication
    // uses a temporary matrix to store the first multiplication:
    for (int np=0; np<Nlatt; ++np) 
      mult(theta_inv, Gd[np].mat, theta, Gd[np].mat);
  }

  // ****************************** OUTPUT ***************************

  // Now, we output EVERYTHING.
  if (VERBOSE) {
    // Human readable:
    printf("# ROTATE: %d\n", ROTATE);
    for (int np=0; np<Nlatt; ++np) {
      printf("# |R| = %8.5lf  Rhat = %8.5lf %8.5lf %8.5lf\n", 
	     Gd[np].Rmagn, Gd[np].Rcart[0], Gd[np].Rcart[1], Gd[np].Rcart[2]);
      printf("# "); print_mat(Gd[np].mat); printf("\n");
    }
  }
  
  sprintf(dump, "%d # [%d %d %d] x [%d %d %d]", Nlatt, 
	  t1_unit[0], t1_unit[1], t1_unit[2],
	  t2_unit[0], t2_unit[1], t2_unit[2]);
  if (ROTATE) sprintf(dump, "%s rotated to t1,t2,n coord", dump);
  write_Dij(stdout, Nlatt, Gd, dump);

  // ************************* GARBAGE COLLECTION ********************
  if (LOOKUP) free_coskR_table(coskR);

  free_kpts(Nkpt, kpt);

  for (int n=0; n<Nkpt; ++n)
    delete[] Gdk[n];
  delete[] Gdk;

  delete[] Dij;
  delete[] Gd;
    
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


// We use this to calculate Gsck for small kpoints.  We return
// (D(k))^-1 - fcut(k) (lambda(k))^-1
// calculated for small k. No more lambda2 part--there's no discontinuity in 1D.
// NOTE: We now use an envelope function on the lambda matrices
// so it requires a reasonable computation of kpoint_gaussian_envelope_1 = 1-fcut.
//
// The algorithm works like this:
// We calculate X(k) which is sum_R (1-1/2(k.R)^2-cos(k.R))D(R) by using
// fourier_smallk().  Then we calculate (lambda(k))^-1 and 
// lambda(k)^-1 lambda2(k) lambda(k)^-1 using eval_lambda_inv and
// eval_lambda2.  We then determine the eigenvalues of X and lambda_inv;
// define b = max(x_i)max(l_i); we make nseries = ln tol / ln b (tol = 1e-8?)
// This value determines if we use a series expansion, or if we directly
// do the inversion.  If nseries < NMAX, we expand out (1-B)^-1 as 
// B+B^2+B^3+... where B = lambda(k)^-1 X(k), add 1-fcut, and
// right multiply everything by lambda^-1.
// If nseries > NMAX, we write A=1-B, and calculate (A^-1 + (1-fcut))lambda^-1.
//
// Whew.  Hopefully, it will all work.

const double SMALLK_TOL = 1e-8;
const double LOG_SMALLK_TOL = log(SMALLK_TOL);
const int SMALLK_NMAX = 10;

void smallk_expand(int Np, point_type* p, double kpt[3], double k_kmax, double pk[9], 
		   double lambda[9][9], const int TYPE, const coskR_table_type *table) 
{
  double Xk[9], lambdak_inv[9];
  double B[9];
  int nseries;

  // 1. Calculate the matrices we need to work with:
  // call according to whether we have a lookup table or not...
  if (table == NULL) fourier_smallk(Np, p, kpt, Xk, TYPE);
  else               fourier_smallk_table(Np, p, kpt, Xk, TYPE, *table);
  eval_lambda_inv(lambda, kpt, lambdak_inv);
  //  eval_lambda2(lambda, lambda2, kpt, lambda2k);
  
  double xi[3], li[3], b;
  eigen(Xk, xi);
  eigen(lambdak_inv, li);
  // the eigenvalues are sorted, so the LARGEST ones are the last;
  // however, that doesn't consider the signs... so xk could be negative.
  // but the li eigenvalues are ALL positive--because we have a
  // stable lattice.
  if (fabs(xi[0]) > fabs(xi[2])) b = fabs(xi[0]);
  else                           b = fabs(xi[2]);
  b *= li[2];
  // Now, we take some logs:
  if (b != 0.0) nseries = lround(LOG_SMALLK_TOL / log(b));
  else nseries = 1;

  // Our B matrix:
  mult(lambdak_inv, Xk, B);
  double fcut_1 = kpoint_gaussian_envelope_1(k_kmax); // 1-fcut, but accurate.
  double fcut = 1-fcut_1;

  double A[9], Ainv[9];
  if (nseries > SMALLK_NMAX) {
    // do the inversion directly:
    for (int d=0; d<9; ++d) A[d] = ident[d] - B[d];  // A = 1-B
    careful_inverse(A, Ainv);
    for (int d=0; d<9; ++d) A[d] = Ainv[d] - ident[d]*fcut; // A = (1-B)^1 - fcut
    mult(A, lambdak_inv, pk);                    // pk = ((1-B)^1 - fcut)lambda^-1
    // our work is done here
  }
  else {
    // else... do the inversion with a power series in B
    double tempA[9];
    for (int d=0; d<9; ++d) tempA[d] = B[d]; // tempA = B
    for (int n=0; n<nseries; ++n) {
      for (int d=0; d<9; ++d) A[d] = ident[d]+tempA[d]; // A = 1+tempA
      mult(B, A, tempA);                            // tempA = B(1+tempA)
    }
    // lastly, we add (fcut-1) and multiply by lambda^-1... 
    for (int d=0; d<9; ++d) tempA[d] += ident[d]*fcut_1;
    mult(tempA, lambdak_inv, pk); 
  }
}

