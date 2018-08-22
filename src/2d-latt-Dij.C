/*
  Program: 2d-latt-Dij.C
  Author:  D. Trinkle
  Date:    September 14, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. dislocation coordinate system
	   3. elastic GF in k-space Fourier series
	   4. second order correction to GF in k-space Fourier series
	   5. dynamical matrix evaluated at a set of points
	   6. a 2d kpoint mesh in our BZ
	   7. a set of points to evaluate the LGF

	   and we compute
	   1. the LGF for each point

	   ==== Note for 2d ====
	   This is the TWO-DIMENSIONAL version of discDij, and so most
	   of the algorithm comes from there, with changes appropriate
	   for a 2-dimensional LGF.  Note: the only thing that makes
	   this "2d" is that we make our LGF periodic along a threading
	   direction t.  You wouldn't believe the complications this
	   produces.
	   ==== Note for 2d ====

	   Changed: for small k points, we calculate the deviation from
	   D(k)^-1 very carefully.  We have to do this, else
	   the error in Gk(n) and G0(n) makes us want to kill ourselves.
	   Seriously: it will induce homi-/sui-cidal feelings when you
	   realize how bad it is.  Grep krealsmall to see when we switch
	   over.

	   It's worth mentioning that the small k expansion is EXACT
	   for ALL k values PROVIDED we are BELOW where the envelope
	   function kicks in... it's just not that efficient, so it's
	   worth our while to only use it when we absolutely have to.

  Param.:  <cell> <disl> <Gk(n)> <G0(n)> <Dij(R)> <kpts> <latt>
           cell:    cell file describing our lattice
	   disl:    dislocation coordinate system
           Gk(n):   Fourier series of elastic GF in k-space
           G0(n):   Fourier series of discontiunity correction
	   Dij(R):  dynamical matrix evaluated at a series of points
	   kpts:    2d kpt mesh (with weights)
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

           ==== disl ====
           t1.1 t1.2 t1.3  # line direction
           m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)
           n1.1 n1.2 n1.3  # mutual perp. vector
           # all three vectors are given in unit cell coord.
           ==== disl ====

	   ==== Gk(n) ====
	   Nmax sxx sxy ... szz # Maximum N value, optional sextic shift
	   n1 R(...) # n value -- real components
	   n1 I(...) # n value -- imag components
	   ...
	   # Read until EOF
	   # NOTE: only EVEN n values should appear!!
	   ==== Gk(n) ====

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

  Algo.:   Stolen (mostly) from discDij and changed accordingly

  Output:  The full-fledged 2d LGF at each point.
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
#include "dislocation.H"
#include "lambda.H" // for evaluation (on the fly) of lambda^-1 and lambda^(2)

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
		   double lambda[9][9], double lambda2[9][9][9],
		   const int TYPE, const coskR_table_type* table);

inline void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
			  double lambda[9][9], double lambda2[9][9][9],
			  const int TYPE, const coskR_table_type &table) 
{ return smallk_expand(Np, p, kpt, pk, lambda, lambda2, TYPE, &table);}
  
inline void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
			  double lambda[9][9], double lambda2[9][9][9],
			  const int TYPE) 
{ return smallk_expand(Np, p, kpt, pk, lambda, lambda2, TYPE, NULL);}



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
const int NUMARGS = 7;
const char* ARGLIST = "<cell> <disl> <Gk(n)> <G0(n)> <Dij(R)> <kpts> <latt>";

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'l', 's', 'r'};
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  disl:    dislocation coordinate system\n\
  Gk(n):   fourier series components of elastic GF in k-space\n\
  G0(n):   fourier series components of discontiunity correction\n\
  Dij(R):  dynamical matrix evaluated at a series of points\n\
  kpts:    symmetrized kpt mesh (with weights)\n\
  latt:    list of lattice points to eval Gdisc\n\
  -l: implement a lookup table for more efficient FT--highly recommended\n\
  -s: add in sextic shift (must be present in Gk(n))\n\
  -r: rotate to (mnt) coordinate system (only for lattice, NOT atomic pos)";

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
==== disl ====\n\
t1.1 t1.2 t1.3  # line direction\n\
m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)\n\
n1.1 n1.2 n1.3  # mutual perp. vector\n\
# all three vectors are given in unit cell coord.\n\
==== disl ====\n\
\n\
==== Gk(n) ====\n\
Nmax sxx sxy ... szz # Maximum N value, optional sextic shift\n\
n1 R(...) # n value -- real components\n\
n1 I(...) # n value -- imag components\n\
...\n\
# Read until EOF\n\
# NOTE: only EVEN n values should appear!!\n\
==== Gk(n) ====\n\
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
  int d, i, k, n; // General counting variables.

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
  int LOOKUP = flagon[0];  // use a lookup table for cos(k.R)--more efficient
  int SHIFT =  flagon[1];  // do we add the sextic shift?
  int ROTATE = flagon[2];  // do we rotate all of our matrices to mnt coord?

  if (TESTING) {
    if (LOOKUP) printf("## Using lookup tables for FT\n");
  }
  
  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char *cell_name = args[0];
  char *disl_name = args[1];
  char *egf_name = args[2];
  char *egf_disc_name = args[3];
  char *Dij_name = args[4];
  char *kpt_name = args[5];
  char *latt_name = args[6];

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

  //++ ==== disl ====
  infile = myopenr(disl_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", disl_name);
    exit(ERROR_NOFILE);
  }

  int t_unit[3], m_unit[3], n_unit[3];
  ERROR = read_disl(infile, t_unit, m_unit, n_unit);
  myclose(infile);
  //-- ==== disl ====

  if (ERROR) {
    fprintf(stderr, "Problem with dislocation: %s\n", disl_name);
    exit(ERROR);
  }
  double t0[3], m0[3], n0[3];
  make_dislcoord(cart, t0, m0, n0, t_unit, m_unit, n_unit, 
                 DISLCOORD_UNIT | DISLCOORD_NOCHANGE);
  if (TESTING) {
    printf("## vect    input                 cartesian\n");
    printf("## t = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           t_unit[0],t_unit[1],t_unit[2], t0[0],t0[1],t0[2]);
    printf("## m = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           m_unit[0],m_unit[1],m_unit[2], m0[0],m0[1],m0[2]);
    printf("## n = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           n_unit[0],n_unit[1],n_unit[2], n0[0],n0[1],n0[2]);
  }

  //++ ==== egf.ft ====
  infile = myopenr(egf_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", egf_name);
    exit(ERROR_NOFILE);
  }
  int Nmax;
  double** gncos, **gnsin;
  
  ERROR = read_2degf(infile, Nmax, gncos, gnsin, dump);  
  myclose(infile);
  //-- ==== egf.ft ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Problem with 2d EGF: %s\n", egf_name);
    exit(ERROR);
  }
  // else, determine the shift value:
  double sextic_shift[9];
  i = sscanf(dump, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     sextic_shift, sextic_shift+1, sextic_shift+2, 
	     sextic_shift+3, sextic_shift+4, sextic_shift+5,
	     sextic_shift+6, sextic_shift+7, sextic_shift+8);
  if (i!=9) {
    // the optional shift was NOT present, so ignore:
    for (d=0; d<9; ++d) sextic_shift[d] = 0;
  }
  if (TESTING) {
    printf("## EGF shift:");
    for (d=0; d<9; ++d) printf(" %.12lf", sextic_shift[d]);
    printf("\n");
    printf("## EGF FT components: %d\n", Nmax*2);
    for (i=0; i<=Nmax; ++i) {
      printf("## cos(%d t):", 2*i);
      print_mat(gncos[i]);
      printf("\n");
      printf("## sin(%d t):", 2*i);
      print_mat(gnsin[i]);
      printf("\n");
    }
  }

  //++ ==== egf.disc.ft ====
  infile = myopenr(egf_disc_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", egf_disc_name);
    exit(ERROR_NOFILE);
  }
  int Nmax0;
  double** g0ncos, **g0nsin;
  
  ERROR = read_2degf(infile, Nmax0, g0ncos, g0nsin);
  myclose(infile);
  //-- ==== egf.ft ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Problem with 2d EGF disc. correction: %s\n", 
	    egf_disc_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## EGF discontinuity correction FT components: %d\n", Nmax0*2);
    for (i=0; i<=Nmax; ++i) {
      printf("## cos(%d t):", 2*i);
      print_mat(g0ncos[i]);
      printf("\n");
      printf("## sin(%d t):", 2*i);
      print_mat(g0nsin[i]);
      printf("\n");
    }
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
  // Check to make sure the kpoint mesh is consistent with our threading
  // direction--i.e., dot(kpt, tlen) = 2*Pi*n 
  double tvect[3], tlen;
  mult_vect(cart, t_unit, tvect);
  tlen = sqrt(dot(tvect, tvect));
  for (k=0; k<Nkpt; ++k) {
    double dottest = fabs(dot(kpt[k], tvect));
    dottest = dottest - 2.0*M_PI*(round(dottest*0.5*M_1_PI));
    if (! zero(dottest)) break;
  }
  if (k != Nkpt) {
    fprintf(stderr, "kpt-set from %s does not appear to be consistent with dislocation coordinate system\n", kpt_name);
    fprintf(stderr, "%d kpoint: dot(kpt,t)= %8.5lf\n", k+1, dot(kpt[k],tvect));
    exit(ERROR_BADFILE);
  }

  if (TESTING) {
    printf("## %d kpoints\n", Nkpt);
    for (i=0; i<Nkpt; ++i)
      printf("## %8.5lf %8.5lf %8.5lf %.5le | m.k= %8.5lf n.k= %8.5lf t.k= %8.5lf\n", 
	     kpt[i][0], kpt[i][1], kpt[i][2], w[i],
	     dot(kpt[i],m0), dot(kpt[i],n0), dot(kpt[i],tvect));
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
  // in the (mn) plane:
  for (n=0; n<Nlatt; ++n)
    Gd[n].Rmagn = hypot(dot(Gd[n].Rcart, m0), dot(Gd[n].Rcart, n0));

  // For ease of use later, let's sort our points:
  sort_points(Gd, Nlatt);
  
  if (TESTING) {
    printf("## lattice points to calc Gd (SORTED):\n");
    printf("## NOTE: magnitude refers to distance in (mn) plane, and coord. are mnt\n");
    for (n=0; n<Nlatt; ++n)
      printf("## %3d %3d %3d  %11.6lf %11.6lf %11.6lf  %11.6lf\n", 
	     Gd[n].Runit[0], Gd[n].Runit[1], Gd[n].Runit[2],
	     dot(Gd[n].Rcart, m0), dot(Gd[n].Rcart, n0),
	     dot(Gd[n].Rcart, t0), Gd[n].Rmagn);
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
  double mink, maxR;
  maxR = Gd[Nlatt-1].Rmagn;
  for (k=0, mink=kmax; k<Nkpt; ++k)
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

  // value for the dc:
  double kmax_0;
  //  kmax_0 = 24*mink;
  kmax_0 = kmax*0.25;
  if (kmax_0 > kmax) kmax_0 = kmax;
  if (TESTING) {
    printf("## ++ Analytic spherical piece integration: discontinuity corr.\n");
    printf("## kmax_0 = %.5lf ; volume fraction = %.5lf\n", 
	   kmax_0, 4*M_PI*kmax_0*kmax_0*kmax_0/(3.*det(cart_rlv)));
  }  

  double ksmall, krealsmall;
  //  ksmall = kmax*0.5;     // the limit below which we don't use series
  ksmall = kmax;     // the limit below which we don't use series
  //  krealsmall = 0.25/maxR; // the limit below which we use expansion of k^2
  krealsmall = kmax*0.1; // the limit below which we really get "clever"
  if (krealsmall > (kpoint_envelope_alpha*kmax))
    krealsmall = (kpoint_envelope_alpha*kmax);
  if (krealsmall > (kpoint_envelope_alpha*kmax_0))
    krealsmall = (kpoint_envelope_alpha*kmax_0);
  
  // We no longer bother with lambda3 and lambda4... because we're
  // ridiculously careful about everything now.  We've got this calculation
  // in a bike-helmet with knee, elbow, and wrist pads.  Not to mention
  // the mylar vest.  It's that careful.
  double lambda[9][9], lambda2[9][9][9];
  make_lambda(Np, Dij, lambda);
  make_lambda2(Np, Dij, lambda2);
  if (TESTING) {
    printf("## ksmall =     %.5lf :below which we don\'t use Fourier series.\n", ksmall);
    printf("## krealsmall = %.5lf :below which we very carefully evaluated Gdk.\n", 
	   krealsmall);
    printf("## lambda: [ab,cd]: --non-zero only\n");
    print_lambda(lambda, "## ");
    printf("## lambda(2): [ab,cdef]: --non-zero only\n");
    print_lambda2(lambda2, "## ");
  }

  // 4 ways we have of dealing with our kpoints:
  const int KREALSMALL = 0;
  const int KSMALL     = 1;
  const int KENV       = 2;
  const int KLARGE     = 3;

  int ktype[Nkpt];
  for (k=0; k<Nkpt; ++k) {
    double* tkpt = kpt[k];
    if ( (! zero( dot(tkpt, t0))) || (tkpt[3] > kmax) ) ktype[k] = KLARGE;
    else 
      if (tkpt[3] > ksmall) ktype[k] = KENV;
      else 
	if (tkpt[3] > krealsmall) ktype[k] = KSMALL;
	else ktype[k] = KREALSMALL;
  }

  int np; // Loop over our points:

  double** Gdk; // the FT of our discrete correction
  Gdk = new double *[Nkpt];
  for (n=0; n<Nkpt; ++n)
    Gdk[n] = new double[9];

  if (TESTING) {
    printf("## Discretization correction FT: Gdk:\n");
  }
  for (k=0; k<Nkpt; ++k) {
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
      for (d=0; d<9; ++d) Gdk[k][d] = 0;
    } 
    else {
      if (ktype[k] == KREALSMALL) {
	// NOTE: this is *exact* for all k where fenv=1... it's just that
	// it's really tuned to work for *small* k values...
	if (LOOKUP) smallk_expand(Np, Dij, tkpt, Gdk[k], lambda, lambda2,
				  KPT_CART, coskR);
	else        smallk_expand(Np, Dij, tkpt, Gdk[k], lambda, lambda2,
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
	  for (d=0; d<9; ++d) fprintf(stderr, " %.15le", Dk[d]);
	  fprintf(stderr, "\nis singular for kpt= %.12lf %.12lf %.12lf.\n",
		  tkpt[0], tkpt[1], tkpt[2]);
	  fprintf(stderr, "det(Dk) = %.15le", det(Dk));
	}

	for (d=0; d<9; ++d)
	  Gdk[k][d] = Dkinv[d];
	if (TESTING) {
	  printf("##   Dk:    "); print_mat(Dk); printf("\n");
	  printf("##   Dk^-1: "); print_mat(Dkinv); printf(" err = %.5le\n",err);
	}
	
	if (ktype[k] != KLARGE) {
	  // 2. Calculate Gsc(k) and Gdisc(k)
	  double fenv = kpoint_envelope(tkpt[3] / kmax);
	  double fenv_0 = kpoint_envelope(tkpt[3] / kmax_0);
	  double Gsck[9], Gdisc[9];
	  if (ktype[k] == KENV) {
	    // Fourier series expansions for each:
	    // theta is the angle in plane, measured relative to m0.
	    double theta = atan2(dot(tkpt, n0), dot(tkpt, m0));
	    double g_scale =  1./(cell_vol * tkpt[3] * tkpt[3]);
	    eval_cossin_expansion(Nmax, theta, gncos, gnsin, Gsck);
	    mult(Gsck, g_scale, Gsck);
	    if (TESTING) {
	      double trueval[9];
	      eval_lambda_inv(lambda, tkpt, trueval);
	      printf("##   EGF (series): "); print_mat(Gsck); printf("\n");
	      printf("##   EGF (true):   "); print_mat(trueval); printf("\n");
	    }
	    
	    g_scale =  1./ cell_vol;
	    eval_cossin_expansion(Nmax0, theta, g0ncos, g0nsin, Gdisc);
	    mult(Gdisc, g_scale, Gdisc);
	    if (TESTING) {
	      double trueval[9];
	      eval_lambda2(lambda, lambda2, tkpt, trueval);
	      printf("##   EGF0 (series): "); print_mat(Gdisc); printf("\n");
	      printf("##   EGF0 (true):   "); print_mat(trueval); printf("\n");
	    }
	  }
	  else {
	    // ktype[k] == KSMALL
	    eval_lambda_inv(lambda, tkpt, Gsck);
	    eval_lambda2(lambda, lambda2, tkpt, Gdisc);
	  }
	  mult(Gsck, fenv, Gsck);
	  mult(Gdisc, fenv_0, Gdisc);
	  if (TESTING) {
	    printf("##   Gsck:  "); print_mat(Gsck); printf("\n");
	    printf("##   Gdisc: "); print_mat(Gdisc); printf("\n");
	  }

	  // 3. Subtract:
	  for (d=0; d<9; ++d)
	    Gdk[k][d] -= (Gsck[d] + Gdisc[d]);
	}
      }
    }
    if (TESTING) {
      printf("##   Gdk:   "); print_mat(Gdk[k]); printf("\n");
    }
  }
  
  // IFT:  We do this in *3* steps for efficiency.
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
	     dot(tGd->Rcart, m0), dot(tGd->Rcart, n0),
	     dot(tGd->Rcart, t0), tGd->Rmagn);
      printf("##   Gd: "); print_mat(tGd->mat); printf("\n");
    }
  }
  
  // IFT.b Invert the EGF piece
  // This is a modification of our gf-disc algorithm.  Also note that
  // we'll be doing the same loop TWICE--once for the EGF, and a second
  // time for the discrete correction.  In a lisp world, this wouldn't
  // happen...
  int Nupper = ( (Nmax > Nmax0) ? Nmax : Nmax0 );
  double **gRncos, **gRnsin; // for handling the inverse FT analytically
  double Fn_val[Nupper+1];
  double gint[9];
  double Rmagn = 0;
  
  ERROR = alloc_2degf(Nupper, gRncos, gRnsin);
  if (ERROR) {
    fprintf(stderr, "Error allocating our real-space Fourier series...\n");
    exit(ERROR);
  }

  if (TESTING) {
    printf("##\n## Inverse FT -- EGF piece:\n");
  }
  for (int np=0; np<Nlatt; ++np) {
    point_type *tGd = Gd + np;
    if (tGd->Rmagn > 0) {
      // only recalculate if we *need* to...
      if (! dcomp(tGd->Rmagn, Rmagn) ) {
	Rmagn = tGd->Rmagn;
	Fn_2d_int(Nmax, kmax*Rmagn, Fn_val, TESTING, 1024);
	// We modify the 0 value twice:
	// 1. we roll our own negative here so we can use (-1)^n later
	// 2. we add in the log(Rmagn) piece which was missing in the integral
	//	Fn_val[0] = log(Rmagn) - Fn_val[0];
	Fn_val[0] = -log(Rmagn) + Fn_val[0];
	if (TESTING) {
	  printf("## int-- Fn_val:\n");
	  for (n=0; n<=Nmax; ++n) {
	    printf("## %d %.12le\n", n*2, Fn_val[n]);
	  }
	}
	double g_scale = 0.5*M_1_PI, n_scale;
	int n_sign;
	for (n=0, n_sign = 1; n<=Nmax; ++n, n_sign = -n_sign) {
	  n_scale = g_scale * n_sign * Fn_val[n];
	  mult(gncos[n], n_scale, gRncos[n]);
	  mult(gnsin[n], n_scale, gRnsin[n]);
	}
      }
      // Get the continuum piece:
      double theta = atan2(dot(tGd->Rcart, n0), dot(tGd->Rcart, m0));
      eval_cossin_expansion(Nmax, theta, gRncos, gRnsin, gint);
    }
    else {
      // For R = 0, the continuum piece is really simple:
      double gscale = 0.5*M_1_PI;
      gscale *= kpoint_envelope_2d_int() + log(kmax); 
      mult(gncos[0], gscale, gint);
    }
    // now divide by |t|:
    mult(gint, 1./tlen, gint);

    if (TESTING) {
      printf("## %3d %3d %3d  %11.6lf %11.6lf %11.6lf  %11.6lf\n", 
	     tGd->Runit[0], tGd->Runit[1], tGd->Runit[2],
	     dot(tGd->Rcart, m0), dot(tGd->Rcart, n0),
	     dot(tGd->Rcart, t0), tGd->Rmagn);
      printf("## discrete FT: "); print_mat(tGd->mat); printf("\n");
      printf("## 2d EGF:      "); print_mat(gint); printf("\n");
    }
    // add correction to the .mat:
    for (d=0; d<9; ++d) tGd->mat[d] += gint[d];
    if (TESTING) {
      printf("## GF sum:      "); print_mat(tGd->mat); printf("\n");
    }
  }


  // IFT.c Invert the discontinuity correction piece
  // Same algorithm as above, but with slight name changes

  if (TESTING) {
    printf("##\n## Inverse FT -- discontinuity correction piece:\n");
  }
  Rmagn = 0;
  for (int np=0; np<Nlatt; ++np) {
    point_type *tGd = Gd + np;
    if (tGd->Rmagn > 0) {
      // only recalculate if we *need* to...
      if (! dcomp(tGd->Rmagn, Rmagn) ) {
	Rmagn = tGd->Rmagn;
	Fn_2d_x2_int(Nmax0, kmax_0*Rmagn, Fn_val, TESTING, 1024);
	if (TESTING) {
	  printf("## int-- Fn_x2_val:\n");
	  for (n=0; n<=Nmax0; ++n) {
	    printf("## %d %.12le\n", n*2, Fn_val[n]);
	  }
	}
	double g_scale = 0.5*M_1_PI/(Rmagn*Rmagn), n_scale;
	int n_sign;
	for (n=0, n_sign = 1; n<=Nmax0; ++n, n_sign = -n_sign) {
	  n_scale = g_scale * n_sign * Fn_val[n];
	  mult(g0ncos[n], n_scale, gRncos[n]);
	  mult(g0nsin[n], n_scale, gRnsin[n]);
	}
      }
      // Get the continuum piece:
      double theta = atan2(dot(tGd->Rcart, n0), dot(tGd->Rcart, m0));
      eval_cossin_expansion(Nmax, theta, gRncos, gRnsin, gint);
    }
    else {
      // For R = 0, the continuum piece is really simple:
      double g0scale = 0.5*M_1_PI;
      g0scale *= kpoint_x2_envelope_2d_int()*(kmax_0*kmax_0); 
      mult(g0ncos[0], g0scale, gint);
    }

    // now divide by |t|:
    mult(gint, 1./tlen, gint);
    if (TESTING) {
      printf("## %3d %3d %3d  %11.6lf %11.6lf %11.6lf  %11.6lf\n", 
	     tGd->Runit[0], tGd->Runit[1], tGd->Runit[2],
	     dot(tGd->Rcart, m0), dot(tGd->Rcart, n0),
	     dot(tGd->Rcart, t0), tGd->Rmagn);
      printf("## discFT+2d EGF: "); print_mat(tGd->mat); printf("\n");
      printf("## 2d discont.:   "); print_mat(gint); printf("\n");
    }
    // add correction to the .mat:
    for (d=0; d<9; ++d) tGd->mat[d] += gint[d];
    if (TESTING) {
      printf("## final GF:      "); print_mat(tGd->mat); printf("\n");
    }
  }

  // garbage collection:
  free_2degf(Nupper, gRncos, gRnsin);

  // SEXTIC.  Add in the shift to get the same results as someone
  // who was using sextic or the integral formalism would.  Just
  // for "backwards compatibility"
  // Note: we have to scale by 1/(2pi*t) to have correct units.
  if (SHIFT) {
    mult(sextic_shift, 0.5*M_1_PI/tlen, sextic_shift);
    if (TESTING) {
      printf("## Adding in sextic shift: ");
      print_mat(sextic_shift);
      printf("\n");
    }
    for (np=0; np<Nlatt; ++np)
      for (d=0; d<9; ++d) Gd[np].mat[d] += sextic_shift[d];
  }
  
  // ROTATE.  Calculate rotation matrix theta = (mnt), and then
  // rotate each matrix to be in those coordinates
  if (ROTATE) {
    double theta[9], theta_inv[9];
    for (d=0; d<3; ++d) {
      theta[  d*3] = m0[d];
      theta[1+d*3] = n0[d];
      theta[2+d*3] = t0[d];
    }
    careful_inverse(theta, theta_inv);
    // NOTE: this theta^-1 G theta -> G works because the multiplication
    // uses a temporary matrix to store the first multiplication:
    for (np=0; np<Nlatt; ++np) 
      mult(theta_inv, Gd[np].mat, theta, Gd[np].mat);
  }

  // ****************************** OUTPUT ***************************

  // Now, we output EVERYTHING.
  if (VERBOSE) {
    // Human readable:
    printf("# SHIFT:  %d\n", SHIFT);
    printf("# ROTATE: %d\n", ROTATE);
    for (np=0; np<Nlatt; ++np) {
      printf("# |R| = %8.5lf  Rhat = %8.5lf %8.5lf %8.5lf\n", 
	     Gd[np].Rmagn, Gd[np].Rcart[0], Gd[np].Rcart[1], Gd[np].Rcart[2]);
      printf("# "); print_mat(Gd[np].mat); printf("\n");
    }
  }
  
  if (SHIFT) 
    if (ROTATE)
      sprintf(dump, "%d %3d %3d %3d # sextic shift applied, rotated to mnt coord", 
	      Nlatt, t_unit[0], t_unit[1], t_unit[2]);
    else
      sprintf(dump, "%d %3d %3d %3d # sextic shift applied", Nlatt, 
	      t_unit[0], t_unit[1], t_unit[2]);
  else
    if (ROTATE)
      sprintf(dump, "%d %3d %3d %3d # rotated to mnt coord", Nlatt, 
	      t_unit[0], t_unit[1], t_unit[2]);
    else
      sprintf(dump, "%d %3d %3d %3d #", Nlatt, 
	      t_unit[0], t_unit[1], t_unit[2]);

  write_Dij(stdout, Nlatt, Gd, dump);

  // ************************* GARBAGE COLLECTION ********************
  if (LOOKUP) free_coskR_table(coskR);

  free_2degf(Nmax, gncos, gnsin);
  free_2degf(Nmax0, g0ncos, g0nsin);
  
  free_kpts(Nkpt, kpt);

  for (n=0; n<Nkpt; ++n)
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
// (D(k))^-1 - (lambda(k))^-1 - (lambda(k))^-1(lambda2(k))(lambda(k))^-1
// calculated for small k.
// NOTE: We DON'T use an envelope function on either of the lambda
// matrices, so this can ONLY be used in kpoint regimes where the
// envelope functions would return 1--which is generally a good idea
// anyway...
//
// The algorithm works like this:
// We calculate X(k) which is sum_R (1-1/2(k.R)^2-cos(k.R))D(R) by using
// fourier_smallk().  Then we calculate (lambda(k))^-1 and 
// lambda(k)^-1 lambda2(k) lambda(k)^-1 using eval_lambda_inv and
// eval_lambda2.  We then determine the eigenvalues of X and lambda_inv;
// define b = max(x_i)max(l_i); we make nseries = ln tol / ln b (tol = 1e-8?)
// This value determines if we use a series expansion, or if we directly
// do the inversion.  If nseries < NMAX, we expand out (1-B)^-1 as 
// B+B^2+B^3+... where B = lambda(k)^-1 X(k).  We subtract the lambda(2)
// piece from B*lambda^-1 before adding it, and left multiply everything
// by lambda^-1.
// If nseries > NMAX, we write A=1-B, and calculate (1-A^-1)lambda^-1.
// From that, we subtract the lambda(2) piece and return it.
//
// Whew.  Hopefully, it will all work.

const double SMALLK_TOL = 1e-8;
const double LOG_SMALLK_TOL = log(SMALLK_TOL);
const int SMALLK_NMAX = 10;

void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
		   double lambda[9][9], double lambda2[9][9][9],
		   const int TYPE, const coskR_table_type *table) 
{
  int d;
  double Xk[9], lambdak_inv[9], lambda2k[9];
  double B[9];
  int nseries;

  // 1. Calculate the matrices we need to work with:
  // call according to whether we have a lookup table or not...
  if (table == NULL) fourier_smallk(Np, p, kpt, Xk, TYPE);
  else               fourier_smallk_table(Np, p, kpt, Xk, TYPE, *table);
  eval_lambda_inv(lambda, kpt, lambdak_inv);
  eval_lambda2(lambda, lambda2, kpt, lambda2k);
  
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

  double A[9], Ainv[9];
  if (nseries > SMALLK_NMAX) {
    // do the inversion directly:
    for (d=0; d<9; ++d) A[d] = ident[d] - B[d];  // A = 1-B
    careful_inverse(A, Ainv);
    for (d=0; d<9; ++d) A[d] = ident[d]-Ainv[d]; // A = 1-(1-B)^1
    mult(A, lambdak_inv, pk);                    // pk = (1-(1-B)^1)lambda^-1
    // finally, subtract off the lambda2k piece:
    for (d=0; d<9; ++d) pk[d] -= lambda2k[d];
    // our work is done here
  }
  else {
    // else... do the inversion with a power series in B
    double tempA[9];
    for (d=0; d<9; ++d) tempA[d] = B[d]; // tempA = B
    for (int n=0; n<nseries; ++n) {
      for (d=0; d<9; ++d) A[d] = ident[d]+tempA[d]; // A = 1+tempA
      mult(B, A, tempA);                            // tempA = B(1+tempA)
    }
    // lastly, we'll multiply by B and lambda^-1... 
    mult(B, tempA, A);          // A = B^2(1+B(1+B(1+...)))
    mult(A, lambdak_inv, tempA); // tempA = B^2(1+B(1+B(1+...)))lambda^-1
    // ... and add B lambda^-1 - lambda2 piece:
    mult(B, lambdak_inv, A);     // A = B lambda^-1
    for (d=0; d<9; ++d) pk[d] = (A[d]-lambda2k[d]) + tempA[d];
  }
}

