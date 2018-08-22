/*
  Program: semicont.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: My attempt at making a semi-continuum GF calculation... eep.
           We now do a *completely* different integration scheme.  We
	   handle the singularity by (almost) doing analytic integration
	   inside a sphere inscribed in the BZ.  The remaining contributions
	   are included by using MC integration over the remaining
	   geometry.  I tend to think this will work, especially as
	   Rc becomes larger than 1/km (the radius of the inscribed
	   sphere).

	   Part of me wonders if maybe this wouldn't all be easier
	   with kubic harmonics... but who am I to judge?

  Param.:  <cell> <G(lm)> <Rcut> <R1> <R2> <R3>
           cell:  cell file describing our lattice
           G(lm): file of spherical harmonic components of elastic GF
	   Rcut:  cutoff value (ta) in the semi-continuum model
	   Ri:    lattice coordinates of point to evaluate
	   MEMORY:   k-point grid to use (??)

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

	   ==== G(lm) ====
	   Lmax d       # Maximum L value, dimensionality of vector
	   l1 m1 R(...) # l,m pair -- real components
	   l1 m1 I(...) # l,m pair -- imag components
	   l2 m2 R(...) 
	   l2 m2 I(...)
	   ...
	   # Read until EOF or l=-1.
	   # NOTE: only EVEN l values should appear!!
	   ==== G(lm) ====

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

// If we are fed our GF in units of A^3/eV, then after we diagonalize
// invert and square root, we multiply by 1/(2Pi)*(ec*Na/10)^(1/2)
// to get the linear frequency nu in THz, once we divide by the
// square root of the mass in amu.  Whew.
const double THz_scale = 15.6333;


//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

#define __SWAP__(x,y) ({temp = x; x = y; y = temp;})
void sort3 (double a[3]) 
{
  double temp;
  if (a[0] > a[1]) __SWAP__(a[0], a[1]);
  if (a[1] > a[2]) __SWAP__(a[1], a[2]);
  if (a[0] > a[1]) __SWAP__(a[0], a[1]);
}
#undef __SWAP__

// Fold into the BZ
void fold_down(double kpt[3], double cart_b[9]);

// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
double max_k_BZ (double cart_b[9]);


const double __RAND_SCALE__ = 1./((double)((1<<31) - 1));
void rand_kpt (double kpt[3], double cart_b[9]) 
{
  kpt[0] = random()*__RAND_SCALE__;
  kpt[1] = random()*__RAND_SCALE__;
  kpt[2] = random()*__RAND_SCALE__;
  fold_down(kpt, cart_b);
}


// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 6;
const char* ARGLIST = "<cell> <G(lm)> <Rcut> <R1> <R2> <R3>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  cell:  cell file describing our lattice\n\
  G(lm): file of spherical harmonic components of elastic GF in k-space\n\
  Rcut:  cutoff value (ta) in the semi-continuum model\n\
  Ri:    lattice coordinates of point to evaluate\n\
  MEMORY:   number of k-points to use in MC integration";

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
";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n; // General counting variables.

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
  char *cell_name, *Ylm_name;
  double Rcut;
  int Runit[3];
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  sscanf(args[2], "%lf", &Rcut);
  sscanf(args[3], "%d", &(Runit[0]));
  sscanf(args[4], "%d", &(Runit[1]));
  sscanf(args[5], "%d", &(Runit[2]));
  
  if (Rcut < 0.) {
    fprintf(stderr, "Rcut less than 0?\n");
    exit(1);
  }

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

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d", &Lmax, &Ndim);
  if ( (Lmax < 0) || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
	    Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    // Turn Lmax into Lmax/2:
    Lmax = Lmax/2;
    // Start allocating...
    RYlm = new double **[Lmax+1];
    IYlm = new double **[Lmax+1];
    for (l=0; l<=Lmax; ++l) {
      RYlm[l] = new double *[4*l+1];
      IYlm[l] = new double *[4*l+1];
      for (m=0; m<=(4*l); ++m) {
	RYlm[l][m] = new double[Ndim];
	IYlm[l][m] = new double[Ndim];
	for (d=0; d<Ndim; ++d) {
	  // Initialize
	  RYlm[l][m][d] = 0.;
	  IYlm[l][m][d] = 0.;
	}
      }
    }
  }

  int no_more;
  int l_index, m_index, l_max_0;
  l_max_0 = 0; // Find the *true* Lmax value...
  for (no_more=0, n=0; (!ERROR) && (!no_more) && (!feof(infile)); ++n) {
    // Read the real elements, and ignore # comments:
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    sscanf(dump, "%d %d", &l, &m);
    if ( l < 0 ) {
      no_more = -1; // Our end of data marker...
      break;
    }
    // Check against:
    // a. l > 2*Lmax
    // b. |m| > l
    // c. l odd
    if ( (l > 2*Lmax) || (abs(m) > l) || ( (l%2) == 1 ) ) {
      fprintf(stderr, "Entry %d has bad lm components = %d %d\n", n+1, l, m);
      ERROR = ERROR_BADFILE;
      break;
    }
    if (l_max_0 < (l/2)) l_max_0 = (l/2);
    m_index = l+m;
    l_index = l/2; // l/2 :)
    
    // Code to pull off each elem one by one...
    char *startp, *endp;
    startp = dump; // Start at the beginning...
    // Need to "prime the pump" by going past first two entries.
    for (k=0; k<2; ++k) {
      strtol(startp, &endp, 10); // use l because we've got two integer elems.
      startp = endp;
    }
    for (d=0; d<Ndim; ++d) {
      RYlm[l_index][m_index][d] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough entries?
    if (d != Ndim) { 
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Didn't read enough entries on real entry %d.\n", n+1);
    }
    
    // Now, repeat with imag piece:
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d", &i, &j);
    if ( (i != l) || (j != m) ) {
      fprintf(stderr,
	      "Entry %d's imag. lm's don't match real lm's = %d %d vs. %d %d\n",
	      n+1, l, m, i, j);
      ERROR = ERROR_BADFILE;
      break;
    }
    
    // Code to pull off each elem one by one...
    startp = dump; // Start at the beginning...
    // Need to "prime the pump" by going past first two entries.
    for (k=0; k<2; ++k) {
      strtol(startp, &endp, 10); // use l because we've got two integer elems.
      startp = endp;
    }
    for (d=0; d<Ndim; ++d) {
      IYlm[l_index][m_index][d] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough entries?
    if (d != Ndim) { 
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Didn't read enough entries on imag entry %d.\n", n+1);
    }
  }
  myclose(infile);
  //-- ==== G(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  // Now, remove the "extra" l elements, if l_max_0 < Lmax.
  if (l_max_0 < Lmax) {
    // Free up the space:
    for (l=(l_max_0+1); l < Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) { // Since l is half it's value...
	delete[] RYlm[l][m];
	delete[] IYlm[l][m];
      }
      delete[] RYlm[l];
      delete[] IYlm[l];
    }
    if (TESTING) {
      printf("## Reset Lmax from %d to %d.\n", Lmax, l_max_0);
    }
    Lmax = l_max_0; // Correct Lmax value
  }

  if (TESTING) {
    // Print out *all* the elements we read in...
    printf("## Ylm components:\n");
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) {
	// Real part
	printf("## %3d %3d R ", 2*l, m-2*l);
	print_mat(RYlm[l][m]); printf("\n");
	// Imag part
	printf("## %3d %3d I ", 2*l, m-2*l);
	print_mat(IYlm[l][m]); printf("\n");
      }
    }
  }

  // VERBOSE: print out information about *size* of each component.
  if (TESTING) {
    printf("## Magnitude of terms:\n");
    double zmagn;
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) {
	zmagn = 0.;
	for (d=0; d<Ndim; ++d)
	  zmagn += RYlm[l][m][d]*RYlm[l][m][d] + 
	    IYlm[l][m][d]*IYlm[l][m][d];
	zmagn = sqrt(zmagn /(double)Ndim);
	printf("## %3d %3d %.8le\n", 2*l, m-2*l, zmagn);
      }
    }
  }


  // ***************************** ANALYSIS **************************
  // We separate our sum into two pieces done numerically.
  // 1 -- sphere insribed inside our BZ, evaluated analytically up
  //      to one piece which is done numerically
  // 2 -- monte carlo integration of remaining BZ

  double Rvect[3], Rmagn;
  double gint[Ndim]; // Our GF evaluated at a point!  
  // Spherical harmonics
  SH_workspace_type * w = allocate_SH_workspace(Lmax);
    
  // Make cartesian coordinates out of our point, and normalize:
  mult_vect(cart, Runit, Rvect);
  Rmagn = sqrt(dot(Rvect, Rvect));
  if (dcomp(Rmagn, 0.)) {
    Rmagn = 0;
    Rvect[0] = 0;
    Rvect[1] = 0;
    Rvect[2] = 0;
  }
  else {
    for (i=0; i<3; ++i) Rvect[i] *= 1./Rmagn;
  }

  // **** Step 1 -- spherical piece
  double gint_cont[Ndim];
  double kmax;
  
  kmax = max_k_BZ(cart_rlv);
  if (TESTING) {
    printf("## ++ Analytic spherical piece integration.\n");
    printf("## kmax = %.5lf ; volume fraction = %.5lf\n", 
	   kmax, 4*M_PI*kmax*kmax*kmax/(3.*det(cart_rlv)));
  }

  if (Rmagn > 0) {
    // Now, we need to do the numerical integration:
    double Fl_val[Lmax+1];
    
    Fn_int (Lmax, kmax*Rmagn, Fl_val, TESTING, 1024);
    
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
    
    double g_scale = 1./(4*M_PI*Rmagn);
    int l_sign;
    double ***RYlm_F, ***IYlm_F;
    RYlm_F = new double **[Lmax+1];
    IYlm_F = new double **[Lmax+1];
    for (l=0, l_sign = 1; l<=Lmax; ++l, l_sign = -l_sign) {
      RYlm_F[l] = new double *[4*l+1];
      IYlm_F[l] = new double *[4*l+1];
      for (m=0; m<=(4*l); ++m) {
	RYlm_F[l][m] = new double[Ndim];
	IYlm_F[l][m] = new double[Ndim];
	for (d=0; d<Ndim; ++d) {
	  // Initialize
	  RYlm_F[l][m][d] = RYlm[l][m][d] * l_sign * g_scale * Fl_val[l];
	  IYlm_F[l][m][d] = IYlm[l][m][d] * l_sign * g_scale * Fl_val[l];
	}
      }
    }
    if (TESTING) {
      // Print out *all* the elements we just calculated
      printf("## Ylm_F components:\n");
      for (l=0; l<=Lmax; ++l) {
	for (m=0; m<=(4*l); ++m) {
	  // Real part
	  printf("## %3d %3d R ", 2*l, m-2*l);
	  print_mat(RYlm_F[l][m]); printf("\n");
	  // Imag part
	  printf("## %3d %3d I ", 2*l, m-2*l);
	  print_mat(IYlm_F[l][m]); printf("\n");
	}
      }
    }
    
    // Get the continuum piece:
    eval_Ylm_expansion_R (Lmax, Rvect, Ndim, RYlm_F, IYlm_F, gint_cont, w);
    
    // garbage collection:
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) {
	delete[] RYlm_F[l][m];
	delete[] IYlm_F[l][m];
      }
      delete[] RYlm_F[l];
      delete[] IYlm_F[l];
    }
    delete[] RYlm_F;
    delete[] IYlm_F;
  }
  else {
    // For R = 0, the continuum piece is really simple:
    double g0_scale, x = kmax*Rcut;
    if ( x < 1e-4 ) {
      // Do a series expansion:
      g0_scale = M_2_SQRTPI*kmax*(1-x*x/6.);
    }
    else       
      g0_scale = erf(x)/Rcut;
    g0_scale *= M_1_PI*M_1_PI * 0.125; // divide by 8Pi^2

    for (d=0; d<Ndim; ++d)
      gint_cont[d] = g0_scale * RYlm[0][0][d];
  }
  
  // **** Step 2 -- Monte carlo integration of remaining pieces

  double gint_mc[Ndim];
  
  for (d=0; d<Ndim; ++d) gint_mc[d] = 0;
  
  // How many MC points to use?
  int Nmc = MEMORY;
  const int rand_seed = 1; // hard-coded for now.
  srandom (rand_seed); // seed it up!
  
  if (MEMORY > 0) {
    double kpt[3];
    double kc[3], kmagn2;  // kc = kpt in cartesian coord.
    double kmax2 = kmax*kmax, Rcut2 = Rcut*Rcut;
    double gtemp[Ndim], cos_kR, exp_kRc;
    double kR;
    int nbin;

    if (TESTING) {
      printf("## Monte Carlo evaluation for %d points.\n", Nmc);
    }

    // New method: we'll bin our results depending on the value of
    // k.R; this will allow us to try to correct some of the possible
    // error coming from adding up oscillating pieces.
    // We look at kpt . Runit, since kpt is in unit coord. and Runit
    // also is; this gives us a number that is multiplied by 2Pi, and
    // then we take the cosine.  Hence, we remove the integer part of
    // the dot product, and use that to determine which bin to put it
    // into.
    
    const int Nbin = 32; // For now, we hard-code this...
    double gint_mc_bin[Nbin][Ndim];
    for (n=0; n<Nbin; ++n) {
      for (d=0; d<Ndim; ++d)
	gint_mc_bin[n][d] = 0.;
    }
    
    for (n=0; n<Nmc; ++n) {
      rand_kpt(kpt, cart_rlv); // produce a random kpoint in BZ
      mult_vect(cart_rlv, kpt, kc);
      kmagn2 = dot(kc, kc);
      if (kmagn2 > kmax2) {
	// We've got a point to add!
	eval_Ylm_expansion_R (Lmax, kc, Ndim, RYlm, IYlm, gtemp, w);
	// k and R are in unit coord
	kR = kpt[0]*Runit[0] + kpt[1]*Runit[1] + kpt[2]*Runit[2];
	kR = kR - (int)kR;
	if (kR < 0) kR += 1.;
	if (kR >= 1) kR -= 1.;
	nbin = (int)(kR * Nbin);
	
	cos_kR = cos (2*M_PI*kR);
	exp_kRc = exp(-kmagn2*Rcut2);      // gaussian / k^2
	double scale = cos_kR * exp_kRc / kmagn2;
	for (d=0; d<Ndim; ++d)
	  gint_mc_bin[nbin][d] += gtemp[d] * scale;
	if (TESTING) {
	  printf("## k = %.5lf %.5lf %.5lf ", kc[0], kc[1], kc[2]);
	  print_mat(gtemp);
	  printf("  bin= %d  cos_kR= %.5lf  exp_kRc= %.5lf  scale= %.5lf\n",
		 nbin, cos_kR, exp_kRc, scale);
	}
      }
    }
    // scale out the number of points we sampled, and add up the
    // bins:
    double gmc_scale = 1./(double)Nmc;
    for (n=0; n<Nbin; ++n)
      for (d=0; d<Ndim; ++d) gint_mc[d] += gmc_scale*gint_mc_bin[n][d];
  }
  
  // **** Final summation:
  if (TESTING) {
    printf("## kmax = %8.5lf  Rcut = %8.5lf  exp(-(kmax * Rcut)^2 = %.4le\n",
	   kmax, Rcut, exp(-(kmax*kmax*Rcut*Rcut)));
  }

  if (TESTING) {
    printf("## Continuum piece of GF:\n");
    printf("## "); print_mat(gint_cont); printf("\n");
  }

  if (TESTING) {
    printf("## Monte carlo piece of GF:\n");
    printf("## "); print_mat(gint_mc); printf("\n");
  }

  for (d=0; d<Ndim; ++d) gint[d] = gint_cont[d] + gint_mc[d];

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.
  if (VERBOSE) {
    // Human readable:
    printf("# |R| = %8.5lf  Rhat = %8.5lf %8.5lf %8.5lf\n", 
	   Rmagn, Rvect[0], Rvect[1], Rvect[2]);
    printf("# "); print_mat(gint); printf("\n");
  }

  printf("1 %d  # Number of points, dimensionality\n", Ndim);
  printf("%d %d %d", Runit[0], Runit[1], Runit[2]);
  for (d=0; d<Ndim; ++d)
    printf(" %.11le", gint[d]);
  printf("\n");

  // ************************* GARBAGE COLLECTION ********************
  free_SH_workspace(w);

  // spherical harmonic components
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(4*l); ++m) {
      delete[] RYlm[l][m];
      delete[] IYlm[l][m];
    }
    delete[] RYlm[l];
    delete[] IYlm[l];
  }
  delete[] RYlm;
  delete[] IYlm;

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


// Fold into the BZ
// our k-point in rlv coordinates -- should be within [-1/2,1/2]^3
// in the beginning
void fold_down(double kpt[3], double cart_b[9]) 
{
  double bb[9];
  int a[3], a0[3], a1[3];
  double q[3], k0[3];
  double kmagn, qmagn;
  int i;
  int better;

  square(cart_b, bb);

  for (i=0; i<3; ++i) {
    kpt[i] = kpt[i] - (int)(kpt[i]);
    if (kpt[i] > 0.5) kpt[i] -= 1;
    if (kpt[i] <= -0.5) kpt[i] += 1;
    // Now, guess our possible shifts
    a0[i] = -2; a1[i] = 2;
    if (kpt[i] != 0) {
      if (kpt[i] > 0) a1[i] = 1;
      if (kpt[i] < 0) a0[i] = -1;
    }
    k0[i] = kpt[i];
  }
  // Magnitudes:
  kmagn = magnsq(bb, k0); // k0 -- our best guess.
  for (a[0]=a0[0]; a[0]<=a1[0]; ++(a[0]))
    for (a[1]=a0[1]; a[1]<=a1[1]; ++(a[1]))
      for (a[2]=a0[2]; a[2]<=a1[2]; ++(a[2])) {
	// Loop through possible shifts:
	for (i=0; i<3; ++i) q[i] = kpt[i] + a[i];
	qmagn = magnsq(bb, q);
	better = (qmagn < kmagn); // Better k-point?
	if (dcomp(qmagn, kmagn)) {
	  // Complicated tests...
	  better = (q[0] > k0[0]);
	  if (dcomp(q[0], k0[0])) {
	    better = (q[1] > k0[1]);
	    if (dcomp(q[1], k0[1]))
	      better = (q[2] > k0[2]);
	  }
	}
	if (better) {
	  kmagn = qmagn;
	  for (i=0; i<3; ++i) k0[i] = q[i];
	}
      }
  // Now we've got our "best" possible k-point:
  for (i=0; i<3; ++i) kpt[i] = k0[i];
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
