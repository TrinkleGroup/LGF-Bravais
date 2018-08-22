/*
  Program: gf-phonon.C
  Author:  D. Trinkle
  Date:    February 20, 2004
  Purpose: Given a spherical harmonic expansion for our elastic Greens
           function in k-space and the discretizaton error at a series
	   of lattice points (that is, our representation of the GF),
	   and a cell file, we go off and calculate the phonons (!) for
	   different k-points.

	   Sounds like a dream come true, eh?

  Param.:  <cell> <Gk(lm)> <Gd(R)> <kpts>
           cell:    cell file describing our lattice
           Gk(lm):  spherical harmonic components of elastic GF in k-space
	   Gd(R):   discrete corrections to GF at series of lattice sites
	   kpts:    list of kpoints to calculate our phonons

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

	   ==== Gd(R) ====
	   N                           # number of points
	   n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG
	   ...
	   nN.1 nN.2 nN.3  Gxx .. Gzz
	   ==== Gd(R) ====

	   ==== kpts ====
	   N k0 <type>     # Number of points, scale factor, C/L for units
	   k1.1 k1.2 k1.3
	   ...
	   kN.1 kN.2 kN.3
	   ==== kpts ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Given our lm spherical harmonic components, and the discrete
           corrections, we calculate our GF at each k-point, and then
	   diagonalize it to get our phonon spectrum.

	   If a (l,m) component pair is *not* given, it is assumed to
	   be zero.  However--we must *always* give both real and imag.
	   parts, even if the imag part is all zero.

  Output:  Just the phonons at each point -- scaling to get meV or THz???
           ==== output ====
           N      # Number of grid points, dimensionality
	   k1.1 k1.2 k1.3  w1.1 w1.2 w1.3  # k-points in Cart, frequencies
	   ...
	   kN.1 kN.2 kN.3  w1.1 w1.2 w1.3  # k-points in Cart, frequencies
           ==== output ====
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include <gsl/gsl_sf_bessel.h>   // Spherical bessel functions
#include <gsl/gsl_integration.h> // Needed for doing integration
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "kpts.H"
#include "eigen.H"
#include "cell.H"  // Do we still use this?
#include "cont-int.H" // For doing the integration over the continuous piece
#include "sphere-harm.H" // Spherical harmonic evaluation


// If we are fed our GF in units of A^3/eV, then after we diagonalize
// invert and square root, we multiply by 1/(2Pi)*(ec*Na/10)^(1/2)
// to get the linear frequency nu in THz, once we divide by the
// square root of the mass in amu.  Whew.
const double THz_scale = 15.6333023;


//****************************** SUBROUTINES ****************************

// Calculate (kk)
void a_mult_a (double a[3], double Cijkl[9][9], double aa[9]) 
{
  int i, j, k, l;
  for (i=0; i<3; ++i) {
    for (j=0; j<i; ++j)
      aa[index(i,j)] = aa[index(j,i)];
    for (   ; j<3; ++j) {
      aa[index(i,j)] = 0.;
      for (k=0; k<3; ++k)
        for (l=0; l<3; ++l)
          aa[index(i,j)] += a[k]*Cijkl[index(k,i)][index(j,l)]*a[l];
    }
  }
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
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <Gk(lm)> <Gd(R)> <kpts>";

const int NFLAGS = 4;
const char USERFLAGLIST[NFLAGS] = {'k', 'n', 'x', 'y'};  // flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice (see note below)\n\
  Gk(lm):  spherical harmonic components of elastic GF in k-space\n\
  Gd(R):   discrete corrections to GF at series of lattice sites\n\
  kpts:    list of kpoints to calculate our phonons\n\
  -y:      use the Ylm expansion for semicontinuum piece\n\
           N.B.: if NOT used, Cij MUST be in eV/A^3\n\
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
==== Gk(lm) ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...) \n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
# NOTE: only EVEN l values should appear!!\n\
==== Gk(lm) ====\n\
\n\
==== Gd(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG\n\
...\n\
nN.1 nN.2 nN.3  Gxx .. Gzz\n\
==== Gd(R) ====\n\
\n\
==== kpts ====\n\
N k0 <type>     # Number of points, scale factor, C/L for units\n\
k1.1 k1.2 k1.3\n\
...\n\
kN.1 kN.2 kN.3\n\
==== kpts ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used
  int KPOINT_CART, KPOINT_NUM; // How do we want to output?
  int XMGRACE;
  
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

  // Flags
  KPOINT_CART = flagon[0];
  KPOINT_NUM = flagon[1];
  XMGRACE = flagon[2];
  if ( (! KPOINT_CART) && (! KPOINT_NUM) )
    KPOINT_NUM = 1;
  int YLMEXPAND = flagon[3];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *Gd_name, *kpts_name;
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  Gd_name = args[2];
  kpts_name = args[3];
  
  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input -- we use it if no YLMEXPAND
  double Cijkl[9][9];  
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

  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);

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


  //++ ==== Gd(R) ====
  int Np;
  double kmax;
  int **Runit;
  double **Gd;
  double **Rcart;
  
  infile = myopenr(Gd_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Gd_name);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %lf", &Np, &kmax);
  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in %s\n", Np, Gd_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    Runit = new int *[Np];
    Gd = new double *[Np];
    Rcart = new double *[Np];
  }
  
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    Runit[n] = new int[3];
    Gd[n] = new double[Ndim];
    sscanf(dump, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	   &(Runit[n][0]), &(Runit[n][1]), &(Runit[n][2]),
	   &(Gd[n][0]), &(Gd[n][1]), &(Gd[n][2]),
	   &(Gd[n][3]), &(Gd[n][4]), &(Gd[n][5]),
	   &(Gd[n][6]), &(Gd[n][7]), &(Gd[n][8]));
    // Cartesian coord:
    Rcart[n] = new double[3];
    mult_vect(cart, Runit[n], Rcart[n]);
  }
  myclose(infile);
  if (n!=Np) {
    fprintf(stderr, "Didn't read enough points in %s\n", Gd_name);
    ERROR = ERROR_BADFILE;
  }
  if (ERROR) exit(ERROR);
  //-- ==== Gd(R) ====

  if (TESTING) {
    printf("## GF discrete correction:\n");
    for (n=0; n<Np; ++n) {
      printf("## %3d %3d %3d ", Runit[n][0], Runit[n][1], Runit[n][2]);
      print_mat(Gd[n]);
      printf("\n");
    }
  }


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

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %lf %c", &Nkpt, &k0, &type_ident);

  if (Nkpt < 1) {
    fprintf(stderr, "No points listed in %s\n", kpts_name);
    ERROR = ERROR_BADFILE;
  }
  if (k0 <= 0.) {
    fprintf(stderr, "Bad scale factor in %s\n", kpts_name);
    ERROR = ERROR_BADFILE;
  }
  type_ident = tolower(type_ident);
  if ( (type_ident != 'l') && (type_ident != 'c') ) {
    fprintf(stderr, "Invalid type character (%c) in %s; should be C or L\n", 
	    type_ident, kpts_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    LATT = (type_ident == 'l'); // convert from rlv to cartesian
  }
  if (ERROR) {
    myclose(infile);
    exit(ERROR);
  }
  else {
    // Allocate!
    kpt = new double *[Nkpt];
    for (i=0; i<Nkpt; ++i) kpt[i] = new double[4];
  }

  // Now, do the readin'!
  for (n=0; (!ERROR) && (n<Nkpt) && (!feof(infile)); ++n) {
    double* kpt_p;
    kpt_p = kpt[n];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", &(kpt_p[0]), &(kpt_p[1]), &(kpt_p[2]));
    for (i=0; i<3; ++i) kpt_p[i] *= k0; // scale!
    // Convert to lattice coord
    double tempvec[3];
    if (! LATT) {
      mult_vect(cart, kpt_p, tempvec);
      for (i=0; i<3; ++i) kpt_p[i] = tempvec[i] / (2*M_PI);
    }
    // Fold into BZ:
    fold_down(kpt_p, cart_rlv);
    // Translate back into cartesian coordinates, and set
    // fourth entry to magnitude:
    mult_vect(cart_rlv, kpt_p, tempvec);
    kpt_p[3] = sqrt(dot(tempvec, tempvec));
    for (i=0; i<3; ++i) kpt_p[i] = tempvec[i];
  }
  // Make sure we read enough points...
  if (n != Nkpt) {
    fprintf(stderr, "Not enough points in %s\n", kpts_name);
    ERROR = ERROR_BADFILE;
  }
  myclose(infile);
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
  // Now, we're gonna compute our GF FT at each k-point.  Note: we're
  // going to divide out the 1/k^2 divergence to do this.

  // Allocate space for our GF FT: G_ab(k_vec):
  double **Gabk;
  Gabk = new double*[Nkpt];
  for (k=0; k<Nkpt; ++k) {
    Gabk[k] = new double[Ndim];
    for (d=0; d<Ndim; ++d) Gabk[k][d] = 0.; // initial value.
  }
  
  // Some workspace:
  // Spherical harmonics
  SH_workspace_type * w = allocate_SH_workspace(Lmax);

  // First piece: continuum piece with envelope function
  for (k=0; k<Nkpt; ++k) {
    // All we do is just evaluate the Ylm expansion
    if ( dcomp(kpt[k][3], 0) || (kpt[k][3] > kmax) ) {
      // Zero magnitude, or beyond our cutoff:
      for (d=0; d<Ndim; ++d) Gabk[k][d] = 0.;
    }
    else {
      // Note carefully: we do this all with (basically) a normalized kpt
      if (YLMEXPAND) {
	eval_Ylm_expansion_R (Lmax, kpt[k], Ndim, RYlm, IYlm, Gabk[k], w);
      }
      else {
	double tkpt[3]; // normalize...
	double g[9];
	for (d=0; d<3; ++d) tkpt[d] = kpt[k][d]*(1./kpt[k][3]);
	a_mult_a(tkpt, Cijkl, g);
	careful_inverse(g, Gabk[k]);
      }
      // Multiply by our envelope function; we also need to scale by vol here:
      double g_scale = kpoint_envelope(kpt[k][3] / kmax) / cell_vol;
      mult(Gabk[k], g_scale, Gabk[k]);
      
    }
  }
  if (TESTING) {
    printf("## Continuum GF:\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## %8.5lf %8.5lf %8.5lf ", kpt[k][0], kpt[k][1], kpt[k][2]);
      print_mat(Gabk[k]);
      printf(" / %8.5lf\n", kpt[k][3]*kpt[k][3]);
    }
  }

  // Second piece: discrete piece
  if (TESTING) {
    printf("## Discrete GF:\n");
  }
  for (k=0; k<Nkpt; ++k) {
    double g_disc[Ndim];
    for (d=0; d<Ndim; ++d) g_disc[d] = 0;

    double cos_kR;
    double kmagn2;
    kmagn2 = kpt[k][3] * kpt[k][3];
    for (n=0; n<Np; ++n) {
      cos_kR = cos(dot(kpt[k], Rcart[n])) * kmagn2;
      for (d=0; d<Ndim; ++d) g_disc[d] += cos_kR * Gd[n][d];
    }
    
    for (d=0; d<Ndim; ++d) Gabk[k][d] += g_disc[d];
    if (TESTING) {
      printf("## %8.5lf %8.5lf %8.5lf ", kpt[k][0], kpt[k][1], kpt[k][2]);
      print_mat(g_disc);
      printf(" / %8.5lf\n", kpt[k][3]*kpt[k][3]);
    }
  }
  
  
  if (TESTING) {
    printf("## **** Final Green's function FT ****\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## G~(%8.5lf %8.5lf %8.5lf) = ",
	     kpt[k][0], kpt[k][1], kpt[k][2]);
      print_mat(Gabk[k]);
      printf(" / %8.5lf\n", kpt[k][3]*kpt[k][3]);
    }
  }

  // *********************** FREQUENCY CALCULATION *******************
  // Compute the frequencies:
  double nu[Nkpt][3], nu_scale;
  nu_scale = THz_scale/sqrt(atomic_mass);
  if (TESTING) {
    printf("## **** Frequency calculation ****\n");
  }
  for (k=0; k<Nkpt; ++k) {
    double lambda[3];
    eigen(Gabk[k], lambda); // Get the eigenvalues.
    if (TESTING) {
      printf("## lambda_%d = %.5lf %.5lf %.5lf\n", k+1, lambda[0],lambda[1],lambda[2]);
    }
    for (d=0; d<3; ++d) {
      if (dcomp(lambda[d], 0.)) {
	if (! dcomp(kpt[k][3], 0.)) {
	  fprintf(stderr, "At kpt %d %.3lf %.3lf %.3lf, got zero inverse frequency without zero k-point magnitude\n",
		  k+1, kpt[k][0], kpt[k][1], kpt[k][2]);
	}
	nu[k][d] = 0.;
	continue;
      }
      if (lambda[d] < 0) 
	nu[k][d] = -nu_scale*kpt[k][3]/sqrt(-lambda[d]);
      else
	nu[k][d] =  nu_scale*kpt[k][3]/sqrt(lambda[d]);
    }
    // Now, we'll order them by size:
    sort3(nu[k]);
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
  free_SH_workspace(w);
  
  // Greens function
  for (n=0; n<Nkpt; ++n) delete[] Gabk[n];
  delete[] Gabk;

  // Grid points
  for (i=0; i<Nkpt; ++i) delete[] kpt[i];
  delete[] kpt;

  // G_disc
  for (n=0; n<Np; ++n) {
    delete[] Runit[n];
    delete[] Gd[n];
    delete[] Rcart[n];
  }
  delete[] Runit;
  delete[] Gd;
  delete[] Rcart;
  
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
