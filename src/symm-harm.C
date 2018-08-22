/*
  Program: symm-harm.C
  Author:  D. Trinkle
  Date:    April 26, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. Ylm expansion of a matrix function

	   ... and we output the fully symmetrized spherical harmonic
	   expansion.  This requires us to compute the Wigner d-matrix
	   for beta=Pi/2.  Using that, we build the full symmetrization
	   matrix operator, which projects out all of the unsymm.
	   pieces.

	   The end result is, for a given Ylm function expansion f^l_mi,
	   the symmetrized version is fsymm^l_mi = SUM_nj S^l_minj f^l_nj.

	   Couldn't be simpler.

	   Previous computation results show that the *upper* limit
	   for l values we can do without incurring reasonable errors
	   is 104.  Above 113, all hell just breaks loose.
	   
  Param.:  <cell> <Ylm>
           cell:    cell file describing our lattice
	   Ylm:     Ylm expansion

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

           ==== Y(lm) ====
           Lmax d       # Maximum L value, dimensionality of vector
           l1 m1 R(...) # l,m pair -- real components
           l1 m1 I(...) # l,m pair -- imag components
           l2 m2 R(...) 
           l2 m2 I(...)
           ...
	   # Must only have EVEN values
           # Read until EOF or l=-1.
           ==== Y(lm) ====


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

// Tolerance for outputing a spherical harmonic component:

const double FFT_TOLER = 1.e-10;

//***************************** SUBROUTINES ****************************

// ** indexing:
// ** n,m->(n+l)*(2l+1) + m+l = i
// ** (i%(2l+1))-l -> m
// ** (i/(2l+1))-l -> n
inline int nm2ind(const int &n, const int &m, const int &l) 
{ return (n+l)*(2*l+1) + m+l; }

// ** mi->(m+l)*Ni+i
inline int mi2ind(const int &m, const int &i, const int &l) 
{ return (m+l)*Ni+i;}

// ** mi,nj->((m+l)*Ni+i)*(2l+1)*Ni + (n+l)*Ni + j
inline int minj2ind(const int &m, const int &i, 
		    const int &n, const int &j, const int &l) 
{ return ((m+l)*Ni+i)*(12*l+6)+(n+l)*Ni+j;}


// neg = (-1)^n
inline int neg(const int &n) 
{ if ( (n%2) == 0 ) return 1;
  else return -1; }

// This is the strange Cm'm = i^(|m'|+m) i^(|m|+m) prefactor Altmann
// and Bradley put in front of the D matrix... really really weird.
inline int Cmm (const int &n, const int &m) 
{
  int Cmm_fact=1;
  if (n>0) Cmm_fact = neg(n);
  if (m>0) Cmm_fact *= neg(m);
  return Cmm_fact;
}

// apply the mirror symmetries for a given l, n (for G):
inline void mirror_sing(double* G_l_nm, const int &n, const int &l) 
{
  int m;
  // vertical mirror  n,m -> n,-m * (-1)^(l+n):
  for (m=1; m<=l; ++m)
    G_l_nm[nm2ind(n,-m, l)] = neg(l+n)*G_l_nm[nm2ind(n,m, l)];
  // rotate 180  n,-m -> -n,m * (-1)^(n+m):
  for (m=-l; m<=l; ++m)
    G_l_nm[nm2ind(-n,m, l)] = neg(n+m)*G_l_nm[nm2ind(n,-m, l)];
}

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
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <Ylm>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {};
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Ylm:     Ylm expansion";

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
==== Y(lm) ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...) \n\
l2 m2 I(...)\n\
...\n\
# Only set up to handle EVEN values\n\
# Read until EOF or l=-1.\n\
==== Y(lm) ====";
  
int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n, p; // General counting variables.

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
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];

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


  //++ ==== Y(lm) ====
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

  const int L_upper_max = 104;
  if (!ERROR && (Lmax > (L_upper_max/2)) ) {
    ERROR = ERROR_BADFILE;
    fprintf(stderr, "Numerical limitations require Lmax less than %d\n",
	    L_upper_max);
  }
  //-- ==== Y(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm, IYlm, "##");
  }

  // ***************************** ANALYSIS **************************
  // Let's get our point group information:
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  double ginv_i_g[MAXop][Ni][Ni];
  // 1. point group operations:
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  for (n=0; n<Nop; ++n) 
    for (i=0; i<Ni; ++i)
      gen_ginv_i_g (i, gcart[n], gcart[inv_index[n]], ginv_i_g[n][i]);

  
  // 2. Wigner D-matrices
  // Remember: l really identifies 2l
  // ** indexing:
  // ** n,m->(n+l)*(2l+1) + m+l = i
  // ** (i%(2l+1))-l -> m
  // ** (i/(2l+1))-l -> n
  double* D_l_nm;
  double* G_l_nm; // repeat the calculation...
  double *SRl_mi_nj, *SIl_mi_nj;

  double alpha, beta, gamma;
  int improp;

  // Need this table for computing D_l_nm_new from G_l_nm:
  double sqrtfact[Lmax*4+1];
  sqrtfact[0] = 1.;
  for (i=1; i<=(4*Lmax); ++i)
    sqrtfact[i] = sqrtfact[i-1]*sqrt(i);
  if (TESTING) {
    printf("## sqrt_factorial:\n");
    for (i=0; i<=(2*Lmax); ++i)
      printf("## %3d  sqrt(i!)= %11.4le\n", i, sqrtfact[i]);
  }

  // ************************ OUTPUT / ANALYSIS **********************
  // We work "in place" to a certain extent--i.e., we output each l,m
  // combination as we compute it, so there's no internal storage.
  double zmax; // based on 00 value:
  for (d=0; d<9; ++d) {
    double z = sqrt(RYlm[0][0][d]*RYlm[0][0][d]+IYlm[0][0][d]*IYlm[0][0][d]);
    if (z > zmax) zmax = z;
  }
  if (zero(zmax)) zmax = 1.;

  printf("%d 9 # Lmax, dim;  l m real(flm_1) .. real(flm_d), then imag\n",
	 Lmax*2);

  // NOTE: we use l_index to be l/2 so things go easier...
  for (int l_index=0; l_index<=Lmax; ++l_index) {
    l = 2*l_index;
    int l2_1 = 2*l+1;
    // 1. Construct the Wigner D-matrices:
    D_l_nm = new double[l2_1*l2_1];
    G_l_nm = new double[l2_1*l2_1];
    if (l==0) {
      D_l_nm[0] = 1.;
      G_l_nm[0] = 1.;
    }
    else {
      // 1/2^l:
      double prod;
      for (i=0, prod=1.; i<l; ++i) prod *= 0.5;
      
      // **** G^l_nm ****
      // Calculate with G's instead:
      for (m=-l; m<=l; ++m)
	G_l_nm[nm2ind(l,m, l)] = 1.;
      // apply mirror symmetries: (reduced from D)
      mirror_sing(G_l_nm, l, l);
      // work our way down for each n...
      // no mirror symmetries?
      for (n=(l-1); n>=0; --n) {
	for (m=l; m>=0; --m)
	  G_l_nm[nm2ind(n,m, l)] = 
	    (n-m+1)*G_l_nm[nm2ind(n+1,m, l)] -
	    (l+m)*G_l_nm[nm2ind(n+1,m-1, l)];
	
	// force G_l_nm[nm2ind(n,0, l)] to be zero if l+n is odd:
	if ( (l%2) != (n%2) )
	  G_l_nm[nm2ind(n,0, l)] = 0.;
	// apply mirror symmetries:
	mirror_sing(G_l_nm, n, l);
      }
      // fix the zero values:
      for (n=-l; n<=l; ++n)
	if ( (l%2) != (abs(n)%2) )
	  G_l_nm[nm2ind(0,n, l)] = 0.;
      // Finally, set our values in D_l_nm_new:
      // prod = 1/2^l
      for (n=-l; n<=l; ++n)
	for (m=-l; m<=l; ++m)
	  // Now we include the Cmn term:
	  D_l_nm[nm2ind(n,m, l)] = prod*G_l_nm[nm2ind(n,m, l)]*Cmm(n,m)
	    *sqrtfact[l+abs(n)]/(sqrtfact[l-abs(m)]*sqrtfact[l+abs(m)]
				 *sqrtfact[l-abs(n)]);
    }
    
    // 2. Construct the symmetrization matrix:
    double Aij;

    d = 36+l*(144+144*l);
    SRl_mi_nj = new double[d];
    SIl_mi_nj = new double[d];
    // init
    for (k=0; k<d; ++k) {
      SRl_mi_nj[k] = 0.;
      SIl_mi_nj[k] = 0.;
    }
    // Let's build that matrix!!
    for (int g=0; g<Nop; ++g) {
      // run through our group operations:
      double *grot = gcart[g];
      improp = euler_rot(grot, alpha, beta, gamma);
      double i_sign = neg(improp*l); // (-1)^l for improper rots
      
      // NOTE: Be VERY careful of the SIGNs we use when adding/subtracting
      // to our matrices:
      if ( zero(beta) ) {
	// beta=0:
	for (m=-l; m<=l; ++m)
	  for (i=0; i<6; ++i)
	    for (j=0; j<6; ++j) {
	      Aij = ginv_i_g[g][j][i];
	      SRl_mi_nj[minj2ind(m,i,m,j, l)] += 
		cos(m*alpha)*i_sign*Aij;
	      SIl_mi_nj[minj2ind(m,i,m,j, l)] -= 
		sin(m*alpha)*i_sign*Aij;
	    }
	continue;
      }
      if ( dcomp(beta, M_PI) ) {
	// beta=Pi:
	i_sign *= neg(l); // include factor of (-1)^l
	for (m=-l; m<=l; ++m)
	  for (i=0; i<6; ++i)
	    for (j=0; j<6; ++j) {
	      Aij = ginv_i_g[g][j][i];
	      SRl_mi_nj[minj2ind(-m,i,m,j, l)] += 
		cos(m*alpha)*i_sign*Aij;
	      SIl_mi_nj[minj2ind(-m,i,m,j, l)] += 
		sin(m*alpha)*i_sign*Aij;
	    }
	continue;
      }
      // else, beta = pi/2
      for (m=-l; m<=l; ++m)
	for (i=0; i<Ni; ++i)
	  for (n=-l; n<=l; ++n)
	    for (j=0; j<Ni; ++j) {
	      Aij = ginv_i_g[g][j][i];
	      SRl_mi_nj[minj2ind(m,i,n,j, l)] +=
		D_l_nm[nm2ind(m,n, l)]*cos(m*alpha+n*gamma)*i_sign*Aij;
	      SIl_mi_nj[minj2ind(m,i,n,j, l)] -=
		D_l_nm[nm2ind(m,n, l)]*sin(m*alpha+n*gamma)*i_sign*Aij;
	    }
    }
    // scale:
    for (k=0; k<d; ++k) {
      SRl_mi_nj[k] *= 1./Nop;
      SIl_mi_nj[k] *= 1./Nop;
    }
    
    // 3. Apply rotations and output:
    // 3.a. Copy into a vector for easy multiplication...
    double *fYlmR = new double[Ni*l2_1];
    double *fYlmI = new double[Ni*l2_1];
    for (m=-l; m<=l; ++m) 
      for (i=0; i<Ni; ++i) {
	// We average in case--somehow--we didn't start off with symm. matrices
	fYlmR[mi2ind(m,i, l)] = 0.5*(RYlm[l_index][l+m][i6tom9[i]]+
				     RYlm[l_index][l+m][trans[i6tom9[i]]]);
	fYlmI[mi2ind(m,i, l)] = 0.5*(IYlm[l_index][l+m][i6tom9[i]]+
				     IYlm[l_index][l+m][trans[i6tom9[i]]]);
      }
    double fYlmR_symm[Ni], fYlmI_symm[Ni];
    // 3.b. Construct the symmetric values, and output:
    for (m=-l; m<=l; ++m) {
      for (i=0; i<Ni; ++i) {
	fYlmR_symm[i] = 0; fYlmI_symm[i] = 0;
	for (n=-l; n<=l; ++n)
	  for (j=0; j<Ni; ++j) {
	    // Real part: RR-II
	    fYlmR_symm[i]+=SRl_mi_nj[minj2ind(m,i,n,j,l)]*fYlmR[mi2ind(n,j,l)]
	      - SIl_mi_nj[minj2ind(m,i,n,j,l)]*fYlmI[mi2ind(n,j,l)];
	    // Imag part: RI+IR
	    fYlmI_symm[i]+=SRl_mi_nj[minj2ind(m,i,n,j,l)]*fYlmI[mi2ind(n,j,l)]
	      + SIl_mi_nj[minj2ind(m,i,n,j,l)]*fYlmR[mi2ind(n,j,l)];
	  }
      }

      // output:
      int allzero = 1;
      for (i=0; (i<Ni) && allzero; ++i) 
	allzero = ( (fabs(fYlmR_symm[i]) < FFT_TOLER*zmax) && 
		    (fabs(fYlmI_symm[i]) < FFT_TOLER*zmax) );
      if (! allzero) {
	printf("%3d %3d", l, m);
	for (d=0; d<9; ++d) printf(" %.12le", fYlmR_symm[m9toi6[d]]);
	printf("\n");
	printf("%3d %3d", l, m);
	for (d=0; d<9; ++d) printf(" %.12le", fYlmI_symm[m9toi6[d]]);
	printf("\n");
      }
    }

    // Garbage collection
    delete[] fYlmR;
    delete[] fYlmI;

    delete[] G_l_nm;
    delete[] D_l_nm;
    
    delete[] SRl_mi_nj;
    delete[] SIl_mi_nj;    
  }
  

  // ************************* GARBAGE COLLECTION ********************
  free_Ylm(Lmax, RYlm, IYlm);
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
