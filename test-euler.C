/*
  Program: test-euler.C
  Author:  D. Trinkle
  Date:    April 26, 2004
  Purpose: This takes in:
           1. unit cell definition

	   and we compute the euler angles for each of the point group
	   operations.
	   ... now it does SO much more.

	   We go ahead and compute the Wigner d-matrices for beta=Pi/2
	   (the only case we really need worry about).  We combine that
	   with the ginv_i_g matrices, and produce a symmetrization
	   matrix S^l_(mi,nj) where m,n=-l..l and i,j are voigt matrices.
	   The end result is, for a given Ylm function expansion f^l_mi,
	   the symmetrized version is fsymm^l_mi = SUM_nj S^l_minj f^l_nj.

	   Couldn't be simpler.

	   We also do lots of tests along the way to make sure we're
	   putting everything together correctly.  The most important
	   is the final test, that shows that S^l is an idempotent operator.
	   This is rather non-obvious, but still true.

	   We test everything up to Lmax.

  Param.:  <cell>
           cell:    cell file describing our lattice

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

// ** indexing:
// ** n,m->(n+l)*(2l+1) + m+l = i
// ** (i%(2l+1))-l -> m
// ** (i/(2l+1))-l -> n
inline int nm2ind(const int &n, const int &m, const int &l) 
{ return (n+l)*(2*l+1) + m+l; }

// ** mi,nj->((n+l)*Ni+i)*
inline int minj2ind(const int &m, const int &i, 
		    const int &n, const int &j, const int &l) 
{ return ((m+l)*Ni+i)*(12*l+6)+(n+l)*Ni+j;}


/*
inline int nm2m(const int &nm, const int &l) 
{ return (nm%(2*l+1))-l; }

inline int nm2n(const int &nm, const int &l) 
{ return (nm/(2*l+1))-l; }
*/

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

// apply all of the mirror symmetries for a given l, n:
inline void mirror(double** D_l_nm, const int &n, const int &l) 
{
  int m;
  // vertical mirror  n,m -> n,-m * (-1)^(l+n):
  for (m=1; m<=l; ++m)
    D_l_nm[l][nm2ind(n,-m, l)] = neg(l+n)*D_l_nm[l][nm2ind(n,m, l)];
  // rotate 180  n,-m -> -n,m * (-1)^(n+m):
  for (m=-l; m<=l; ++m)
    D_l_nm[l][nm2ind(-n,m, l)] = neg(n+m)*D_l_nm[l][nm2ind(n,-m, l)];
  // diagonal mirror  n,m  ->  m,n * (-1)^(n+m):
  //             and  n,-m -> -m,n * (-1)^(n+m):
  for (m=(-n+1); m<n; ++m) {
    D_l_nm[l][nm2ind(m,n, l)] = neg(n+m)*D_l_nm[l][nm2ind(n,m, l)];
    D_l_nm[l][nm2ind(-m,n, l)] = neg(n+m)*D_l_nm[l][nm2ind(n,-m, l)];
  }
}

// apply the mirror symmetries for a given l, n (for G):
inline void mirror_sing(double** G_l_nm, const int &n, const int &l) 
{
  int m;
  // vertical mirror  n,m -> n,-m * (-1)^(l+n):
  for (m=1; m<=l; ++m)
    G_l_nm[l][nm2ind(n,-m, l)] = neg(l+n)*G_l_nm[l][nm2ind(n,m, l)];
  // rotate 180  n,-m -> -n,m * (-1)^(n+m):
  for (m=-l; m<=l; ++m)
    G_l_nm[l][nm2ind(-n,m, l)] = neg(n+m)*G_l_nm[l][nm2ind(n,-m, l)];
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
const int NUMARGS = 1;
const char* ARGLIST = "<cell>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice";

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
==== cell ====";

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
  char *cell_name;
  
  // Let's pull off the args:
  cell_name = args[0];

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
  // ** indexing:
  // ** n,m->(n+l)*(2l+1) + m+l = i
  // ** (i%(2l+1))-l -> m
  // ** (i/(2l+1))-l -> n

  // HEAVILY rewritten... 
  // following Altmann+Bradley, Phil.Trans. A225, 193-8 (1963)
  //  const int Lmax = 128;
  const int Lmax = 8;
  double* D_l_nm[Lmax+1];
  double* G_l_nm[Lmax+1]; // repeat the calculation...

  // Need this table for computing D_l_nm_new from G_l_nm:
  double sqrtfact[Lmax*2+1];
  sqrtfact[0] = 1.;
  for (i=1; i<=(2*Lmax); ++i)
    sqrtfact[i] = sqrtfact[i-1]*sqrt(i);
  for (i=0; i<=(2*Lmax); ++i)
    printf("%3d  sqrt(i!)= %11.4le\n", i, sqrtfact[i]);

  int l2_1;
  for (l=0; l<=Lmax; ++l) {
    l2_1 = 2*l+1;
    D_l_nm[l] = new double[l2_1*l2_1];
    G_l_nm[l] = new double[l2_1*l2_1];
    if (l==0) {
      D_l_nm[l][0] = 1.;
      G_l_nm[l][0] = 1.;
      continue;
    }

    // 1/2^l:
    double prod;
    for (i=0, prod=1.; i<l; ++i) prod *= 0.5;

    // **** G^l_nm ****
    // Calculate with G's instead:
    for (m=-l; m<=l; ++m)
      G_l_nm[l][nm2ind(l,m, l)] = 1.;
    // apply mirror symmetries: (reduced from D)
    mirror_sing(G_l_nm, l, l);
    // work our way down for each n...
    // no mirror symmetries?
    for (n=(l-1); n>=0; --n) {
      for (m=l; m>=0; --m)
	G_l_nm[l][nm2ind(n,m, l)] = 
	  (n-m+1)*G_l_nm[l][nm2ind(n+1,m, l)] -
	  (l+m)*G_l_nm[l][nm2ind(n+1,m-1, l)];
      
      // force G_l_nm[l][nm2ind(n,0, l)] to be zero if l+n is odd:
      if ( (l%2) != (n%2) )
	G_l_nm[l][nm2ind(n,0, l)] = 0.;
      // apply mirror symmetries:
      mirror_sing(G_l_nm, n, l);
    }
    // fix the zero values:
    for (n=-l; n<=l; ++n)
      if ( (l%2) != (abs(n)%2) )
	G_l_nm[l][nm2ind(0,n, l)] = 0.;
    // Finally, set our values in D_l_nm_new:
    // prod = 1/2^l
    for (n=-l; n<=l; ++n)
      for (m=-l; m<=l; ++m)
	// Now we include the Cmn term:
	D_l_nm[l][nm2ind(n,m, l)] = prod*G_l_nm[l][nm2ind(n,m, l)]*Cmm(n,m)
	  *sqrtfact[l+abs(n)]/(sqrtfact[l-abs(m)]*sqrtfact[l+abs(m)]
			       *sqrtfact[l-abs(n)]);
  }


  
  // ****************************** OUTPUT ***************************
  double alpha, beta, gamma;
  int improp;

  for (n=0; n<Nop; ++n) {
    double *grot = gcart[n];
    improp = euler_rot(grot, alpha, beta, gamma);
    printf("       [ %4.1lf  %4.1lf  %4.1lf ]\n", grot[0],grot[1],grot[2]);
    printf("g_%2d = [ %4.1lf  %4.1lf  %4.1lf ]\n",n+1,grot[3],grot[4],grot[5]);
    printf("       [ %4.1lf  %4.1lf  %4.1lf ]\n", grot[6],grot[7],grot[8]);
    printf("  i=%d  alpha = %9.5lfpi  beta = %9.5lfpi  gamma = %9.5lfpi\n\n",
	   improp, alpha/M_PI, beta/M_PI, gamma/M_PI);
  }

  for (l=0; l<=Lmax; ++l) {
    printf("D^(%d):\n", l);
    printf("m =   ");
    for (m=l; m>=-l; --m) printf("%7d      ", m);
    printf("\n");
    for (n=l; n>=-l; --n) {
      if (n==0) printf("n = %2d", n);
      else      printf("    %2d", n);
      for (m=l; m>=-l; --m)
	//	printf(" %7.4lf ", D_l_nm[l][nm2ind(n,m,l)]);
	// Now with the Cmm factor:
	printf(" %11.4le ", D_l_nm[l][nm2ind(n,m,l)]);
      printf("\n");
    }
    // Normalization test:
    printf("normerr:");
    double max_normerr = 0.;
    double normerr;
    for (n=l; n>=-l; --n) {
      double sum=0;
      for (m=-l; m<=l; ++m) 
	sum += D_l_nm[l][nm2ind(n,m,l)]*D_l_nm[l][nm2ind(n,m,l)];
      normerr = fabs(1.-sqrt(sum));
      printf(" n=%d: %6.4le", n, normerr);
      if (normerr > max_normerr) max_normerr = normerr;
    }
    printf("\nmax_normerr (%d): %6.4le\n", l, max_normerr);
    printf("\n");
  }

  printf("\n====\n");
  
  // Now, let's try our hand at testing this expansion:
  double x[3] = {0.8, 0.6, 0.4};
  double xmagn;
  xmagn = 1./sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); // normalize
  for (i=0; i<3; ++i) x[i] *= xmagn;
  double cost, phi; // base values
  calc_angle(x, cost, phi);
  for (l=0; l<=Lmax; ++l) {
    for (m=-l; m<=l; ++m) {
      // run through our l,m values
      printf("l= %3d  m= %4d\n", l, m);
      for (int g=0; g<Nop; ++g) {
	// run through our group operations:
	double *grot = gcart[g];
	improp = euler_rot(grot, alpha, beta, gamma);

	// rotate x first, and evaluate:
	double xrot[3];
	mult_vect(grot, x, xrot);
	double cost_rot, phi_rot;
	calc_angle(xrot, cost_rot, phi_rot);
	double YlmR, YlmI, Plm;
	Plm = gsl_sf_legendre_sphPlm(l, abs(m), cost_rot);
	YlmR = Plm * cos(m*phi_rot);
	YlmI = Plm * sin(m*phi_rot);

	printf("  %9.6lf %9.6lf %9.6lf | %9.6lf %9.6lf %9.6lf\n",
	       x[0], x[1], x[2], xrot[0], xrot[1], xrot[2]);
	printf("    i= %d  alpha= %5.2lfpi  beta= %4.2lfpi  gamma= %5.2lfpi\n",
	       neg(improp), alpha/M_PI, beta/M_PI, gamma/M_PI);
	printf("    Y%3d%4d(gx)= %.12le %+.12le i\n", l,m, YlmR, YlmI);

	// next, let's do the expansion:
	double i_sign = neg(improp*l); // (-1)^l for improper rots
	double Plm_t, YlmR_t, YlmI_t;
	if ( zero(beta) || dcomp(beta, M_PI) ) {
	  Plm_t = gsl_sf_legendre_sphPlm(l, abs(m), cost) * i_sign;
	  if (zero(beta)) {
	    YlmR_t = Plm_t * cos(m*(phi-alpha));
	    YlmI_t = Plm_t * sin(m*(phi-alpha));
	  }
	  else {
	    YlmR_t = Plm_t * cos(-m*(phi-alpha)) * neg(l);
	    YlmI_t = Plm_t * sin(-m*(phi-alpha)) * neg(l);
	  }
	}
	else {
	  // beta = pi/2
	  YlmR_t = 0.; YlmI_t = 0.;
	  for (n=-l; n<=l; ++n) {
	    //	    lambda = m*alpha+n*gamma;
	    // now the Cmm piece is already in D_l_nm
	    Plm_t = gsl_sf_legendre_sphPlm(l, abs(n), cost) * i_sign
	      * D_l_nm[l][nm2ind(n,m, l)];
	    YlmR_t += Plm_t*cos(n*phi-n*alpha-m*gamma);
	    YlmI_t += Plm_t*sin(n*phi-n*alpha-m*gamma);
	  }
	}
	printf("    -> transform= %.12le %+.12le i\n", YlmR_t, YlmI_t);
	if ( (! dcomp(YlmR, YlmR_t)) || (! dcomp(YlmI, YlmI_t)) ) {
	  printf("    MISMATCH\n");
	}
	
	printf("\n");
      }
    }
  }

  printf("\n====\n\n");
  double *SRl_mi_nj, *SIl_mi_nj;
  double Aij;
  for (l=0; l<=Lmax; ++l) {
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
	      //	      Aij = ginv_i_g[g][i][j];
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
	      //	      Aij = ginv_i_g[g][i][j];
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
	      //	      Aij = ginv_i_g[g][i][j];
	      SRl_mi_nj[minj2ind(m,i,n,j, l)] +=
		D_l_nm[l][nm2ind(m,n, l)]*cos(m*alpha+n*gamma)*i_sign*Aij;
	      //		D_l_nm[l][nm2ind(n,m, l)]*cos(m*alpha+n*gamma)*i_sign*Aij;
	      SIl_mi_nj[minj2ind(m,i,n,j, l)] -=
		D_l_nm[l][nm2ind(m,n, l)]*sin(m*alpha+n*gamma)*i_sign*Aij;
	      //		D_l_nm[l][nm2ind(n,m, l)]*sin(m*alpha+n*gamma)*i_sign*Aij;
	    }
    }
    // scale:
    for (k=0; k<d; ++k) {
      SRl_mi_nj[k] *= 1./Nop;
      SIl_mi_nj[k] *= 1./Nop;
    }
    
    // Output... eep.
    int printed;
    printf("l= %d\n", l);
    for (m=-l; m<=l; ++m)
      for (i=0; i<Ni; ++i) {
	printf("  F(%3d,%3d).%s= ", l,m,v_coord[i]);
	printed=0;
	for (n=-l; n<=l; ++n)
	  for (j=0; j<Ni; ++j) {
	    double SRl = SRl_mi_nj[minj2ind(m,i,n,j, l)];
	    double SIl = SIl_mi_nj[minj2ind(m,i,n,j, l)];
	    if (! zero(SRl)) {
	      if (printed) printf("                 ");
	      else printed = 1;
	      printf("%+.4le", SRl);
	    }
	    if (! zero(SIl)) {
	      if (printed) printf("                 ");
	      else printed = 1;
	      printf("%+.4le*I", SIl);
	    }
	    if ( (! zero(SRl)) || (! zero(SIl)) )
	      printf(" F%d,%d.%s\n", l,n,v_coord[j]);
	  }
	if (! printed) printf("0\n");
      }
    printf("\n");

    // Need to test idempotency of the Sl matrix by calculating
    // Sl^2 - Sl (or Sl*(Sl-1) ) and comparing to zero.

    // We run through each mi,nj point.  The value
    // (Sl*Sl)_(mi,nj) = SUM_pk Sl_mi,pk Sl_pk,nj 
    // which needs to equal Sl_mi,nj.  We do this by initializing our
    // sum with -Sl_mi,nj, and then adding the elements--the end
    // result should be zero for all mi,nj.
    printf("IDEMPOTENCY test:\n");
    double sumSR, sumSI;
    for (m=-l; m<=l; ++m)
      for (i=0; i<Ni; ++i) 
	for (n=-l; n<=l; ++n)
	  for (j=0; j<Ni; ++j) {
	    printf("Testing S^(%d)_(%3d.%s,%3d.%s)... ", l, m,v_coord[i],
		   n,v_coord[j]);
	    sumSR = -SRl_mi_nj[minj2ind(m,i,n,j, l)];
	    sumSI = -SIl_mi_nj[minj2ind(m,i,n,j, l)];
	    for (p=-l; p<=l; ++p)
	      for (k=0; k<Ni; ++k) {
		// Real part: (RR-II)
		sumSR += SRl_mi_nj[minj2ind(m,i,p,k, l)]*
		  SRl_mi_nj[minj2ind(p,k,n,j, l)] 
		  - SIl_mi_nj[minj2ind(m,i,p,k, l)]*
		  SIl_mi_nj[minj2ind(p,k,n,j, l)];
		// Imag part: (RI+IR)
		sumSI += SRl_mi_nj[minj2ind(m,i,p,k, l)]*
		  SIl_mi_nj[minj2ind(p,k,n,j, l)] 
		  + SIl_mi_nj[minj2ind(m,i,p,k, l)]*
		  SRl_mi_nj[minj2ind(p,k,n,j, l)];
	      }
	    if (!zero(sumSR)) printf("%+.8le", sumSR);
	    if (!zero(sumSI)) printf("%+.8le*I", sumSI);
	    if ( (! zero(sumSR)) || (! zero(sumSI)) )
	      printf(" --FAILED\n");
	    else printf("\n");
	  }

    delete[] SRl_mi_nj;
    delete[] SIl_mi_nj;    
  }
  

  // ************************* GARBAGE COLLECTION ********************
  for (l=0; l<=Lmax; ++l) {
    //    delete[] D_l_nm[l];
    delete[] G_l_nm[l];
    delete[] D_l_nm[l];
  }
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
