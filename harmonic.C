/*
  Program: harmonic.C
  Author:  D. Trinkle
  Date:    February 12, 2004
  Purpose: Calculate the spherical harmonic components for a vector 
           function defined on a (theta, phi) grid, where phi is gridded
	   uniformly and theta is gridded by gaussian quadratures.

  Param.:  <gauss.J> <grid.IJ>
           gauss.J:  file of gaussian quadrature points and weights
           grid.IJ:  grid points and the function evaluated at those points

	   ==== grid.IJ ====
	   N d       # Number of points = I*J, dimensionality of vector
	   l1 m1 n1  f1_1 .. f1_d # Direction cosines (l^2 + m^2 + n^2 = 1)
	   l2 m2 n2  f2_1 .. f2_d
	   ...
	   lN mN nN  fN_1 .. fN_d
	   ==== grid.IJ ====

  Flags:   MEMORY:  Used to force a given Lmax value
	   VERBOSE: ??
	   TESTING: ??

  Algo.:   Given our function at each point, we calculate the lm spherical
           harmonic components.

           We're *very* generous about the exact order of the points read
	   in from grid.IJ; as we read each point, we calc. what theta_j
	   and phi_i it corresponds to, and put it in the correct place.
	   As a sanity check, we also keep track of which points have
	   been hit, and make sure we mapped exactly one point to each.

	   Call this fd_tp[d=0..d-1][theta_j][phi_i].  We transform in
	   two steps.  First, define Lmax = [I/2] - 1 (maximum l value
	   we'll have).

	   1. FFT the phi_i component to fd_tm[d][theta_j][m], where
	   m = 0..(I-1).  We talk about -m as I-m.  This is a *complex*
	   array, stored as packed reals.  We do this with the GSL
	   FFT routines.

	   2. Gaussian quadrature the theta_j component to get
	   fd_lm[d][l][m], where l = 0..Lmax.  Again, we end up with
	   a *complex* array, stored as packed reals.

  Output:  Finally, it's a simple matter of outputing all of the
	   components:

	   ==== output ====
	   Lmax dim  # maximum l value in run
	   0  0 R(f00_1) .. R(f00_d)  # Real components...
	   0  0 I(f00_1) .. I(f00_d)  # Real components...
	   1 -1 R(f1-1_1) ..
	   1 -1 I(f1-1_1) ..
	   1  0 R(f10_1) ..
	   1  0 I(f10_1) ..
	   ...
	   Lmax Lmax R(fLL_1) ..
	   Lmax Lmax I(fLL_1) ..  # Last value.
	   ==== output ====

	   Changed so that we *only* output a line if it is not
	   identically zero (both real and imag).

	   There are a total of (Lmax+1)^2 components to output.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_fft_complex.h> // FFT routines
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include "io.H"   // All of our "read in file", etc.

// Tolerance for outputing a spherical harmonic component:

const double FFT_TOLER = 1.e-11;

//****************************** SUBROUTINES ****************************

// Needed for our packed complex arrays:
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

inline double zmagn (const double x, const double y) 
{
  return sqrt(x*x+y*y);
}

inline int test_close(double x, double y) 
{
  return (fabs(x-y) < 1.e-5);
}

int binary_search (double n, double *xg, int N) 
{
  int lo, hi, mid;
  double mid_val;
  lo = 0;
  hi = N-1;
  mid = (lo+hi)/2; mid_val = xg[mid];
  while ( (!test_close(n, mid_val)) && (lo <= hi) ) {
    if (n < mid_val) {
      hi = mid-1;
      mid = (lo+hi)/2;
    }
    else {
      lo = mid+1;
      mid = (lo+hi)/2;
    }
    // Sanity check... not sure we need this?
    if (mid < 0) mid = 0;
    if (mid >= N) mid = N-1;
    mid_val = xg[mid];
  }
  if (lo > hi) // didn't find it...
    return -1;
  else
    return mid;
}
  

// Returns what ival (for phi) and jval (for theta) the point
// matches; returns 0 if a "reasonable" match is found; else
// returns a non-zero number (i.e., ERROR)
// Takes in direction cosines lmn...
int check_ij (double lmn[3], int Nphi, int Ntheta, double *xg,
	      int &ival, int &jval) 
{
  // First find lmn[2] in xg list (quicky binary search)
  jval = binary_search(lmn[2], xg, Ntheta);
  // Now, guess what value of phi we have...
  double phi;
  phi = atan2(lmn[1], lmn[0]);
  if (phi < 0) phi += 2*M_PI; // Make it between 0 and 2Pi
  ival = lround(Nphi*phi/(2*M_PI)); // Should get a value between 0..Nphi-1
  if (! test_close(phi, ival*2*M_PI/(double)Nphi) )
    ival = -1; // we failed... *sigh*
  // Sanity check on the range.
  if (ival < 0) ival += Nphi;
  if (ival >= Nphi) ival -= Nphi;
  // Now, just return a failure value :)
  return ((ival == -1) || (jval == -1));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<gauss.J> <grid.IJ>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'e'}; // Would be the flag characters.

const char* ARGEXPL = 
"  gauss.J:  file of gaussian quadrature points and weights\n\
  grid.IJ:  grid points and the function evaluated at those points\n\
  -e: output even l values only (enforce symmetry)";

const char* FILEEXPL =
"==== grid.IJ ====\n\
N d       # Number of points = I*J, dimensionality of vector\n\
l1 m1 n1  f1_1 .. f1_d # Direction cosines (l^2 + m^2 + n^2 = 1)\n\
l2 m2 n2  f2_1 .. f2_d\n\
...\n\
lN mN nN  fN_1 .. fN_d\n\
==== grid.IJ ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = -1;  // Lmax, fed to us by the user.
  int EVEN = 0;

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
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  EVEN = flagon[0];

  // ****************************** INPUT ****************************
  // Command line parameters:
  char *gauss_name, *grid_name;
  // Let's pull off the args:
  gauss_name = args[0];
  grid_name = args[1];

  char dump[512];
  FILE* infile;


  //++ ==== gauss ====
  int Ntheta; //  Number of theta (quadrature) points
  double *xg, *wg; // Gaussian quadrature points and weights.

  infile = myopenr(gauss_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", gauss_name);
    exit(ERROR_NOFILE);
  }

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Ntheta);
  if (Ntheta < 1) {
    fprintf(stderr, "No points listed in %s\n", gauss_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    xg = new double[Ntheta];
    wg = new double[Ntheta];
  }

  for (n=0; (!ERROR) && (n<Ntheta) && (!feof(infile)); ++n) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf", &(xg[n]), &(wg[n]));
  }
  // Make sure we read enough points...
  if (n != Ntheta) {
    fprintf(stderr, "Not enough points in %s\n", gauss_name);
    ERROR = ERROR_BADFILE;
  }
  myclose(infile);
  //-- ==== gauss ====

  if (ERROR) exit(ERROR); // time to bail out yet?
  

  //++ ==== grid.IJ ====

  int Ngrid; // Total number of grid points
  int Nphi; // Number of phi points
  int Ndim; // Dimensionality
  double*** fd_tp; // Our data! fd_tp[d=0..Ndim-1][j=0..Ntheta-1][i=0..Nphi-1]
  
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d", &Ngrid, &Ndim);

  Nphi = Ngrid/Ntheta;
  if (Ngrid < 1) {
    fprintf(stderr, "No points listed in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  if ( (Nphi*Ntheta) != Ngrid) {
    fprintf(stderr, "Ngrid = %d not divisible by Ntheta = %d ; incompatible files?\n",
	    Ngrid, Ntheta);
    ERROR = ERROR_BADFILE;
  }
  if ( Nphi < 2) {
    fprintf(stderr, "Nphi must be at least 2.\n");
    ERROR = ERROR_BADFILE;
  }
  if (Ndim < 1) {
    fprintf(stderr, "No dimensionality parameter in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  if (ERROR) {
    myclose(infile);
    exit(ERROR);
  }

  // Now, start allocating...
  fd_tp = new double **[Ndim];
  for (d=0; d<Ndim; ++d) {
    fd_tp[d] = new double *[Ntheta];
    for (j=0; j<Ntheta; ++j)
      fd_tp[d][j] = new double[Nphi];
  }
  
  int **tp; // Our theta/phi grid checklist :)
  tp = new int *[Ntheta];
  for (j=0; j<Ntheta; ++j) {
    tp[j] = new int[Nphi];
    for (i=0; i<Nphi; ++i) tp[j][i] = 0; // Haven't hit this point yet.
  }
  
  // Now, do the readin'!
  for (n=0; (!ERROR) && (n<Ngrid) && (!feof(infile)); ++n) {
    double lmn[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", &(lmn[0]), &(lmn[1]), &(lmn[2]));
    // Normalize
    double magn;
    magn = sqrt(lmn[0]*lmn[0] + lmn[1]*lmn[1] + lmn[2]*lmn[2]);
    if (dcomp(magn, 0.)) {
      fprintf(stderr, "Point %d on grid (%.5lf $.5lf %.5lf) has zero magnitude.\n", 
	      n+1, lmn[0], lmn[1], lmn[2]);
      ERROR = ERROR_BADFILE;
      break;
    }
    else {
      for (i=0; i<3; ++i) lmn[i] *= 1./magn;
    }
    ERROR = check_ij(lmn, Nphi, Ntheta, xg, i, j);
    if (ERROR) {
      fprintf(stderr, "Couldn't put %d point on grid (%.5lf %.5lf %.5lf)\n",
	      n+1, lmn[0], lmn[1], lmn[2]);
      ERROR = ERROR_BADFILE;
      break;
    }
    ERROR = (tp[j][i] != 0); // See if we've already read this point 
    if (ERROR) {
      fprintf(stderr, "Point %d on grid (%.5lf %.5lf %.5lf) already read?\n",
	      n+1, lmn[0], lmn[1], lmn[2]);
      break;
    }
    // Else, we can read in the matrix elements and put 'em in.
    tp[j][i] = 1;
    // Code to pull off each elem one by one...
    char *startp, *endp;
    startp = dump; // Start at the beginning...
    // Need to "prime the pump" by going past first three entries.
    for (k=0; k<3; ++k) {
      strtod(startp, &endp);
      startp = endp;
    }
    for (d=0; d<Ndim; ++d) {
      fd_tp[d][j][i] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough C_ij's?
    if (d != Ndim) { 
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Didn't read enough entries on line %d.\n", n+1);
    }
  }
  // Make sure we read enough points...
  if (n != Ngrid) {
    fprintf(stderr, "Not enough points in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  myclose(infile);
  // Garbage collect tp:
  for (j=0; j<Ntheta; ++j) delete[] tp[j];
  delete[] tp;
  //-- ==== grid.IJ ====

  // ***************************** ANALYSIS **************************
  // We do our spherical harmonic transform in two steps: first
  // FFT phi_i, the gaussian quadrature on cos(theta_j).  We need to
  // calc our Lmax as well. (Maybe we'll use -m for this?)

  int Lmax;
  Lmax = Nphi/2 - 1;
  // Did our user give us a "reasonable" value?
  if ( (MEMORY >= 0) && (MEMORY < Lmax) )
    Lmax = MEMORY;

  //++ ==== STEP 1 ====
  
  // allocate our FFT elements, and copy data in.
  double ***fd_tm;
  fd_tm = new double **[Ndim];
  for (d=0; d<Ndim; ++d) {
    fd_tm[d] = new double *[Ntheta];
    for (j=0; j<Ntheta; ++j) {
      fd_tm[d][j] = new double[Nphi*2]; // packed array!
      // Go ahead and copy fd_tp into it...
      for (i=0; i<Nphi; ++i) {
	REAL(fd_tm[d][j], i) = fd_tp[d][j][i];
	IMAG(fd_tm[d][j], i) = 0.;
      }
    }
  }
  // Now, we'll use the GSL routines to do our FFT's on phi.
  // OSC: older version of GSL FFT routines; no workspace, etc.
#ifndef OSC
  gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(Nphi);
  gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(Nphi);
#else
  gsl_fft_wavetable_complex *wavetable = gsl_fft_complex_alloc(Nphi); 
#endif

  if (TESTING) {
    printf("## FFT info:\n");
    for (size_t iwav=0; iwav<wavetable->nf; ++iwav) {
      printf ("## factor %d: %d\n", (int)iwav, (int)(wavetable->factor[iwav]));
    }
  }
  
  double Nphi_scale = 2*M_PI/(double)Nphi; // Include 2Pi factor here.
  for (d=0; d<Ndim; ++d) {
    for (j=0; j<Ntheta; ++j) {
#ifndef OSC
      gsl_fft_complex_forward (fd_tm[d][j], 1, Nphi, 
			       wavetable, workspace);
#else
      gsl_fft_complex_forward (fd_tm[d][j], 1, Nphi, wavetable);
#endif
      // Now, scale all by Nphi:
      for (i=0; i<(2*Nphi); ++i)
	fd_tm[d][j][i] *= Nphi_scale;
    }
  }
  
  // Garbage collection...
#ifndef OSC
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
#else
  gsl_fft_complex_free(wavetable);
#endif

  //-- ==== STEP 1 ====
  
  //++ ==== STEP 2 ====
  // Now, to do our gaussian quadrature.
  double ***fd_lm; // Our lm elements.
  int l, m;
  int Mmax = 2*Lmax+1;
  // Note: to get the m^th element in fd_lm, add l value to m.
  //       so [0] -> -m, [1] -> -m+1, ... [l] -> 0, ... [2*l] -> m
  // BUT we have to access using REAL and IMAG, because of the packing.
  fd_lm = new double**[Ndim];
  for (d=0; d<Ndim; ++d) {
    fd_lm[d] = new double *[Lmax+1];
    for (l=0; l<=Lmax; ++l) {
      fd_lm[d][l] = new double[Mmax*2];
    }
  }

  // Go through and calculate our spherical harmonic elements at x_j:
  double ***sphP;
  // Allocate first...
  sphP = new double **[Lmax+1];
  for (l=0; l<=Lmax; ++l) {
    sphP[l] = new double *[l+1]; // Only calc for m=0..l
    for (m=0; m<=l; ++m) {
      sphP[l][m] = new double[Ntheta];
    }
  }
  // Now fill in the elements: (use gsl_sf_legendre_sphPlm_array)
  double* result_temp;
  result_temp = new double[Lmax+1];
  for (m=0; m<=Lmax; ++m) {
    for (j=0; j<Ntheta; ++j) {
      gsl_sf_legendre_sphPlm_array(Lmax, m, xg[j], result_temp);
      // Now, put in the appropriate spots:
      for (l=m; l<=Lmax; ++l)
	sphP[l][m][j] = result_temp[l-m]; // result_temp goes from m to Lmax
    }
  }
  delete[] result_temp; // Garbage collection...
  
  // Now we're ready to do our quadrature:
  int m_index, m_j;
  for (d=0; d<Ndim; ++d) {
    for (l=0; l<=Lmax; ++l) {
      for (m=-l; m<=l; ++m) {
	m_index = m+l;
	if (m>=0)
	  m_j = m;
	else
	  m_j = Nphi + m; // negative indices go at the end of our FFT.
	REAL(fd_lm[d][l], m_index) = 0.;
	IMAG(fd_lm[d][l], m_index) = 0.;
	for (j=0; j<Ntheta; ++j) {
	  REAL(fd_lm[d][l], m_index) +=
	    wg[j]*REAL(fd_tm[d][j], m_j)*sphP[l][abs(m)][j];
	  IMAG(fd_lm[d][l], m_index) +=
	    wg[j]*IMAG(fd_tm[d][j], m_j)*sphP[l][abs(m)][j];
	}
      }
    }
  }

  // Garbage collect sphP:
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=l; ++m) 
      delete[] sphP[l][m];
    delete[] sphP[l];
  }
  delete[] sphP;
  //-- ==== STEP 2 ====

  // Final step... we're gonna only want to output values that are
  // above zero; so we'll find out the largest magnitude component
  // for each dimension:

  double zmax=0, zval;
  for (l=0; l<=Lmax; ++l)
    for (m_index=0; m_index<=(2*l); ++m_index) 
      for (d=0; d<Ndim; ++d) {
	zval = zmagn(REAL(fd_lm[d][l], m_index), IMAG(fd_lm[d][l], m_index));
	if (zval > zmax) zmax = zval;
      }
  
  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.

  printf("%d %d # Lmax, dim;  l m real(flm_1) .. real(flm_d), then imag\n",
	 Lmax, Ndim);
  for (l=0; l<=Lmax; ++l) {
    if (EVEN && ((l%2)==1)) continue;
    for (m=-l; m<=l; ++m) {
      m_index = m+l;
      // Check to see that we've got something non-zero:
      // We'll need to do this condition *correctly* at some point;
      // we should be comparing to the largest expansion element, or
      // something...

      int is_zero;
      for (is_zero=-1, d=0; (d<Ndim) && (is_zero); ++d)
	is_zero = (zmagn(REAL(fd_lm[d][l], m_index), 
			 IMAG(fd_lm[d][l], m_index)) <= FFT_TOLER*zmax);
      if (is_zero) continue;

      // Real part:
      printf("%3d %3d", l, m);
      for (d=0; d<Ndim; ++d)
	printf(" %.12le", REAL(fd_lm[d][l], m_index));
      printf("\n");
      // Imag part:
      printf("%3d %3d", l, m);
      for (d=0; d<Ndim; ++d)
	printf(" %.12le", IMAG(fd_lm[d][l], m_index));
      printf("\n");
    }
  }


  // ************************* GARBAGE COLLECTION ********************

  // Our spherical harmonic elements:
  for (d=0; d<Ndim; ++d) {
    for (l=0; l<=Lmax; ++l)
      delete[] fd_lm[d][l];
    delete[] fd_lm[d];
  }
  delete[] fd_lm;
  
  // Our FFT elements
  for (d=0; d<Ndim; ++d) {
    for (j=0; j<Ntheta; ++j) delete[] fd_tm[d][j];
    delete[] fd_tm[d];
  }
  delete[] fd_tm;

  // Our input data.
  for (d=0; d<Ndim; ++d) {
    for (j=0; j<Ntheta; ++j)
      delete[] fd_tp[d][j];
    delete[] fd_tp[d];
  }
  delete[] fd_tp;

  // Gaussian quad points:
  delete[] xg;
  delete[] wg;

  return 0;
}


  /* 
  // Our sanity test -- reoutput everything.
  printf("%d %d\n", Ntheta*Nphi, Ndim);
  for (j=0; j<Ntheta; ++j) {
    double cost, sint;
    cost = xg[j];
    sint = sqrt(1-cost*cost);
    for (i=0; i<Nphi; ++i) {
      double phi, x, y;
      phi = 2*M_PI*i/(double)Nphi;
      x = sint*cos(phi);
      y = sint*sin(phi);
      printf("%15.12lf %15.12lf %15.12lf ", x, y, cost);
      for (d=0; d<Ndim; ++d)
	printf(" %.8le", fd_tp[d][j][i]);
      printf("\n");
    }
  }
  */

