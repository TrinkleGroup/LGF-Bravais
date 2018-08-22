/*
  Program: harmonic-inter.C
  Author:  D. Trinkle
  Date:    February 12, 2004
  Purpose: Given a spherical harmonic expansion for a vector function,
           evaluate the function at a set of points on the unit sphere.

  Param.:  <Ylm> <grid>
           Ylm:  file of spherical harmonic components
           grid: grid points to evaluate function

	   ==== Ylm ====
	   Lmax d       # Maximum L value, dimensionality of vector
	   l1 m1 R(...) # l,m pair -- real components
	   l1 m1 I(...) # l,m pair -- imag components
	   l2 m2 R(...)
	   l2 m2 I(...)
	   ...
	   # Read until EOF or l=-1.
	   ==== Ylm ====

  Flags:   MEMORY:  not used
	   VERBOSE: print out information about size of components
	   TESTING: repeat what Ylm values were read in

  Algo.:   Given our lm spherical harmonic components, we calc our
           function at each grid point given.

	   If a (l,m) component pair is *not* given, it is assumed to
	   be zero.  However--we must *always* give both real and imag.
	   parts, even if the imag part is all zero.

  Output:  Just the function values at each point (doubled due to real
           and imaginary parts)
           ==== output ====
           N d  # Number of grid points, dimensionality
	   l1 m1 n1  R(f1_1) .. R(f1_d)
	   l1 m1 n1  I(f1_1) .. I(f1_d)
	   ...
	   lN mN nN  R(fN_1) .. R(fN_d)
	   lN mN nN  I(fN_1) .. I(fN_d)
           ==== output ====
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include "io.H"   // All of our "read in file", etc.

//****************************** SUBROUTINES ****************************

inline void calc_angle(double lmn[3], double &cost, double &phi) 
{
  cost = lmn[2];
  phi = atan2(lmn[1], lmn[0]);
}
 

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<Ylm> <grid>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'r','i'}; // Would be the flag characters.

const char* ARGEXPL = 
"  Ylm:  file of spherical harmonic components\n\
  grid: grid points to evaluate function\n\
\n\
  -r:  output only real components\n\
  -i:  output only imaginary\n\
  (note: -r -i outputs both, which is the default)";

const char* FILEEXPL =
"==== Ylm ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...)\n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
==== Ylm ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used
  int REALOUT = 1;  // Which do we want to output?
  int IMAGOUT = 1;

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

  // Flags--note: this is backwards from normally how we do it.
  // Here, if -r is chosen, we set IMAGOUT to 0, and vice versa
  // for -i.  Then, if *both* are set, we turn 'em back on :)
  if (flagon[0]) IMAGOUT = 0;
  if (flagon[1]) REALOUT = 0;
  if (IMAGOUT == REALOUT) IMAGOUT = REALOUT = 1;

  // ****************************** INPUT ****************************
  // Command line parameters:
  char *Ylm_name, *grid_name;
  // Let's pull off the args:
  Ylm_name = args[0];
  grid_name = args[1];

  char dump[512];
  FILE* infile;


  //++ ==== Ylm ====
  int Lmax;  // Maximum l value
  int Ndim;  // Dimensionality
  double ***RYlm, ***IYlm; // Separate our real and imag. components

  infile = myopenr(Ylm_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Ylm_name);
    exit(ERROR_NOFILE);
  }

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d", &Lmax, &Ndim);
  if ( (Lmax < 0) || (Ndim < 1) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
	    Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    // Start allocating...
    RYlm = new double **[Lmax+1];
    IYlm = new double **[Lmax+1];
    for (l=0; l<=Lmax; ++l) {
      RYlm[l] = new double *[2*l+1];
      IYlm[l] = new double *[2*l+1];
      for (m=0; m<=(2*l); ++m) {
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
  int m_index, l_max_0;
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
    if ( (l > Lmax) || (abs(m) > l) ) {
      fprintf(stderr, "Entry %d has bad lm components = %d %d\n", n+1, l, m);
      ERROR = ERROR_BADFILE;
      break;
    }
    if (l_max_0 < l) l_max_0 = l;
    m_index = l+m;
    
    // Code to pull off each elem one by one...
    char *startp, *endp;
    startp = dump; // Start at the beginning...
    // Need to "prime the pump" by going past first two entries.
    for (k=0; k<2; ++k) {
      strtol(startp, &endp, 10); // use l because we've got two integer elems.
      startp = endp;
    }
    for (d=0; d<Ndim; ++d) {
      RYlm[l][m_index][d] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough C_ij's?
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
      IYlm[l][m_index][d] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough C_ij's?
    if (d != Ndim) { 
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Didn't read enough entries on imag entry %d.\n", n+1);
    }
  }
  myclose(infile);
  //-- ==== Ylm ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  // Now, remove the "extra" l elements, if l_max_0 < Lmax.
  if (l_max_0 < Lmax) {
    // Free up the space:
    for (l=(l_max_0+1); l < Lmax; ++l) {
      for (m=0; m<=(2*l); ++m) {
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
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(2*l); ++m) {
	// Real part
	printf("## %3d %3d R", l, m-l);
	for (d=0; d<Ndim; ++d)
	  printf(" %11.4le", RYlm[l][m][d]);
	printf("\n");
	// Imag part
	printf("## %3d %3d I", l, m-l);
	for (d=0; d<Ndim; ++d)
	  printf(" %11.4le", IYlm[l][m][d]);
	printf("\n");
      }
    }
  }

  // VERBOSE: print out information about *size* of each component.
  if (VERBOSE) {
    double zmagn;
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(2*l); ++m) {
	zmagn = 0.;
	for (d=0; d<Ndim; ++d)
	  zmagn += RYlm[l][m][d]*RYlm[l][m][d] + 
	    IYlm[l][m][d]*IYlm[l][m][d];
	zmagn = sqrt(zmagn /(double)Ndim);
	printf("# %3d %3d %.8le\n", l, m-l, zmagn);
      }
    }
  }

  //++ ==== grid ====

  int Ngrid;   // Total number of grid points
  double** lmn; // Our points
  
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Ngrid);

  if (Ngrid < 1) {
    fprintf(stderr, "No points listed in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    lmn = new double *[Ngrid];
    for (i=0; i<Ngrid; ++i) lmn[i] = new double[3];
  }
			
  if (ERROR) {
    myclose(infile);
    exit(ERROR);
  }
  // Now, do the readin'!
  for (n=0; (!ERROR) && (n<Ngrid) && (!feof(infile)); ++n) {
    double* lmn_p;
    lmn_p = lmn[n];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", &(lmn_p[0]), &(lmn_p[1]), &(lmn_p[2]));
    // Normalize
    double magn;
    magn = sqrt(lmn_p[0]*lmn_p[0] + lmn_p[1]*lmn_p[1] + lmn_p[2]*lmn_p[2]);
    if (dcomp(magn, 0.)) {
      fprintf(stderr, "Point %d on grid (%.5lf $.5lf %.5lf) has zero magnitude.\n", 
	      n+1, lmn_p[0], lmn_p[1], lmn_p[2]);
      ERROR = ERROR_BADFILE;
      break;
    }
    else {
      for (i=0; i<3; ++i) lmn_p[i] *= 1./magn;
    }
  }
  // Make sure we read enough points...
  if (n != Ngrid) {
    fprintf(stderr, "Not enough points in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  myclose(infile);
  //-- ==== grid ====

  // ***************************** ANALYSIS **************************

  // Go through and calculate our spherical harmonic elements at each
  // grid point.
  // Note: RsphP[l][-m] = RsphP[l][m], and IsphP[l][-m] = -IsphP[l][m]
  //       so we only calc for m=0..l
  double ***RsphP, ***IsphP;
  // Allocate first...
  RsphP = new double **[Lmax];
  IsphP = new double **[Lmax];
  for (l=0; l<=Lmax; ++l) {
    RsphP[l] = new double *[l+1]; // Only calc for m=0..l
    IsphP[l] = new double *[l+1];
    for (m=0; m<=l; ++m) {
      RsphP[l][m] = new double[Ngrid];
      IsphP[l][m] = new double[Ngrid];
    }
  }

  // Now fill in the elements: (use gsl_sf_legendre_sphPlm_array)
  double* result_temp;
  result_temp = new double[Lmax+1];
  for (m=0; m<=Lmax; ++m) {
    for (n=0; n<Ngrid; ++n) {
      double phi, cosm, sinm;
      double cost; // cos(theta) for lmn
      calc_angle(lmn[n], cost, phi);
      // exp(I*m*phi):
      sinm = sin(m*phi);
      cosm = cos(m*phi);
      // Plm(cost):
      gsl_sf_legendre_sphPlm_array(Lmax, m, cost, result_temp);
      // Now, put in the appropriate spots:
      for (l=m; l<=Lmax; ++l) {
	// result_temp goes from m to Lmax
	RsphP[l][m][n] = cosm*result_temp[l-m];
	IsphP[l][m][n] = sinm*result_temp[l-m];
      }
    }
  }
  delete[] result_temp; // Garbage collection...
  
  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.

  printf("%d %d # Ngrid, dim; x y z  ", Ngrid, Ndim);
  if (REALOUT) printf("real");
  if (IMAGOUT) printf("imag");
  printf("\n");
  int mplus, mneg;
  for (n=0; n<Ngrid; ++n) {
    double Rval[Ndim], Ival[Ndim];
    for (d=0; d<Ndim; ++d) {
      Rval[d] = 0.;
      Ival[d] = 0.;
      for (l=0; l<=Lmax; ++l) {
	// m = 0:
	// note: xYlm[l][m_index], m_index = l+m, so for m=0, we use l.
	m_index = l;
	Rval[d] += RYlm[l][m_index][d]*RsphP[l][0][n];
	Ival[d] += IYlm[l][m_index][d]*RsphP[l][0][n];
	// Next, rest of m terms:
	for (m=1; m<=l; ++m) {
	  mplus = m+l;
	  mneg = -m+l;
	  Rval[d] += (RYlm[l][mplus][d] + RYlm[l][mneg][d])*RsphP[l][m][n]
	    - (IYlm[l][mplus][d] - IYlm[l][mneg][d])*IsphP[l][m][n];
	  Ival[d] += (RYlm[l][mplus][d] - RYlm[l][mneg][d])*IsphP[l][m][n]
	    + (IYlm[l][mplus][d] + IYlm[l][mneg][d])*RsphP[l][m][n];
	}
      }
    }
    // Now, output!
    // Real part:
    if (REALOUT) {
      printf("%15.12lf %15.12lf %15.12lf ", lmn[n][0], lmn[n][1], lmn[n][2]);
      for (d=0; d<Ndim; ++d)
	printf(" %.12le", Rval[d]);
      printf("\n");
    }
    // Imag part:
    if (IMAGOUT) {
      printf("%15.12lf %15.12lf %15.12lf ", lmn[n][0], lmn[n][1], lmn[n][2]);
      for (d=0; d<Ndim; ++d)
	printf(" %.12le", Ival[d]);
      printf("\n");
    }
  }

  // ************************* GARBAGE COLLECTION ********************

  // Garbage collect sphP:
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=l; ++m) {
      delete[] RsphP[l][m];
      delete[] IsphP[l][m];
    }
    delete[] RsphP[l];
    delete[] IsphP[l];
  }
  delete[] RsphP;
  delete[] IsphP;

  // Grid points
  for (i=0; i<Ngrid; ++i) delete[] lmn[i];
  delete[] lmn;

  // spherical harmonic components
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(2*l); ++m) {
      delete[] RYlm[l][m];
      delete[] IYlm[l][m];
    }
    delete[] RYlm[l];
    delete[] IYlm[l];
  }
  delete[] RYlm;
  delete[] IYlm;

  return 0;
}
