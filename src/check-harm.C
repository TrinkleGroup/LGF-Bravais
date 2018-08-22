/*
  Program: check-harm.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Reads in a Ylm expansion and a grid of points where it's been
           evaluated and computes the maximum error and average error on
	   the grid of points.

  Param.:  <f(lm)> <f(k)>
	   f(lm):  spherical harmonic expansion f(k)
	   f(k):   function evaluated at a set of k-points

	   NOTE: currently only supported even expansions.

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our points and our spherical harmonic expansion

  Output:  
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_eigen.h>  // for determining stability
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "sphere-harm.H"

//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<f(lm)> <f(k)>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  f(lm):  spherical harmonic expansion f(k)\n\
  f(k):   function evaluated at a set of k-points";

const char* FILEEXPL =
"==== f(lm) ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...) \n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
# NOTE: only EVEN l values should appear!!\n\
==== f(lm) ====\n\
\n\
==== f(k) ====\n\
N d                         # number of points, dimensionality\n\
k1.x k1.y k1.z  f1 .. fd    # cart. coord, f.d\n\
...\n\
kN.x kN.y kN.z  f1 .. fd\n\
==== Dij(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (i=0; i<NFLAGS; ++i) flagon[i] = 0;
  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }
  
  // flags:

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *Ylm_name, *grid_name;
  
  // Let's pull off the args:
  Ylm_name = args[0];
  grid_name = args[1];

  //++ ==== f(lm) ====
  // NOTE: this is different than other versions, since we're *forced*
  // to only have EVEN l values.  So [l] is really 2*l.
  int Lmax;  // Maximum l value / 2
  int Ndim;  // Dimensionality
  double ***RYlm, ***IYlm; // Separate our real and imag. components

  infile = myopenr(Ylm_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Ylm_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Ylm_even(infile, Lmax, Ndim, RYlm, IYlm);
  myclose(infile);

  if ( ERROR ) {
    fprintf(stderr, "Bad Lmax (%d) value in %s\n", 
            Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== f(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm, IYlm, "##");
  }

  //++ ==== grid ====
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }

  int Ngrid, Ndimg;
  double** lmn_list, **mat;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d", &Ngrid, &Ndimg);
  if (Ngrid <= 0) {
    fprintf(stderr, "No points listed in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  if (Ndimg != Ndim) {
    fprintf(stderr, "Ndim= %d in %s but %d in %s\n", Ndimg, grid_name,
	    Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  if (!ERROR) {
    lmn_list = new double* [Ngrid];
    mat = new double* [Ngrid];
    for (n=0; (!ERROR) && (n<Ngrid) && (!feof(infile)); ++n) {
      double magn;
      lmn_list[n] = new double[3];
      mat[n] = new double[Ndim];
      fgets(dump, sizeof(dump), infile);
      sscanf(dump, "%lf %lf %lf", 
	     &(lmn_list[n][0]), &(lmn_list[n][1]), &(lmn_list[n][2]) );
      char *startp, *endp;
      startp = dump; // Start at the beginning...
      // Need to "prime the pump" by going past first three entries.
      for (k=0; k<3; ++k) {
	strtod(startp, &endp); // three doubles to pass up...
	startp = endp;
      }
      for (d=0; d<Ndim; ++d) {
	mat[n][d] = strtod(startp, &endp);
	// check for conversion:
	if (startp == endp) break;
	startp = endp;
      }
      // TEST: Did we read enough entries?
      if (d != Ndim) { 
	ERROR = ERROR_BADFILE;
	fprintf(stderr, "Didn't read enough func values in %s at entry %d.\n",
		grid_name, n+1);
      }

      // Normalize.
      magn = sqrt(dot(lmn_list[n],lmn_list[n]));
      if (dcomp(magn, 0.)) {
	fprintf(stderr, 
		"%s at line %d: %.5le %.5le %.5le with magn = %.5le\n",
		grid_name, i+1,
		lmn_list[0][0], lmn_list[0][1], lmn_list[0][2], magn);
	ERROR = ERROR_BADFILE;
      }
      else {
	magn = 1./magn;
	for (i=0; i<3; ++i) lmn_list[n][i] *= magn;
      }
    }
    // Make sure we read enough points...
    if (n != Ngrid) ERROR = ERROR_BADFILE;
  }
  myclose(infile);
  //-- ==== grid ====
  if (ERROR) {
    fprintf(stderr, "Problem with %s\n", grid_name);
    exit(ERROR);
  }


  // ***************************** ANALYSIS **************************
  // Spherical harmonics
  SH_workspace_type * work = allocate_SH_workspace(Lmax);
  double thiserr, maxerr, err;
  int worst;
  
  maxerr=0;
  err=0;
  worst=0;
  
  for (n=0; n<Ngrid; ++n) {
    double *lmn = lmn_list[n];
    double funcYlm[Ndim];
    eval_Ylm_expansion_R(Lmax, lmn, Ndim, RYlm, IYlm, funcYlm, work);

    for (thiserr=0, d=0; d<Ndim; ++d)
      thiserr += fabs(mat[n][d]-funcYlm[d]);
    if (thiserr>maxerr) {
      maxerr=thiserr;
      worst=n;
    }
    err += thiserr;
  }
  err *= 1./Ngrid; // scale by number of points

  // ******************************* OUTPUT **************************
  printf("%.5le  :accumulated average error\n", err);
  printf("%.5le  :maximum error at %11.8lf %11.8lf %11.8lf\n",
	 maxerr, lmn_list[worst][0], lmn_list[worst][1], lmn_list[worst][2]);

  // ************************* GARBAGE COLLECTION ********************
  free_SH_workspace(work);
  free_Ylm(Lmax, RYlm, IYlm);
  for (n=0; n<Ngrid; ++n) {
    delete[] lmn_list[n];
    delete[] mat[n];
  }
  delete[] lmn_list;
  delete[] mat;

  return 0;
}
