/*
  Program: func.C
  Author:  D. Trinkle
  Date:    February 9, 2004
  Purpose: Calc. functions on a sphere.

  Param.:  <grid>
           grid:  "grid" of direction cosines for which to compute Gij(R)

	   ==== grid ====
	   N         # Number of points
	   l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)
	   l2 m2 n2
	   ...
	   lN mN nN
	   ==== grid ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: not used

  Algo.:   N/A

  Output:  Function evaluated at each point.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.

//****************************** SUBROUTINES ****************************

inline double dot(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<grid>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  grid:  \"grid\" of direction cosines for which to compute";

const char* FILEEXPL =
"==== grid ====\n\
N         # Number of points\n\
l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)\n\
...\n\
lN mN nN\n\
==== infile ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used

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
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // ****************************** INPUT ****************************
  // Command line parameters:
  char *infile_name;
  // Let's pull off the args:
  infile_name = args[0];

  char dump[512];
  FILE* infile;

  //++ ==== grid ====
  infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", infile_name);
    exit(ERROR_NOFILE);
  }

  int Npoints;
  double** lmn_list;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Npoints);
  if (Npoints <= 0) {
    fprintf(stderr, "No points listed in %s\n", infile_name);
    ERROR = ERROR_BADFILE;
  }
  else lmn_list = new double* [Npoints];
  for (n=0; (!ERROR) && (n<Npoints) && (!feof(infile)); ++n) {
    double magn;
    lmn_list[n] = new double[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", 
	   &(lmn_list[n][0]), &(lmn_list[n][1]), &(lmn_list[n][2]) );
    // Normalize.
    magn = sqrt(dot(lmn_list[n],lmn_list[n]));
    if (dcomp(magn, 0.)) {
      fprintf(stderr, 
	      "%s at line %d: %.5le %.5le %.5le with magn = %.5le\n",
	      infile_name, i+1,
	      lmn_list[0][0], lmn_list[0][1], lmn_list[0][2], magn);
      ERROR = ERROR_BADFILE;
    }
    else {
      magn = 1./magn;
      for (i=0; i<3; ++i) lmn_list[n][i] *= magn;
    }
  }
  // Make sure we read enough points...
  if (n != Npoints) ERROR = ERROR_BADFILE;
  myclose(infile);
  //-- ==== grid ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }
  

  // ***************************** ANALYSIS **************************
  int Ndim = 4;
  double fint[Ndim];
  
  // Spit out our first bit of info:
  printf("%d %d # Number of points and dimensionality\n",
	 Npoints, Ndim);
  // For each point...
  for (n=0; n<Npoints; ++n) {
    double *lmn;
    lmn = lmn_list[n]; // Our specific direction.
    
    // Eval three simple functions:
    fint[0] = 1.;
    fint[1] = lmn[0];
    fint[2] = lmn[1];
    fint[3] = lmn[2];

    // ****************************** OUTPUT ***************************
    // Spit back out the lmn point, and include the gij matrix values:
    printf("%15.12lf %15.12lf %15.12lf", lmn[0], lmn[1], lmn[2]);
    for (i=0; i<Ndim; ++i) printf(" %.8le", fint[i]);
    printf("\n");
  }
  
  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Npoints; ++n)
    delete lmn_list[n];
  delete[] lmn_list;
  
  return 0;
}
