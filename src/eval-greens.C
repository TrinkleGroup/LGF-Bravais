/*
  Program: eval-greens.C
  Author:  D. Trinkle
  Date:    February 9, 2004
  Purpose: Given an elastic greens function (calculated by elastic-greens)
           defined on a triangular spherical grid (from sphere), we
	   wish to linearly interpolate our results for an arbitrary
	   point on a sphere.

  Param.:  <GF-grid> <points>
           GF-grid:  triangular grid of direction cosines with elastic GF
	   points:   set of points (direction cosines) to eval GF
	   
	   ==== GF-grid ====
	   N_GF      # Total number of grid points
	   l1 m1 n1  g[xx xy xz  yx yy yz  zx zy zz]
	   ...
	   lN mN nN  g[...]
	   ==== GF-grid ====

	   ==== points ====
	   N         # Number of points
	   l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)
	   l2 m2 n2
	   ...
	   lN mN nN
	   ==== points ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   We read in our GF-grid output from elastic-greens.  The *only*
           sanity check we do on this file is to make sure that N_GF
	   corresponds to a triangular grid (i.e., there exists an int
	   N such that N_GF = (N+1)(N+2)/2).  We are *assuming* that
	   the grid goes along x slowest, then from y to z.  This means
	   our entries goes as (N+1)i - i(i-1)/2 + j.

	   In the future it will probably behoove us to check this and
	   possibly rotate the grid as appropriate.

	   We read in each point in our points file; for each, we do
	   linear interpolation on our triangular grid.  This is somewhat
	   non-trivial, so it will be explained:

	   For a point, we calculate l2, m2, and n2 (the direction cosines
	   squared).  We also check that they're normalized.  Now, given
	   N (the triangular grid param), we calc i0 = floor(N l2),
	   j0 = floor(N m2), and k0 = floor(N n2).  There are three
	   possibilities:

	   1. i0+j0+k0 = N.  In this rare case, we're bang on a grid
	      point, so use it.

	   2. i0+j0+k0 = N-1.  We're in an "up" triangle:
	      define di = N l2 - i0, and so on.  Then
	      f = (100)+(010)+(001) + di((100)-(001)) + dj((010)-(100))
	                            + dk((001)-(010))
	      where (001) = f(i0,j0,k0+1), and so on.
				   
           3. i0+j0+k0 = N-2.  We're in a "down" triangle.
	      Shift i0=i0+1, etc.  Define di = i0 - N l2, and so on.  Then
	      f = (100)+(010)+(001) + di((100)-(001)) + dj((010)-(100))
	                            + dk((001)-(010))
	      where (001) = f(i0,j0,k0-1), and so on.
	   
	   The final point to remember is that to access the (i,j,k)
	   entry in our grid, go to (N+1)i - i(i-1)/2 + j.

  Output:  I figure for now, I'll just output g_ij for each direction
	   point, and that should suffice.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"


/*================================= main ==================================*/

inline double dot(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


// Our equation for the entry index:
inline int entry (int N, int i, int j) 
{
  return (N+1)*i - (i*i-i)/2 + j;
}


// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<GF-grid> <points>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  GF-grid:  triangular grid of direction cosines with elastic GF\n\
  points:   set of points (direction cosines) to eval GF";

const char* FILEEXPL =
"==== GF-grid ====\n\
N_GF      # Total number of grid points\n\
l1 m1 n1  g[xx xy xz  yx yy yz  zx zy zz]\n\
...\n\
lN mN nN  g[...]\n\
==== GF-grid ====\n\
\n\
==== points ====\n\
N         # Number of points\n\
l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)\n\
l2 m2 n2\n\
...\n\
lN mN nN\n\
==== points ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 16384; // 2^14, default (gets turned into Nsteps for int.)

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

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

  int Nsteps;
  
  // ****************************** INPUT ****************************
  // Command line parameters:
  char *grid_name, *infile_name;
  // Let's pull off the args:
  grid_name = args[0];
  infile_name = args[1];

  char dump[512];
  FILE* infile;

  //++ ==== GF-grid ====
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }
  // Read in our grid file
  int Ngrid, Ntria;
  double** grid_list;
  double** Gij_list;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Ngrid);
  // Check that Ngrid makes sense:
  Ntria = (int)(sqrt(2*Ngrid + 0.5) - 1.4); // include a slight rounding...
  if ( (Ntria+1)*(Ntria+2) != (Ngrid*2) ) {
    fprintf(stderr, "Ngrid = %d, not a recognized grid value (tried %d).\n", 
	    Ngrid, Ntria);
    myclose(infile);
    exit(ERROR_BADFILE);
  }
  
  if (Ngrid <= 0) {
    fprintf(stderr, "No points listed in %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    grid_list = new double* [Ngrid]; 
    Gij_list = new double* [Ngrid];
  }
  for (n=0; (!ERROR) && (n<Ngrid) && (!feof(infile)); ++n) {
    grid_list[n] = new double[3];
    Gij_list[n] = new double[9];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	   &(grid_list[n][0]), &(grid_list[n][1]), &(grid_list[n][2]),
	   &(Gij_list[n][0]), &(Gij_list[n][1]), &(Gij_list[n][2]), 
	   &(Gij_list[n][3]), &(Gij_list[n][4]), &(Gij_list[n][5]), 
	   &(Gij_list[n][6]), &(Gij_list[n][7]), &(Gij_list[n][8]));
    // In the future, might add sanity checks here...
  }
  // Make sure we read enough points...
  if (n != Ngrid) ERROR = ERROR_BADFILE;

  myclose(infile);
  //-- ==== GF-grid ====


  //++ ==== points ====
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
  //-- ==== points ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }
  

  // ***************************** ANALYSIS **************************

  // Spit out our first bit of info:
  printf("%d # Number of points: l m n  g[xx xy xz  yx yy yz  zx zy zz]\n",
	 Npoints);
  // For each point...
  for (n=0; n<Npoints; ++n) {
    double *lmn;
    double gint[9]; // Our interpolated Gij
    lmn = lmn_list[n]; // Our specific direction.

    // Now, do our interpolation.
    int i0[3];    // i0, j0, k0
    double di[3]; // di, dj, dk
    double l2; // Holds our direction cosine squared
    
    for (i=0; i<3; ++i) {
      l2 = lmn[i]*lmn[i];
      i0[i] = (int)(Ntria * l2);
      di[i] = Ntria*l2 - i0[i];
    }
    if (TESTING) {
      printf("## Point    %.7lf %.7lf %.7lf\n", lmn[0], lmn[1], lmn[2]);
      printf("## Squared: %.7lf %.7lf %.7lf\n", 
	     lmn[0]*lmn[0], lmn[1]*lmn[1], lmn[2]*lmn[2]);
      printf("## i0 = %d  j0 = %d  k0 = %d\n", i0[0], i0[1], i0[2]);
      printf("## di = %.7lf  dj = %.7lf  dk = %.7lf\n", di[0], di[1], di[2]);
    }
      
    // Now check our three cases:
    int Nsum;
    Nsum = i0[0] + i0[1] + i0[2];
    if (Nsum == Ntria) {
      if (TESTING) printf("## Bang on\n");
      for (i=0; i<9; ++i)
	gint[i] = Gij_list[entry(Ntria, i0[0], i0[1])][i];
    }
    else {
      int index[3]; // indices for 100 (i), 010 (j), 001 (k)
      int base_entry;
      if (Nsum == (Ntria-1)) {
	// Up triangle
	if (TESTING) printf("## Up triangle\n");
	base_entry = entry(Ntria, i0[0], i0[1]);
	index[0] = base_entry+Ntria+1-i0[0];
	index[1] = base_entry+1;
	index[2] = base_entry;
      }
      if (Nsum == (Ntria-2)) {
	// Down triangle
	if (TESTING) printf("## Down triangle\n");
	for (i=0; i<3; ++i) { // Redefine i0 and di:
	  ++(i0[i]);
	  di[i] = 1 - di[i];
	}
	base_entry = entry(Ntria, i0[0], i0[1]);
	index[0] = base_entry-Ntria-2+i0[0];
	index[1] = base_entry-1;
	index[2] = base_entry;
      }
      if (Nsum < (Ntria - 2)) {
	// Scream bloody murder
	fprintf(stderr, "Nsum = %d\n", Nsum);
      }
      if (TESTING) {
	// Output the three points we interpolate
	for (i=0; i<3; ++i) {
	  switch (i) {
	  case 0: printf("## 100 "); break;
	  case 1: printf("## 010 "); break;
	  case 2: printf("## 001 "); break;
	  }
	  printf("(%.5lf %.5lf %.5lf) =",
		 grid_list[index[i]][0], 
		 grid_list[index[i]][1], 
		 grid_list[index[i]][2]);
	  for (j=0; j<9; ++j)
	    printf(" %.3le", Gij_list[index[i]][j]);
	  printf("\n");
	}
      }
      // Now, interpolate!
      for (j=0; j<9; ++j) {
	gint[j] = 0;
	for (i=0; i<3; ++i)
	  gint[j] += di[i]*Gij_list[index[i]][j];
      }
    }

    // ****************************** OUTPUT ***************************
    // Human readable (sorta) first:
    
    if (VERBOSE) {
      // What, exactly, is "verbose" output...?
    }
    // Spit back out the lmn point, and include the gij matrix values:
    printf("%15.12lf %15.12lf %15.12lf", lmn[0], lmn[1], lmn[2]);
    for (i=0; i<3; ++i) {
      printf(" ");
      for (j=0; j<3; ++j) printf(" %.8le", gint[3*i+j]);
    }
    printf("\n");
  }
  
  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Ngrid; ++n) {
    delete[] Gij_list[n];
    delete[] grid_list[n];
  }
  delete[] Gij_list;
  delete[] grid_list;
    
  for (n=0; n<Npoints; ++n) 
    delete[] lmn_list[n];
  delete[] lmn_list;

  return 0;
}
