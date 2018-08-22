#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"

// Make (and output) a (phi, theta) angular grid in cartesian coordinates

// We do the entire sphere now, and read in our grid on theta from a
// file containing gaussian quadrature points.  Our phi grid is still
// uniform, however.

// Output routine:
// Designed to convert from spherical to cartesian coord:
// x = r*sin(theta)*cos(phi)
// y = r*sin(theta)*sin(phi)
// z = r*cos(theta)
inline void print_rpt(const double r, const double phi, const double theta) 
{
  double x, y, z;
  x = r*sin(theta)*cos(phi);
  y = r*sin(theta)*sin(phi);
  z = r*cos(theta);
  printf("%18.15lf %18.15lf %18.15lf\n", x, y, z);
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<Nphi> <gauss.N>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  Outputs a (theta,phi) grid that's uniform in phi and uses gaussian\n\
  quadrature points for theta = acos(x_i).\n\
\n\
  Nphi:  number of phi points from 0..2*Pi\n\
  gauss.N: points for a gaussian quadrature on [-1,1] (N first, then entries)";


int main (int argc, char **argv) {
  int i, j, n;

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // unused

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
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  int N, Nphi, Ntheta;
  double dphi;
  double phi, theta, r;
  double *xg; // Our gaussian quadrature points
  char* gauss_filename;

  sscanf(argv[1], "%d", &Nphi);
  if (Nphi<=0) Nphi = 1;
  gauss_filename = argv[2];
  
  // Now read in the gaussian quadrature points:
  char dump[512];
  FILE *infile;
  
  infile = myopenr(gauss_filename);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", gauss_filename);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Ntheta);
  if (Ntheta < 1) {
    fprintf(stderr, "Gaussian file %s doesn't contain any points?\n",
	    gauss_filename);
    exit(ERROR_BADFILE);
  }
  // Get them points!
  xg = new double[Ntheta];
  for (n=0; (!ERROR) && (n<Ntheta) && (!feof(infile)); ++n) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf", &(xg[n]));
  }
  if (n != Ntheta) ERROR = ERROR_BADFILE;
  myclose(infile);
  
  if (ERROR) {
    fprintf(stderr, "Not enough points in %s?\n", gauss_filename);
    delete[] xg;
    exit(ERROR);
  }
  
  // Now, do the real work.
  r = 1; // Just a sphere... not very interesting...

  dphi = 2*M_PI/(double)Nphi;  
  printf("%d  # %d x %d grid on (theta, phi); theta major order\n", 
	 Nphi*Ntheta, Ntheta, Nphi);
  for (i=0; i<Ntheta; ++i) {
    theta = acos(xg[i]); // Get theta value...
    for (j=0; j<Nphi; ++j) {
      phi = dphi*j;
      print_rpt(r, phi, theta);
    }
  }

  // Garbage collection:
  delete[] xg;
  
  return 0;
}
