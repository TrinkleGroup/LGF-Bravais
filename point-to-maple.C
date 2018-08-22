#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"

// Make (and output) a (phi, theta) angular grid in maple-friendly
// format for plotting and what-not.

// We assume that what's fed to us is an NxN grid; we read it in,
// and output something maple can handle.  Basically, we will set the
// radius at each point to the Gij value there, and spit out the
// 6 values: Gxx, Gyy, Gzz, Gxy, Gyz, Gzx.

/*================================= main ==================================*/

inline void print_xyzr(const double x, const double y, const double z,
		       const double r) 
{
  printf("   [%.12lf, %.12lf, %.12lf]", x*r, y*r, z*r);
}

// Order of values to print, and their names:
const int Nval = 6;
const char valname[Nval][3] = {"xx", "yy", "zz", "xy", "yz", "zx"};
const int valindex[Nval] = {0, 4, 8, 1, 5, 2};


// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<grid>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  grid:  square angular grid on the sphere with Gij matrix\n";

const char* FILEEXPL =
"==== grid ====\n\
N         # Total number of grid points -- must be square\n\
l1 m1 n1  g[xx xy xz  yx yy yz  zx zy zz]\n\
...\n\
lN mN nN  g[...]\n\
==== grid ====\n";

int main (int argc, char **argv) {
  int i, j, k, n;

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

  // ****************************** INPUT ****************************
  // Command line parameters:
  char *grid_name;
  // Let's pull off the args:
  grid_name = args[0];

  char dump[512];
  FILE* infile;

  //++ ==== GF-grid ====
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }
  // Read in our grid file
  int Ngrid, Ns;
  double** grid_list;
  double** Gij_list;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Ngrid);
  // Check that Ngrid makes sense:
  Ns = (int)(sqrt(Ngrid) + 0.1); // include a slight rounding...
  if ( (Ns*Ns) != (Ngrid) ) {
    fprintf(stderr, "Ngrid = %d, not a recognized grid value (tried %d).\n", 
            Ngrid, Ns);
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

  // Now that we've got all the data, let's pull it apart.
  int val;
  for (val = 0; val < Nval; ++val) {
    printf("printf(\"# read in plot_%s ; %d x %d grid\\n\");\n", 
	   valname[val], Ns, Ns);
    printf("plot_%s:=[\n", valname[val]);
    for (i=0; i<Ns; ++i) {
      // Print out a "comment" line
      printf("  [ # %.5lf %.5lf %.5lf -> %.5lf %.5lf %.5lf\n",
	     grid_list[i*Ns][0], grid_list[i*Ns][1], grid_list[i*Ns][2], 
	     grid_list[(i+1)*Ns-1][0], grid_list[(i+1)*Ns-1][1], 
	     grid_list[(i+1)*Ns-1][2]);
      for (j=0; j<Ns; ++j) {
	n = i*Ns + j;
	print_xyzr(grid_list[n][0], grid_list[n][1], grid_list[n][2], 
		   Gij_list[n][valindex[val]]);
	if (j != (Ns-1))
	  printf(",\n");
	else
	  printf(" ]");
      }
      if (i != (Ns-1))
	printf(",\n");
      else
	printf("\n]:\n"); // colon keeps maple from printing when it reads.
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

  return 0;
}
