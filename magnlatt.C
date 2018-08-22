/*
  Program: magnlatt.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Utility code to read in a lattice function and output magnitudes

  Param.:  <cell> <latt)>
           cell: the usual.
           latt: lattice functions

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and calculate.

  Output:  We output the points in unit coord.  (Switch for lattice coord?)
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "cell.H"
#include "Dij.H"


//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <latt>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell: file describing our unit cell\n\
  latt: lattice function";

const char* FILEEXPL =
"==== latt(R) ====\n\
N kmax                      # number of points, kmax\n\
n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij\n\
...\n\
nN.1 nN.2 nN.3  Lxx .. lzz\n\
==== latt(R) ====\n";

int main ( int argc, char **argv ) 
{
  int i, d, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used

  char* args[argc];
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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters: 
  char *cell_name, *latt_name;
  
  cell_name = args[0];
  latt_name = args[1];

  //++ ==== cell ====
  double cart[9], cart2[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input -- we read it, but don't use it.
  int Natoms = NO_ATOMS;
  double** u_atoms = NULL;
  double atomic_mass;

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


  //++ ==== latt(R) ====
  int Np, Ndim = 9;
  point_type *p; // set of points

  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij (infile, cart, Np, p, READ_DIJ_MAT);
  myclose(infile);
  
  if (ERROR) {
    delete p;
    exit(ERROR);
  }
  //-- ==== latt(R) ====

  // ***************************** ANALYSIS **************************

  // ****************************** OUTPUT ***************************
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    double sum;
    for(sum=0, d=0; d<9; ++d) sum += tp->mat[d]*tp->mat[d];
    sum = sqrt(sum/9);
    printf("%.8le %.8le\n", tp->Rmagn, sum);
  }

  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);
  
  return ERROR;
}
