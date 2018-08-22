/*
  Program: rotlatt.C
  Author:  D. Trinkle
  Date:    August 24, 2004
  Purpose: Utility code to read in a lattice function and rotate it.

  Param.:  <latt(R)> <rot>
           latt(R): lattice function
	   rot:     rotation matrix

           ==== latt(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij
           ...
           nN.1 nN.2 nN.3  Lxx .. Lzz
           ==== latt(R) ====


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
#include "matrix.H"

//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<latt(R)> <rot>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  latt(R): lattice function\n\
  rot:     rotation matrix";

const char* FILEEXPL =
"==== latt(R) ====\n\
N kmax                      # number of points, kmax\n\
n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij\n\
...\n\
nN.1 nN.2 nN.3  Lxx .. lzz\n\
==== latt(R) ====\n\
\n\
rotation matrix: Rxx Rxy Rxz ... Rzz\n\
\n\
We calculate R*L*R^T\n";

int main ( int argc, char **argv ) 
{
  int i, d, np; // General counting variables.

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters: 
  char *latt_name, *rot_name;
  
  // Let's pull off the args:
  latt_name = args[0];
  rot_name = args[1];

  //++ ==== latt(R) ====
  int Np, Ndim = 9;
  point_type *p; // set of points

  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  double cart[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // really just a dummy
       
  ERROR = read_Dij (infile, cart, Np, p, READ_DIJ_MAT);
  myclose(infile);
  if (ERROR) {delete p; exit(ERROR);}
  //-- ==== latt(R) ====

  //++ ==== rot ====
  double rot[9];
  infile = myopenr(rot_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", rot_name);
    exit(ERROR_NOFILE);
  }
  // Read the matrix
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &(rot[0]), &(rot[1]), &(rot[2]), 
         &(rot[3]), &(rot[4]), &(rot[5]), 
         &(rot[6]), &(rot[7]), &(rot[8]));
  myclose(infile);
  //-- ==== rot ====

  // ***************************** ANALYSIS **************************
  double rot_trans[9];
  transpose(rot, rot_trans);
  for (np=0; np<Np; ++np)
    mult(rot, p[np].mat, rot_trans, p[np].mat);
  
  // ****************************** OUTPUT ***************************
  write_Dij(stdout, Np, p);
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;

  return 0;
}
