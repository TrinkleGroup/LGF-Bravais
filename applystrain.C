/*
  Program: applystrain.C
  Author:  D. Trinkle
  Date:    October 13, 2004
  Purpose: Read in:
           1. Set of positions (cartesian)
	   2. Strain direction (matrix)
	   3. Strain magnitude -- to scale the matrix

	   ... and we displace the atoms according to the strain matrix.
	   We do make the displacements sum to 0.

  Param.:  <pos-for> <strain matrix> <magn>
           pos-for: atomic positions and forces in cartesian coord.
	   strain:  9 index version of our strain matrix
	   magn:    scale factor for strain matrix

           ==== pos-for ====
           N [t1 t2 t3] # number of points, optional threading direction
           r1.x r1.y r1.z  # position in cart. coord
           ...
           rN.x rN.y rN.z
           ==== pos-for ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   

  Output:  We output the new positions in cartesian coord. from the update.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "Dij.H"
#include "supercell.H" // indexing routines for easy lookup


//***************************** SUBROUTINES ****************************

// structure to hold atomic positions, forces, and displacements:
typedef struct 
{
  double R[3], f[3], u[3];
} posfor_type;

int read_posfor (FILE* infile, int& Np, posfor_type* &p,
		 int t_unit[3], int THREAD, int READ_FORCES);

inline void setmax (int &a, const int &b) 
{if (a<b) a=b;}



/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<pos-for> <strain matrix> <magn>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  pos-for: atomic positions and forces in cartesian coord.\n\
  matrix:  matrix describing the strain; scaled by...\n\
  magn:    scale factor for strain";

const char* FILEEXPL =
"==== pos-for ====\n\
N [t1 t2 t3] # number of points, optional threading direction\n\
r1.x r1.y r1.z f1.x f1.y f1.z  # position, force in cart. coord\n\
...\n\
rN.x rN.y rN.z fN.x fN.y fN.z\n\
==== pos-for ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used

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

  // flags

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char* posfor_name = args[0];
  char* matrix_name = args[1];
  double strain_magn;
  sscanf(args[2], "%lf", &strain_magn);

  //++ ==== posfor ====
  int Np; 
  posfor_type* p;

  infile = myopenr(posfor_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", posfor_name);
    exit(ERROR_NOFILE);
  }
  // no forces or threading direction that we care about
  ERROR = read_posfor(infile, Np, p, NULL, 0, 0);
  myclose(infile);
  //-- ==== posfor ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", posfor_name);
    exit(ERROR);
  }

  //++ ==== matrix ====
  double matrix[9];

  infile = myopenr(matrix_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", matrix_name);
    exit(ERROR_NOFILE);
  }
  nextnoncomment(infile, dump, sizeof(dump));
  ERROR = (sscanf(dump, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		  matrix, matrix+1, matrix+2, matrix+3, matrix+4,
		  matrix+5, matrix+6, matrix+7, matrix+8) != 9);
  myclose(infile);
  //-- ==== LGF ====

  if (ERROR) {
    delete p; 
    fprintf(stderr, "Error with %s\n", matrix_name);
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  mult(matrix, strain_magn, matrix);
  double usum[3] = {0,0,0};
  for (np=0; np<Np; ++np) {
    posfor_type *tp = p+np;
    mult_vect(matrix, tp->R, tp->u);
    for (d=0; d<3; ++d) usum[d] += tp->u[d];
  }
  for (d=0; d<3; ++d) usum[d] *= -1./(double)Np;
  for (np=0; np<Np; ++np) {
    posfor_type *tp = p+np;
    for (d=0; d<3; ++d) tp->u[d] += usum[d];
  }
  
  // ****************************** OUTPUT ***************************
  printf("%d\n", Np);
  for (np=0; np<Np; ++np) {
    posfor_type* tp = p+np;
    printf("%18.12lf %18.12lf %18.12lf\n", tp->R[0]+tp->u[0],
	   tp->R[1]+tp->u[1], tp->R[2]+tp->u[2]);
  }
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;

  return 0;
}


int read_posfor (FILE* infile, int& Np, posfor_type* &p,
		 int t_unit[3], int THREAD, int READ_FORCES) 
{
  int n, i, icount;
  int ERROR = 0;
  char dump[512];

  nextnoncomment(infile, dump, sizeof(dump));
  if (THREAD) {
    i = sscanf(dump, "%d %d %d %d", &Np, t_unit, t_unit+1, t_unit+2);
    if ( (i != 4) || zero(t_unit) ) {
      fprintf(stderr, "Asked for threading, but bad threading direction in read_posfor\n");
      return ERROR_BADFILE;
    }
  }
  else
    sscanf(dump, "%d", &Np);

  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in read_posfor\n", Np);
    return ERROR_BADFILE;
  }
  //    if (p != NULL) delete[] p;
  p = new posfor_type[Np];
  if (p == NULL) {
    fprintf(stderr, "Error allocating memory...?\n");
    return ERROR_MEMORY;
  }

  if (READ_FORCES) icount = 6;
  else             icount = 3;
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    posfor_type *tp = p + n; // this point
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) break;
    // Parse it.
    if (READ_FORCES)
      i = sscanf(dump, "%lf %lf %lf %lf %lf %lf", 
		 &(tp->R[0]), &(tp->R[1]), &(tp->R[2]),
		 &(tp->f[0]), &(tp->f[1]), &(tp->f[2]));
    else
      i = sscanf(dump, "%lf %lf %lf", 
		 &(tp->R[0]), &(tp->R[1]), &(tp->R[2]));
    if (i != icount) ERROR = ERROR_BADFILE;
    for (int d=0; d<3; ++d) tp->u[d] = 0;
    if (! READ_FORCES)
      for (int d=0; d<3; ++d) tp->f[d] = 0;
  }
  if (ERROR) {
    fprintf(stderr, "Bad read line in read_posfor\n");
  }
  if (n!=Np) {
    fprintf(stderr, "Not enough points in read_posfor\n");
    ERROR = ERROR_BADFILE;
  }
  return ERROR;
}
