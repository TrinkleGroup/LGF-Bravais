/*
  Program: applystrain-contcar.C
  Author:  D. Trinkle
  Date:    October 13, 2004
  Purpose: Read in:
           1. CONTCAR file
	   2. Strain direction (matrix)
	   3. Strain magnitude -- to scale the matrix

	   ... and we modify the supercell according to the strain matrix.

  Param.:  <strain matrix> <magn> <CONTCAR> 
	   strain:  9 index version of our strain matrix
	   magn:    scale factor for strain matrix
           CONTCAR: output from VASP

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   

  Output:  Output new POSCAR file
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
const char* ARGLIST = "<strain matrix> <magn> <CONTCAR>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  matrix:  matrix describing the strain; scaled by...\n\
  magn:    scale factor for strain\n\
  CONTCAR: position file from VASP";

const char* FILEEXPL =
"";

int main ( int argc, char **argv ) 
{
  int d, i, n; // General counting variables.

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
  char* matrix_name = args[0];
  double strain_magn;
  sscanf(args[1], "%lf", &strain_magn);
  char* contcar_name = args[2];

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
  //-- ==== matrix ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", matrix_name);
    exit(ERROR);
  }
  mult(matrix, strain_magn, matrix);

  //++ ==== contcar ====
  int Np;
  double super[9], a0;
  
  infile = myopenr(contcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", contcar_name);
    exit(ERROR_NOFILE);
  }
  
  // comment line
  nextnoncomment(infile, dump, sizeof(dump));
  char* strp;
  strp = strchr(dump, '\n');
  if (strp != NULL) *strp = '\0'; // null terminate that sucka.
  printf("%s %s strain %.9lf\n", dump, matrix_name, strain_magn);
  
  // a0
  nextnoncomment(infile, dump, sizeof(dump));
  printf("%s", dump);
  sscanf(dump, "%lf", &a0);

  for (d=0; d<3; ++d) {
    // lattice:
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%lf %lf %lf", super+d, super+d+3, super+d+6);
  }
  // Now, strain our lattice...
  double str_super[9];
  for (d=0; d<9; ++d) matrix[d] += ident[d];
  mult(matrix, super, str_super);
  for (d=0; d<3; ++d)
    printf("  %20.15lf  %20.15lf  %20.15lf\n", 
	   str_super[d], str_super[d+3], str_super[d+6]);

  // Number of atoms:
  nextnoncomment(infile, dump, sizeof(dump));
  printf("%s", dump);
  char *startp, *endp = dump;
  Np=0;
  do {
    startp = endp;
    Np += strtol(startp, &endp, 10);
  } while (startp != endp);
  if (Np < 1) {
    fprintf(stderr, "Need at least one atom.\n");
    exit(ERROR_BADFILE);
  }
  // go till we get the atoms:
  nextnoncomment(infile, dump, sizeof(dump));
  if ( (dump[0] == 's') || (dump[0] == 'S') ) {
    printf("%s", dump);
    nextnoncomment(infile, dump, sizeof(dump));
  }

  // now, determine if we have cartesian or direct coord:
  int DIRECT = ( (dump[0] == 'D') || (dump[0] == 'd') );
  printf("Direct\n"); // FORCE this, but handle cartesian coord, too.
  double super_inv[9];
  char *tags;
  mult(super, a0, super); // scale by a0
  careful_inverse(super, super_inv);
  for (n=0; n<Np; ++n) {
    double r[3], u[3];
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) break;
    // now, skip past the first three digits:
    char* endp;
    for (tags=dump, d=0; d<3; ++d, tags=endp) strtod(tags, &endp);
    sscanf(dump, "%lf %lf %lf", r, r+1, r+2);
    if (!DIRECT) {
      mult_vect(super_inv, r, u);
      for (d=0; d<3; ++d) {
	if (u[d] >= 1.0) u[d] -= 1.;
	if (u[d] <  0.0) u[d] += 1.;
      }
    }
    else
      for (d=0; d<3; ++d) u[d] = r[d];
    printf("%18.12lf %18.12lf %18.12lf %s", u[0], u[1], u[2], tags);
  }
  myclose(infile);
  ERROR = (n<Np);
  //-- ==== contcar ====
  
  // ************************* GARBAGE COLLECTION ********************

  return 0;
}
