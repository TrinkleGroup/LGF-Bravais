/*
  Program: posfor2cont.C
  Author:  D. Trinkle
  Date:    October 9, 2004
  Purpose: Read in:
           1. Set of cartesian positions
	   2. POSCAR file
	   
	   and we output a brand-spanking new POSCAR file that had the
	   same header, and uses the same selective dynamics.

  Param.:  <pos-for> <contcar>
	   pos-for: atomic positions and forces in cartesian coord.
           contcar: original CONTCAR/POSCAR file

           ==== pos-for ====
           N [t1 t2 t3] # number of points, optional threading direction
           r1.x r1.y r1.z f1.x f1.y f1.z  # position, force in cart. coord
           ...
           rN.x rN.y rN.z fN.x fN.y fN.z
           ==== pos-for ====

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   

  Output:  
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

//***************************** SUBROUTINES ****************************

// structure to hold atomic positions, forces, and displacements:
typedef struct 
{
  double R[3], f[3], u[3];
} posfor_type;

int read_posfor (FILE* infile, int& Np, posfor_type* &p,
		 int t_unit[3], int THREAD, int READ_FORCES);


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<pos-for> <contcar>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  pos-for: atomic positions in cartesian coord\n\
  contcar: VASP position input file";

const char* FILEEXPL =
"==== pos-for ====\n\
N [t1 t2 t3] # number of points, optional threading direction\n\
r1.x r1.y r1.z # position, force in cart. coord\n\
...\n\
rN.x rN.y rN.z\n\
==== pos-for ====\n";

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
  char* posfor_name = args[0];
  char* contcar_name = args[1];


  //++ ==== posfor ====
  int Np, t_unit[3]; // optional threading direction
  posfor_type* p;

  infile = myopenr(posfor_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", posfor_name);
    exit(ERROR_NOFILE);
  }
  // don't bother reading forces if we're supposed to generate forces
  ERROR = read_posfor(infile, Np, p, t_unit, 0, 0);
  myclose(infile);
  //-- ==== posfor ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", posfor_name);
    exit(ERROR);
  }

  //++ ==== contcar ====
  double super[9], a0;
  
  infile = myopenr(contcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", contcar_name);
    exit(ERROR_NOFILE);
  }
  
  // comment line
  nextnoncomment(infile, dump, sizeof(dump));
  printf("%s", dump);
  // a0
  nextnoncomment(infile, dump, sizeof(dump));
  printf("%s", dump);
  sscanf(dump, "%lf", &a0);
  for (d=0; d<3; ++d) {
    // lattice:
    nextnoncomment(infile, dump, sizeof(dump));
    printf("%s", dump);
    sscanf(dump, "%lf %lf %lf", super+d, super+d+3, super+d+6);
  }
  mult(super, a0, super);
  // Number of atoms:
  nextnoncomment(infile, dump, sizeof(dump));
  printf("%s", dump);
  // go till we get the atoms:
  nextnoncomment(infile, dump, sizeof(dump));
  printf("%s", dump);
  if ( (dump[0] == 's') || (dump[0] == 'S') ) {
    nextnoncomment(infile, dump, sizeof(dump));
    printf("%s", dump);
  }
  // now, determine if we have cartesian or direct coord:
  int DIRECT = ( (dump[0] == 'D') || (dump[0] == 'd') );
  double super_inv[9];
  char* tags;
  careful_inverse(super, super_inv);
  for (n=0; n<Np; ++n) {
    posfor_type *tp = p + n;
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) break;
    // now, skip past the first three digits:
    char* endp;
    for (tags=dump, d=0; d<3; ++d, tags=endp) strtod(tags, &endp);
    if (DIRECT) {
      mult_vect(super_inv, tp->R, tp->u);
      for (d=0; d<3; ++d) {
	if (tp->u[d] < 0.) tp->u[d] += 1.;
	if (tp->u[d] >= 1.) tp->u[d] -= 1.;
      }
      printf("%18.12lf %18.12lf %18.12lf %s", 
	     tp->u[0], tp->u[1], tp->u[2], tags);
    } else {
      printf("%18.12lf %18.12lf %18.12lf %s", 
	     tp->R[0], tp->R[1], tp->R[2], tags);
    }
  }

  // ***************************** ANALYSIS **************************
  
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
