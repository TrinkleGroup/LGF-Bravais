/*
  Program: make-slab.C
  Author:  D. Trinkle
  Date:    October 6, 2004
  Purpose: Make a slab of atoms out to a given cutoff (in the grand
           tradition of the make-* family of programs...)

  Param.:  <cell> <disl> <Rcut>
           cell:   cell file describing our lattice
	   disl:   dislocation coordinate system
	   Rcut:   cutoff

           ==== cell ====
           a0                            # Scale factor for unit cell
           a1.x a1.y a1.z                # Cartesian coord of unit cell
           a2.x a2.y a2.z
           a3.x a3.y a3.z
           crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
           Natoms                        # Number of atoms in first unit cell
           u1.1 u1.2 u1.3                # Atom locations, in direct coord.
           ...
           uN.1 uN.2 uN.3
           ==== cell ====

           ==== disl ====
           t1.1 t1.2 t1.3  # line direction
           m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)
           n1.1 n1.2 n1.3  # mutual perp. vector
           # all three vectors are given in unit cell coord.
           ==== disl ====

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and pass off to slab.H routines.

  Output:  We output the points in unit coord.  (Switch for lattice coord?)
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
#include "slab.H" 


//***************************** SUBROUTINES ****************************


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <disl> <Rcut>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'c'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  disl:   dislocation coordinate system\n\
  Rcut:   cutoff for sphere\n\
  -c  output in cartesian coord (default is unit cell coord)";

const char* FILEEXPL =
"==== cell ====\n\
a0                            # Scale factor for unit cell\n\
a1.x a1.y a1.z                # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.\n\
Natoms                        # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3                # Atom locations, in direct coord.\n\
...\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== disl ====\n\
t1.1 t1.2 t1.3  # line direction\n\
m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)\n\
n1.1 n1.2 n1.3  # mutual perp. vector\n\
# all three vectors are given in unit cell coord.\n\
==== disl ====\n";

int main ( int argc, char **argv ) 
{
  int i, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used
  int CARTOUT = 0;  // output in cartesian coord?

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
  CARTOUT = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  
  // Let's pull off the args:
  char* cell_name = args[0];
  char* disl_name = args[1];
  double Rcut;
  sscanf(args[2], "%lf", &Rcut);

  if (Rcut <= 0) {
    fprintf(stderr, "Rcut (%.5lf) less than 0?\n", Rcut);
    exit(1);
  }

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input -- we read it, but don't use it.
  int Natoms = NO_ATOMS;
  double** u_atoms = NULL;
  double atomic_mass;

  //++ ==== cell ====
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

  //++ ==== disl ====
  int t_unit[3], m_unit[3], n_unit[3];
  
  infile = myopenr(disl_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", disl_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_disl(infile, t_unit, m_unit, n_unit);
  myclose(infile);
  //-- ==== disl ====
  if (ERROR) {
    fprintf(stderr, "Bad dislocation coordinate system in %s\n", disl_name);
    exit(ERROR);
  }
  if (TESTING) {
    double t_vect[3], m_vect[3], n_vect[3];
    make_dislcoord(cart, t_vect, m_vect, n_vect, t_unit, m_unit, n_unit,
		   DISLCOORD_UNIT + DISLCOORD_NOCHANGE);
    printf("## vect    input              cartesian\n");
    printf("## t = (%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
           t_unit[0],t_unit[1],t_unit[2], t_vect[0],t_vect[1],t_vect[2]);
    printf("## m = (%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
           m_unit[0],m_unit[1],m_unit[2], m_vect[0],m_vect[1],m_vect[2]);
    printf("## n = (%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
           n_unit[0],n_unit[1],n_unit[2], n_vect[0],n_vect[1],n_vect[2]);
  }

  // ***************************** ANALYSIS **************************
  point_type* p;
  int Np;

  // This is where the magic happens:
  ERROR = construct_slab(cart, Rcut, t_unit, m_unit, n_unit, Np, p);

  if (ERROR) {
    fprintf(stderr, "Error encountered in construct_slab.\n");
    exit(ERROR);
  }

  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  printf("%d %3d %3d %3d # Number of points, threading direc. (unit coord)\n",
	 Np, t_unit[0], t_unit[1], t_unit[2]);
  for (n=0; n<Np; ++n) {
    point_type* tp = p+n;
    printf("%3d %3d %3d", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    if (CARTOUT) {
      printf(" %18.12lf %18.12lf %18.12lf  %18.12lf\n", 
	     tp->Rcart[0], tp->Rcart[1], tp->Rcart[2], tp->Rmagn);
    }
    else
      printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  delete[] p;

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
