/*
  Program: make-ball.C
  Author:  D. Trinkle
  Date:    February 9, 2004
  Purpose: Make a ball of atoms out to a given cutoff.

  Param.:  <cell> <Rcut>
           cell:   cell file describing our lattice
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

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and pass off to ball.H routines.

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
#include "ball.H" 


//***************************** SUBROUTINES ****************************


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <Rcut>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'c'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
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
==== cell ====\n";

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
  char *cell_name;
  double Rcut;
  
  // Let's pull off the args:
  cell_name = args[0];
  sscanf(args[1], "%lf", &Rcut);

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


  // ***************************** ANALYSIS **************************
  int **Runit;
  int Np;

  // This is where the magic happens:
  atom_sphere_unit(cart, Rcut, Np, Runit);

  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  printf("%d  # Number of lattice points\n", Np);
  for (n=0; n<Np; ++n) {
    printf("%3d %3d %3d", Runit[n][0], Runit[n][1], Runit[n][2]);
    if (CARTOUT) {
      double Rvect[3];
      mult_vect(cart, Runit[n], Rvect);
      printf(" %18.12lf %18.12lf %18.12lf  %18.12lf\n", 
	     Rvect[0], Rvect[1], Rvect[2], sqrt(dot(Rvect, Rvect)));
    }
    else
      printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Np; ++n)
    delete[] Runit[n];
  delete[] Runit;

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
