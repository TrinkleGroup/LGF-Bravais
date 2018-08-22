#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.H"
#include "dcomp.H"
#include "cell.H"
#include "io.H"
#include "Dij.H"
#include <gsl/gsl_rng.h>
#include "expr-programs.H" // this is the code to read in and execute our func.


// A fun little procedure... takes in a set of positions, and outputs
// random matrices with an envelope function based on the distance
// from the origin.
// The envelope function is now read in as an "assembler" procedure

// We determine what to pass based on the number of arguments;
// SINGLE means pass the magnitude; 3D means pass all three components.
// NONE means it doesn't need any parameters (very rare)
const int ARGTYPE_NONE = 0;
const int ARGTYPE_SINGLE = 1;
const int ARGTYPE_3D = 3;


// ***************************** SUBROUTINES *************************
// alpha between 0 and 1; x between 0 and 1.  At x=alpha, we transition
// to a smooth cubic down to zero.
inline double smooth_cutoff (double x, double alpha) 
{
  if (x <= alpha) return 1;
  if (x >= 1.0) return 0;
  double y = (x-alpha)/(1.-alpha);
  return 1 + y*y*(-3 + 2*y);
}


// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <grid> <asm>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'c'};
// flag characters.

const char* ARGEXPL = 
"  cell:  cell file (-h for format)\n\
  grid:  positions to evaluate, in lattice coordinates\n\
  asm:   file with \"assembler\" version of function to evaluate (from compiler)\n\
  -c     use a smooth cutoff to zero\n\
  -m <seed>  for random number generator\n\
\n\
  NOTE: assembly program must take either 0, 1 or 3 arguments\n\
  0 passed no arguments\n\
  1 passes the magnitude of the point\n\
  3 passes the cartesian coordinates of the point\n";

const char* FILEEXPL =
"==== cell ====\n\
a0                               # Scale factor for unit cell\n\
a1.x a1.y a1.z                   # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.\n\
Natoms                           # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3                   # Atom locations, in direct coord.\n\
...\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== grid ====\n\
N  # number of points\n\
i1.1 i1.2 i1.3  # position in lattice coordinates\n\
...\n\
==== grid ====\n";

int main ( int argc, char **argv ) 
{
  int i, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Used as our seed...
  int SMOOTHCUT = 0;

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
      // Note: we don't print out the elastic const. info... mainly
      // because we just ignore all that stuff anyway.
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags:

  SMOOTHCUT = flagon[0];

#ifdef MACOSX
  if (MEMORY == 0)
    srandomdev(); // initialize from /dev/random (not reproducible...)
  else
#endif
    srandom(MEMORY);

  char *cell_name, *grid_name, *asm_name;
  cell_name = args[0];
  grid_name = args[1];
  asm_name =  args[2];
  
  FILE* infile;
  //  char dump[512];
  
  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  int Natoms = NO_ATOMS;
  double** u_atoms = NULL;

  //++ ==== cell ====
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms);
  myclose(infile);
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }
  if (TESTING)
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
  //-- ==== cell ====

  double cart2[9];
  square(cart, cart2);

  //++ ==== grid ====
  infile = myopenr(grid_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", grid_name);
    exit(ERROR_NOFILE);
  }

  int Np;
  point_type *p;
  
  ERROR = read_Dij(infile, cart, Np, p, READ_DIJ_NOMAT);
  myclose(infile);
  //-- ==== grid ====
  if (ERROR) {
    fprintf(stderr, "Problem with %s\n", grid_name);
    ERROR = ERROR_BADFILE;
  }

  //++ ==== asm ====
  infile = myopenr(asm_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", asm_name);
    exit(ERROR_NOFILE);
  }
  program_pointer_type* mach = NULL;
  
  // Just go ahead and assemble the code
  read_mach(infile, mach);
  myclose(infile);

  //-- ==== asm ====
  int ARGTYPE = mach->Narg;
  double* arglist = NULL;
  ERROR = (ARGTYPE != ARGTYPE_NONE) &&
    (ARGTYPE != ARGTYPE_SINGLE) &&
    (ARGTYPE != ARGTYPE_3D);
  
  if (ERROR) {
    fprintf(stderr, "Your assembler program must have 0, 1, or 3 arguments.\n");
    exit(ERROR);
  }
  // We no longer bother allocating this... but in the future, we
  // may need to change this:
  // if (ARGTYPE > 0) arglist = new double[ARGTYPE];
  

  // *************************** ANALYSIS ****************************
  double maxmagn=0, alpha=0.25; // may want to change this later...
  for (n=0; n<Np; ++n)
    if (p[n].Rmagn > maxmagn) maxmagn = p[n].Rmagn;

  double magn;
  double scale = 1.;
  for (n=0; n<Np; ++n) {
    point_type* tp = p+n;
    magn = tp->Rmagn;
    // set arglist to point to the correct place:
    switch (ARGTYPE) {
    case ARGTYPE_NONE: break;
    case ARGTYPE_SINGLE: 
      arglist = &magn;
      break;
    case ARGTYPE_3D:
      arglist = tp->Rcart;
      break;
    }
    if ( zero(magn) )
      for (i=0; i<9; ++i) tp->mat[i] = 0.;
    else {
      if (SMOOTHCUT) scale = smooth_cutoff(magn/maxmagn, alpha);
      for (i=0; i<9; ++i)
	tp->mat[i] = evalprogram(mach, arglist) * scale;
    }
  }

  // **************************** OUTPUT *****************************

  write_Dij(stdout, Np, p);


  // ************************* GARBAGE COLLECTION ********************  
  // no garbage collection on arglist, since we don't declare memory:
  //  if (arglist != NULL) delete[] arglist;
  delete[] p;
  free_program_pointer(mach);
  free_cell(Cmn_list, u_atoms, Natoms);
  delete[] Cmn_list;

  return 0;
}
