/*
  Program: flatten.C
  Author:  D. Trinkle
  Date:    October 6, 2004
  Purpose: Utility code to read in a lattice function and flatten it
           along a given direction.

  Param.:  <cell> <latt(R)> <t1> <t2> <t3>
           cell:    cell file (see below for format)
           latt(R): lattice function
	   t.i:     threading direction in unit cell coord

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
#include "supercell.H"
#include "shell.H"

//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 5;
const char* ARGLIST = "<cell> <latt(R)> <t1> <t2> <t3>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:    cell file (see below for format)\n\
  latt(R): lattice function\n\
  t.i:     threading direction in unit cell coord.";

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
==== latt(R) ====\n\
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
  // Let's pull off the args:
  char* cell_name = args[0];
  char* latt_name = args[1];
  int t_unit[3];
  sscanf(args[2], "%d", t_unit+0);
  sscanf(args[3], "%d", t_unit+1);
  sscanf(args[4], "%d", t_unit+2);

  if ( (t_unit[0] == 0) &&  (t_unit[1] == 0) &&  (t_unit[2] == 0)) {
    fprintf(stderr, "Bad threading direction--must not be 0.\n");
    exit(ERROR_BADFILE);
  }

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

  //++ ==== latt(R) ====
  int Np;
  point_type *p; // set of points

  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij (infile, cart, Np, p, READ_DIJ_MAT);
  myclose(infile);
  if (ERROR) {delete p; exit(ERROR);}
  //-- ==== latt(R) ====

  // ***************************** ANALYSIS **************************
  // 0. Setup to flatten
  double tdotvect[3], tvect[3];
  double cart2[9];
  square(cart, cart2);
  double tlen2 = innerprod(t_unit, cart2, t_unit);
  mult_vect(cart2, t_unit, tdotvect);
  for (i=0; i<3; ++i) tdotvect[i] *= 1./tlen2;
  mult_vect(cart, t_unit, tvect);
  
  // 1. Run through all of the elements and flatten positions
  int maxn=0, alpha;
  double x;
  for (n=0; n<Np; ++n) {
    point_type* tp = p + n;
    x = dot(tp->Runit, tdotvect);
    alpha = lround(x);
    if ( (x - alpha) < -TOLER) --alpha;
    for (i=0; i<3; ++i) {
      tp->Runit[i] += -alpha*t_unit[i];
      if (abs(tp->Runit[i]) > maxn) maxn = abs(tp->Runit[i]);
    }
    mult_vect(cart, tp->Runit, tp->Rcart);
    x = dot(tp->Rcart, tvect) / tlen2;
    double r[3] = {tp->Rcart[0]-x*tvect[0], tp->Rcart[1]-x*tvect[1], 
		   tp->Rcart[2]-x*tvect[2]};
    tp->Rmagn = sqrt(dot(r,r));
  }
  
  // 2. Flatten lattice function
  int*** index;
  make_index(maxn, index);
  point_type* fp = new point_type[Np]; // flattened list...
  int Nflat = 0;
  
  for (n=0; n<Np; ++n) {
    point_type* tp = p + n;
    i = get_index(maxn, index, tp->Runit);
    if (i == -1) {
      // New point:
      copy_point(*tp, fp[Nflat]);
      set_index(maxn, index, tp->Runit, Nflat);
      ++Nflat;
    }
    else {
      // Old point:
      for (d=0; d<9; ++d) fp[i].mat[d] += tp->mat[d];
    }
  }

  // 3. Sort
  sort_points(fp, Nflat);
  
  // ****************************** OUTPUT ***************************
  sprintf(dump, "%d %3d %3d %3d  # number of points, threading direction",
	  Nflat, t_unit[0], t_unit[1], t_unit[2]);
  write_Dij(stdout, Nflat, fp, dump);
  
  // ************************* GARBAGE COLLECTION ********************
  free_index(maxn, index);
  delete[] p;
  delete[] fp;
  
  free_cell(Cmn_list, u_atoms, Natoms);

  return 0;
}
