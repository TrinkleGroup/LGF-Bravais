/*
  Program: difflatt.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Utility code to read in two lattice functions and determine
           the difference between them.  We have a few ways to output this;
	   we define the difference as:

	   abserr = sqrt( sum_ij (lat1_ij - lat2_ij)^2 )
	   relerr = sqrt( sum_ij (lat1^-1*(lat1-lat2))_ij^2 )
	   
	   1. absolute error as a function of distance
	   2. relative error as a function of distance

  Param.:  <latt.1(R)> <latt.2(R)>
           latt.i(R): lattice functions

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


//***************************** SUBROUTINES ****************************

int __COMPARE_LATT_TYPE__ (const void* a_elem, const void* b_elem) 
{
  point_type *x = (point_type*)a_elem;
  point_type *y = (point_type*)b_elem;
  int a, b, comp=0, d;
  for (d=0; (d<3) && (!comp); ++d) {
    a = abs(x->Runit[d]);  b = abs(y->Runit[d]);
    if (a<b) comp = -1;
    else
      if (a>b) comp = 1;
      else {
        a = x->Runit[d];  b = y->Runit[d];
	if (a<b) comp = -1;
	else
	  if (a>b) comp = 1;
	  else comp = 0;
      }
  }
  return comp;
}
  
int min_point(int Nlatt, point_type **tp, int* use) 
{
  int i, j;
  int comp, best;
  
  // initialize "best" by skipping past the NULL pointers.
  for (best=0; (best<Nlatt) && (tp[best] == NULL); ++best) 
    use[best] = 0;
  for (i=best; i<Nlatt; ++i) {
    if (tp[i] == NULL) 
      use[i] = 0;
    else {
      comp = __COMPARE_LATT_TYPE__((void*)(tp[i]), (void*)(tp[best]));
      if (comp == 0) {
	use[i] = 1; // equal points:
      }
      else {
	if (comp < 0) {
	  // new best point:
	  for (j=i-1; j>=0; --j) use[j] = 0; // zero out old points
	  best = i;
	  use[i] = 1;
	}
	else
	  use[i] = 0; // worse than best point
      }
    }
  }
  if (best==Nlatt) return -1;
  else return best;
}


inline void print_mat (double mat[9]) 
{
  printf("%.3le %.3le %.3le %.3le %.3le %.3le %.3le %.3le %.3le",
	 mat[0], mat[1], mat[2], mat[3], mat[4], mat[5], 
	 mat[6], mat[7], mat[8]);
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <latt.1(R)> <latt.2(R)";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'r'}; // flag characters.

const char* ARGEXPL = 
"  cell:      cell\n\
  latt.i(R): lattice functions\n\
  -r         output relative error instead of absolute error\n\
  absolute error = sqrt(sum_ij (lat1-lat2)_ij^2 )\n\
  relative error = sqrt(sum_ij lat1^-1*(lat1-lat2)_ij ^2 )";

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
  int RELERROR = 0;

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

  // flags
  RELERROR = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters: 
  const int Nlatt = 2;
  char *cell_name;
  char *latt_name[Nlatt];

  cell_name = args[0];
  // Let's pull off the args:
  for (i=0; i<Nlatt; ++i) 
    latt_name[i] = args[i+1];

  double cart[9], cart2[9];
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
    printf("# Flag settings:\n");
    if (RELERROR)
      printf("# --relative error not absolute\n");
  }
  square(cart, cart2);
  //-- ==== cell ====

  //++ ==== latt(R) ====
  int Np[Nlatt], Ndim = 9;
  point_type *p[Nlatt]; // set of points

  for (i=0; i<Nlatt; ++i) {
    infile = myopenr(latt_name[i]);
    if (infile == NULL) {
      fprintf(stderr, "Couldn't open %s for reading.\n", latt_name[i]);
      exit(ERROR_NOFILE);
    }
    ERROR = read_Dij (infile, cart, Np[i], p[i], READ_DIJ_MAT);
    myclose(infile);
  }

  // To facilitate the combination of all of our lattice functions,
  // we first sort each set:
  for (i=0; i<Nlatt; ++i)
    qsort((void*)(p[i]), Np[i], sizeof(point_type), __COMPARE_LATT_TYPE__);

  if (Np[0] != Np[1]) {
    fprintf(stderr, "Two lattice functions must have the same number of points.\n");
    ERROR = ERROR_BADFILE;
  }
  else {
    for (n=0; (n<Np[0]) && (!ERROR); ++n)
      ERROR = ! equal_vect(p[0][n].Runit, p[1][n].Runit);
    if (ERROR) {
      fprintf(stderr, "Two lattice functions must have equal points.\n");
    }
  }
    
  if (ERROR) {
    for (i=0; i<Nlatt; ++i) delete p[i];
    exit(ERROR);
  }
  //-- ==== latt(R) ====

  // ***************************** ANALYSIS **************************
  if (TESTING) {
    for (i=0; i<Nlatt; ++i) {
      printf("## ==== lattice function %d ====\n", i);
      for (n=0; n<Np[i]; ++n) {
	printf("# %3d %3d %3d ", p[i][n].Runit[0], p[i][n].Runit[1], 
	       p[i][n].Runit[2]);
	print_mat(p[i][n].mat);
	printf("\n");
      }
      printf("\n");
    }
  }
  
  int Nsum = Np[0];
  double Rmagn[Nsum], error[Nsum];
  double diff[9];
  double temp[9], invlatt[9], scale;

  for (n=0; n<Nsum; ++n) {
    Rmagn[n] = p[0][n].Rmagn;
    for (d=0; d<9; ++d) 
      diff[d] = p[0][n].mat[d] - p[1][n].mat[d];
    if (RELERROR) {
      for (scale=0, d=0; d<9; ++d) scale += p[0][n].mat[d]*p[0][n].mat[d];
      scale = sqrt(scale/3.);
      if (! dcomp(scale, 0))
	for (d=0; d<9; ++d) diff[d] = diff[d]/scale;
      else
	for (d=0; d<9; ++d) diff[d] = 0;
    }
    error[n] = 0;
    for (d=0; d<9; ++d) error[n] += diff[d]*diff[d];
    error[n] = sqrt(error[n]/3);
  }
  
  // ****************************** OUTPUT ***************************
  for (n=0; n<Nsum; ++n)
    if ( (!RELERROR) || (! zero(error[n])) )
      printf("%.5lf %.8le\n", Rmagn[n], error[n]);
  
  // ************************* GARBAGE COLLECTION ********************
  for (i=0; i<Nlatt; ++i)
    delete[] p[i];
  
  return ERROR;
}
