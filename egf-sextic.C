/*
  Program: egf-sextic.C
  Author:  D. Trinkle
  Date:    August 30, 2004
  Purpose: Given our cell, our sextic roots (from sextic) and a set of
           atomic positions, we compute the 2d EGF at each point.  We do
	   NOT rotate our GF, so it's in cartesian coordinates.

  Param.:  <cell> <sextic> <latt>
           cell:    cell file (see below for format)
           sextic:  dislocation coord, and sextic roots
	   latt:    set of atomic positions (unit cell coord)

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
	   
	   ==== sextic ====
	   t1.x t1.y t1.z t1.1 t1.2 t1.3 # line direction (cart and unit)
	   m1.x m1.y m1.z # dislocation cut vector (perp to t, in slip plane)
	   n1.x n1.y n1.z # mutual perp. vector
	   p1.r p1.i
	   Re(A1.x) Im(A1.x) Re(A1.y) Im(A1.y) Re(A1.z) Im(A1.z)
	   Re(L1.x) Im(L1.x) Re(L1.y) Im(L1.y) Re(L1.z) Im(L1.z)
	   p2.r p2.i
	   Re(A2.x) Im(A2.x) Re(A2.y) Im(A2.y) Re(A2.z) Im(A2.z)
	   Re(L2.x) Im(L2.x) Re(L2.y) Im(L2.y) Re(L2.z) Im(L2.z)
	   p3.r p3.i
	   Re(A3.x) Im(A3.x) Re(A3.y) Im(A3.y) Re(A3.z) Im(A3.z)
	   Re(L3.x) Im(L3.x) Re(L3.y) Im(L3.y) Re(L3.z) Im(L3.z)
	   ==== sextic ====

	   ==== latt ====
	   Np
	   u1.1 u1.2 u1.3
	   ...
	   uN.1 uN.2 uN.3
	   ==== latt ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   We take in our sextic roots and eigenvectors, and then
           evaluate the 2d EGF.

  Output:  EGF at each point, except R=0.
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
#include "elastic.H"
#include "cell.H"
#include "Dij.H"

//****************************** SUBROUTINES ****************************

inline void print_mat (double a[9]) 
{
  printf("xx=%10.3le yy=%10.3le zz=%10.3le xy=%10.3le yz=%10.3le zx=%10.3le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}

inline void print_vect (double a[3]) 
{
  printf(" %.5lf %.5lf %.5lf", a[0], a[1], a[2]);
}

inline void print_vect (int a[3]) 
{
  printf(" %d %d %d", a[0], a[1], a[2]);
}


inline double log(double mx, double nx, double pr, double pi) 
{
  double r = hypot(mx + nx*pr, nx*pi);
  return log(r);
}

// We "cast" pi to >0 and handle the branch cut as carefully as we
// can.  IN THEORY we shouldn't need this... but it's here for safety's
// sake.
inline double arg(double mx, double nx, double pr, double pi) 
{
  double x = mx + nx*pr;
  double y = nx*fabs(pi);
  double theta;
  if (zero(y)) {
    if (x>0) theta = 0;
    else     theta = M_PI;
  }
  else {
    theta = atan2(y, x);
    if (theta < 0) theta += 2*M_PI;
  }
  if (dcomp(theta, 2*M_PI)) theta = 0;
  if (pi > 0)
    return theta;
  else
    return -theta;
}




/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <sextic> <latt>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:    cell file (see below for format)\n\
  sextic:  dislocation coord, and sextic roots\n\
  latt:    set of atomic positions (unit cell coord)";

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
==== sextic ====\n\
t1.x t1.y t1.z  t1.1 t1.2 t1.3 # line direction (cart and unit)\n\
m1.x m1.y m1.z # dislocation cut vector (perp to t, in slip plane)\n\
n1.x n1.y n1.z # mutual perp. vector\n\
p1.r p1.i\n\
Re(A1.x) Im(A1.x) Re(A1.y) Im(A1.y) Re(A1.z) Im(A1.z)\n\
Re(L1.x) Im(L1.x) Re(L1.y) Im(L1.y) Re(L1.z) Im(L1.z)\n\
p2.r p2.i\n\
Re(A2.x) Im(A2.x) Re(A2.y) Im(A2.y) Re(A2.z) Im(A2.z)\n\
Re(L2.x) Im(L2.x) Re(L2.y) Im(L2.y) Re(L2.z) Im(L2.z)\n\
p3.r p3.i\n\
Re(A3.x) Im(A3.x) Re(A3.y) Im(A3.y) Re(A3.z) Im(A3.z)\n\
Re(L3.x) Im(L3.x) Re(L3.y) Im(L3.y) Re(L3.z) Im(L3.z)\n\
==== sextic ====\n\
\n\
==== latt ====\n\
Np\n\
u1.1 u1.2 u1.3\n\
...\n\
uN.1 uN.2 uN.3\n\
==== latt ====\n";

int main ( int argc, char **argv ) 
{
  int a, i, j, k, d; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used

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
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "Crystal classes:\n%s\n", CRYSTAL_CLASS);
      fprintf(stderr, "\nElastic constants ordering:\n");
      for (k=0; k<NCLASSES; ++k) {
        fprintf(stderr, "  Class %2d (%d):", k, class_len[k]);
        for (i=0; i<class_len[k]; ++i)
          fprintf(stderr, " C_%2d", class_Cij[k][i]);
        fprintf(stderr, "\n");
      }
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  
  // ****************************** INPUT ****************************
  // Command line parameters:
  char *cell_name, *sextic_name, *latt_name;
  // Let's pull off the args:
  cell_name = args[0];
  sextic_name = args[1];
  latt_name = args[2];

  char dump[512];
  FILE* infile;

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

  //++ ==== sextic ====
  infile = myopenr(sextic_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", sextic_name);
    exit(ERROR_NOFILE);
  }

  double t_vect[3], m_vect[3], n_vect[3];
  int t_unit[3];
  double pr[3], pi[3];
  double AR_ai[3][3], AI_ai[3][3];
  // we don't bother reading in L

  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%lf %lf %lf %d %d %d", t_vect, t_vect+1, t_vect+2,
	 t_unit, t_unit+1, t_unit+2);
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%lf %lf %lf", m_vect, m_vect+1, m_vect+2);
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%lf %lf %lf", n_vect, n_vect+1, n_vect+2);

  // sextic roots and eigenvectors:
  for (a=0; a<3; ++a) {
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%lf %lf", pr+a, pi+a);
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%lf %lf %lf %lf %lf %lf", AR_ai[a], AI_ai[a],
	   AR_ai[a]+1, AI_ai[a]+1, AR_ai[a]+2, AI_ai[a]+2);
    nextnoncomment(infile, dump, sizeof(dump));
    // L line--ignored
  }

  if ( (t_unit[0] == 0) &&  (t_unit[0] == 0) &&  (t_unit[0] == 0) ) {
    fprintf(stderr, "t=000 in unit coord from %s\n", sextic_name);
    ERROR = ERROR_BADFILE;
  }

  myclose(infile);
  //-- ==== sextic ====

  if (TESTING) {
    printf("## threading direction t:");
    print_vect(t_vect); print_vect(t_unit);
    printf("\n");
    printf("## sextic roots\n");
    for (a=0; a<3; ++a) {
      printf("## p%d= %+.15lf%+.15lfI\n", a+1, pr[a], pi[a]);
      printf("## A%d=", a+1);
      for (i=0; i<3; ++i)
	printf(" %+.4lf%+.4lfI", AR_ai[a][i], AI_ai[a][i]);
      printf("\n");
    }
  }

  //++ ==== latt ====
  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  int Np;
  point_type *latt;
  
  ERROR = read_Dij(infile, cart, Np, latt, READ_DIJ_NOMAT);

  myclose(infile);
  //-- ==== latt ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }
  
  // ***************************** ANALYSIS **************************
  // We need to make our Asq matrices: Ar and Ai; there's one for
  // each alpha:
  double Ar[3][9], Ai[3][9];
  
  for (a=0; a<3; ++a) {
    for (d=0; d<9; ++d) {
      Ar[a][d] = 0;
      Ai[a][d] = 0;
    }
    for (i=0; i<3; ++i)
      for (j=0; j<3; ++j) {
	d = 3*i+j;
	Ar[a][d] = AR_ai[a][i]*AR_ai[a][j] - AI_ai[a][i]*AI_ai[a][j];
	Ai[a][d] = AI_ai[a][i]*AR_ai[a][j] + AI_ai[a][i]*AR_ai[a][j];
      }
  }

  // scale by -1/Pi:
  double scale = -1./M_PI;
  for (a=0; a<3; ++a)
    for (d=0; d<9; ++d) {
      Ar[a][d] *= scale;
      Ai[a][d] *= scale;
    }
  
  // check our most basic sum rule:
  for (d=0; d<9; ++d) {
    double sum=0;
    for (a=0; a<3; ++a) sum += Ar[a][d];
    if (! zero(sum)) {
      fprintf(stderr, "ERROR: the eigenvectors do not obey the sum rule\n");
    }
    ERROR = ERROR_BADFILE;
  }

  // Now, go point by point and calculate!
  for (int n=0; n<Np; ++n) {
    point_type *tp = latt+n;
    double nx, mx;
    mx = dot(m_vect, tp->Rcart);
    nx = dot(n_vect, tp->Rcart);
    for (d=0; d<9; ++d) tp->mat[d] = 0;
    
    if ( (! zero(mx)) || (! zero(nx)) ) {
      for (a=0; a<3; ++a) {
	double ln_x, arg_x;
	ln_x = log(mx, nx, pr[a], pi[a]);
	arg_x = arg(mx, nx, pr[a], pi[a]);
	for (d=0; d<9; ++d)
	  tp->mat[d] += Ai[a][d]*ln_x + Ar[a][d]*arg_x;
      }
    }
    if (TESTING) {
      printf("## Runit= %3d %3d %3d  Rcart= %12.8lf %12.8lf %12.8lf\n", 
	     tp->Runit[0], tp->Runit[1], tp->Runit[2],
	     tp->Rcart[0], tp->Rcart[1], tp->Rcart[2]);
      printf("## m.x= %.5lf n.x= %.5lf\n## G:", mx, nx);
      print_mat(tp->mat);
      printf("\n##\n");
    }
  }

  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  if (VERBOSE) {
    // What, exactly, is "verbose" output...?
  }

  // we roll our own first line, so as to include information
  // about the threading direction
  sprintf(dump, "%d 9  %d %d %d  # Np dim  t in unit coord", Np,
	  t_unit[0], t_unit[1], t_unit[2]);
  write_Dij(stdout, Np, latt, dump);

  // ************************* GARBAGE COLLECTION ********************
  delete[] latt;
  free_cell(Cmn_list, u_atoms, Natoms);
  delete[] Cmn_list;

  return ERROR;
}
