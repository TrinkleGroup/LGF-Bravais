/*
  Program: 2d-egf-test.C
  Author:  D. Trinkle
  Date:    September 9, 2004
  Purpose: Interpolate using the 2d-egf results around a circle.

  Param.:  <cell> <infile> <2d.EGF>
           cell:     cell file (see below for format)
           infile:   input file (see below for format)
	   2d.EGF:   FT of 2d EGF

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
	   
	   ==== infile ====
	   t1 t2 t3     # dislocation line direction (unit cell coord.)
	   m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)
	   n1 n2 n3
	   ==== infile ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
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
#include "elastic.H"
#include "cell.H"
#include "integrate.H"
#include "dislocation.H"

//****************************** SUBROUTINES ****************************

inline void m_theta(double m[3], double n[3], double theta, double mt[3]) 
{
  mt[0] = cos(theta)*m[0] + sin(theta)*n[0];
  mt[1] = cos(theta)*m[1] + sin(theta)*n[1];
  mt[2] = cos(theta)*m[2] + sin(theta)*n[2];
}

inline void n_theta(double m[3], double n[3], double theta, double nt[3]) 
{
  nt[0] = -sin(theta)*m[0] + cos(theta)*n[0];
  nt[1] = -sin(theta)*m[1] + cos(theta)*n[1];
  nt[2] = -sin(theta)*m[2] + cos(theta)*n[2];
}


//void rotate (const double* &G, const double* x0, 
//	     const double* &x1, const double* &x2,
//	     double* Gii) 
void rotate (double G[9], const double x0[3], const double x1[3], const double x2[3],
	     double Gii[6]) 
{
  Gii[0] = innerprod(x0, G, x0);
  Gii[1] = innerprod(x1, G, x1);
  Gii[2] = innerprod(x2, G, x2);

  Gii[3] = innerprod(x0, G, x1);
  Gii[4] = innerprod(x1, G, x2);
  Gii[5] = innerprod(x2, G, x0);
}

inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <infile> <2d.EGF>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'r', 's'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:     cell file (-h for format)\n\
  infile:   input file (-h for format)\n\
  2d.EGF:   2d EGF FT, from 2d-egf-FT\n\
  -r               output rotated matrix, not cartesian matrix\n\
  -s               rotated matrix, where the components are *also* rotated";

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
==== infile ====\n\
t1 t2 t3     # dislocation line direction (unit cell coord.)\n\
m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)\n\
n1 n2 n3     # perpendicular direction\n\
==== infile ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, d; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 16384; // 2^14, default (gets turned into Nsteps for int.)

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
      printf("Input file format:\n%s\n", FILEEXPL);
      printf("Crystal classes:\n%s\n", CRYSTAL_CLASS);
      printf("\nElastic constants ordering:\n");
      for (k=0; k<NCLASSES; ++k) {
        printf("  Class %2d (%d):", k, class_len[k]);
        for (i=0; i<class_len[k]; ++i)
          printf(" C_%2d", class_Cij[k][i]);
        printf("\n");
      }
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags:
  int ROTATE = flagon[0];
  int SPIN = flagon[1];
  if (SPIN) ROTATE = 1;
  int Nsteps = MEMORY;
  
  // We're going to use the number of steps according to our preferred
  // amount of memory allocation.
  if (Nsteps < 4) {
    fprintf(stderr, "Nsteps (MEMORY = %d) must be 4 or larger.\n", Nsteps);
    exit(ERROR_BADFILE);
  }  
  
  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  int Natoms;
  double** u_atoms;

  // Let's pull off the args:
  char* cell_name = args[0];
  char* disl_name = args[1];
  char* egf_name = args[2];
  
  //++ ==== cell ====
  // First, read in the cell.
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(ERROR_NOFILE);
  }
  Natoms = 0;
  ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms);
  myclose(infile);
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    fprintf(stderr, "Problem with cell: %s\n", cell_name);
    exit(ERROR);
  }

  if (TESTING)
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
  //-- ==== cell ====

  //++ ==== disl ====
  infile = myopenr(disl_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", disl_name);
    exit(ERROR_NOFILE);
  }

  int t_unit[3], m_unit[3], n_unit[3];
  ERROR = read_disl(infile, t_unit, m_unit, n_unit);
  myclose(infile);
  //-- ==== disl ====

  if (ERROR) {
    fprintf(stderr, "Problem with dislocation: %s\n", disl_name);
    exit(ERROR);
  }

  //++ ==== egf.ft ====
  infile = myopenr(egf_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", egf_name);
    exit(ERROR_NOFILE);
  }
  int N;
  double** gncos, **gnsin;
  
  ERROR = read_2degf(infile, N, gncos, gnsin);
  myclose(infile);
  //-- ==== egf.ft ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Problem with 2d EGF: %s\n", egf_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## FT components: %d\n", N*2);
    for (i=0; i<=N; ++i) {
      printf("## cos(%d t):", 2*i);
      print_mat(gncos[i]);
      printf("\n");
      printf("## sin(%d t):", 2*i);
      print_mat(gnsin[i]);
      printf("\n");
    }
  }

  // ***************************** ANALYSIS **************************
  // 0. Construct our coordinate system
  double t0[3], m0[3], n0[3];

  make_dislcoord(cart, t0, m0, n0, t_unit, m_unit, n_unit, 
		 DISLCOORD_UNIT | DISLCOORD_NOCHANGE);

  if (TESTING) {
    printf("## vect    input                 cartesian\n");
    printf("## t = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           t_unit[0],t_unit[1],t_unit[2], t0[0],t0[1],t0[2]);
    printf("## m = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           m_unit[0],m_unit[1],m_unit[2], m0[0],m0[1],m0[2]);
    printf("## n = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           n_unit[0],n_unit[1],n_unit[2], n0[0],n0[1],n0[2]);
  }

  // 1. Scale EGF to get the real-space fourier series:
  int isign = -1;
  double scale;
  for (i=1; i<=N; ++i, isign = -isign) {
    scale = isign / (4.*M_PI*i);
    mult(gncos[i], scale, gncos[i]);
    mult(gnsin[i], scale, gnsin[i]);
  }
  scale = -1/(2.*M_PI);
  mult(gncos[0], scale, gncos[0]);

  // ****************************** OUTPUT ***************************
  // Now some evaluating of integrals
  double theta;
  double dtheta = 2*M_PI / Nsteps;
  double Gw[9], Gii[6];
  double mt[3], nt[3];

  if (!SPIN) {
    m_theta(m0, n0, 0, mt);
    n_theta(m0, n0, 0, nt);
  }
  for (k=0; k<=Nsteps; ++k) {
    theta = k*dtheta;
    for (d=0; d<9; ++d) Gw[d] = 0;
    for (i=1; i<=N; ++i) {
      for (d=0; d<9; ++d) Gw[d] += cos(2*i*theta)*gncos[i][d];
      for (d=0; d<9; ++d) Gw[d] += sin(2*i*theta)*gnsin[i][d];
    }
    if (SPIN) {
      m_theta(m0, n0, theta, mt);
      n_theta(m0, n0, theta, nt);
    }
    if (ROTATE) rotate(Gw, mt, nt, t0, Gii);
    else        rotate(Gw, ident, ident+3, ident+6, Gii);
    printf("%.8lf", theta);
    for (d=0; d<6; ++d) printf(" %.8le", Gii[d]);
    printf("\n");
  }

  printf("&\n");

  if (ROTATE) rotate(gncos[0], m0, n0, t0, Gii);
  else        rotate(gncos[0], ident, ident+3, ident+6, Gii);
  printf("0");
  for (d=0; d<6; ++d) printf(" %.8le", Gii[d]);
  printf("\n");


  // ************************* GARBAGE COLLECTION ********************
  for (i=0; i<=N; ++i) {
    delete[] gncos[i];
    delete[] gnsin[i];
  }
  delete[] gncos;
  delete[] gnsin;

  free_cell(Cmn_list, u_atoms, Natoms);

  return 0;
}
