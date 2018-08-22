/*
  Program: 2d-fourth-FT.C
  Author:  D. Trinkle
  Date:    September 11, 2004
  Purpose: Read in a dynamical matrix, and calculate:

                       1
	   [ab,cd] = - - SUM D   x x
                       2  x   ab  c d

                          1
	   [ab,cdef] = - -- SUM D   x x x x
                         24  x   ab  c d e f

	   From these, we construct the matrices:

	   lambda (k) = SUM [ab,cd]k k
	                 cd         c d

                 (2)
	   lambda  (k) =  SUM [ab,cdef]k k k k 
	                 cdef           c d e f

	   where k is normalized.  We compute the matrix:
                      -1       (2)               -1
	   [lambda(k)]  [lambda  (k)] [lambda(k)]

	   for each k point in our mesh around a circle.
	   This has units of A^3/eV (i.e., same as our GF).  We write
	   this as a Fourier series.

  Param.:  <cell> <Dij(R)> <dislcoord>
           cell:   cell file describing our lattice
	   Dij(R): dynamical matrix
	   disl:   dislocation coordinate system

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

           ==== Dij(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij
           ...
           nN.1 nN.2 nN.3  Dxx .. Dzz
           ==== Dij(R) ====

           ==== disl ====
           t1.1 t1.2 t1.3  # line direction
           m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)
           n1.1 n1.2 n1.3  # mutual perp. vector
           # all three vectors are given in unit cell coord.
           ==== disl ====

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
#include <gsl/gsl_fft_complex.h> // GSL FFT routines.
#include <gsl/gsl_integration.h> // GSL integration routines (for shift)
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "voigt.H"
#include "Dij.H"
#include "sphere-harm.H"
#include "lambda.H"
#include "dislocation.H"

const double FFT_TOLER = 1e-9; // tolerance for when we truncate our series

//***************************** SUBROUTINES ****************************

// Needed for our packed complex arrays:
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

inline double zmagn (const double x, const double y) 
{
  return sqrt(x*x+y*y);
}

inline void print_val (const double &x) 
{
  printf(" %+.12le", x);
  //   if (zero(x)) printf(" 0.0");
  //   else         printf(" %+.12le", x);
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
const char* ARGLIST = "<cell> <Dij(R)> <disl>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'e'}; // flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  Dij(R): dynamical matrix\n\
  disl:   dislocation coordinate system\n\
  -e      only output even components";

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
==== Dij(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij\n\
...\n\
nN.1 nN.2 nN.3  Dxx .. Dzz\n\
==== Dij(R) ====\n\
\n\
==== disl ====\n\
t1.1 t1.2 t1.3  # line direction\n\
m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)\n\
n1.1 n1.2 n1.3  # mutual perp. vector\n\
# all three vectors are given in unit cell coord.\n\
==== disl ====\n";

int main ( int argc, char **argv ) 
{
  int d, i; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 128;   // not used

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (i=0; i<NFLAGS; ++i) flagon[i] = 0;
  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  if (MEMORY < 2) {
    fprintf(stderr, "Bad MEMORY value--needs to be >=1\n");
    ERROR = ERROR_BADFILE;
  }
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
  int EVEN = flagon[0];
  
  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char *cell_name = args[0];
  char *Dij_name = args[1];
  char *disl_name = args[2];


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
  square(cart, cart2);
  //-- ==== cell ====

  //++ ==== Dij(R) ====
  int Np;
  point_type *p; // set of points

  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, p, READ_DIJ_MAT);
  if (ERROR) {delete p; exit(ERROR);}
  //-- ==== Dij(R) ====


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

  // 1. Use Voigt notation to generate the symmetrized elements:
  double lambda[9][9], lambda2[9][9][9];

  make_lambda(Np, p, lambda);
  make_lambda2(Np, p, lambda2);
  
  if (TESTING) {
    printf("## lambda: [ab,cd]: --non-zero only\n");
    print_lambda(lambda, "## ");
    printf("## lambda(2): [ab,cdef]: --non-zero only\n");
    print_lambda2(lambda2, "## ");
  }

  // 2. Now build our function to fourier transform (packed array)
  // MEMORY is our size
  int Nelem = MEMORY;
  double dtheta = 2.*M_PI/(double)Nelem, theta;
  double** FFT_data;
  double Gdc[9], ktheta[3];
  double gf_scale = det(cart); // scale to have the same units as EGF
  FFT_data = new double*[9];
  for (d=0; d<9; ++d) FFT_data[d] = new double[Nelem*2];
  
  for (i=0; i<Nelem; ++i) {
    theta = dtheta * i;
    for (d=0; d<3; ++d) ktheta[d] = m0[d]*cos(theta)+n0[d]*sin(theta);
    eval_lambda2(lambda, lambda2, ktheta, Gdc);
    mult(Gdc, gf_scale, Gdc);
    for (d=0; d<9; ++d) {
      REAL(FFT_data[d], i) = Gdc[d];
      IMAG(FFT_data[d], i) = 0.;
    }
    if (TESTING) {
      printf("## i=%3d theta=%.5lf  Gdc:", i, theta);
      print_mat(Gdc);
      printf("\n");
    }
  }
  
  // 3. FFT our data.
  gsl_fft_complex_wavetable *wavetable=gsl_fft_complex_wavetable_alloc(Nelem);
  gsl_fft_complex_workspace *workspace=gsl_fft_complex_workspace_alloc(Nelem);

  for (d=0; d<9; ++d) {
    gsl_fft_complex_forward (FFT_data[d], 1, Nelem, wavetable, workspace);
    // scale
    for (i=0; i<(2*Nelem); ++i) FFT_data[d][i] *= 1./(double)Nelem;
  }
  // garbage collection
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);

  // 4. Tag the elements we want to output:
  int out[Nelem/2+1];
  double zmax=0, zval;
  for (i=0; i<=(Nelem/2); ++i) {
    zval = 0;
    for (d=0; d<9; ++d) 
      zval += gsl_pow_2(REAL(FFT_data[d],i))
        + gsl_pow_2(IMAG(FFT_data[d],i));
    zval = sqrt(zval);
    if (zval > zmax) zmax = zval;
  }
  // now, let's tag 'em:
  for (i=0; i<=(Nelem/2); ++i) {
    zval = 0;
    for (d=0; d<9; ++d) 
      zval += gsl_pow_2(REAL(FFT_data[d],i))
        + gsl_pow_2(IMAG(FFT_data[d],i));
    zval = sqrt(zval);
    out[i] = (zval > FFT_TOLER*zmax);
  }

  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  if (VERBOSE) {
  }

  // output...
  printf("%d  # fourth-order\n", Nelem/2);
  // if we try to use an odd number of elements, this will do it correctly:

  for (i=0; i<((Nelem/2) - (1 - (Nelem%2))); ++i)
    if (out[i])
      if ( (!EVEN) || ((i%2)==0) ) {
        printf("%3d", i);
        for (d=0; d<9; ++d) print_val(REAL(FFT_data[d], i));
        printf("\n");
        printf("%3d", i);
        for (d=0; d<9; ++d) print_val(IMAG(FFT_data[d], i));
        printf("\n");
        if (i != 0) {
          printf("%3d", -i);
          for (d=0; d<9; ++d) print_val(REAL(FFT_data[d], Nelem-i));
          printf("\n");
          printf("%3d", -i);
          for (d=0; d<9; ++d) print_val(IMAG(FFT_data[d], Nelem-i));
          printf("\n");
        }
      }
  if ( (Nelem%2) == 0)
    if (out[Nelem/2]) {
      i = Nelem/2;
      printf("%3d", i);
      for (d=0; d<9; ++d) print_val(0.5*REAL(FFT_data[d], i));
      printf("\n");
      printf("%3d", i);
      for (d=0; d<9; ++d) print_val(0.5*IMAG(FFT_data[d], i));
      printf("\n");
      
      printf("%3d", -i);
      for (d=0; d<9; ++d) print_val(0.5*REAL(FFT_data[d], i));
      printf("\n");
      printf("%3d", -i);
      for (d=0; d<9; ++d) print_val(0.5*IMAG(FFT_data[d], i));
      printf("\n");
    }
    
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
