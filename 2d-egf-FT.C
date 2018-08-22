/*
  Program: 2d-egf-FT.C
  Author:  D. Trinkle
  Date:    September 8, 2004
  Purpose: Calculate the Fourier expansion of the 2d EGF.

	   Note: Cij are read in in GPa, and converted to eV/A^3 unless
	   told not to do so.  I know--it doesn't make a lot of sense.

  Param.:  <cell> <disl coord>
           cell:  cell file (see below for format)
           disl:  dislocation coordinates in unit cell coord.

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

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   We solve the sextic problem using a 6th order eigensolver,
           followed by root polishing (?).  The eigenvectors are 
	   normalized using Stroh's extended sextic definition.

  Output:  The dislocation coordinates, sextic roots, and eigenvectors.
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
#include "elastic.H"
#include "cell.H"
#include "integrate.H"

// 1 GPa = 1e9 J/m^3 = 10/1.60217733 meV/A^3 = 6.241506363 meV/A^3
// or 0.00624 eV/A^3.
// This makes the units of gij(R) A^3 / eV; multiplied by a force
// in eV/A and divided by the distance in A, gives the displacement
// in A.
const double meV_GPa = 10/1.60217733;  // = 10/(ec in coulumbs x 10^19)
const double eV_GPa = 0.01/1.60217733;  // = 0.01/(ec in coulumbs x 10^19)

const double FFT_TOLER = 1e-9; // tolerance for when we truncate our series

//****************************** SUBROUTINES ****************************

// Needed for our packed complex arrays:
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

inline double zmagn (const double x, const double y) 
{
  return sqrt(x*x+y*y);
}


void calc_shift(double Cijkl[9][9], double n_vect[3], double m_vect[3], 
		double shift[9]);


// Calculate (aa) and (ab):
void a_mult_b (double a[3], double b[3], double Cijkl[9][9],
               double ab[9]) 
{
  int i, j, k, l;
  for (i=0; i<3; ++i)
    for (j=0; j<3; ++j) {
      ab[index(i,j)] = 0.;
      for (k=0; k<3; ++k)
        for (l=0; l<3; ++l)
          ab[index(i,j)] += a[k]*Cijkl[index(k,i)][index(j,l)]*b[l];
    }
}

void a_mult_a (double a[3], double Cijkl[9][9], double aa[9]) 
{
  int i, j, k, l;
  for (i=0; i<3; ++i) {
    for (j=0; j<i; ++j)
      aa[index(i,j)] = aa[index(j,i)];
    for (   ; j<3; ++j) {
      aa[index(i,j)] = 0.;
      for (k=0; k<3; ++k)
	for (l=0; l<3; ++l)
	  aa[index(i,j)] += a[k]*Cijkl[index(k,i)][index(j,l)]*a[l];
    }
  }
}


inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}

void print_mat (double a[9], const char* pre) 
{
  int i, j;
  for (i=0; i<3; ++i) {
    printf("%s", pre);
    for (j=0; j<3; ++j)
      printf("  %10.5lf", a[index(i,j)]);
    printf("\n");
  }
}

inline void print_val (const double &x) 
{
  printf(" %+.12le", x);
  //   if (zero(x)) printf(" 0.0");
  //   else         printf(" %+.12le", x);
}

void maple_print_mat (double a[9], const char* name) 
{
  printf("%s := matrix(3,3, [", name);
  for (int d=0; d<8; ++d) { 
    print_val(a[d]); printf(",");
  }
  print_val(a[8]);
  printf("]);\n");
}

void make_unit (double cart[9], double v[3], int u[3]);



/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <disl coord>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'g', 'e'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:  cell file (see below for format)\n\
  disl:  dislocation coordinates in unit cell coord.\n\
  -g               don't convert from GPa to eV/A^3\n\
  -e               only output even values (force symmetrization)\n\
  MEMORY value determines size of FFT grid";

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
==== disl ====\n\
t1.1 t1.2 t1.3  # line direction\n\
m1.1 m1.2 m1.3  # dislocation cut vector (perp to t, in slip plane)\n\
n1.1 n1.2 n1.3  # mutual perp. vector\n\
# all three vectors are given in unit cell coord.\n\
==== disl ====\n";

int main ( int argc, char **argv ) 
{
  int i, k, d; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 128; // Our size for FFT

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

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
  int EV_CONV = !(flagon[0]);  // Convert from GPa to eV/A^3?
  int EVEN = flagon[1];
  
  // ****************************** INPUT ****************************
  // Command line parameters:
  char *cell_name, *disl_name;
  // Let's pull off the args:
  cell_name = args[0];
  disl_name = args[1];

  char dump[512];
  FILE* infile;

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  double Cijkl[9][9];
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

  if (EV_CONV) {
    // Convert Cmn_list from GPa to eV/ang^3:
    for (int n=0; n<class_len[crystal]; ++n)
      Cmn_list[n] *= eV_GPa;
    if (TESTING) {
      printf("## Converted elastic constants:\n");
      verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    }
  }
  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);


  //++ ==== disl ====
  infile = myopenr(disl_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", disl_name);
    exit(ERROR_NOFILE);
  }

  int t_unit[3], m_unit[3], n_unit[3];

  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", t_unit, t_unit+1, t_unit+2);
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", m_unit, m_unit+1, m_unit+2);
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", n_unit, n_unit+1, n_unit+2);

  if ( (t_unit[0]==0) && (t_unit[1]==0) && (t_unit[2]==0) ) {
    fprintf(stderr, "t vector is zero\n");
    ERROR = ERROR_BADFILE;
  }
  if ( (m_unit[0]==0) && (m_unit[1]==0) && (m_unit[2]==0) ) {
    fprintf(stderr, "m vector is zero\n");
    ERROR = ERROR_BADFILE;
  }
  if ( (n_unit[0]==0) && (n_unit[1]==0) && (n_unit[2]==0) ) {
    fprintf(stderr, "n vector is zero\n");
    ERROR = ERROR_BADFILE;
  }

  myclose(infile);
  //-- ==== disl ====
  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  // 0. Construct our coordinate system
  double t_vect[3], m_vect[3], n_vect[3];

  mult_vect(cart, t_unit, t_vect);
  mult_vect(cart, m_unit, m_vect);
  mult_vect(cart, n_unit, n_vect);

  // Normalize:
  double v_magn;
  v_magn = sqrt(dot(t_vect, t_vect));
  for (d=0; d<3; ++d) t_vect[d] /= v_magn;
  // make sure m is normal to t:
  v_magn = dot(m_vect, t_vect);
  for (d=0; d<3; ++d) m_vect[d] -= v_magn * t_vect[d];
  // normalize
  v_magn = sqrt(dot(m_vect, m_vect));
  if (zero(v_magn)) {
    fprintf(stderr, "t and m are parallel...?  t_unit=(%d,%d,%d)  m_unit=(%d,%d,%d)\n",
	    t_unit[0],t_unit[1],t_unit[2], m_unit[0],m_unit[1],m_unit[2]);
    ERROR = ERROR_BADFILE;
    exit(ERROR);
  }
  else 
    for (d=0; d<3; ++d) m_vect[d] /= v_magn;
  // make n as the cross product of t and m:
  crossprod(t_vect, m_vect, n_vect);
  
  int t_new[3], m_new[3], n_new[3];
  make_unit(cart, t_vect, t_new);
  make_unit(cart, m_vect, m_new);
  make_unit(cart, n_vect, n_new);

  if (TESTING) {
    printf("## vect    input        corrected             cartesian\n");
    printf("## t = (%3d %3d %3d)->(%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
	   t_unit[0],t_unit[1],t_unit[2], t_new[0],t_new[1],t_new[2], 
	   t_vect[0],t_vect[1],t_vect[2]);
    printf("## m = (%3d %3d %3d)->(%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
	   m_unit[0],m_unit[1],m_unit[2], m_new[0],m_new[1],m_new[2],
	   m_vect[0],m_vect[1],m_vect[2]);
    printf("## n = (%3d %3d %3d)->(%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
	   n_unit[0],n_unit[1],n_unit[2], n_new[0],n_new[1],n_new[2], 
	   n_vect[0],n_vect[1],n_vect[2]);
  }


  // 1. Now build our function to fourier transform (packed array)
  // MEMORY is our size
  int Nelem = MEMORY;
  double dtheta = 2.*M_PI/(double)Nelem, theta;
  double** FFT_data;
  double kCk[9], EGF[9], ktheta[3];
  double gf_scale = 1.;  // used 1./det(cart) before... not appropriate?
  FFT_data = new double*[9];
  for (d=0; d<9; ++d) FFT_data[d] = new double[Nelem*2];
  
  for (i=0; i<Nelem; ++i) {
    theta = dtheta * i;
    for (d=0; d<3; ++d) ktheta[d] = m_vect[d]*cos(theta)+n_vect[d]*sin(theta);
    a_mult_a(ktheta, Cijkl, kCk);
    careful_inverse(kCk, EGF); // invert
    mult(EGF, gf_scale, EGF);  // scale
    for (d=0; d<9; ++d) {
      REAL(FFT_data[d], i) = EGF[d];
      IMAG(FFT_data[d], i) = 0.;
    }
    if (TESTING) {
      printf("## i=%3d theta=%.5lf  EGF:", i, theta);
      print_mat(EGF);
      printf("\n");
    }
  }
  
  // 2. FFT our data.
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

  // 3. Tag the elements we want to output:
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

  // 4. Calculate our shift
  double shift[9] = {0,0,0, 0,0,0, 0,0,0};
  
  calc_shift(Cijkl, n_vect, m_vect, shift);

  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  if (VERBOSE) {
  }

  // output...
  printf("%d", Nelem/2);
  for (d=0; d<9; ++d) print_val(shift[d]);
  printf(" # FFT elements, shift\n");
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
  for (d=0; d<9; ++d) delete[] FFT_data[d];
  delete[] FFT_data;
  free_cell(Cmn_list, u_atoms, Natoms);

  return 0;
}

typedef struct 
{
  double Cijkl[9][9];
  double* n_vect;
  double* m_vect;
  int d;
} shift_param_type;

/*
// OLD STYLE SHIFT CALCULATION... fails numerically in some cases
// such as Al fcc edge dislocation (t=112)
// Not really sure why...

double shift_int (double x, void *param) 
{
  double *n0 = ((shift_param_type *)param)->n_vect;
  double *m0 = ((shift_param_type *)param)->m_vect;
  int d = ((shift_param_type *)param)->d;
  double nv[3], kCk[9], EGF[9];
  
  nv[0] = -sin(x)*m0[0] + cos(x)*n0[0];
  nv[1] = -sin(x)*m0[1] + cos(x)*n0[1];
  nv[2] = -sin(x)*m0[2] + cos(x)*n0[2];
  a_mult_a(nv, ((shift_param_type *)param)->Cijkl, kCk);
  careful_inverse(kCk, EGF);

  return EGF[d] * log(2*fabs(sin(x)));
}


void calc_shift(double Cijkl[9][9], double n_vect[3], double m_vect[3], 
		double shift[9]) 
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1024);
  double result, error;
  gsl_function F;
  shift_param_type p;
  double pts[2] = {0, M_PI};

  F.function = &shift_int;  
  F.params = (void *)(&p);
  for (int i=0; i<9; ++i)
    for (int j=0; j<9; ++j)
      p.Cijkl[i][j] = Cijkl[i][j];
  p.n_vect = n_vect;
  p.m_vect = m_vect;

  for (int d=0; d<9; ++d) {
    p.d = d;
    gsl_integration_qagp (&F, pts, 2, 0, 1e-9, 1024, w, &result, &error);
    //    gsl_integration_qags (&F, pts[0], pts[1], 0, 1e-5, 1024, w, &result, &error);
    shift[d] = 0.5*M_1_PI*M_1_PI * result;
  }

  gsl_integration_workspace_free(w);
}
*/


double shift_int (double x, void *param) 
{
  double *n0 = ((shift_param_type *)param)->n_vect;
  double *m0 = ((shift_param_type *)param)->m_vect;
  int d = ((shift_param_type *)param)->d;
  double nv[3], kCk[9], EGF[9];
  
  //  double theta = asin(0.5*x); // should give something between 0 and Pi/2
  double sint = 0.5*x;
  double cost = 0.5*sqrt(4-x*x);
  double sum;
  
  nv[0] = -sint*m0[0] + cost*n0[0];
  nv[1] = -sint*m0[1] + cost*n0[1];
  nv[2] = -sint*m0[2] + cost*n0[2];
  a_mult_a(nv, ((shift_param_type *)param)->Cijkl, kCk);
  careful_inverse(kCk, EGF);
  sum = EGF[d];

  // now do Pi-theta: sin->sin, cos->-cos
  nv[0] = -sint*m0[0] - cost*n0[0];
  nv[1] = -sint*m0[1] - cost*n0[1];
  nv[2] = -sint*m0[2] - cost*n0[2];
  a_mult_a(nv, ((shift_param_type *)param)->Cijkl, kCk);
  careful_inverse(kCk, EGF);
  sum -= EGF[d];

  return sum/sqrt(2+x);
}


void calc_shift(double Cijkl[9][9], double n_vect[3], double m_vect[3], 
		double shift[9]) 
{

  gsl_integration_qaws_table * t = 
    gsl_integration_qaws_table_alloc (0, -0.5, 1, 0);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1024);
  double result, error;
  gsl_function F;
  shift_param_type p;

  F.function = &shift_int;  
  F.params = (void *)(&p);
  for (int i=0; i<9; ++i)
    for (int j=0; j<9; ++j)
      p.Cijkl[i][j] = Cijkl[i][j];
  p.n_vect = n_vect;
  p.m_vect = m_vect;

  for (int d=0; d<9; ++d) {
    p.d = d;
    gsl_integration_qaws (&F, 0, 2, t, 1e-9, 1e-9, 1024, w, &result, &error);
    shift[d] = 0.5*M_1_PI*M_1_PI * result;
  }
  
  // garbage collection:
  gsl_integration_qaws_table_free(t);
  gsl_integration_workspace_free(w);
}


inline double min_magn (double a[3]) 
{
  double min=1e30;
  for (int i=0; i<3; ++i) 
    if (! zero(a[i]))
      if (fabs(a[i]) < min) min = fabs(a[i]);
  return min;
}

const double DIV_TOL = 1e-3;

double divisor (double a[3]) 
{
  int m, i, fail;
  double am, am_i;

  m=0; fail=1;
  while (fail) {
    ++m;
    fail=0;
    for (i=0; i<3; ++i) {
      am = a[i] * m;
      am_i = round(am);
      if ( fabs(am-am_i) > DIV_TOL) fail = 1;
    }
  }
  return (double)m;
}

void make_unit (double cart[9], double v[3], int u[3]) 
{
  double cart_inv[9], scale;
  double uinit[3];
  int i;
  
  careful_inverse(cart, cart_inv);
  mult_vect(cart_inv, v, uinit);
  
  scale = min_magn(uinit);
  for (i=0; i<3; ++i) uinit[i] *= 1./scale;
  scale = divisor(uinit);
  for (i=0; i<3; ++i) u[i] = lround(uinit[i]*scale);
}
