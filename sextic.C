/*
  Program: sextic.C
  Author:  D. Trinkle
  Date:    August 27, 2004
  Purpose: Calculate the sextic roots and eigenvectors for a given
           crystal relative to a set of axes.  The roots, eigenvectors,
	   and axes are output to be used by further routines.

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
#include <vecLib/clapack.h>  // Lapack routines, from vecLib
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

//****************************** SUBROUTINES ****************************

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
  if (zero(x)) printf("0");
  else         printf("%.8le", x);
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
const char USERFLAGLIST[NFLAGS] = {'g', 'o'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:  cell file (see below for format)\n\
  disl:  dislocation coordinates in unit cell coord.\n\
  -g               don't convert from GPa to eV/A^3\n\
  -o               print maple friendly output";

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
  int a, b, i, j, k, d; // General counting variables.

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
  int EV_CONV = !(flagon[0]);  // Convert from GPa to eV/A^3?
  int MAPLE = flagon[1];
  
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

  
  // 1. Construct our sextic problem to solve as an eigenproblem.
  double* Nmat = new double[6*6];
  double nn[9], nn_inv[9];
  double nm[9], mn[9], mm[9];

  a_mult_a(n_vect, Cijkl, nn);
  a_mult_a(m_vect, Cijkl, mm);
  careful_inverse(nn, nn_inv);
  a_mult_b(n_vect, m_vect, Cijkl, nm);
  a_mult_b(m_vect, n_vect, Cijkl, mn); // is (mn) = (nm)^T ?

  double nn_nm[9], mn_nn[9], mn_nn_nm[9];
  mult(nn_inv, nm, nn_nm);
  mult(mn, nn_nm, mn_nn_nm);
  for (d=0; d<9; ++d) mn_nn_nm[d] -= mm[d];
  mult(mn, nn_inv, mn_nn);

  if (TESTING) {
    printf("## (nn):\n"); print_mat(nn, "## ");
    printf("## (mm):\n"); print_mat(mm, "## ");
    printf("## (nn)^-1:\n"); print_mat(nn_inv, "## ");
  }
  
  // finally, put together Nmat:
  for (i=0; i<3; ++i)
    for (j=0; j<3; ++j) {
      Nmat[(i+0)*6 + j+0] = nn_nm[i*3+j];
      Nmat[(i+0)*6 + j+3] = nn_inv[i*3+j];
      Nmat[(i+3)*6 + j+0] = mn_nn_nm[i*3+j];
      Nmat[(i+3)*6 + j+3] = mn_nn[i*3+j];
    }
  // negate...
  for (d=0; d<36; ++d) Nmat[d] = -Nmat[d];
  
  if (TESTING) {
    printf("## Nmat:\n");
    for (i=0; i<6; ++i) {
      printf("##");
      for (j=0; j<6; ++j) printf(" %12.8lf", Nmat[i*6+j]);
      printf("\n");
    }
  }

  // 2. calculate eigenvectors and values
  char jobvl = 'N'; // calc left eigenvectors (see below)
  char jobvr = 'V'; // no right eigenvectors
  long int N = 6;        // order of our matrix
  long int lda = 6;      // leading dimension
  double* wr = new double[N]; // real and imag comp of eigenvalues
  double* wi = new double[N];
  double* vl = new double[N*N]; // left eigenvectors
  long int ldvl = N;            // left eigenvector leading dim
  double* vr = new double[N*N]; // right eigenvectors
  long int ldvr = N;            // right eigenvector leading dim
  long int lwork = N*N;                  // workspace size
  double* work = new double[lwork]; // workspace
  long int info;

  // NOW: the lapack routine uses fortran's braindead definition for
  // column vs. row.  Which means amat is the "transpose" of Nmat...
  double* amat = new double[N*N];
  for (d=0; d<(N*N); ++d) 
    amat[(d/N)*N + d%N] = Nmat[(d%N)*N + d/N];

  // call that routine!
  dgeev_(&jobvl, &jobvr, &N, amat, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
	 work, &lwork, &info);
  delete[] amat;
  delete[] vl;

  // 2.b. Store eigenvectors, etc.:
  double eig_r[N][N], eig_i[N][N];
  double pr[N], pi[N];
  // A_alpha,k = eig_r[alpha=0..5][k=0..2] + I eig_i[alpha][k=0..2]
  // L_alpha,k = eig_r[alpha=0..5][k+3]    + I eig_i[alpha][k+3]
  for (a=0; a<3; ++a) {
    pr[a  ] = wr[2*a];
    pi[a  ] = wi[2*a];
    pr[a+3] = wr[2*a+1];
    pi[a+3] = wi[2*a+1];
    for (k=0; k<N; ++k) {
      eig_r[a  ][k] =  vr[N* 2*a    + k];
      eig_i[a  ][k] =  vr[N*(2*a+1) + k];
      eig_r[a+3][k] =  vr[N* 2*a    + k];
      eig_i[a+3][k] = -vr[N*(2*a+1) + k];
    }
  }
  if ( zero(pr[0]) && (zero(pr[1])) && (zero(pr[2])) ) {
    fprintf(stderr, 
	    "ACK--you've chosen the isotropic case... this is not handled correctly.\n");
  }
      
  for (a=0; a<N; ++a) {
    // Do the normalization:
    double scale_r=0, scale_i=0;
    for (k=0; k<3; ++k) {
      scale_r += 2*(eig_r[a][k]*eig_r[a][k+3] - eig_i[a][k]*eig_i[a][k+3]);
      scale_i += 2*(eig_r[a][k]*eig_i[a][k+3] + eig_i[a][k]*eig_r[a][k+3]);
    }
    if (zero(scale_r) && zero(scale_i)) {
      scale_r = 1;
      scale_i = 0;
    }
    // Now, this is supposed to be 1... so we need to divide
    // by the sqrt of this.
    double theta = atan2(scale_i, scale_r);
    double magn = hypot(scale_r, scale_i);
    magn = 1./sqrt(magn);
    theta *= -0.5;
    scale_r = magn * cos(theta);
    scale_i = magn * sin(theta);
    for (k=0; k<N; ++k) {
      double temp = eig_r[a][k];
      eig_r[a][k] = scale_r*eig_r[a][k] - scale_i*eig_i[a][k];
      eig_i[a][k] = scale_i*temp + scale_r*eig_i[a][k];
    }
  }

  if (TESTING) {
    printf("## dgeev summary:\n");
    printf("## info= %ld\n", info);
    printf("## eigenvalue pairs:\n");
    for (a=0; a<3; ++a)
      printf("## %.15lf%+.15lfI  %.15lf%+.15lfI\n",
	     pr[a], pi[a], pr[a+3], pi[a+3]);
    printf("## eigenvector testing:\n");
  }
  if (TESTING || MAPLE) {
    for (a=0; a<N; ++a) {
      // First, construct the matrix problem for A:
      double tpr = pr[a];
      double tpi = pi[a];
      double CmatR[9], CmatI[9];
      for (d=0; d<9; ++d) {
	CmatR[d] = mm[d] + (tpr*tpr-tpi*tpi)*nn[d] + tpr*(mn[d]+nm[d]);
	CmatI[d] = 2*tpr*tpi*nn[d] + tpi*(mn[d]+nm[d]);
      }

      double Nei_r[N], Nei_i[N];
      // Now make the Ar and Ai matrices:
      double Ar[9], Ai[9];
      for (int i=0; i<3; ++i) 
	for (int j=0; j<3; ++j) {
	  d = 3*i+j;
	  Ar[d] = eig_r[a][i]*eig_r[a][j] - eig_i[a][i]*eig_i[a][j];
	  Ai[d] = eig_r[a][i]*eig_i[a][j] + eig_i[a][i]*eig_r[a][j];
	}

      for (j=0; j<N; ++j) {
	Nei_r[j] = 0;
	Nei_i[j] = 0;
	for (k=0; k<N; ++k) {
	  Nei_r[j] += Nmat[N*j+k]*eig_r[a][k];
	  Nei_i[j] += Nmat[N*j+k]*eig_i[a][k];
	}
	// Now, subtract off the eigenvalue * eig
	Nei_r[j] -= (tpr*eig_r[a][j] - tpi*eig_i[a][j]);
	Nei_i[j] -= (tpi*eig_r[a][j] + tpr*eig_i[a][j]);
      }
      double rms=0;
      for (j=0; j<N; ++j)
	rms += Nei_r[j]*Nei_r[j] + Nei_i[j]*Nei_i[j];
      rms = sqrt(rms/N);
      
      // print out how well we did:
      if (TESTING) {
	printf("## p_%d= %.8lf%+.8lfI\n", a+1, tpr, tpi);
	printf("## Cmat=\n");
	for (int i=0; i<3; ++i) {
	  printf("##");
	  for (int j=0; j<3; ++j) printf(" %.4lf%+.4lfI", 
					 CmatR[3*i+j], CmatI[3*i+j]);
	  printf("\n");
	}
	printf("## L = ( )*A:\n");
	for (int i=0; i<3; ++i) {
	  printf("##");
	  for (int j=0; j<3; ++j) printf(" %.4lf%+.4lfI", 
					 -nm[3*i+j]-tpr*nn[3*i+j], 
					 -tpi*nn[3*i+j]);
	  printf("\n");
	}      
	printf("## A_%d= ", a+1);
	for (j=0; j<3; ++j)
	  printf("  %10.6lf%+10.6lfI", eig_r[a][j], eig_i[a][j]);
	printf("\n## L_%d= ", a+1);
	for (j=0; j<3; ++j)
	  printf("  %10.6lf%+10.6lfI", eig_r[a][j+3], eig_i[a][j+3]);
	printf("\n## |N*AL - p*AL|= %.5le\n", rms); // remainder
	printf("## A_%dr and A_%di:\n", a+1, a+1);
	for (int i=0; i<3; ++i) {
	  printf("##");
	  for (int j=0; j<3; ++j) printf(" %7.4lf", Ar[3*i+j]);
	  printf("  ");
	  for (int j=0; j<3; ++j) printf(" %7.4lf", Ai[3*i+j]);
	  printf("\n");
	}
	printf("##\n");
      }
      if (MAPLE) {
	// Do maple input:
	printf("p[%d] := %.8le %+.8le * I;\n", a+1, tpr, tpi);
	sprintf(dump, "AR[%d]", a+1);
	maple_print_mat(Ar, dump);
	sprintf(dump, "AI[%d]", a+1);
	maple_print_mat(Ai, dump);
      }
    }
  }
  if (TESTING) {
      
    // Normalization testing:
    printf("## sum rules:\n");
    printf("## sum_i A_a,i L_b,i + L_a,i A_b,i = del_ab\n");
    double sumr, sumi;
    for (a=0; a<N; ++a) {
      printf("##");
      for (b=0; b<N; ++b) {
	sumr = 0;
	sumi = 0;
	for (i=0; i<3; ++i) {
	  sumr += eig_r[a][i]*eig_r[b][i+3] - eig_i[a][i]*eig_i[b][i+3];
	  sumr += eig_r[b][i]*eig_r[a][i+3] - eig_i[b][i]*eig_i[a][i+3];

	  sumi += eig_r[a][i]*eig_i[b][i+3] + eig_i[a][i]*eig_r[b][i+3];
	  sumi += eig_r[b][i]*eig_i[a][i+3] + eig_i[b][i]*eig_r[a][i+3];
	}
	if (zero(sumr)) sumr = 0;
	if (zero(sumi)) sumi = 0;
	printf(" %+9.2le%+9.2leI", sumr, sumi);
      }
      printf("\n");
    }
    printf("## sum_a A_a,i A_a,j = 0\n");
    for (i=0; i<3; ++i) {
      printf("##");
      for (j=0; j<3; ++j) {
	sumr = 0;
	sumi = 0;
	for (a=0; a<N; ++a) {
	  sumr += eig_r[a][i]*eig_r[a][j] - eig_i[a][i]*eig_i[a][j];
	  sumi += eig_r[a][i]*eig_i[a][j] + eig_i[a][i]*eig_r[a][j];
	}
	if (zero(sumr)) sumr = 0;
	if (zero(sumi)) sumi = 0;
	printf(" %+9.2le%+9.2leI", sumr, sumi);
      }
      printf("\n");
    }
    printf("## sum_a L_a,i L_a,j = 0\n");
    for (i=0; i<3; ++i) {
      printf("##");
      for (j=0; j<3; ++j) {
	sumr = 0;
	sumi = 0;
	for (a=0; a<N; ++a) {
	  sumr += eig_r[a][i+3]*eig_r[a][j+3] - eig_i[a][i+3]*eig_i[a][j+3];
	  sumi += eig_r[a][i+3]*eig_i[a][j+3] + eig_i[a][i+3]*eig_r[a][j+3];
	}
	if (zero(sumr)) sumr = 0;
	if (zero(sumi)) sumi = 0;
	printf(" %+9.2le%+9.2leI", sumr, sumi);
      }
      printf("\n");
    }
    printf("## sum_a A_a,i L_a,j = del_ij\n");
    for (i=0; i<3; ++i) {
      printf("##");
      for (j=0; j<3; ++j) {
	sumr = 0;
	sumi = 0;
	for (a=0; a<N; ++a) {
	  sumr += eig_r[a][i]*eig_r[a][j+3] - eig_i[a][i]*eig_i[a][j+3];
	  sumi += eig_r[a][i]*eig_i[a][j+3] + eig_i[a][i]*eig_r[a][j+3];
	}
	if (zero(sumr)) sumr = 0;
	if (zero(sumi)) sumi = 0;
	printf(" %+9.2le%+9.2leI", sumr, sumi);
      }
      printf("\n");
    }

  }
  
  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  if (VERBOSE) {
    // What, exactly, is "verbose" output...?
  }

  printf("%18.15lf %18.15lf %18.15lf %2d %2d %2d # t (unit cart)\n",
	 t_vect[0],t_vect[1],t_vect[2], t_new[0],t_new[1],t_new[2]);
  printf("%18.15lf %18.15lf %18.15lf %2d %2d %2d # m (unit cart)\n",
	 m_vect[0],m_vect[1],m_vect[2], m_new[0],m_new[1],m_new[2]);
  printf("%18.15lf %18.15lf %18.15lf %2d %2d %2d # n (unit cart)\n",
	 n_vect[0],n_vect[1],n_vect[2], n_new[0],n_new[1],n_new[2]);
  for (a=0; a<3; ++a) {
    printf("%.15lf %.15lf  # Re(p%d)  Im(p%d)\n", pr[a], pi[a], a+1, a+1);
    for (i=0; i<3; ++i)
      printf(" %.15le %.15le", eig_r[a][i], eig_i[a][i]);
    printf(" # Re(A%d,1) Im(A%d,1) ...\n", a+1, a+1);
    for (i=0; i<3; ++i)
      printf(" %.15le %.15le", eig_r[a][i+3], eig_i[a][i+3]);
    printf(" # Re(L%d,1) Im(L%d,1) ...\n", a+1, a+1);
  }

  // ************************* GARBAGE COLLECTION ********************
  delete[] wr;
  delete[] wi;
  delete[] vr;
  delete[] work;
  
  delete[] Nmat;
  free_cell(Cmn_list, u_atoms, Natoms);

  return 0;
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
