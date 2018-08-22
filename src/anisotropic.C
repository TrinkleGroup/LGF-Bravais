/*
  Program: anisotropic.C
  Author:  D. Trinkle
  Date:    August 14, 2003
  Purpose: Calculate the anisotropic elastic solution for a general
           dislocation given:
	   1. dislocation line vector t = |t| (t1, t2, t3)
	   2. burgers vector          b =     (b1, b2, b3)
	   3. dislocation cut vector  m =     (m1, m2, m3) (normalized)
	   4. elastic constants Cmn, crystal class c

	   Fixed to use the correct slip plane definition (important
	   for edge and mixed dislocations); that is, the slip plane
	   vector normal is:

	     n0 = t x b

	   unless it's zero; then n0 = t x m0.  Note: m0 is to be 
	   perpendicular to t and in the plane of n0.

  Param.:  <atomname> <cell> <infile>
           atomname: appended to each line of xyz files
           cell:     cell file (see below for format)
           infile:   input file (see below for format)

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
	   b1 b2 b3 bd  # burgers vector (unit cell coord.)/bd
	   m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)
	   c1 c2 c3 cd  # center of dislocation in unit cell ([c1 c2 c3]/cd)
	   c1' c2' c3'  # center of dislocation (shifts are added)
	   ==== infile ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   Read in everything, and just go *nuts* a calculatin'.

           First, we make sure that m0 is perp. to t and to b, and 
	   normalized.  We also construct n0 = t x m0.

	   We then define the vectors m(theta) and n(theta) as:

	   m(theta) =  m0*cos(theta) + n0*sin(theta)
	   n(theta) = -m0*sin(theta) + n0*cos(theta)

	   and the matrices (ab)_ij as

	   (ab)_ij = sum  a_k C_ikjl b_l
	              kl

           We have to do four integrals and store two of them as
	   functions of theta, namely, the two constant matrices:

	            1   Pi      -1
	   S_ij = - -  Int  (nn)  (nm)   dtheta
	            Pi  0       ik    kj

	            1     Pi                    -1
	   B_ij = -----  Int  (mm)   - (mn) (nn)  (nm)   dtheta
	          4Pi^2   0       ij      ik    kl    lj

	     (Note: B_ij = B_ji)

  	   and the two matrices as a function of theta:

	                    theta     -1
	   N_ij(theta) = 4Pi Int  (nn)   dtheta
	                      0       ij

	                theta     -1
	   L_ij(theta) = Int  (nn)  (nm)   dtheta
	                  0       ik    kj

	     (Note: N_ij = N_ji)

	   Also, due to the theta periodicity, we only have to evaluate
	   these from 0..Pi, since 

	     N_ij(theta+Pi) = N_ij(theta) + N_ij(Pi)

	   and similarly for L_ij.

	   Also, S_ij = -1/Pi * L_ij(Pi), so we have only three integrals
	   to do.

	   *Then*, to turn these all into our displacement using the 
	   equation:

	   u (|x|, theta) = [-S  ln |x| + N  B   + L  S  ] b  / 2 Pi
	    i                  is          ik ks    ik ks   s

	   Voila! (whew)

	   We integrate using a stepping scheme based on Simpson's
	   extended rule; basically, to integrate f(x) from 0 to x,
	   we calculate f(x) at a grid of points, and if F(N-1) is
	   the integral to x-h, and f[i] = f(x-ih), then:

   	     F(N) = F(N-1) + h*SUM(i=0..3, int_weight[i]*f[i])

	   which works very well.  To get the first two points, we 
	   actually have to evaluate f at -h and -2h, and then
	   start with F(0) = 0.  It works (go figure).

	   We do 16384 integration steps (2^14... woohoo!) to make
	   sure that we have something reasonable :)

  Output:  If we're verbose, we'll output the theta dependence of u_i, and
           also the ln |x| prefactor.

	   For fun, and profit, we can output the energy prefactor:

	   E = b_i B_ij b_j  : self-energy prefactor per length

	   The real meat of the code, though, is to actually put these
	   displacements to work by outputting the XYZ files for a
	   cylindrical slab material.  We do this by adding in the 
	   displacements... not too hard.
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

// 1 GPa = 1e9 J/m^3 = 10/1.60217733 meV/A^3 = 6.241506363 meV/A^3
// or 0.00624 eV/A^3.
// This makes the units of gij(R) A^3 / eV; multiplied by a force
// in eV/A and divided by the distance in A, gives the displacement
// in A.
const double meV_GPa = 10/1.60217733;  // = 10/(ec in coulumbs x 10^19)
const double eV_GPa = 0.01/1.60217733;  // = 0.01/(ec in coulumbs x 10^19)

//****************************** SUBROUTINES ****************************

inline double dot(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


void m_theta(double theta, double m[3], double n[3], double mt[3]) 
{
  mt[0] = cos(theta)*m[0] + sin(theta)*n[0];
  mt[1] = cos(theta)*m[1] + sin(theta)*n[1];
  mt[2] = cos(theta)*m[2] + sin(theta)*n[2];
}

void n_theta(double theta, double m[3], double n[3], double nt[3]) 
{
  nt[0] = -sin(theta)*m[0] + cos(theta)*n[0];
  nt[1] = -sin(theta)*m[1] + cos(theta)*n[1];
  nt[2] = -sin(theta)*m[2] + cos(theta)*n[2];
}

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

void print_mat (double a[9]) 
{
  int i, j;
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j)
      printf("  %10.5lf", a[index(i,j)]);
    printf("\n");
  }
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


void make_unit (double cart[9], double v[3], int u[3]);



/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <infile>";

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'g', 'r', 's'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:     cell file (-h for format)\n\
  infile:   input file (-h for format)\n\
  -g               don't convert from GPa to eV/A^3\n\
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
  int EV_CONV = !(flagon[0]);  // Convert from GPa to eV/A^3?
  int ROTATE = flagon[1];
  int SPIN = flagon[2];
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
  char *cell_name, *disl_name;
  FILE* infile;

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  double Cijkl[9][9];
  int Natoms;
  double** u_atoms;

  // Let's pull off the args:
  cell_name = args[0];
  disl_name = args[1];
  
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
  double t0[3], m0[3], n0[3];

  mult_vect(cart, t_unit, t0);
  mult_vect(cart, m_unit, m0);
  mult_vect(cart, n_unit, n0);

  // Normalize:
  double v_magn;
  v_magn = sqrt(dot(t0, t0));
  for (d=0; d<3; ++d) t0[d] /= v_magn;
  // make sure m is normal to t:
  v_magn = dot(m0, t0);
  for (d=0; d<3; ++d) m0[d] -= v_magn * t0[d];
  // normalize
  v_magn = sqrt(dot(m0, m0));
  if (zero(v_magn)) {
    fprintf(stderr, "t and m are parallel...?  t_unit=(%d,%d,%d)  m_unit=(%d,%d,%d)\n",
            t_unit[0],t_unit[1],t_unit[2], m_unit[0],m_unit[1],m_unit[2]);
    ERROR = ERROR_BADFILE;
    exit(ERROR);
  }
  else 
    for (d=0; d<3; ++d) m0[d] /= v_magn;
  // make n as the cross product of t and m:
  crossprod(t0, m0, n0);
  
  int t_new[3], m_new[3], n_new[3];
  make_unit(cart, t0, t_new);
  make_unit(cart, m0, m_new);
  make_unit(cart, n0, n_new);

  if (TESTING) {
    printf("## vect    input        corrected             cartesian\n");
    printf("## t = (%3d %3d %3d)->(%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
           t_unit[0],t_unit[1],t_unit[2], t_new[0],t_new[1],t_new[2], 
           t0[0],t0[1],t0[2]);
    printf("## m = (%3d %3d %3d)->(%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
           m_unit[0],m_unit[1],m_unit[2], m_new[0],m_new[1],m_new[2],
           m0[0],m0[1],m0[2]);
    printf("## n = (%3d %3d %3d)->(%3d %3d %3d) %12.8lf %12.8lf %12.8lf\n", 
           n_unit[0],n_unit[1],n_unit[2], n_new[0],n_new[1],n_new[2], 
           n0[0],n0[1],n0[2]);
  }


  // Now some evaluating of integrals :)
  double theta;
  double dtheta;
  dtheta = M_PI / Nsteps;
  double mt[3], nt[3];
  double nnt[9], mmt[9], nmt[9], mnt[9], nnti[9];

  // We have to integrate three functions.
  double **Nint, **Lint;
  double Bint[9];

  Nint = new double*[Nsteps+1];
  Lint = new double*[Nsteps+1];
  for (i=0; i<=Nsteps; ++i) {
    Nint[i] = new double[9];
    Lint[i] = new double[9];
  }
  
  // Function evaluations, stored for integration purposes.
  double nn_old[4][9], nnnm_old[4][9], mnnnnm_old[4][9];

  // First, prime the integration pump:
  for (k=1; k<=3; ++k) {
    theta = -(k-1)*dtheta;
    // Eval (nn), (nm), (mn), (mm), and (nn)^-1
    m_theta(theta, m0, n0, mt);
    n_theta(theta, m0, n0, nt);
    a_mult_a(mt, Cijkl, mmt);
    a_mult_a(nt, Cijkl, nnt);
    a_mult_b(nt, mt, Cijkl, nmt);
    a_mult_b(mt, nt, Cijkl, mnt);
    careful_inverse(nnt, nnti);
    // Now, put into the function evaluations:
    for (i=0; i<9; ++i) nn_old[k][i] = nnti[i];
    mult(nnti, nmt, nnnm_old[k]);
    mult(mnt, nnnm_old[k], mnnnnm_old[k]);
    for (i=0; i<9; ++i)
      mnnnnm_old[k][i] = mmt[i] - mnnnnm_old[k][i];
    // And HERE's where we'd integrate, if we wanted to... :)
  }

  // Now we've got EVERYTHING, let's integrate!
  // theta = 0 is easy...
  for (i=0; i<9; ++i) {
    Nint[0][i] = 0.;
    Lint[0][i] = 0.;
    Bint[i] = 0.;
  }

  for (k=1; k<=Nsteps; ++k) {
    theta = k*dtheta;
    // Eval (nn), (nm), (mn), (mm), and (nn)^-1
    m_theta(theta, m0, n0, mt);
    n_theta(theta, m0, n0, nt);
    a_mult_a(mt, Cijkl, mmt);
    a_mult_a(nt, Cijkl, nnt);
    a_mult_b(nt, mt, Cijkl, nmt);
    a_mult_b(mt, nt, Cijkl, mnt);
    careful_inverse(nnt, nnti);

    // Now, put into the function evaluations:
    for (i=0; i<9; ++i) nn_old[0][i] = nnti[i];
    mult(nnti, nmt, nnnm_old[0]);
    mult(mnt, nnnm_old[0], mnnnnm_old[0]);
    for (i=0; i<9; ++i)
      mnnnnm_old[0][i] = mmt[i] - mnnnnm_old[0][i];
    // Now, we can integrate!
    for (i=0; i<9; ++i) {
      Nint[k][i] = Nint[k-1][i];
      Lint[k][i] = Lint[k-1][i];
      for (j=0; j<4; ++j) {
	Nint[k][i] += dtheta*int_weight[j]*nn_old[j][i];
	Lint[k][i] += dtheta*int_weight[j]*nnnm_old[j][i];
	Bint[i]    += dtheta*int_weight[j]*mnnnnm_old[j][i];
      }
    }
    // Now, we slide down all of our "old" values:
    for (j=3; j>0; --j)
      for (i=0; i<9; ++i) {
	nn_old[j][i] = nn_old[j-1][i];
	nnnm_old[j][i] = nnnm_old[j-1][i];
	mnnnnm_old[j][i] = mnnnnm_old[j-1][i];
      }
    // And do it all again!
  }
  // Finally, define S and Q, and scale everything appropriately.
  double Sint[9], Qint[9];

  for (k=0; k<=Nsteps; ++k)
    for (i=0; i<9; ++i) 
      Nint[k][i] *= (4.*M_PI);
  
  for (i=0; i<9; ++i) {
    Sint[i] = -Lint[Nsteps][i] * M_1_PI;
    Qint[i] = -0.25*M_1_PI*M_1_PI * Nint[Nsteps][i];
    Bint[i] *= 0.25*M_1_PI*M_1_PI;
  }      
  
  // ****************************** OUTPUT ***************************

  // Human readable (sorta) first:
  double LQ[9], NS[9], ST[9], Gw[9], Gw0[9];
  double Gii[6];

  if (!SPIN) {
    m_theta(0, m0, n0, mt);
    n_theta(0, m0, n0, nt);
  }
  
  transpose(Sint, ST);
  for (k=0; k<=Nsteps; ++k) {
    theta = k*dtheta;
    mult(Lint[k], Qint, LQ);
    mult(Nint[k], ST, NS);
    for (d=0; d<9; ++d) Gw[d] = -LQ[d] - 0.25*M_1_PI*NS[d];
    mult(Gw, 0.5*M_1_PI, Gw);

    //++ totally changing the definition here...
    //    for (d=0; d<9; ++d) Gw[d] = Nint[k][d];
    //-- totally changing the definition here...

    if (SPIN) {
      m_theta(theta, m0, n0, mt);
      n_theta(theta, m0, n0, nt);
    }
    if (ROTATE) rotate(Gw, mt, nt, t0, Gii);
    else        rotate(Gw, ident, ident+3, ident+6, Gii);
    printf("%.8lf", theta);
    for (d=0; d<6; ++d) printf(" %.8le", Gii[d]);
    printf("\n");
  }
  for (d=0; d<9; ++d) Gw0[d] = Gw[d];
  for (k=1; k<=Nsteps; ++k) {
    theta = (k+Nsteps)*dtheta;
    mult(Lint[k], Qint, LQ);
    mult(Nint[k], ST, NS);
    for (d=0; d<9; ++d) Gw[d] = Gw0[d] - LQ[d] - 0.25*M_1_PI*NS[d];
    mult(Gw, 0.5*M_1_PI, Gw);

    //++ totally changing the definition here...
    //    for (d=0; d<9; ++d) Gw[d] = Gw0[d] + Nint[k][d];
    //-- totally changing the definition here...

    if (SPIN) {
      m_theta(theta, m0, n0, mt);
      n_theta(theta, m0, n0, nt);
    }
    if (ROTATE) rotate(Gw, mt, nt, t0, Gii);
    else        rotate(Gw, ident, ident+3, ident+6, Gii);
    printf("%.8lf", theta);
    for (d=0; d<6; ++d) printf(" %.8le", Gii[d]);
    printf("\n");
  }
  printf("&\n");

  mult(Qint, 0.5*M_1_PI, Qint);
  if (ROTATE) rotate(Qint, m0, n0, t0, Gii);
  else        rotate(Qint, ident, ident+3, ident+6, Gii);
  printf("0");
  for (d=0; d<6; ++d) printf(" %.8le", Gii[d]);
  printf("\n");


  // ************************* GARBAGE COLLECTION ********************
  for (i=0; i<=Nsteps; ++i) {
    delete[] Nint[i];
    delete[] Lint[i];
  }
  delete[] Nint;
  delete[] Lint;

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
