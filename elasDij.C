/*
  Program: elasDij.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a dynamical matrix, and calculate the elastic constants
           using Born and Huang's expression for the case of a single
	   atom Bravais lattice:

                        1
	   [ab,cd] = - --- SUM D   x x
                       2Vc  x   ab  c d

	   then

	   C_abcd = [ac,bd] + [bc,ad] - [cd,ab]

           There are also 15 equivalence relations:

	   [ab,cd] = [cd,ab]

	   that must be satisfied if the dynamical matrix obeys rotational
	   invariance _and_ if the lattice is stress-free.  For a lattice
	   obeying *cubic* symmetry, I *believe* these equivalence
	   relations are trivially obeyed if Dij has cubic symmetry.

  Param.:  <cell> <Dij(R)>
           cell:   cell file describing our lattice
	   Dij(R): dynamical matrix

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
#include <gsl/gsl_eigen.h>  // for determining stability
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "eigen.H" // the "other" eigenroutines--for anisotropy
#include "cell.H"
#include "voigt.H"
#include "pointgroup.H"
#include "Dij.H"
#include "lambda.H" // just use the lambda_k routine


//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <Dij(R)>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'s', 'a'}; // flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  Dij(R): dynamical matrix\n\
  -s      calculate stability (output 1 if stable, 0 if not)\n\
  -a      determine the anisotropy ratio (MEMORY determines grid to use)";

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
==== Dij(R) ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 8;   // grid to use for checking anisotropy

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

  // flags
  int STABILITY = flagon[0];
  int ANISOTROPY = flagon[1];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Dij_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Dij_name = args[1];

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
  int Np, Ndim = 9;
  point_type *p; // set of points

  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, p);
  if (ERROR) {delete p; exit(ERROR);}
  //-- ==== Dij(R) ====


  // ***************************** ANALYSIS **************************
  double bracket[9][9];
  double vol;
  
  vol = 1./det(cart);
  // Use Voigt notation to generate the symmetrized elements:
  int ab, cd;
  int a, b, c, d;
  for (cd=0; cd<6; ++cd) {
    double brack[9];
    k2ij(cd, c, d);
    j = index(c,d);
    i = index(d,c);
    for (ab=0; ab<9; ++ab) brack[ab] = 0;
    for (np=0; np<Np; ++np) 
      for (ab=0; ab<9; ++ab)
	brack[ab] += (p[np].mat[ab]) * (p[np].Rcart[c]) * (p[np].Rcart[d]);
    for (ab=0; ab<9; ++ab) bracket[ab][j] = -0.5*vol*brack[ab];
    if (i!=j)
      for (ab=0; ab<9; ++ab) bracket[ab][i] = bracket[ab][j];
  }
  
  if (TESTING) {
    for (a=0; a<3; ++a)
      for (b=0; b<3; ++b) {
	ab = index(a,b);
	for (c=0; c<3; ++c)
	  for (d=0; d<3; ++d) {
	    cd = index(c,d);
	    printf("## [%c%c,%c%c] = %.12le\n",
		   coord[a], coord[b], coord[c], coord[d], bracket[ab][cd]);
	  }
      }
  }
  
  // First, test the invariance relations:
  int inv_test[6][6];
  for (ab=0; ab<6; ++ab) {
    k2ij(ab, a, b);
    i = index(a,b);
    for (cd=0; cd<ab; ++cd) inv_test[ab][cd]=inv_test[cd][ab];
    inv_test[ab][ab] = 1;
    for (cd=(ab+1); cd<6; ++cd) {
      k2ij(cd, c, d);
      j = index(c,d);
      inv_test[ab][cd] = dcomp(bracket[i][j], bracket[j][i]);
      if (! inv_test[ab][cd]) {
	fprintf(stderr, "Failed invariance test [%c%c,%c%c]=[%c%c,%c%c]\n",
		coord[a], coord[b], coord[c], coord[d],
		coord[c], coord[d], coord[a], coord[b]);
	fprintf(stderr, "%.12le != %.12le\n", bracket[i][j], bracket[j][i]);
      }
    }
  }
  
  if (VERBOSE) {
    printf("# Invariance tests:\n");
    for (ab=0; ab<6; ++ab) {
      k2ij(ab, a, b);
      i = index(a,b);
      for (cd=(ab+1); cd<6; ++cd) {
	k2ij(cd, c, d);
	j = index(c,d);
	printf("# [%c%c,%c%c] = %.12le  [%c%c,%c%c] = %.12le",
	       coord[a], coord[b], coord[c], coord[d], bracket[i][j],
	       coord[c], coord[d], coord[a], coord[b], bracket[j][i]);
	if (inv_test[ab][cd]) printf("\n");
	else                  printf(" --FAILED!\n");
      }
    }
  }
  

  // Now, calculate our elastic constants:
  double Cij[6][6];
  for (ab=0; ab<6; ++ab) {
    k2ij(ab, a, b);
    for (cd=0; cd<ab; ++cd) 
      Cij[ab][cd] = Cij[cd][ab];
    for (cd=ab; cd<6; ++cd) {
      k2ij(cd, c, d);
      Cij[ab][cd] = bracket[index(a,c)][index(b,d)]
	+ bracket[index(b,c)][index(a,d)] - bracket[index(c,d)][index(a,b)];
    }
  }

  if (VERBOSE) {
    printf("# Elastic constants:\n");
    for (ab=0; ab<6; ++ab) {
      k2ij(ab, a, b);
      for (cd=0; cd<6; ++cd) {
	k2ij(cd, c, d);
	printf("# C_%c%c%c%c = %.12le\n",
	       coord[a], coord[b], coord[c], coord[d], Cij[ab][cd]);
      }
    }
  }
  
  // Now, let's test our Cij against the declared crystal class; then
  // we can output the Cij line correctly
  double Cabcd[9][9], Cabcd_class[9][9];
  // Note: Cmn_list has already been declared, so we'll simply
  // overwrite the values.
  // First, populate Cabcd accordingly:
  for (i=0; i<9; ++i) {
    ab = ij2k(i/3, i%3);
    for (j=0; j<9; ++j) {
      cd = ij2k(j/3, j%3);
      Cabcd[i][j] = Cij[ab][cd];
    }
  }
  // Next, populate Cmn_list accordingly:
  for (i=0; i<class_len[crystal]; ++i) {
    j = class_Cij[crystal][i];
    ab = (j/10) - 1;
    cd = (j%10) - 1;
    Cmn_list[i] = Cij[ab][cd];
  }
  // Now, make Cabcd according to Cmn_list:
  make_Cijkl(crystal, Cmn_list, Cabcd_class);
  // Finally, test:
  for (i=0; (i<9) && (!ERROR); ++i)
    for (j=0; (j<9) && (!ERROR); ++j) 
      ERROR = !(dcomp(Cabcd[i][j], Cabcd_class[i][j]));
  if (ERROR) {
    fprintf(stderr, "Error--crystal class %d not correct?\n", crystal);
    fprintf(stderr, "Disagreement at C_%c%c%c%c:\n",
	    coord[i/3], coord[i%3], coord[j/3], coord[j%3]);
    fprintf(stderr, "Calculated %.12lf, should be %.12lf\n",
	    Cabcd[i][j], Cabcd_class[i][j]);
  }

  // Output crystal class and Cij info
  printf("%d", crystal);
  for (i=0; i<class_len[crystal]; ++i) printf(" %.12le", Cmn_list[i]);
  printf(" # class");
  for (i=0; i<class_len[crystal]; ++i) printf(" C%2d", class_Cij[crystal][i]);
  printf("\n");

  // Calculate stability?
  if (STABILITY && (!ERROR)) {
    // We do this by ensuring that Cij is positive definite, by calculating
    // the eigenvalues and making sure they're all >0.
    double data[36];
    for (ab=0; ab<6; ++ab)
      for (cd=0; cd<6; ++cd)
	data[ab*6+cd] = Cij[ab][cd];
    gsl_matrix_view m = gsl_matrix_view_array(data, 6, 6);
    gsl_vector *eval = gsl_vector_alloc (6);
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(6);
    gsl_eigen_symm(&m.matrix, eval, w);
    gsl_eigen_symm_free(w);
    int stable=-1;
    for (i=0; (i<6) && (stable); ++i)
      stable = (gsl_vector_get(eval,i) > 0.);

    if (stable) printf("1  # stable elastic constants\n");
    else        printf("0  # unstable elastic constants\n");
  }

  // Calculate approximate anisotropic ratio?
  if (ANISOTROPY && (!ERROR)) {
    // remember: bracket == lambda
    int l2, m2, n2;
    double kpt[3], invN = 1./MEMORY;
    double lambdak[9], eig[3];
    double mineig, maxeig;
    mineig = 1e300; maxeig = -1e300;
    for (l2=0; (l2<=MEMORY) && (mineig > 0.); ++l2) {
      kpt[0] = sqrt(l2*invN);
      for (m2=0; (m2<=(MEMORY-l2)) && (mineig > 0.); ++m2) {
	kpt[1] = sqrt(m2*invN);
	n2 = MEMORY-l2-m2;
	kpt[2] = sqrt(n2*invN);
	lambda_k(kpt, bracket, lambdak);
	eigen(lambdak, eig); // these are already sorted, too.
	if (eig[0] < mineig) mineig = eig[0];
	if (eig[2] > maxeig) maxeig = eig[2];
      }
    }
    double ani = 0.;
    if (! zero(mineig)) ani = maxeig/mineig;
    printf("%.5lf # anisotropy ratio= %.8le / %.8le\n", ani,
	   maxeig, mineig);
    if (mineig > 0.) printf("1  # stable elastic phonons\n");
    else             printf("0  # unstable elastic phonons\n");
  }
  

  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
