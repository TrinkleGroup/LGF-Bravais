/*
  Program: super-poscar.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. a periodic supercell with points and a GF
	   3. a test force in eV/A

	   and we output 
	   1. POSCAR input file with all of the positions (displaced) and
	   the applied forces

  Param.:  <cell> <period> <fx> <fy> <fz> <Rcut>
           cell:    cell file describing our lattice
           period:  supercell and periodic GF
	   fxyz:    test force, in eV/A
	   Rcut:    cutoff force freezing atoms in A

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

	   ==== period ====
	   n11 n12 n13 ... n33  # supercell definition
	   N  # number of points
	   n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG
	   ...
	   nN.1 nN.2 nN.3  Gxx .. Gzz
	   ==== period ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Currently just testing the gridding of k-points

  Output:  A correct POSCAR file.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "cont-int.H" // For doing the integration over the continuous piece
#include "sphere-harm.H" // Spherical harmonic evaluation
#include "supercell.H" // supercell 

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 6;
const char* ARGLIST = "<cell> <period> <fx> <fy> <fz> <Rcut>";


const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'f'}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  period:  supercell definition, points, and periodic GF\n\
  fxyz:    test force in eV/A\n\
  Rcut:    cutoff force freezing atoms in A (=0 -- no cutoff)\n\
  -f       fold down into supercell";

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
==== period ====\n\
n11 n12 n13 ... n33  # supercell definition\n\
N  # number of points\n\
n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG\n\
...\n\
nN.1 nN.2 nN.3  Gxx .. Gzz\n\
==== period ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter
  int FOLDDOWN;

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
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  FOLDDOWN = (flagon[0]); // default: no folddown.  

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *period_name;
  double testf[3];
  double Rcut;
  
  // Let's pull off the args:
  cell_name = args[0];
  period_name = args[1];
  // Get the test force:
  for (i=0; i<3; ++i)
    sscanf(args[i+2], "%lf", &(testf[i]));
  sscanf(args[5], "%lf", &Rcut);
  
  double cart[9];
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
  }
  //-- ==== cell ====

  //++ ==== period ====
  int super_n[9];
  int Nsuper;
  int Ndim = 9;
  int **Runit;
  double **Gperiod;
  
  infile = myopenr(period_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", period_name);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d %d %d %d %d %d %d",
	 &(super_n[0]), &(super_n[1]), &(super_n[2]), 
	 &(super_n[3]), &(super_n[4]), &(super_n[5]), 
	 &(super_n[6]), &(super_n[7]), &(super_n[8]));

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Nsuper);
  if (Nsuper != abs(det(super_n))) {
    fprintf(stderr, "Nsuper value (%d) doesn't match supercell?\n", Nsuper);
    ERROR = ERROR_BADFILE;
  }
  if (Nsuper < 1) {
    fprintf(stderr, "Bad Nsuper (%d) value in %s\n", Nsuper, period_name);
    ERROR = ERROR_BADFILE;
  }
  if (!ERROR) {
    Runit = new int *[Nsuper];
    Gperiod = new double *[Nsuper];
  }
  
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Nsuper); ++n) {
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    Runit[n] = new int[3];
    Gperiod[n] = new double[Ndim];
    sscanf(dump, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	   &(Runit[n][0]), &(Runit[n][1]), &(Runit[n][2]),
	   &(Gperiod[n][0]), &(Gperiod[n][1]), &(Gperiod[n][2]),
	   &(Gperiod[n][3]), &(Gperiod[n][4]), &(Gperiod[n][5]),
	   &(Gperiod[n][6]), &(Gperiod[n][7]), &(Gperiod[n][8]));
  }
  myclose(infile);
  if (n!=Nsuper) {
    fprintf(stderr, "Didn't read enough points in %s\n", period_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== Gperiod ====

  double cartA[9];
  mult(cart, super_n, cartA);

  if (TESTING) {
    printf("## Supercell: %d atoms\n", Nsuper);
    for (i=0; i<3; ++i) {
      printf("## A%d =", i+1);
      for (j=0; j<3; ++j)
	printf(" %9.5lf", cartA[i+3*j]);
      printf("\n");
    }
    printf("## Periodic GF:\n");
    for (np=0; np<Nsuper; ++np) {
      printf("## %3d %3d %3d  ", Runit[np][0],  Runit[np][1],  Runit[np][2]);
      print_mat(Gperiod[np]);
      printf("\n");
    }
  }

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  // Need to construct
  // 1. positions in supercell (they're given in unit cell direct coord)
  // 2. displacements in supercell (they're given in cartesian coord)
  // 3. list of which atoms are frozen
  double Rpos[Nsuper][3];
  int inv_n[Ndim];
  inverse(super_n, inv_n);
  double cartA_inv[Ndim];
  inverse(cartA, cartA_inv);
  mult(cartA_inv, 1/det(cartA), cartA_inv);

  double Rcut2 = Rcut*Rcut;
  int cutoff;
  cutoff = (Rcut2 > 0.0);
  int frozen[Nsuper];
  double cartA_sq[9];
  square(cartA, cartA_sq);

  for (np=0; np<Nsuper; ++np) {
    // position in supercell
    int rlatt[3];
    mult_vect(inv_n, Runit[np], rlatt);
    for (k=0; k<3; ++k) Rpos[np][k] = rlatt[k] * (1./Nsuper);
    // frozen?
    if (magnsq(cartA_sq, Rpos[np]) < Rcut2)  frozen[np] = 0;
    else                                     frozen[np] = -1;
    // Now, fold back into supercell:
    if (FOLDDOWN) {
      for (k=0; k<3; ++k) {
	if (Rpos[np][k] < 0.)  ++(Rpos[np][k]);
	if (Rpos[np][k] >= 1.) --(Rpos[np][k]);
      }
    }
    // displacement
    double disp_cart[3], disp_unit[3];
    mult_vect(Gperiod[np], testf, disp_cart);
    mult_vect(cartA_inv, disp_cart, disp_unit);

    // Put together:
    for (k=0; k<3; ++k) Rpos[np][k] += disp_unit[k];
  }

  // ****************************** OUTPUT ***************************

  if (VERBOSE) {
    // create an xyz version
    printf("%d\n", Nsuper);
    printf("created with %s %s  f0 = %.5le %.5le %.5le\n", 
	 cell_name, period_name, testf[0], testf[1], testf[2]);
    char* atomname = "Mo"; // may want to extract this instead...
    for (np=0; np<Nsuper; ++np) {
      double Rcart[3];
      mult_vect(cartA, Rpos[np], Rcart);
      printf("%s %18.12lf %18.12lf %18.12lf\n", atomname,
	     Rcart[0], Rcart[1], Rcart[2]);
    }
  }

  // Roll your own POSCAR file:
  printf("<sysname>: created with %s %s  f0 = %.5le %.5le %.5le", 
	 cell_name, period_name, testf[0], testf[1], testf[2]);
  if (cutoff)
    printf(" Rcut = %.3lf\n", Rcut);
  else
    printf("\n");
  printf("    1.000\n");
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j)
      printf(" %16.12lf", cartA[i+3*j]);
    printf("\n");
  }
  printf("    %3d\n", Nsuper);
  printf("external force\n");
  if (cutoff) printf("selective dynamics\n");
  printf("direct\n");
  for (np=0; np<Nsuper; ++np) {
    printf("%15.12lf %15.12lf %15.12lf", 
	   Rpos[np][0], Rpos[np][1], Rpos[np][2]);
    if (! cutoff) printf("\n");
    else {
      if (frozen[np]) printf(" .f. .f. .f.\n");
      else            printf(" .t. .t. .t.\n");
    }
  }
  printf("cartesian\n"); // output in cartesian!!
  for (np=0; np<Nsuper; ++np) {
    double force[3];
    if ( (Runit[np][0]==0) && (Runit[np][1]==0) && (Runit[np][2]==0) )
      for (i=0; i<3; ++i)
	force[i] = testf[i] * (1. - 1./Nsuper);
    else
      for (i=0; i<3; ++i)
	force[i] = testf[i] * (   - 1./Nsuper);
    printf("%15.12lf %15.12lf %15.12lf\n", 
	   force[0], force[1], force[2]);
  }

  // ************************* GARBAGE COLLECTION ********************
  for (np=0; np<Nsuper; ++np) {
    delete[] Runit[np];
    delete[] Gperiod[np];
  }
  delete[] Runit;
  delete[] Gperiod;
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
