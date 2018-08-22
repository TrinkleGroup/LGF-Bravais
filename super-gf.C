/*
  Program: super-gf.C
  Author:  D. Trinkle
  Date:    August 24, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. the dynamical matrix defined at a set of points in a supercell
	   3. a periodic supercell

	   and we compute
	   1. periodic version of the lattice GF evaluated at points
	      inside of the supercell only

  Param.:  <cell> <Dij(R))> <supercell>
           cell:    cell file describing our lattice
	   Dij(R):  dynamical matrix function evaluated at a set of points
	   supercell: supercell geometry

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
	   N kmax                      # number of points, kmax
	   n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij
	   ...
	   nN.1 nN.2 nN.3  Lxx .. Lzz
	   ==== Dij(R) ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We generate points in our supercell, and k equivalent to Gamma.
           We FT Dij(R) at the k points, then invert the matrices to get
	   G(k), and IFT at the points in our supercell.

  Output:  Going to output a series of lattice points inside the supercell, 
           and the periodic GF evaluated at those points.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"      // cell info
#include "supercell.H" // supercell filling
#include "Dij.H"       // point_type
#include "fourier.H"   // fourier and IFT

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
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <Dij(R)> <supercell>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:      cell file describing our lattice\n\
  Dij(R):    dynamical matrix evaluated at points on supercell\n\
  supercell: supercell geometry";

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
n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Lij\n\
...\n\
nN.1 nN.2 nN.3  Dxx .. Dzz\n\
==== Dij(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, m, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *latt_name, *super_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  latt_name = args[1];
  super_name = args[2];
  
  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
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

  //++ ==== latt(R) ====
  int Np;
  point_type *latt;
  
  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, latt, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== latt(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", latt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (n=0; n<Np; ++n) {
      point_type *tp = latt + n;
      printf("## %3d %3d %3d ", 
             tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      print_mat(tp->mat);
      printf("\n");
    }
  }


  //++ ==== supercell ====
  int super_n[9];
  infile = myopenr(super_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", super_name);
    exit(ERROR_NOFILE);
  }
  // For now, just read the matrix (later on, read the cartesian coord?)
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d %d %d %d %d %d %d",
	 &(super_n[0]), &(super_n[1]), &(super_n[2]), 
	 &(super_n[3]), &(super_n[4]), &(super_n[5]), 
	 &(super_n[6]), &(super_n[7]), &(super_n[8]));
  myclose(infile);
  
  //-- ==== supercell ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }


  // ***************************** ANALYSIS **************************
  // 1. Make our supercell points
  point_type *period;
  int Nsuper;
  fill_super(cart, super_n, Nsuper, period);

  if (TESTING) {
    printf("## Nsuper= %d\n", Nsuper);
    for (n=0; n<Nsuper; ++n) {
      printf("## %4d: %3d%3d%3d  %15.8lf %15.8lf %15.8lf\n", n+1,
	     period[n].Runit[0], period[n].Runit[1], period[n].Runit[2],
	     period[n].Rcart[0], period[n].Rcart[1], period[n].Rcart[2]);
    }
  }


  // 2. Time to make k = 0 (mod[n]) and the weights
  // 2.a. First, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[9], cartB[9];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

  // 2.b. Make our set of gamma equivalent kpoints:
  int** klist;
  int Nkpt;
  int transn[9]; // transpose of super_n
  transpose(super_n, transn);
  int inv_transn[9];
  Nkpt = inverse(transn, inv_transn);
  mult(cartb, inv_transn, cartB);
  mult(cartB, 1./Nkpt, cartB);

  // Now use fill_super to make the points
  fill_super(cartB, transn, Nkpt, klist);
  // Finally, write in terms of rlv and cartesian coordinates:
  // (throw away gamma itself):
  //  double klist_unit[Nkpt-1][3], klist_cart[Nkpt-1][4]; // [3] = magn
  double** klist_unit = new double*[Nkpt-1];
  double** klist_cart = new double*[Nkpt-1];
  for (k=0; k<(Nkpt-1); ++k) {
    klist_unit[k] = new double[3];
    klist_cart[k] = new double[4];
    int temp[3];
    mult_vect(inv_transn, klist[k+1], temp);
    for (i=0; i<3; ++i) klist_unit[k][i] = temp[i]*(1./Nkpt);
    // Now, make cartesian:
    mult_vect(cartb, klist_unit[k], klist_cart[k]);
    // magnitude:
    klist_cart[k][3] = sqrt(dot(klist_cart[k], klist_cart[k]));
  }
  // Garbage collection...
  free_super(Nkpt, klist);
  --Nkpt;

  if (TESTING) {
    printf("## %d kpts: k = 0(mod [n])\n", Nkpt);
    for (k=0; k<Nkpt; ++k) {
      printf("## %4d: %8.5lf %8.5lf %8.5lf | %8.5lf %8.5lf %8.5lf  %8.5lf\n",
	     k+1,
	     klist_unit[k][0], klist_unit[k][1], klist_unit[k][2],
	     klist_cart[k][0], klist_cart[k][1], klist_cart[k][2],
	     klist_cart[k][3]);
    }
  }

  // 3. Fourier transform Dij, and invert to get Gk
  double** Gk = new double *[Nkpt];
  coskR_table_type coskR; // FT and IFT table lookup
  init_coskR_table(Nkpt, klist_unit, KPT_LATT, cart, coskR);

  if (TESTING) printf("## Fourier transform:\n");
  for (k=0; k<Nkpt; ++k) {
    double Dk[9];
    fourier_table(Np, latt, klist_unit[k], Dk, KPT_LATT, coskR);
    
    Gk[k] = new double[9];
    careful_inverse(Dk, Gk[k]);
    if (TESTING) {
      printf("## kpt: %8.5lf %8.5lf %8.5lf\n");
      printf("##   Dk:"); print_mat(Dk); printf("\n");
      printf("##   Gk:"); print_mat(Gk[k]); printf("\n");
    }
  }

  // 4. Inverse FT Gk to fill in period
  double scale = (double)Nkpt/(double)Nsuper;
  for (n=0; n<Nsuper; ++n) {
    // use the *simpler* version of inv_fourier_table:
    inv_fourier_table(period[n], Nkpt, Gk, klist_unit, KPT_LATT, coskR);
    // scale by Nkpt/Nsuper (because inv_fourier divides by number of points)
    mult(period[n].mat, scale, period[n].mat);
  }

  // ****************************** OUTPUT ***************************
  // roll our own first line:
  sprintf(dump, "%d 9 # super= %d %d %d  %d %d %d  %d %d %d",
	  Nsuper, super_n[0], super_n[1], super_n[2], 
	  super_n[3], super_n[4], super_n[5], 
	  super_n[6], super_n[7], super_n[8]);
  write_Dij(stdout, Nsuper, period, dump);

  // ************************* GARBAGE COLLECTION ********************
  free_coskR_table(coskR);
  for (k=0; k<Nkpt; ++k) {
    delete[] klist_unit[k];
    delete[] klist_cart[k];
    delete[] Gk[k];
  }
  delete[] klist_unit;
  delete[] klist_cart;
  delete[] Gk;
  delete[] period;
  delete[] latt;
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
