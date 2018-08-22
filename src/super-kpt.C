/*
  Program: super-kpt.C
  Author:  D. Trinkle
  Date:    October 15, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. a periodic supercell

	   and we compute
	   1. the folding matrices for kpoints in the supercell BZ

  Param.:  <cell> <supercell>
           cell:    cell file describing our lattice
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

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   

  Output:  
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "supercell.H" // supercell 
#include "Dij.H"
#include "voigt.H"
#include "pointgroup.H"
#include "shell.H"
#include "fourier.H"
#include "kpts.H" // fold_down

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

// "positive" mod... a (mod m) returns >=0 regardless of sign of a or m.
inline int pmod (const int &a, const int& m) 
{
  int y = (a)%(abs(m));
  if (y>=0) return y;
  else      return y+abs(m);
}


// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}



void insert_kpt(int kfold[][3], int* weight, int &Nfold, int kpt[3]) 
{
  int n;
  for (n=0; (n<Nfold) && (! equal_vect(kfold[n], kpt)); ++n) ;
  if (n==Nfold) {
    for (int d=0; d<3; ++d) kfold[n][d] = kpt[d];
    weight[n] = 1;
    ++Nfold;
  }
  else
    ++(weight[n]);
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <supercell>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:      cell file describing our lattice\n\
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
==== cell ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, k, n; // General counting variables.

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
  // Let's pull off the args:
  char* cell_name = args[0];
  char* super_name = args[1];
  
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
  // 0. Time to make k = 0 (mod[n]) and the weights
  // 0.a. First, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[9], cartB[9];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

  // 0.b. Make our set of gamma equivalent kpoints:
  int** klist;
  int Nkpt;
  int transn[9]; // transpose of super_n
  transpose(super_n, transn);
  int inv_transn[9];
  Nkpt = inverse(transn, inv_transn);
  if (Nkpt < 0) {
    mult(inv_transn, -1, inv_transn);
    Nkpt = -Nkpt;
  }
  mult(cartb, inv_transn, cartB);
  mult(cartB, 1./Nkpt, cartB);

  // Now use fill_super to make the points
  fill_super(cartB, transn, Nkpt, klist);


  // 1. Pointgroup info
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  // 1.a Point group operations:
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  // 1.b Point group operations on our supercell:
  int g_eta[MAXop][9];
  for (n=0; n<Nop; ++n) {
    int temp_mat[9];
    i = inverse(gunit[n], temp_mat);
    mult(temp_mat, i, temp_mat);
    self_transpose(temp_mat);
    // temp_mat = gunit^T^-1
    mult(transn, temp_mat, inv_transn, g_eta[n]);
    for (d=0; d<9; ++d) g_eta[n][d] = pmod(g_eta[n][d], Nkpt);
  }
  
  if (VERBOSE) {
    for (n=0; n<Nop; ++n) {
      printf("# g%2d:", n+1);
      for (i=0; i<3; ++i) {
	if (i!=0) printf("#     ");
	for (int j=0; j<3; ++j)
	  printf(" %3d", gunit[n][i+3*j]);
	printf("\n");
      }
      printf("# gn~:");
      for (i=0; i<3; ++i) {
	if (i!=0) printf("#     ");
	for (int j=0; j<3; ++j)
	  printf(" %3d", g_eta[n][i+3*j]);
	printf("\n");
      }
      printf("#\n");
    }
  }

  // 2. Scale out our kpoints
  int kfold[Nkpt*Nop][3], Nfold=0;
  int weight[Nkpt*Nop];
  int zero_kpt[3] = {0,0,0};
  int krot[3];
  
  insert_kpt(kfold, weight, Nfold, zero_kpt);
  weight[0] = 0;

  for (n=0; n<Nop; ++n)
    if (zero(g_eta[n])) 
      weight[0] += Nkpt;
    else {
      for (k=0; k<Nkpt; ++k) {
	mult_vect(g_eta[n], klist[k], krot);
	for (d=0; d<3; ++d) krot[d] = pmod(krot[d], Nkpt);
	insert_kpt(kfold, weight, Nfold, krot);
      }
    }

  // ****************************** OUTPUT ***************************
  double Nweight = 1./(double)(Nkpt*Nop);
  double latscale = 1./(double)Nkpt;
  printf("%d 1.0 L # Np, scale, lattice coord; supercell=", Nfold);
  for (d=0; d<9; ++d) printf(" %d", super_n[d]);
  printf("\n");
  for (k=0; k<Nfold; ++k)
    printf("%18.15lf %18.15lf %18.15lf %.15le\n", 
	   kfold[k][0]*latscale,
	   kfold[k][1]*latscale,
	   kfold[k][2]*latscale,
	   Nweight * weight[k]);

  // ************************* GARBAGE COLLECTION ********************
  free_super(Nkpt, klist);
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
