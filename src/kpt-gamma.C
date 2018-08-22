/*
  Program: kpt-gamma.C
  Author:  D. Trinkle
  Date:    Sept. 9, 2004
  Purpose: Remember repus?  Yeah--this is it, all over again.  Except
           it's not as advanced.  It generates the kpt mesh of points
	   equivalent to gamma in a given supercell.

  Param.:  <cell> <supercell>
           cell:       cell file describing our lattice
	   supercell:  supercell geometry

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
  Output:  kpt mesh, baby.
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
#include "kpts.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <supercell>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'g'}; // determine how we handle k=0
// flag characters.

const char* ARGEXPL = 
"  cell:       cell file describing our lattice\n\
  supercell:  supercell geometry\n\
  -g          exclude gamma from list (does NOT affect weights)";

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
\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, n; // General counting variables.

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
  int NOGAMMA = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char* cell_name = args[0];
  char* super_name = args[1];
  
  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; 
  double Cijkl[9][9];
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
  //-- ==== cell ====
  if (TESTING) {
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    printf("# atomic mass = %.3lf\n", atomic_mass);
  }
  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);


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
    fprintf(stderr, "Error encountered with file %s\n", super_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Supercell definition:\n");
    for (i=0; i<3; ++i) {
      printf("## A%d =", i+1);
      for (j=0; j<3; ++j) printf(" %+3d a%d", super_n[j*3+i], j+1);
      printf("\n");
    }
    printf("## %d atoms\n", Natoms*det(super_n));
  }

  // Let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[9], cartB[9];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

  // ***************************** ANALYSIS **************************
  // 1. Time to make k = 0 (mod[n]) and the weights

  // 1.a. Make our set of gamma equivalent kpoints:
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
  // but *keep* gamma:
  double klist_unit[Nkpt][3], klist_cart[Nkpt][4]; // [3] = magn
  for (n=0; n<Nkpt; ++n) {
    int temp[3];
    mult_vect(inv_transn, klist[n], temp);
    for (i=0; i<3; ++i) klist_unit[n][i] = temp[i]*(1./Nkpt);
    // Now, make cartesian:
    mult_vect(cartb, klist_unit[n], klist_cart[n]);
    // magnitude:
    klist_cart[n][3] = sqrt(dot(klist_cart[n], klist_cart[n]));
  }
  // Garbage collection...
  free_super(Nkpt, klist);

  // 2.c. Calculate the weights for each kpoint
  double w[Nkpt];
  for (k=0; k<Nkpt; ++k) w[k] = 1./Nkpt;

  if (TESTING) {
    printf("## %d kpoints equivalent to gamma:\n", Nkpt);
    printf("## equal weights = %.8le\n", w[0]);
    for (n=0; n<Nkpt; ++n) {
      printf("## %3d : %8.5lf %8.5lf %8.5lf | %8.5lf %8.5lf %8.5lf | %8.5lf\n", 
             n+1,
             klist_unit[n][0], klist_unit[n][1], klist_unit[n][2], 
             klist_cart[n][0], klist_cart[n][1], klist_cart[n][2],
             klist_cart[n][3]);
    }
  }


  // ****************************** OUTPUT ***************************
  int nstart;
  if (NOGAMMA) nstart=1;
  else         nstart=0;
  printf("%d 1.0 L   # number of kpoints, scaling factor, lattice coord\n",
         Nkpt-nstart);
  for (n=nstart; n<Nkpt; ++n)
    printf("%18.15lf %18.15lf %18.15lf  %.14le\n", 
           klist_unit[n][0], klist_unit[n][1], klist_unit[n][2], w[n]);

  // ************************* GARBAGE COLLECTION ********************

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
