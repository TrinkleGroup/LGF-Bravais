/*
  Program: egf-FT.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: Calculate the FT of the anisotropic elastic greens function
           for a given set of k directions given the elastic constants
	   Cmn, crystal class c

	   Note: Cij are read in in GPa, and converted to eV/A^3.

  Param.:  <cell> <grid>
           cell:  cell file (see below for format)
           grid:  "grid" of direction cosines for which to compute Gij(k)

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
	   
	   ==== grid ====
	   N         # Number of points
	   l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)
	   l2 m2 n2
	   ...
	   lN mN nN
	   ==== grid ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   Read in everything, and calculate.  All we're gonna output
           is g~(k)_ab = (k_c C_acdb k_d)^-1, where k is normalized

  Output:  Pretty self-explanatory
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

// Calculate (kk)
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


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <grid>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'g','c'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:  cell file (see below for format)\n\
  grid:  \"grid\" of direction cosines for which to compute Gij(k)\n\
\n\
  -g               don't convert from GPa to eV/A^3\n\
  -c               \"careful\" inverse (more accurate)";

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
==== grid ====\n\
N         # Number of points\n\
l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)\n\
...\n\
lN mN nN\n\
==== infile ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, n; // General counting variables.

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

  int EV_CONV = !(flagon[0]);  // Convert from GPa to eV/A^3?
  int CAREFUL = (flagon[1]);   // use careful inversion
  
  // ****************************** INPUT ****************************
  // Command line parameters:
  char *cell_name, *infile_name;
  // Let's pull off the args:
  cell_name = args[0];
  infile_name = args[1];

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
    for (n=0; n<class_len[crystal]; ++n)
      Cmn_list[n] *= eV_GPa;
    if (TESTING) {
      printf("## Converted elastic constants:\n");
      verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    }
  }
  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);


  //++ ==== grid ====
  infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", infile_name);
    exit(ERROR_NOFILE);
  }

  int Npoints;
  double** lmn_list;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Npoints);
  if (Npoints <= 0) {
    fprintf(stderr, "No points listed in %s\n", infile_name);
    ERROR = ERROR_BADFILE;
  }
  else lmn_list = new double* [Npoints];
  for (n=0; (!ERROR) && (n<Npoints) && (!feof(infile)); ++n) {
    double magn;
    lmn_list[n] = new double[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", 
	   &(lmn_list[n][0]), &(lmn_list[n][1]), &(lmn_list[n][2]) );
    // Normalize.
    magn = sqrt(dot(lmn_list[n],lmn_list[n]));
    if (dcomp(magn, 0.)) {
      fprintf(stderr, 
	      "%s at line %d: %.5le %.5le %.5le with magn = %.5le\n",
	      infile_name, i+1,
	      lmn_list[0][0], lmn_list[0][1], lmn_list[0][2], magn);
      ERROR = ERROR_BADFILE;
    }
    else {
      magn = 1./magn;
      for (i=0; i<3; ++i) lmn_list[n][i] *= magn;
    }
  }
  // Make sure we read enough points...
  if (n != Npoints) ERROR = ERROR_BADFILE;
  myclose(infile);
  //-- ==== grid ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }
  

  // ***************************** ANALYSIS **************************

  // Spit out our first bit of info:
  printf("%d 9 # Number of points: l m n  g[xx xy xz  yx yy yz  zx zy zz]",
	 Npoints);
  if (EV_CONV)
    printf(" (A^3/eV)\n");
  else
    printf(" (--)\n");
  // For each point...
  for (n=0; n<Npoints; ++n) {
    double *lmn;
    double g[9], ginv[9], gdet;
    lmn = lmn_list[n]; // Our specific direction.
    
    a_mult_a(lmn, Cijkl, g);
    if (CAREFUL) {
      double err;
      careful_inverse(g, ginv, err);
      if (TESTING) {
	printf("## estimated error in inverse: %.5le\n", err);
      }
    }
    else {
      gdet = symm_inverse(g, ginv);
      mult(ginv, (1./gdet), ginv); // scale properly
    }
    
    // ****************************** OUTPUT ***************************
    // Human readable (sorta) first:
    
    if (VERBOSE) {
      // What, exactly, is "verbose" output...?
    }
    // Spit back out the lmn point, and include the gij matrix values:
    printf("%15.12lf %15.12lf %15.12lf", lmn[0], lmn[1], lmn[2]);
    for (i=0; i<9; ++i) printf(" %.15le", ginv[i]);
    printf("\n");
  }
  
  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Npoints; ++n)
    delete lmn_list[n];
  delete[] lmn_list;
  
  free_cell(Cmn_list, u_atoms, Natoms);
  delete[] Cmn_list;

  return 0;
}
