/*
  Program: symmetrize.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a lattice GF and symmetrize it.  For now, we generate
           all of the point group operations on the fly.
           

  Param.:  <cell> <GL(R)>
           cell:   cell file describing our lattice
	   GL(R):  lattice green function

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

           ==== GL(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG
           ...
           nN.1 nN.2 nN.3  Gxx .. Gzz
           ==== GL(R) ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and pass off to ball.H routines.

  Output:  We output the points in unit coord.  (Switch for lattice coord?)
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "ball.H" 
#include "pointgroup.H"
#include "Dij.H"
#include "shell.H"


//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <GL(R)>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'s','z'}; // flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  GL(R):  lattice greens function\n\
  -s  list of points is already sorted by magnitude\n\
  -z  if (0,0,0) is in file, set it to -sum_(r>0) GL(R); for Dij";

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
==== GL(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GF\n\
...\n\
nN.1 nN.2 nN.3  Gxx .. Gzz\n\
==== GL(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used
  int SORTED = 0;   // assume input is already sorted;
  int FIXZERO = 0;  // fix the origin -- for dynamical matrices only  

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

  // flags:
  SORTED = flagon[0];
  FIXZERO = flagon[1];
  
  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *GL_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  GL_name = args[1];

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
  if (TESTING) {
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    printf("# atomic mass = %.3lf\n", atomic_mass);
    printf("# Flag settings:\n");
    if (SORTED)
      printf("#  -s: assuming sorted.\n");
    if (FIXZERO)
      printf("#  -z: fixing the value at the origin (if present).\n");
  }
  square(cart, cart2);
  //-- ==== cell ====

  //++ ==== GL(R) ====
  int Np, Ndim = 9;
  point_type *p; // set of points

  infile = myopenr(GL_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", GL_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, p, dump); // get the first line
  if (ERROR) {delete p; exit(ERROR);}
  //-- ==== GL(R) ====
  // peel off a newline on dump, if it's there
  i = strlen(dump);
  if (dump[i-1] == '\n') dump[i-1] = '\0';
  

  // ***************************** ANALYSIS **************************
  // First, generate our point group operations:
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;

  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  if (TESTING) {
    printf("## Nops = %d\n", Nop);
    for (n=0; n<Nop; ++n) {
      printf("## %2d: %2d %2d %2d  %2d %2d %2d  %2d %2d %2d  inv= %d\n", n+1, 
	     gunit[n][0], gunit[n][1], gunit[n][2], 
	     gunit[n][3], gunit[n][4], gunit[n][5], 
	     gunit[n][6], gunit[n][7], gunit[n][8], inv_index[n]+1);
      printf("##     %6.3lf %6.3lf %6.3lf  %6.3lf %6.3lf %6.3lf  %6.3lf %6.3lf %6.3lf\n",
	     gcart[n][0], gcart[n][1], gcart[n][2], 
	     gcart[n][3], gcart[n][4], gcart[n][5], 
	     gcart[n][6], gcart[n][7], gcart[n][8]);
    }
  }


  // Now, sort the list of points
  sort_points(p, Np);

  //++ ==== generate the symmetric shells ====
  shell_type* shell_list;
  int Nshell;
  
  ERROR = gen_shell_list(p, Np, gunit, Nop, shell_list, Nshell);

  if (TESTING) {
    printf("## Nshells = %d\n", Nshell);
    for (n=0; n<Nshell; ++n) {
      shell_type *ts = shell_list + n; // this shell
      printf("## shell %d -- %d entries\n", n, ts->Nelem);
      for (i=0; i<(ts->Nelem); ++i) {
        np = ts->elem[i];
        printf("##     %3d: %3d %3d %3d\n", np,
               p[np].Runit[0], p[np].Runit[1], p[np].Runit[2]);
      }
      printf("##   symmetry info:\n##  ");
      for (i=0; i<Nop; ++i) {
	printf(" g[%d]*R0=%d", i+1, ts->gR0[i]);
      }
      printf("\n");
    }
  }
  
  if (ERROR) {exit(ERROR);}

  //-- ==== generate the symmetric shells ====

  //++ ==== symmetrize all entries ====
  // FIRST step: make sure every element is it's own transpose:
  double gtrans[Ndim];
  for (np=0; np<Np; ++np) {
    point_type *tp = p + np;
    transpose(tp->mat, gtrans);
    for (d=0; d<Ndim; ++d) tp->mat[d] = 0.5*(tp->mat[d] + gtrans[d]);
  }

  // Note: we now do all of this in place, by working one shell at a time
  double Gsymm[Ndim];     // temporary holding place for our symmetrized GL
  double gscale = 1./Nop;
  
  if (TESTING) {
    printf("## Npoints = %d\n", Np);
  }

  for (i=0; i<Nshell; ++i) {
    shell_type *ts = shell_list + i; // this shell...
    
    double grot[Ndim], temp[Ndim];
    // Start by calculating the symmetrized version for Runit[0]
    for (d=0; d<Ndim; ++d) Gsymm[d] = 0;
    for (n=0; n<Nop; ++n) {
      // calculate the rotated GF value: g^-1 G(gR) g
      // gR = indexed by ts->gR0[n]
      mult(gcart[inv_index[n]], p[ts->gR0[n]].mat, temp);
      mult(temp, gcart[n], grot);
      // add contribution:
      for (d=0; d<Ndim; ++d) Gsymm[d] += grot[d];
    }
    // Scale by 1/Nop:
    for (d=0; d<Ndim; ++d) Gsymm[d] *= gscale;

    // Now, put each one in it's place:
    for (j=0; j<ts->Nelem; ++j) {
      // we need to find which group element transforms j to Runit[0]
      np = ts->elem[j];
      for (n=0; np != ts->gR0[n]; ++n) ;
      // Now, g[n].R0 = R[np], so GL[np] = g[n] GL[R0] g[n]^-1...
      mult(gcart[n], Gsymm, temp);
      mult(temp, gcart[inv_index[n]], p[np].mat);
    }
  }

  // Fix zero?
  if (FIXZERO) {
    point_type *tp = p;
    if ( (tp->Runit[0]==0)&&(tp->Runit[1]==0)&&(tp->Runit[2]==0) ) {
      // We've got zero, so let's sum up everybody else:
      double gsum[Ndim];
      for (d=0; d<Ndim; ++d) gsum[d] = 0.;
      for (np=1; np<Np; ++np)
	for (d=0; d<Ndim; ++d) gsum[d] += p[np].mat[d];
      for (d=0; d<Ndim; ++d) tp->mat[d] = -gsum[d];
    }
    else {
    fprintf(stderr, 
	    "You selected -0 (fix zero) but zero isn't in the list...?\n");
    fprintf(stderr, "the smallest value was %d %d %d\n", 
	    tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    }
  }

  // ****************************** OUTPUT ***************************

  write_Dij(stdout, Np, p, dump);

  // ************************* GARBAGE COLLECTION ********************
  free_shell(shell_list, Nshell);
  delete[] p;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
