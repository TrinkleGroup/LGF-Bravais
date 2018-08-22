/*
  Program: errorDij.C
  Author:  D. Trinkle
  Date:    April 14, 2004
  Purpose: This takes in:
           1. unit cell definition
           2. dynamical matrix evaluated at a set of points

           and we compute
	   1. the approximate relative error and
	   2. the quadratic error estimate

  Param.:  <cell> <Dij(R)> <Rmax>
           cell:    cell file describing our lattice
           Dij(R):  dynamical matrix evaluated at a series of points
	   Rmax:    maximum value to eval. the errors

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
           VERBOSE: not used
           TESTING: usual screen diahrea

  Algo.:   just add it up.

  Output:
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "pointgroup.H"
#include "Dij.H"

// ***************************** SUBROUTINES *************************

inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <Dij(R)> <Rmax>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {};
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Dij(R):  dynamical matrix evaluated at a series of points\n\
  Rmax:    maximum value of R to evaluate (in lattice coord)";

int main ( int argc, char **argv ) 
{
  int i, j, k, n, d; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter
  int SIN = 0;

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
  char *cell_name, *Dij_name;
  double Rmax;
  
  // Let's pull off the args:
  cell_name = args[0];
  Dij_name = args[1];
  sscanf(args[2], "%lf", &Rmax);

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
  }
  square(cart, cart2);
  //-- ==== cell ====

  //++ ==== Dij(R) ====
  int Np;
  point_type* Dij;
  
  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Dij(infile, cart, Np, Dij, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== Dij(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", Dij_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (n=0; n<Np; ++n) {
      printf("## %3d %3d %3d ", 
             Dij[n].Runit[0], Dij[n].Runit[1], Dij[n].Runit[2]);
      print_mat(Dij[n].mat);
      printf("\n");
    }
  }

  // We make our R's go along the a1 axis, upto Rmax.
  int Nmax;
  Nmax = (int)(Rmax/sqrt(cart[0]*cart[0]+cart[3]*cart[3]+cart[6]*cart[6]) + 1);
  int R[3] = {0,0,0};

  double D0[9], relD[9], Rmagn_2, Rmagn_1;
  // find the zero point...
  for (n=0; (n<Np) && (! equal_vect(Dij[n].Runit, R)); ++n) ;
  for (d=0; d<9; ++d) D0[d] = Dij[n].mat[d];
  double invD0[9], scale, temp[9];
  scale = inverse(D0, invD0);
  mult(scale, invD0, invD0);  
  
  for (R[0]=1; R[0]<=Nmax; ++(R[0])) {
    Rmagn_2 = 1./magnsq(cart2, R);
    Rmagn_1 = sqrt(Rmagn_2);
    // Let's sum it up.
    for (d=0; d<9; ++d) relD[d] = 0;

    for (n=0; n<Np; ++n) {
      point_type *tp = Dij+n;
      double delG = -1; // the one constant piece involved...
      double xdotR_R2, xdivR2; // (x.R)/ R^2 and (x/R)^2
      xdotR_R2 = innerprod(tp->Runit, cart2, R) * Rmagn_2;
      xdivR2 = tp->Rmagn * tp->Rmagn * Rmagn_2;
      if (tp->Runit[0] != R[0]) // note: need to change if we do different R's
	delG += 0.5/sqrt(1 - 2.*xdotR_R2 + xdivR2);
      if (tp->Runit[0] != -R[0]) // note: need to change if we do different R's
	delG += 0.5/sqrt(1 + 2.*xdotR_R2 + xdivR2);
      delG -= 1.5*xdotR_R2*xdotR_R2 - 0.5*xdivR2;
      
      // Just do the scaling to put it in...
      mult(delG, tp->mat, temp); // now, we've multiplied delG * Dij(x)
      for (d=0; d<9; ++d) relD[d] += temp[d];
    }
    // scale by invD0:
    mult(relD, invD0, temp);

    // get the RMS for the error estimate:
    for (d=0, scale=0; d<9; ++d) scale += relD[d]*relD[d];
    scale = sqrt(scale);
    
    printf("%.5lf %.8le\n", 1./Rmagn_1, scale);
  }
  
  printf("&\n");

  int NCij = (int)(Nmax + 1);
  double Cij[NCij+1]; // contributions to each shell...
  double dC;
  double Cij0 = 0;
  for (n=0; n<=NCij; ++n) Cij[n] = 0;
  
  for (n=0; n<Np; ++n) {
    point_type *tp = Dij + n;
    for (dC=0, d=0; d<9; ++d) dC += tp->mat[d] * tp->mat[d];
    dC = sqrt(dC) * tp->Rmagn * tp->Rmagn;
    for (i=NCij; i>(tp->Rmagn); --i) Cij[i] += dC;
    Cij0 += dC; // always add to this one...
  }

  // print error estimate from quadratic piece:
  double quaderr[Nmax+1];
  for (scale=0, d=0; d<9; ++d) scale += D0[d]*D0[d];
  scale = 1./sqrt(scale);
  R[1]=0;  R[2]=0;
  for (R[0]=1; R[0]<=Nmax; ++(R[0])) {
    double Rmagn = sqrt(magnsq(cart2, R));
    quaderr[R[0]] = fabs(Cij0-Cij[R[0]])*scale/(Rmagn*Rmagn);
    printf("%.5lf %.8le\n", Rmagn, quaderr[R[0]]);
  }
  
  return 0;
}
