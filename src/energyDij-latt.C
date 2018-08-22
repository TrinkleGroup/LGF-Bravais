/*
  Program: energyDij-latt.C
  Author:  D. Trinkle
  Date:    December 15, 2004
  Purpose: Read in:
	   1. Cell info
           2. Set of positions (original latt positions and cartesian coord)
	   3. Dynamical matrix
	   4. Supercell info

	   and we apply Dij to the positions (or dynamical matrix), to compute
	   the harmonic energy.  We do this in Fourier space to make it more
	   correct for "large" displacements.

	   We've recently given up on trying to calc. the original latt
	   position on the fly--it's too hard to do correctly for large
	   temperatures (which is what we care about).  Drift is partly
	   to blame, I believe.

  Param.:  <cell> <pos> <Dij> <super>
           cell:    cell file describing our lattice
	   pos:     atomic positions and forces in cartesian coord.
	   Dij:     dynamical matrix
	   super:   supercell definition

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

           ==== pos ====
           N                             # number of points
           R1.1 R1.2 R1.3 r1.x r1.y r1.z # original latt position and cart.
           ...
           RN.1 RN.2 RN.3 rN.x rN.y rN.z
           ==== pos ====

	   ==== Dij ====
	   N                          # number of points
	   u1.1 u1.2 u1.3 Dxx .. Dzz  # position (unit) and Dij
	   ...
	   uN.1 uN.2 uN.3 Dxx .. Dzz
	   ==== Dij ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   

  Output:  We output the new positions in cartesian coord. from the update.
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
#include "cell.H"
#include "Dij.H"
#include "supercell.H" // indexing routines for easy lookup
#include "fourier.H"


//***************************** SUBROUTINES ****************************

// structure to hold atomic positions, forces, and displacements:
typedef struct 
{
  int R[3];
  double r[3], u[3];
} posfor_type;

inline void setmax (int &a, const int &b) 
{if (a<b) a=b;}


inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <pos> <Dij> <super>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  pos:     atomic positions and forces in cartesian coord.\n\
  Dij:     dynamical matrix\n\
  super:   supercell definition";

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
==== pos ====\n\
N                             # number of points\n\
R1.1 R1.2 R1.3 r1.x r1.y r1.z # original latt position and cart.\n\
...\n\
RN.1 RN.2 RN.3 rN.x rN.y rN.z\n\
==== pos ====\n\
\n\
==== Dij ====\n\
N                          # number of points\n\
u1.1 u1.2 u1.3 Dxx .. Dzz  # position (unit) and Dij\n\
...\n\
uN.1 uN.2 uN.3 Dxx .. Dzz\n\
==== Dij ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  // Let's pull off the args:
  char* cell_name = args[0];
  char* pos_name = args[1];
  char* Dij_name = args[2];
  char* super_name = args[3];

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
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%lf", &atomic_mass);
  myclose(infile);
  //-- ==== cell ====
  
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

  //++ ==== pos ====
  int Np;
  posfor_type* p;

  infile = myopenr(pos_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", pos_name);
    exit(ERROR_NOFILE);
  }
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d", &Np);

  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in read_posfor\n", Np);
    return ERROR_BADFILE;
  }
  //    if (p != NULL) delete[] p;
  p = new posfor_type[Np];
  if (p == NULL) {
    fprintf(stderr, "Error allocating memory...?\n");
    exit(ERROR_MEMORY);
  }
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    posfor_type *tp = p + n; // this point
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) {ERROR = ERROR_BADFILE; break;}
    sscanf(dump, "%d %d %d %lf %lf %lf", 
	   &(tp->R[0]), &(tp->R[1]), &(tp->R[2]),
	   &(tp->u[0]), &(tp->u[1]), &(tp->u[2]));
    mult_vect(cart, tp->R, tp->r);
  }
  myclose(infile);
  //-- ==== pos ====

  if (ERROR) {
    fprintf(stderr, "Error with %s\n", pos_name);
    exit(ERROR);
  }

  //++ ==== Dij ====
  int Nlatt;
  point_type *Dij; // set of points

  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Nlatt, Dij);
  myclose(infile);
  //-- ==== Dij ====

  if (ERROR) {
    delete Dij; 
    fprintf(stderr, "Error with %s\n", Dij_name);
    exit(ERROR);
  }
  
  //++ ==== super ====
  int super[9];
  
  infile = myopenr(super_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", super_name);
    exit(ERROR_NOFILE);
  }
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d %d %d %d %d %d %d", super, super+1, super+2,
	 super+3, super+4, super+5, super+6, super+7, super+8);
  myclose(infile);
  //-- ==== super ====


  // ***************************** ANALYSIS **************************
  // 0. Make our set of gamma equivalent kpoints:
  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[9], cartB[9];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

  // Make our set of gamma equivalent kpoints:
  int** klist;
  int Nkpt;
  int transn[9]; // transpose of super
  transpose(super, transn);
  int inv_transn[9];
  Nkpt = inverse(transn, inv_transn);
  mult(cartb, inv_transn, cartB);
  mult(cartB, 1./Nkpt, cartB);

  // Now use fill_super to make the points
  fill_super(cartB, transn, Nkpt, klist);
  // Finally, write in terms of rlv and cartesian coordinates:
  // take out gamma:
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

  if (TESTING) {
    printf("## %d kpoints equivalent to gamma:\n", Nkpt);
    for (n=0; n<Nkpt; ++n) {
      printf("## %3d : %8.5lf %8.5lf %8.5lf | %8.5lf %8.5lf %8.5lf | %8.5lf\n", 
	     n+1,
	     klist_unit[n][0], klist_unit[n][1], klist_unit[n][2], 
	     klist_cart[n][0], klist_cart[n][1], klist_cart[n][2],
	     klist_cart[n][3]);
    }
  }  

  // 1. Now we iterate over all of our positions and forces...
  double du[3];
  double utotal[3] = {0,0,0};
  double cart_super[9], cart_super_inv[9];
  mult(cart, super, cart_super);
  careful_inverse(cart_super, cart_super_inv);

  // determine u (displacement)
  if (TESTING) {
    printf("## %d atoms: (NOTE: COM shift NOT removed)\n", Np);
  }
  for (np=0; np<Np; ++np) {
    posfor_type *tp = p + np;
    if (TESTING) {
      printf("## x%d =  %8.5lf %8.5lf %8.5lf\n",np,tp->u[0],tp->u[1],tp->u[2]);
    }
    for (d=0; d<3; ++d) tp->u[d] -= tp->r[d];
    mult_vect(cart_super_inv, tp->u, du);
    for (d=0; d<3; ++d) {
      if (du[d] < -0.5) du[d] += 1.;
      if (du[d] >= 0.5) du[d] -= 1.;
    }
    mult_vect(cart_super, du, tp->u);
    for (d=0; d<3; ++d) utotal[d] += tp->u[d];
    
    if (TESTING) {
      printf("##  u%d = %8.5lf %8.5lf %8.5lf",np,tp->u[0],tp->u[1],tp->u[2]);
      printf("  %8.5lf\n", sqrt(dot(tp->u, tp->u)));
      printf("##  R%d = %8.5lf %8.5lf %8.5lf\n",np,tp->r[0],tp->r[1],tp->r[2]);
    }
  }
  // COM shift removal:
  double rms = 0; // calculate rms deviation too
  for (d=0; d<3; ++d) utotal[d] *= -1./(double)Np;
  for (np=0; np<Np; ++np) {
    for (d=0; d<3; ++d) p[np].u[d] += utotal[d];
    rms += dot(p[np].u, p[np].u);
  }
  rms = sqrt(rms/(double)Np);
  
  // 1.b. Check consistency by calculating sum_R cos(k.R) for all k.
  //      should be equal to N delta(k)
  if (TESTING) printf("## sum test:\n");
  for (n=0; n<Nkpt; ++n) {
    double *tk = klist_cart[n];
    double sumk = 0;
    for (np=0; np<Np; ++np) sumk += cos(dot(p[np].r, tk));
    sumk *= 1./(double)Np;
    ERROR = ( zero(tk[3]) && (! dcomp(sumk, 1.)) ) ||
      ( (! zero(tk[3])) && (! zero(sumk)) ) ;
    if (TESTING)
      printf("## k%03d = %8.5lf %8.5lf %8.5lf, sum= %.4le\n",
	      n+1, tk[0], tk[1], tk[2], sumk);
    if (ERROR) {
      fprintf(stderr, "kptsum: %d = %8.5lf %8.5lf %8.5lf, sum= %.7lf\n",
	      n+1, tk[0], tk[1], tk[2], sumk);
      break;
    }
  }
  if (ERROR) {
    fprintf(stderr, "Supercell inconsistent with atomic positions...?\n");
    exit(1);
  }

  // 2. Fourier transform it all, and sum as we go:
  double W = 0;
  double uRk[3], uIk[3], Dk[9];
  
  for (n=0; n<Nkpt; ++n) {
    double *tk = klist_cart[n];
    // uk:
    for (d=0; d<3; ++d) {uRk[d] = 0; uIk[d] = 0;}
    for (np=0; np<Np; ++np) {
      posfor_type *tp = p+np;
      double coskR = cos(dot(tp->r, tk));
      double sinkR = sin(dot(tp->r, tk));
      for (d=0; d<3; ++d) 
	{uRk[d] += coskR*tp->u[d]; uIk[d] += sinkR*tp->u[d];}
    }
    for (d=0; d<3; ++d) {uRk[d] *= 1./(double)Np; uIk[d] *= 1./(double)Np;}
    // Dk:
    for (d=0; d<9; ++d) Dk[d] = 0;
    for (np=0; np<Nlatt; ++np) {
      point_type *tp = Dij+np;
      double coskR = cos(dot(tp->Rcart, tk));
      for (d=0; d<9; ++d) Dk[d] += coskR*tp->mat[d];
    }
    W += magnsq(Dk, uRk) + magnsq(Dk, uIk);
    if (TESTING) {
      printf("## k%d = %8.5lf %8.5lf %8.5lf\n", n+1, tk[0], tk[1], tk[2]);
      printf("##   u(k)= %8.5lf %8.5lf %8.5lf  +I %8.5lf %8.5lf %8.5lf\n",
	     uRk[0], uRk[1], uRk[2], uIk[0], uIk[1], uIk[2]);
      printf("##   D(k) ");
      print_mat(Dk);
      printf("\n");
      printf("##   uR.D.uR= %.5le\n", magnsq(Dk, uRk));
      printf("##   uI.D.uI= %.5le\n##\n", magnsq(Dk, uIk));
    }
  }
  W *= 0.5*(double)Np; // make into an energy per cell
  
  // ****************************** OUTPUT ***************************
  printf("%12le\n", W);
  if (VERBOSE) printf("# rms= %12le\n", rms);
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  delete[] Dij;
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}

