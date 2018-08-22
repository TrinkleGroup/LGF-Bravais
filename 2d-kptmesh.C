/*
  Program: 2d-kptmesh.C
  Author:  D. Trinkle
  Date:    September 9, 2004
  Purpose: Generate a kptmesh in a 2d plane of the BZ

  Param.:  <cell> <infile> <M> <N>
           cell:     cell file (see below for format)
           infile:   input file (see below for format)
	   M,N:      divisions in the BZ directions

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
	   m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)
	   n1 n2 n3
	   ==== infile ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:  
  Output:  
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
#include "dislocation.H"
#include "kpts.H"
#include "fourier.H"  // included just to get our hands on the gcd algorithm

//****************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <infile> <M> <N>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'s'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:     cell file (-h for format)\n\
  infile:   input file (-h for format)\n\
  M,N:      divisions in the BZ directions\n\
  -s        apply a shift off of gamma in the plane";

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
  int i, k, n, d; // General counting variables.

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
  int SHIFT = flagon[0];
  
  
  // ****************************** INPUT ****************************
  //  char dump[512];
  FILE* infile;

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  int Natoms;
  double** u_atoms;

  // Let's pull off the args:
  char* cell_name = args[0];
  char* disl_name = args[1];
  int M; sscanf(args[2], "%d", &M);
  int N; sscanf(args[3], "%d", &N);
  int T; // number of divisions for t direction; calculated later.
  
  if ( (M<1) || (N<1) ) {
    fprintf(stderr, "Error: both M and N must be >0\n");
    exit(ERROR_BADFILE);
  }
  
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
    fprintf(stderr, "Problem with cell: %s\n", cell_name);
    exit(ERROR);
  }

  if (TESTING)
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
  //-- ==== cell ====

  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cart_rlv[9];
  double cell_vol;
  cell_vol = inverse(cart, cart_rlv);        // invert
  self_transpose(cart_rlv);                  // transpose in place
  mult(cart_rlv, 2*M_PI/cell_vol, cart_rlv); // scale

  if (TESTING) {
    printf("## Reciprocal lattice:\n");
    printf("## b0 = %16.12lf %16.12lf %16.12lf\n", 
           cart_rlv[0], cart_rlv[3], cart_rlv[6]);
    printf("## b1 = %16.12lf %16.12lf %16.12lf\n", 
           cart_rlv[1], cart_rlv[4], cart_rlv[7]);
    printf("## b2 = %16.12lf %16.12lf %16.12lf\n", 
           cart_rlv[2], cart_rlv[5], cart_rlv[8]);
  }

  //++ ==== disl ====
  infile = myopenr(disl_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", disl_name);
    exit(ERROR_NOFILE);
  }

  int t_unit[3], m_unit[3], n_unit[3];
  ERROR = read_disl(infile, t_unit, m_unit, n_unit);
  myclose(infile);
  //-- ==== disl ====

  if (ERROR) {
    fprintf(stderr, "Problem with dislocation: %s\n", disl_name);
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  // 0. Construct our coordinate system
  double t0[3], m0[3], n0[3];
  int disl[9], disl_inv[9], disl_det;

  // 0.a. go ahead and overwrite our unit coord:
  make_dislcoord(cart, t0, m0, n0, t_unit, m_unit, n_unit, 
		 DISLCOORD_UNIT | DISLCOORD_CHANGE);

  // 0.b. make our supercell
  for (d=0; d<3; ++d) {
    disl[3*d+0] = m_unit[d];
    disl[3*d+1] = n_unit[d];
    disl[3*d+2] = t_unit[d];
  }
  disl_det = inverse(disl, disl_inv);
  self_transpose(disl_inv);
  // Now, we need to determine our "multipliers":
  int scale; // the last one will be the gcd for the t direction
  for (i=0; i<3; ++i) {
    scale = gcd (disl_inv[i], disl_inv[i+3], disl_inv[i+6]);
    for (d=0; d<3; ++d) disl_inv[i+3*d] /= scale;
  }
  T = disl_det/scale; // I believe this is correct.

  if (TESTING) {
    printf("## vect    input                 cartesian\n");
    printf("## t = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           t_unit[0],t_unit[1],t_unit[2], t0[0],t0[1],t0[2]);
    printf("## m = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           m_unit[0],m_unit[1],m_unit[2], m0[0],m0[1],m0[2]);
    printf("## n = (%3d %3d %3d)-> %12.8lf %12.8lf %12.8lf\n", 
           n_unit[0],n_unit[1],n_unit[2], n0[0],n0[1],n0[2]);
    printf("## kpt directions:\n");
    double kptcell[9];
    mult(cart_rlv, disl_inv, kptcell);
    for (i=0; i<3; ++i) {
      printf("## k%d:", i+1);
      for (d=0; d<3; ++d) printf(" %3d", disl_inv[i+3*d]);
      printf(" | ");
      for (d=0; d<3; ++d) printf(" %8.5lf", kptcell[i+3*d]);
      printf("  %8.5lf\n", 
	     sqrt(kptcell[i]*kptcell[i]+kptcell[i+3]*kptcell[i+3]+
		  kptcell[i+6]*kptcell[i+6]));
    }
    printf("## Mesh: %d x %d x %d\n", M, N, T);
  }


  // 1. Let's make us some kpoints.
  // This code is swiped from kptmesh, with some minor modifications
  // For one, we no longer bother doing the symmetry stuff
  int ni[3], Ni[3] = {M, N, T};
  int Nkpt;
  double** kpt; // Our points--unit cell coord, but LAST index is magn.
  Nkpt = 8*(Ni[0]*Ni[1]*Ni[2]); // Maximum number of points...
  // For now, we only allocation the list of *pointers*...
  kpt = new double *[Nkpt];

  double q0[3], q[3], qfold[3], dN[3];
  double qmagn;
  double dq;
  //  for (d=0; d<3; ++d) dN[d] = 0.5/Ni[d];
  dN[0] = 0.5/Ni[0];
  dN[1] = 0.5/Ni[1];
  dN[2] = 1.0/Ni[2];
  if (SHIFT) dq=0.5; else dq=0.;
  
    if (TESTING) {
    printf("## kpoint search...\n");
  }
  Nkpt = 0;
  for (ni[0]=-Ni[0]+1; ni[0] <= Ni[0]; ++(ni[0])) {
    q0[0] = dN[0] * (ni[0]+dq);
    for (ni[1]=-Ni[1]+1; ni[1] <= Ni[1]; ++(ni[1])) {
      q0[1] = dN[1] * (ni[1]+dq);
      for (ni[2]=-Ni[2]+1; ni[2] <= Ni[2]; ++(ni[2])) {
        q0[2] = dN[2] * (ni[2]);
	mult_vect(disl_inv, q0, q); // added to put us on them planes
        qmagn = fold_down (q, cart_rlv, qfold);
        if (TESTING) {
          printf("## (%3d%3d%3d) q = %8.5lf %8.5lf %8.5lf -> %8.5lf %8.5lf %8.5lf  %8.5lf\n",
                 ni[0], ni[1], ni[2], 
		 q[0], q[1], q[2], qfold[0], qfold[1], qfold[2], qmagn);
        }
        // Now, we need to see if the point is in our list.
        // 1. do a binary search on kptinfo:
        int lo, hi;
        if (Nkpt > 0) {
          int probe, otherlo;
          lo = -1;
          hi = Nkpt;
          // get the floor first
          while ( (hi-lo) > 1) {
            probe = (hi+lo)/2;
            if (kpt[probe][3] < (qmagn-TOLER)) lo = probe;
            else                               hi = probe;
          }
          hi = Nkpt;
          otherlo = lo;
          while ( (hi - otherlo) > 1) {
            probe = (hi+otherlo)/2;
            if (kpt[probe][3] > (qmagn+TOLER)) hi = probe;
            else                               otherlo = probe;
          }
          if (lo < 0) lo = 0;
          if (kpt[lo][3]<(qmagn-TOLER)) ++lo;
          if (hi >= Nkpt) hi = Nkpt-1;
          if (kpt[hi][3]>(qmagn+TOLER)) --hi;
        } else {
          lo = 0;
          hi = -1; // We didn't find it; we also insert it at hi+1
        }
        // So now, we should have all of the kpoints in the range from
        // lo..hi.  So now:
        // kpt[0]..kpt[lo-1] < qmagn
        // kpt[lo]..kpt[hi] == qmagn (within TOLER)
        // kpt[hi+1]..kpt[Nkpt-1] > qmagn
        int found;
        found = ( (hi-lo) >= 0 ); // need at least one point to bother...
        if (found) {
          found = 0;
          // See if it's really there...
          for (n=lo; (n<=hi) && (!found); ++n)
	    found = equal_vect(qfold, kpt[n]);
          --n; // now, n indexes the one we matched
        }
        // Now, if we HAVEN'T found it, we need to add it:
        if (!found) {
          if (TESTING) {
            printf("##  -- New point!\n");
          }
          // first, move everything up one:
          for (n=(Nkpt-1); n>hi; --n)
            kpt[n+1] = kpt[n];
          // Now, insert our point:
          n = hi+1;
          kpt[n] = new double[4];
          for (d=0; d<3; ++d) kpt[n][d] = qfold[d];
          kpt[n][3] = qmagn;
          ++Nkpt;
        } else {
          if (TESTING) {
            printf("##  -- matches: %8.5lf %8.5lf %8.5lf  %8.5lf\n", 
                   kpt[n][0], kpt[n][1], kpt[n][2], kpt[n][3]);
          }
        }
        /*
        // EXTREME VERBOSITY WARNING
        if (TESTING) {
          printf("##   Nkpt = %d\n", Nkpt);
          for (n=0; n<Nkpt; ++n)
            printf("##   %3d %8.5lf %8.5lf %8.5lf  %8.5lf\n",
                   n, kpt[n][0], kpt[n][1], kpt[n][2], kpt[n][3]);
        }
        */
      }
    }
  }
  
  if (VERBOSE) {
    printf("# %d final kpoints (no weights yet)\n", Nkpt);
    for (i=0; i<Nkpt; ++i) {
      printf("# %15.12lf %15.12lf %15.12lf  %.5lf\n",
             kpt[i][0], kpt[i][1], kpt[i][2], kpt[i][3]);
    }
  }

  // ****************************** OUTPUT ***************************
  double w = 1./(double)Nkpt; // equal weights...
  
  // Now, we output EVERYTHING.
  printf("%d 1.0 L # Nkpt, scale, lattice coord: %dx%d mesh along [%d %d %d]\n",
         Nkpt, M, N, t_unit[0], t_unit[1], t_unit[2]);
  for (n=0; n<Nkpt; ++n)
    printf("%18.15lf %18.15lf %18.15lf  %.14le\n", 
           kpt[n][0], kpt[n][1], kpt[n][2], w);

  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Nkpt; ++n) delete[] kpt[n];
  delete[] kpt;

  free_cell(Cmn_list, u_atoms, Natoms);

  return 0;
}
