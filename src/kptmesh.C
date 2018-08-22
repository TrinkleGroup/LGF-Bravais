/*
  Program: kptmesh.C
  Author:  D. Trinkle
  Date:    April 11, 2004
  Purpose: Construct a symmetrized, folded down MP kpt mesh.

  Param.:  <cell> <N1> <N2> <N3>
           cell:    cell file describing our lattice
	   Ni:      MP mesh (includes the origin)

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

  Algo.:   Just generate, folddown, and calculate weights.  No big whoop.
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
#include "kpts.H"
#include "pointgroup.H"


//****************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <N1> <N2> <N3>";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'s', 'u'};  // flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Ni:      MP mesh (includes the origin)\n\
  -s: shift off of gamma\n\
  -u: unsymmetrized (turn off symmetry)";

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
==== cell ====";

int main ( int argc, char **argv ) 
{
  int d, i, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used
  int SHIFT = 0;    // do we shift off of gamma?
  int SYMMETRY = 1; // symmetry?
  
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

  // Flags
  SHIFT = flagon[0];
  SYMMETRY = ! flagon[1];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name;
  int Ni[3];
  // Let's pull off the args:
  cell_name = args[0];
  sscanf(args[1], "%d", &(Ni[0]));
  sscanf(args[2], "%d", &(Ni[1]));
  sscanf(args[3], "%d", &(Ni[2]));

  if ( (Ni[0] < 0) && (Ni[1] < 0) && (Ni[2] < 0) ) {
    fprintf(stderr, "Bad mesh: %d x %d x %d\n", Ni[0], Ni[1], Ni[2]);
    exit(-1);
  }
  
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

  // ***************************** ANALYSIS **************************
  // First, generate our point group operations:
  // NOTE: we do this on the reciprocal lattice!!
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;

  if (SYMMETRY) {
    Nop = gen_pointgroup(cart_rlv, gcart, gunit, inv_index);
    ERROR = (Nop < 1);
    if (ERROR) {
      fprintf(stderr, "Somehow... gen_pointgroup failed.  This should never happen...\n");
      exit(ERROR);
    }
  } else {
    Nop = 1;
    for (d=0; d<9; ++d) gcart[0][d] = ident[d];
    for (d=0; d<9; ++d) gunit[0][d] = ident_i[d];
    inv_index[0] = 0;
  }
  
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


  // Now, let's go make them kpoints!
  int Nkpt;   // Total number of grid points
  
  // We store the UNIT coordinates of our kpoints, but the last index
  // is the magnitude, because we keep them sorted by magnitude to make
  // checking a lot easier.  kpt[3] = kmagn
  double** kpt;          // Our points
  double** gk_list;      // g.k for each group operation

  // Allocation
  if (SHIFT)
    MEMORY = 64*(Ni[0]*Ni[1]*Ni[2]); // Maximum number of points...
  else
    MEMORY = 8*(Ni[0]*Ni[1]*Ni[2]); // Maximum number of points...
  // max out at 8.8e07 (2GB memory...)
  if (MEMORY > 88000000) MEMORY = 88000000;
  Nkpt = MEMORY; // Maximum number of points...
  // For now, we only allocation the list of *pointers*...
  kpt = new double *[Nkpt];
  gk_list = new double *[Nkpt];
  
  int ni[3];
  double q[3], qfold[3], dN[3];
  double qmagn;
  double dq;
  for (d=0; d<3; ++d) dN[d] = 0.5/Ni[d];
  if (SHIFT) dq=0.5; else dq=0.;
  
  if (TESTING) {
    printf("## kpoint search...\n");
  }
  Nkpt = 0;
  for (ni[0]=-Ni[0]+1; (ni[0] <= Ni[0]) && !ERROR; ++(ni[0])) {
    q[0] = dN[0] * (ni[0]+dq);
    for (ni[1]=-Ni[1]+1; (ni[1] <= Ni[1]) && !ERROR; ++(ni[1])) {
      q[1] = dN[1] * (ni[1]+dq);
      for (ni[2]=-Ni[2]+1; (ni[2] <= Ni[2]) && !ERROR; ++(ni[2])) {
	q[2] = dN[2] * (ni[2]+dq);
	qmagn = fold_down (q, cart_rlv, qfold);
	if (TESTING) {
	  printf("## q = %8.5lf %8.5lf %8.5lf -> %8.5lf %8.5lf %8.5lf  %8.5lf\n",
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
	  for (n=lo; (n<=hi) && (!found); ++n) {
	    for (i=0; (i<Nop) && (!found); ++i)
	      found = equal_vect(qfold, gk_list[n] + 3*i);
	  }
	  --n; // now, n indexes the one we matched
	}
	// Now, if we HAVEN'T found it, we need to add it:
	if (!found) {
	  if (Nkpt == MEMORY) {
	    // not enough memory preallocated... bail out.
	    fprintf(stderr, "Fatal error: insufficient memory?\n");
	    ERROR = 1;
	    continue;
	  }
	  if (TESTING) {
	    printf("##  -- New point!\n");
	  }
	  // first, move everything up one:
	  for (n=(Nkpt-1); n>hi; --n) {
	    kpt[n+1] = kpt[n];
	    gk_list[n+1] = gk_list[n];
	  }
	  // Now, insert our point:
	  n = hi+1;
	  kpt[n] = new double[4];
	  for (d=0; d<3; ++d) kpt[n][d] = qfold[d];
	  kpt[n][3] = qmagn;
	  gk_list[n] = new double[Nop*3];
	  for (i=0; i<Nop; ++i) {
	    mult_vect(gunit[i], qfold, gk_list[n]+i*3);
	    // double qrot[3];
	    // mult_vect(gunit[i], qfold, qrot);
	    // fold_down(qrot, cart_rlv, gk_list[n]+i*3);
	  }
	  ++Nkpt;
          if (VERBOSE && (Nkpt % 10 == 0)) {
            fprintf(stderr, "%8d\n", Nkpt);
          }
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
  
  if (TESTING) {
    printf("# %d final kpoints (no weights yet)\n", Nkpt);
    for (i=0; i<Nkpt; ++i) {
      printf("# %15.12lf %15.12lf %15.12lf  %.5lf\n",
	     kpt[i][0], kpt[i][1], kpt[i][2], kpt[i][3]);
    }
  }
  
  // Now we go through and do the weights and "prettify" the kpoint listings
  if (VERBOSE) {
    fprintf(stderr, "# Weight calculation:\n");
  }
  double *w = new double[Nkpt], sum=0;
  double kdiff[3], kint[3];
  for (n=0; n<Nkpt; ++n) {
    double *tkpt = kpt[n];
    double *gk = gk_list[n];
    int count = 0;
    for (i=0; i<Nop; ++i) {
      for (d=0; d<3; ++d) {
	kdiff[d] = tkpt[d] - gk[3*i+d];
	//	kint[d] = (int)(kdiff[d]);
	kint[d] = round(kdiff[d]);
      }
      if (equal_vect(kdiff, kint)) ++count;
    }
    if (VERBOSE) {
      fprintf(stderr, "# %8.5lf %8.5lf %8.5lf   %8.5lf %d\n",
	      tkpt[0], tkpt[1], tkpt[2], tkpt[3], count);
    }
    w[n] = 1./count;
    sum += w[n];

    // prettify:
    int better;
    for (i=0; i<Nop; ++i) {
      better = (tkpt[0] < gk[3*i]);
      if (dcomp(tkpt[0], gk[3*i])) {
	better = (tkpt[1] < gk[3*i+1]);
	if (dcomp(tkpt[1], gk[3*i+1])) {
	  better = (tkpt[2] < gk[3*i+2]);
	}
      }
      if (better)
	for (d=0; d<3; ++d) tkpt[d] = gk[3*i+d];
    }
  }
  // scale 'em
  sum = 1./sum;
  for (n=0; n<Nkpt; ++n) w[n] *= sum;

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.
  if (!ERROR) {
    printf("%d 1.0 L   # number of kpoints, scaling factor, lattice coord\n",
	   Nkpt);
    for (n=0; n<Nkpt; ++n)
      printf("%18.15lf %18.15lf %18.15lf  %.14le\n", 
	     kpt[n][0], kpt[n][1], kpt[n][2], w[n]);
  }
    
  // ************************* GARBAGE COLLECTION ********************

  free_cell(Cmn_list, u_atoms, 0);

  delete[] w;
  for (n=0; n<Nkpt; ++n) {
    delete[] kpt[n];
    delete[] gk_list[n];
  }
  delete[] kpt;
  delete[] gk_list;

  return ERROR;
}
