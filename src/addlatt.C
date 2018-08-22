/*
  Program: addlatt.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Utility code to read in two lattice functions and add 'em

  Param.:  <latt.1(R)> <latt.2(R)>
           latt.i(R): lattice functions

           ==== latt(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij
           ...
           nN.1 nN.2 nN.3  Lxx .. Lzz
           ==== latt(R) ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and calculate.

  Output:  We output the points in unit coord.  (Switch for lattice coord?)
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "cell.H"
#include "Dij.H"


//***************************** SUBROUTINES ****************************

int __COMPARE_LATT_TYPE__ (const void* a_elem, const void* b_elem) 
{
  point_type *x = (point_type*)a_elem;
  point_type *y = (point_type*)b_elem;
  int a, b, comp=0, d;
  if (dcomp(x->Rmagn, y->Rmagn)) {
    for (d=0; (d<3) && (!comp); ++d) {
      a = abs(x->Runit[d]);  b = abs(y->Runit[d]);
      if (a<b) comp = -1;
      else
	if (a>b) comp = 1;
	else {
	  a = x->Runit[d];  b = y->Runit[d];
	  if (a<b) comp = -1;
	  else
	    if (a>b) comp = 1;
	    else comp = 0;
	}
    }
  }
  else {
    if (x->Rmagn < y->Rmagn) comp = -1;
    else                     comp = 1;
  }
  return comp;
}
  
int min_point(int Nlatt, point_type **tp, int* use) 
{
  int i, j;
  int comp, best;
  
  // initialize "best" by skipping past the NULL pointers.
  for (best=0; (best<Nlatt) && (tp[best] == NULL); ++best) 
    use[best] = 0;
  for (i=best; i<Nlatt; ++i) {
    if (tp[i] == NULL) 
      use[i] = 0;
    else {
      comp = __COMPARE_LATT_TYPE__((void*)(tp[i]), (void*)(tp[best]));
      if (comp == 0) {
	use[i] = 1; // equal points:
      }
      else {
	if (comp < 0) {
	  // new best point:
	  for (j=i-1; j>=0; --j) use[j] = 0; // zero out old points
	  best = i;
	  use[i] = 1;
	}
	else
	  use[i] = 0; // worse than best point
      }
    }
  }
  if (best==Nlatt) return -1;
  else return best;
}


inline void print_mat (double mat[9]) 
{
  printf("%.3le %.3le %.3le %.3le %.3le %.3le %.3le %.3le %.3le",
	 mat[0], mat[1], mat[2], mat[3], mat[4], mat[5], 
	 mat[6], mat[7], mat[8]);
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const char* ARGLIST = "<latt.i(R)> ...";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  latt.i(R): lattice functions";

const char* FILEEXPL =
"==== latt(R) ====\n\
N kmax                      # number of points, kmax\n\
n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij\n\
...\n\
nN.1 nN.2 nN.3  Lxx .. lzz\n\
==== latt(R) ====\n";

int main ( int argc, char **argv ) 
{
  int i, d, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used

  int NUMARGS = -1; // variable number of arguments
  char* args[argc];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (i=0; i<NFLAGS; ++i) flagon[i] = 0;
  // Read our commandline.
  ERROR = parse_commandline_var(argc, argv, NUMARGS, args,
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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters: 
  int Nlatt = NUMARGS;
  char *latt_name[Nlatt];
  
  // Let's pull off the args:
  for (i=0; i<Nlatt; ++i) 
    latt_name[i] = args[i];

  //++ ==== latt(R) ====
  int Np[Nlatt], Ndim = 9;
  point_type *p[Nlatt]; // set of points

  for (i=0; i<Nlatt; ++i) {
    infile = myopenr(latt_name[i]);
    if (infile == NULL) {
      fprintf(stderr, "Couldn't open %s for reading.\n", latt_name[i]);
      exit(ERROR_NOFILE);
    }
    double cart[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // really just a dummy
       
    ERROR = read_Dij (infile, cart, Np[i], p[i], READ_DIJ_MAT);
    myclose(infile);
  }
  
  if (ERROR) {
    for (i=0; i<Nlatt; ++i) delete p[i];
    exit(ERROR);
  }
  //-- ==== latt(R) ====

  // ***************************** ANALYSIS **************************
  // To facilitate the combination of all of our lattice functions,
  // we first sort each set:
  for (i=0; i<Nlatt; ++i)
    qsort((void*)(p[i]), Np[i], sizeof(point_type), __COMPARE_LATT_TYPE__);

  if (TESTING) {
    for (i=0; i<Nlatt; ++i) {
      printf("## ==== lattice function %d ====\n", i);
      for (n=0; n<Np[i]; ++n) {
	printf("# %3d %3d %3d ", p[i][n].Runit[0], p[i][n].Runit[1], 
	       p[i][n].Runit[2]);
	print_mat(p[i][n].mat);
	printf("\n");
      }
      printf("\n");
    }
  }
  
  point_type *sum;
  int Nsum=0;
  for (i=0; i<Nlatt; ++i) Nsum += Np[i];
  sum = new point_type[Nsum];

  int np[Nlatt];  // counter for each...
  point_type **tp = 
    new point_type *[Nlatt]; // a pointer to the current point in each
  int *use = new int[Nlatt]; // index which points we need to use
  int more = 1;   // do we have more points to do?

  // initialize:
  for (i=0; i<Nlatt; ++i) {
    np[i] = 0;
    if (np[i]<Np[i])
      tp[i] = p[i];
    else
      tp[i] = NULL;
  }
  Nsum = 0;
  while (more) {
    n = min_point(Nlatt, tp, use);
    if (n<0) {
      fprintf(stderr, "Weird... some error in min_point; shouldn't happen.\n");
      for (i=0; i<Nlatt; ++i)
	fprintf(stderr, "np[%d] = %d  Np[%d] = %d\n",
		i, np[i], i, Np[i]);
      ERROR = -1;
      more = 0;
    }
    else {
      for (d=0; d<3; ++d)
	sum[Nsum].Runit[d] = tp[n]->Runit[d];
      for (d=0; d<9; ++d) 
	sum[Nsum].mat[d] = 0.;
      for (i=0; i<Nlatt; ++i) {
	if (use[i]) {
	  for (d=0; d<9; ++d)
	    sum[Nsum].mat[d] += tp[i]->mat[d];
	  // increment counter:
	  ++(np[i]);
	  if (np[i]<Np[i]) ++(tp[i]);
	  else             tp[i] = NULL;
	}
      }
      ++Nsum;
      // if we have any more points...
      for (more=0,i=0; (i<Nlatt) && (!more); ++i)
	more = (np[i]<Np[i]);
    }
  }
  
  // ****************************** OUTPUT ***************************
  write_Dij(stdout, Nsum, sum);
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] tp;
  delete[] use;
  for (i=0; i<Nlatt; ++i)
    delete[] p[i];
  delete[] sum;
  
  return ERROR;
}
