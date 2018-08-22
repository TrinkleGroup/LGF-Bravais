/*
  Program: test-tricubic.C
  Author:  D. Trinkle
  Date:    June 3, 2004
  Purpose: I'm using this to test out some of the pieces of the algorithm
           for doing a tricubic spline.  It kinda looks like a little
	   nightmare... hopefully it won't be so bad when it's all said
	   and done.

  Param.:  as needed...

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Testing...

  Output:  Whatever I need to output.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "pointgroup.H"  // Do we still use this?
#include "Dij.H"
#include "kpts.H"
#include "shell.H"
#include "tricubic.H"

//****************************** STRUCTURES ****************************

// our function to fit:
inline double f0(double x) 
{
  return 0.5*(1-cos(x));
}


//***************************** SUBROUTINES ****************************

// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  //  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
  printf("xx= %15.8le yy= %15.8le zz= %15.8le  xy= %15.8le yz= %15.8le zx= %15.8le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<N0> <N1> <N2>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  N  number of points along a given direction for spline\n";

const char* FILEEXPL =
"";

int main ( int argc, char **argv ) 
{
  int d, i, j, k; // General counting variables.

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
  int Ni[3];

  // Let's pull off the args:
  sscanf(args[0], "%d", Ni);
  sscanf(args[0], "%d", Ni+1);
  sscanf(args[0], "%d", Ni+2);
  
  if ( (Ni[0] < 1) && (Ni[1] < 1) && (Ni[2] < 1) ) {
    fprintf(stderr, "To get anything useful, Ni>0\n");
    exit(1);
  }

  // ***************************** ANALYSIS **************************

  // Make a simple function in 3d using Ni as a grid size:
  int N=Ni[0]*Ni[1]*Ni[2];
  int*** index, n[3], ind;
  int Ndim=1;
  double** func = new double*[Ndim];
  for (d=0; d<Ndim; ++d) func[d] = new double[N];

  double alpha[3]={2*M_PI/(double)Ni[0], 2*M_PI/(double)Ni[1],
		   2*M_PI/(double)Ni[2]};
  
  
  index = new int**[Ni[0]];
  for (ind=0, n[0]=0; n[0]<Ni[0]; ++(n[0])) {
    index[n[0]] = new int*[Ni[1]];
    for (n[1]=0; n[1]<Ni[1]; ++(n[1])) {
      index[n[0]][n[1]] = new int[Ni[1]];
      for (n[2]=0; n[2]<Ni[2]; ++(n[2]), ++ind) {
	//	func[ind] = f0(alpha*n[0])*f0(alpha*n[1])*f0(alpha*n[2]);
	for (d=0; d<Ndim; ++d) {
	  func[d][ind] = f0(alpha[2-d]*n[2-d]); // run z,y,x
	  //	  func[d][ind] = f0(alpha[0]*n[0]) + f0(2*alpha[1]*n[1])
	  //	    + f0(3*alpha[2]*n[2]);
	}
	index[n[0]][n[1]][n[2]] = ind;
      }
    }
  }
  tricubic_type* t;
  calc_tricubic(Ni[0], Ni[1], Ni[2], Ndim, index, func, t);

  // ****************************** OUTPUT ***************************

  const double dz = 0.001;
  
  for (double z=0; z<1.000001; z += dz) {
    double x[3] = {z,z,z};
    printf("%.8lf", z);
    for (d=0; d<Ndim; ++d)
      printf(" %.15lf", eval_tricubic(x, t[d]));
    // test derivatives as well:
    for (d=0; d<Ndim; ++d)
      printf(" %.15lf", eval_tricubic(x, 1,0,0, t[d]));
    for (d=0; d<Ndim; ++d)
      printf(" %.15lf", eval_tricubic(x, 0,1,0, t[d]));
    for (d=0; d<Ndim; ++d)
      printf(" %.15lf", eval_tricubic(x, 0,0,1, t[d]));
    printf("\n");
  }

  if (TESTING) {
    for (d=0; d<3; ++d) {
      double* HinvDel;
      calc_HinvDel(Ni[d], HinvDel);

      printf("# HinvDel (Ni=%d):\n", Ni[d]);
      printf("#     ");
      for (j=0; j<Ni[d]; ++j) printf("%5d   ", j);
      printf("\n");
      for (i=0; i<Ni[d]; ++i) {
	printf("# %2d |", i);
	for (j=0; j<Ni[d]; ++j) printf("%8.4lf", HinvDel[i*Ni[d]+j]);
	printf(" |\n");
      }
      
      delete[] HinvDel;
    }
  }
  

  if (VERBOSE) {
    write_tricubic(stdout, t, Ndim);
  }
  

  // ************************* GARBAGE COLLECTION ********************
  for (n[0]=0; n[0]<Ni[0]; ++(n[0])) {
    for (n[1]=0; n[1]<Ni[1]; ++(n[1])) 
      delete[] index[n[0]][n[1]];
    delete[] index[n[0]];
  }
  delete[] index;
  for (d=0; d<Ndim; ++d) delete[] func[d];
  delete[] func;
  
  for (d=0; d<Ndim; ++d)
    free_tricubic(t[d]);
  delete[] t;

  return 0;
}
