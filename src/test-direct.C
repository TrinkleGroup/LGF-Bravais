/*
  Program: test-direct.C
  Author:  D. Trinkle
  Date:    June 3, 2004
  Purpose: I'm using this to test out some of the pieces of the algorithm
           for constraining a direct force expansion.

  Param.:  as needed...

  Flags:   MEMORY:  not used
           VERBOSE: not used
           TESTING: usual screen diahrea

  Algo.:   Testing... We're gonna try our hands at GSL's interface to
           BLAS; at some point, we'll probably have to bite the bullet
	   and really *do* BLAS and LAPACK.  Might also need to rewrite
	   the code in a more efficient / accurate manner as needed.

  Output:  Whatever I need to output.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include <gsl/gsl_blas.h>   // BLAS baby!
#include <gsl/gsl_eigen.h>  // eigensolver (replaced by LAPACK in the future)
#include <gsl/gsl_linalg.h> // LU decomp solver
#include "constraint.H"

//****************************** STRUCTURES ****************************

//***************************** SUBROUTINES ****************************

double func_eval (double x) 
{
  //  return 0.5*(1-cos(2*M_PI*x));
  return -(cos(2*M_PI*x) + 0.25*cos(4*M_PI*x) + 1./36.*cos(6*M_PI*x));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<N>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'w'};  // flag characters.

const char* ARGEXPL = 
"  N  number of points along a given direction for spline\n\
  -w weight the midpoint doubly";

const char* FILEEXPL =
"";

int main ( int argc, char **argv ) 
{
  int d, i, j; // General counting variables.

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
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  int WEIGHT = flagon[0];

  // ****************************** INPUT ****************************
  int Ni;
  int Ngrid;
  
  sscanf(args[0], "%d", &Ngrid);
  if (Ngrid<1) {
    fprintf(stderr, "Need at least one point\n");
    exit(1);
  }

  Ni = (Ngrid/2)+1; // a little weird, I know... but it's true.

  // ***************************** ANALYSIS **************************
  // 0. Construct our function to fit, including the second derivative
  // at the origin
  double f0;
  double fpp0; // our second derivative at 0 SCALED by dx^2 (deriv. wrt "i")
  double* k = new double[Ni];
  double* func = new double[Ni];
  double* weight = new double[Ni];
  // our usual function...
  f0 = func_eval(0); // subtracted off, just in case...  
  for (i=0; i<Ni; ++i) {
    k[i] = i/(double)Ngrid;
    func[i] = func_eval(k[i]) - f0;
  }
  // calculated numerically...
  double dx=0.0001;
  fpp0 = (func_eval(dx)+func_eval(-dx) - 2*f0)/(dx*dx);
  for (i=0; i<Ni; ++i) weight[i] = 1;
  if (WEIGHT && ((Ngrid%2) == 0) )
    // only if we have an even grid...
    weight[Ni-1] = 2; // weight the "edge" point more

  // 1. Compute our alpha and alpha' matrices and our a and a' vectors
  // These describe what degrees of freedom we have, and how our variables
  // determine the values and derivatives at each lattice site from 0 to N-1.
  // We have N degrees of freedom, as well.
  int Nfree = Ni;
  double* alpha = new double[Ni*Nfree];
  double* R = new double[Nfree];
  gsl_matrix_view alphaview = gsl_matrix_view_array(alpha, Ni, Nfree);
  
  for (j=0; j<Nfree; ++j)
    R[j] = 2*M_PI*j;

  // initialize to zero
  for (i=0; i<(Ni*Nfree); ++i)
    alpha[i]=0;

  for (i=0; i<Ni; ++i)
    for (j=0; j<Nfree; ++j)
      alpha[i*Nfree+j] = cos(k[i]*R[j]);

  if (TESTING) {
    printf("## Ni= %d  Nfree= %d\n", Ni, Nfree);
    printf("## f''0 = %.5le\n", fpp0);
    for (i=0; i<Ni; ++i)
      printf("## k0_%d = %8.5lf f0_%d = %+.5le\n", i, k[i], i, func[i]);
    printf("##\n");
    for (i=0; i<Ni; ++i) {
      printf("## f_%d = f(k= %.5lf)", i, k[i]);
      for (j=0; j<Nfree; ++j)
	if (! zero(alpha[i*Nfree+j])) {
	  printf(" %+.5le D(%3.1lfpi)", alpha[i*Nfree+j], R[j]/M_PI);
	}
      printf("\n");
    }
  }

  // 3. Construct our constraint matrices and vectors.
  // NOTE: New definition in use here--the constraint is
  // now Cf = c'
  int Nconst = 3; // 3 in 1d, 6 in 2d, 10 in 3d
  double* Cmat = new double[Nconst*Nfree];
  double* cvect = new double[Nconst];

  // initialize
  for (i=0; i<(Ni*Nfree); ++i) Cmat[i]=0;
  for (i=0; i<Ni; ++i) cvect[i]=0;

  // 0: D(0) = 0
  for (j=0; j<Nfree; ++j) Cmat[j] = 1;
  cvect[0] = 0;
  // 1: D'(0) = 0 -- not a real constraint in our case,
  // since we're using cosines.  But for the general case... (not
  // sure it matters there, either)
  for (j=0; j<Nfree; ++j) Cmat[Nfree+j] = 0;
  //  if (EVEN) Cmat[Nfree+(Nfree/2)] = 0;
  cvect[1] = 0;
  // 2: D''(0) = fpp0
  for (j=0; j<Nfree; ++j) 
    Cmat[2*Nfree+j] = -R[j]*R[j];
  cvect[2] = fpp0;

  if (TESTING) {
    printf("## Cmat (%d x %d) + cvect (%d):\n", Nconst, Nfree, Nconst);
    for (i=0; i<Nconst; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %6.2lf", Cmat[i*Nfree+j]);
      printf(" |   | %6.2lf |\n", cvect[i]);
    }
  }

  // Solve using constrain_fit:
  double* FC=new double[Nfree*Nconst];
  double* FE=new double[Nfree*Ni];
  constrain_fit(Nfree, Nconst, Ni, Cmat, alpha, weight, FC, FE, TESTING);

  double* fvec=new double[Nfree];
  for (i=0; i<Nfree; ++i) {
    fvec[i]=0;
    for (j=0; j<Nconst; ++j) fvec[i] += FC[i*Nconst+j]*cvect[j];
    for (j=0; j<Ni; ++j)  fvec[i] += FE[i*Ni+j]*func[j];
  }
  delete[] FC;
  delete[] FE;

  // Now, instead, solve the problem without any constraints (that
  // is, invert and solve alpha*f = func
  double* fvec_unc = new double[Nfree];
  gsl_vector_view funcview = gsl_vector_view_array(func, Ni);
  gsl_vector_view fvec_uncview = gsl_vector_view_array(fvec_unc, Ni);
  int s;
  gsl_permutation * p = gsl_permutation_alloc (Nfree);
  gsl_linalg_LU_decomp (&alphaview.matrix, p, &s);
  double det_alpha;
  det_alpha = gsl_linalg_LU_det(&alphaview.matrix, s);
  if (fabs(det_alpha) < 1e-3) {
    fprintf(stderr, "Something horrible happened... don\'t rightly know how.\n");
  }
  else {
    gsl_linalg_LU_solve (&alphaview.matrix, p, 
			 &funcview.vector, &fvec_uncview.vector);
  }

  // GC:
  gsl_permutation_free(p);

  // ****************************** OUTPUT ***************************
  // turn back into values and derivatives at each point

  if (VERBOSE) {
    printf("# fvec constrained  unconstrained:\n");
    for (j=0; j<Nfree; ++j)
      printf("# D(%+8.5lfpi) = %+.5le  %+.5le\n", 
	     R[j]/M_PI, fvec[j], fvec_unc[j]);
  }

  double finter; // interpolation value
  double finter_unc;
  for (double x=0; x<=1.0000001; x += 0.001) {
    finter = 0;
    finter_unc = 0;
    for (j=0; j<Nfree; ++j) {
      finter += fvec[j] * cos(x*R[j]);
      finter_unc += fvec_unc[j] * cos(x*R[j]);
    }
    printf("%.5lf %.8le %.8le %.8le\n", x, func_eval(x)-f0, finter, finter_unc);
  }

  // ************************* GARBAGE COLLECTION ********************
  delete[] fvec;
  delete[] fvec_unc;
  
  delete[] Cmat;
  delete[] cvect;
  
  delete[] alpha;
  delete[] R;

  delete[] func;
  delete[] k;
  delete[] weight;
  
  return 0;
}
