/*
  Program: test-disc.C
  Author:  D. Trinkle
  Date:    April 18, 2004
  Purpose: My effort at testing some of the aspects involved in computing
           the discretization correction directly from the dynamical
	   matrix.  We use the simple 1D example, where

	          / -1: |x| = a
	   D(x) = |  2:  x  = 0
                  \  0: |x| > a

	   We choose a=1 to make everything unitless; then k goes from
	   -Pi to Pi.  The elastic GF is the FT of 1/k^2: 2Pi/r, the
	   dynamical matrix in k-space is 2(1-cos(k)), and so the
	   lattice GF is the FT of 0.5/(1-cos(k)), whatever *that* is.

	   We evaluate the discretization correction, which is
	   the FT of [D(k)]^-1 - GE*fenv(k), where fenv(k) is an envelope
	   function that takes k smoothly to 0 at Pi.

  Param.:  

  Flags:   MEMORY:  not used
           VERBOSE: ??
           TESTING: output practically everything as we do it.

  Algo.:   

  Output:  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "io.H"

const double kpoint_envelope_alpha = 0.5;
const double kpoint_envelope_alpha_1 = 1./(1.-kpoint_envelope_alpha);
inline double kpoint_envelope (double x) 
{
  double y = 1.-(1.-x)*kpoint_envelope_alpha_1;
  if (x <= kpoint_envelope_alpha) return 1;
  if (x >= 1.0) return 0;
  return 1 + y*y*(-3 + 2*y);
  // return 1 + x*x*(-3 + 2*x);
}

const double kmax_1 = 1./M_PI; // 1/kmax

// fixed to handle behavior near 0.
double func_noenv(double k, void * params) 
{
  if (fabs(k) < 1e-3) return 1./12. + 1./240.*k*k;
  return 0.5/(1.-cos(k)) - 1/(k*k);
}

double func_env(double k, void * params) 
{
  if (fabs(k) < 1e-3) return 1./12. + 1./240.*k*k;
  return 0.5/(1.-cos(k)) - kpoint_envelope(fabs(k)*kmax_1)/(k*k);
}


// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<Nk> <NR>";

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'e', 'c', 'i'}; // flag characters.

const char* ARGEXPL = 
"  Nk: number of kpoints in integration\n\
  NR: number of R points to test\n\
  -e  no envelope function\n\
  -c  use the linearized correction term\n\
  -i  compare with exact integration\n\
  MEMORY parameter: determines number of subdivisions in oscillatory int.";

int main (int argc, char **argv) 
{
  int i, j, R;

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// used for numerical intergration only
  int NOENV = 0;
  int CORREC = 0;
  int EXACTINT = 0;
    
  char* args[argc];
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
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  NOENV = flagon[0];
  CORREC = flagon[1];
  EXACTINT = flagon[2];

  // ***************************** INPUT *****************************
  int Nk, NR;
  sscanf(args[0], "%d", &Nk);
  sscanf(args[1], "%d", &NR);
  
  // integration...
  gsl_integration_workspace *work;
  gsl_integration_qawo_table *t;
  gsl_function fint;
  if (EXACTINT) {
    work = gsl_integration_workspace_alloc (MEMORY);
    t = gsl_integration_qawo_table_alloc(0, 2.*M_PI, GSL_INTEG_COSINE,
					 MEMORY);
    if (NOENV) fint.function = &func_noenv;
    else       fint.function = &func_env;
    fint.params = NULL;
  }
  

  double Gd, Dk, Gek;
  double k, w;
  double Nk_1 = 1./(double)Nk;
  double dk = M_PI*Nk_1;
  printf("Nk = %d   NR = %d\n", Nk, NR);
  for (R=0; R<=NR; ++R) {
    Gd = 0;
    double Gdp = 0, Gdn = 0;
    double correc = 1; // default value
    if ( (CORREC) && (R!=0) ) correc = sin(dk*R)/(dk*R);
    for (i=0; i<=Nk; ++i) {
      k = dk*i;
      // calculate the weight; all points except 0 and Pi get
      // the same weight--these points are weighted by 1/2 due
      // to symmetry.
      if ( (i==Nk) || (i==0) ) w = 0.5*Nk_1;
      else                     w = Nk_1;
      
      double coskR = w*correc*cos(k*R);
      double func;
      Dk = 0.5/(1.-cos(k));
      Gek = 1./(k*k);
      if (NOENV)
	func = coskR*func_noenv(k, NULL);
      else
	func = coskR*func_env(k, NULL);
      if (func >= 0) Gdp += func;
      else           Gdn += func;
    }
    Gd = Gdp+Gdn;
    printf("%d %.8le", R, Gd);

    if (EXACTINT) {
      // do the exact integration numerically...
      gsl_integration_qawo_table_set(t, (double)R, 2.*M_PI, GSL_INTEG_COSINE);
      double res, err;
      ERROR = gsl_integration_qawo(&fint, -M_PI, 0, 2e-5, MEMORY, work, t,
			     &res, &err);
      if (ERROR) {
	fprintf(stderr, "Error in computing the exact integral; try with larger MEMORY value.\n");
	res = 0;
      }
      res *= 1./(2.*M_PI);
      printf(" %.8le %.8le\n", res, fabs((Gd-res)/res));
    }
    else {
      printf("\n");
    }
  }

  // Garbage collection
  if (EXACTINT) {
    gsl_integration_workspace_free(work);
    gsl_integration_qawo_table_free(t);
  }
  
  return 0;
}

    
