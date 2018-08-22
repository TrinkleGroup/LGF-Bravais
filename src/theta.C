/*
  Program: theta.C
  Author:  D. Trinkle
  Date:    2006 June 29
  Purpose: Integrate the inverse fourier transform for sgn(x) on a lattice

  Param.:  <x>
           x: value

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   direct computation

  Output:  mean and standard deviation of weakest source
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include <gsl/gsl_integration.h>  // integration routine

//****************************** STRUCTURES ****************************

//***************************** SUBROUTINES ****************************

const double A0 = 1;

typedef struct 
{
  double x;
  double u2, u4, uthresh;
} FT_type;

double FT_func (double u, void * VOIDparams) 
{
  FT_type *param = (FT_type *) VOIDparams;
  if (u< (param->uthresh) ) 
    return param->x*(2 + u*u*(param->u2 + u*u*param->u4));
  else
    return sin(param->x * u) * cos(0.5*u)/sin(0.5*u);
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<x>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {};
// flag characters.

const char* ARGEXPL = 
"  x: value";

const char* FILEEXPL =
"\n";

int main ( int argc, char **argv ) 
{
  int i; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 25;  // Our default stepping parameter

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
  // Command line parameters:
  double x;
  
  sscanf(args[0], "%lf", &x);

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  
  // Integration
  const int Nwork = 1024;
  gsl_integration_workspace * work = gsl_integration_workspace_alloc (Nwork);
  double result, error;

  gsl_function F;
  FT_type param;
  F.function = &FT_func;
  F.params = (void *)(&param);

  double xa2 = (x*x)/(A0*A0);
  param.x = x/A0;
  param.u2 = - (1 + xa2*2)/6;
  param.u4 = - (1 + xa2*(10*-xa2*6))/360;
  param.uthresh = (1 + xa2*(-7 + xa2*(21 + xa2*6)))/15120;
  param.uthresh = exp( log(1e-8/fabs(param.uthresh))/6. );


  double result0;
  result = 0;
  double x0, x1 = 0;
  double dx = M_PI*A0/x;
  do {
    x0 = x1;
    x1 = x0 + dx;
    if (x1>M_PI) x1 = M_PI;
    gsl_integration_qag(&F, x0, x1, 1e-7, 1e-7, Nwork, GSL_INTEG_GAUSS61,
			work, &result0, &error); 
    result0 *= 1./M_PI;
    if (TESTING) 
      printf("### int(%.8le..%.8le)= %.12le +- %.12le\n", x0, x1, 
	     result0, error);
    result += result0;
  } while (x1 < M_PI);
  
  printf("%.8le %.12lf\n", x, result);
  
  // ************************* GARBAGE COLLECTION ********************

  gsl_integration_workspace_free(work);
  return 0;
}
