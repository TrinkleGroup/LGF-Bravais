// This is a test-bed for our numeric evaluation of the integral:
//
// 2/Pi * int( j_2l (x) * fenv(x/a), x=0..a) = F(a,b)
// 
// numerically.  j_2l is the 2l-th order spherical bessel function.
//
// We also get to do it with cylindrical bessel functions J_2n(x),
// and using different powers of x in the integral. 
//
// both a will always be finite, though both can be large, since
// a = r k_m, where r is the distance from the origin.

// See cont-int for all the gory details of what we call when.

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>   // Bessel functions
#include <gsl/gsl_integration.h> // Numerical integration
#include <math.h>
#include "io.H"
#include "cont-int.H"

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<xmax> <Lmax>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'c'}; // Would be the flag characters.

const char* ARGEXPL = 
"  Calculates F_2l (xmax) = 2/Pi int(j_2l(x)*fenv(x/xmax), x=0..xmax) for\n\
  multiple l values from 0 to Lmax\n\
\n\
  xmax:  integration endpoint\n\
  Lmax:  maximum l value to use (goes up to 2Lmax)\n\
  -c     do with cylindrical Bessel functions instead (2d)";

int main (int argc, char **argv) {
  int i, l;

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 16;  // number of steps

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
                            VERBOSE, TESTING, MEMORY, 
                            NFLAGS, USERFLAGLIST, flagon);
  if (MEMORY < 1) {
    fprintf(stderr, "-m value must be at least 1.\n");
    ERROR = ERROR_BADFILE;
  }
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  int CYLINDRICAL = flagon[0];
  
  double xmax;
  int Lmax;
  
  sscanf(args[0], "%lf", &xmax);
  sscanf(args[1], "%d", &Lmax);

  double dx;
  dx = xmax/MEMORY;

  double fval[Lmax+1], f2val[Lmax+1];

  /*
    // Commented out by D. Trinkle... I don't remember what all of this
    // is for...
    double fval_series[Lmax+1], f2val_series[Lmax+1];
    
    // log_2n1_fact[n] = ln (1/(2*(2n)+1)!!)
    double log_2n1_fact[Lmax+5], pow2n[Lmax+5];
    
    log_2n1_fact[0] = 0.;
    pow2n[0] = 1.;
    for (l=1; l<=(Lmax+4); ++l) {
    //    log_2n1_fact[l] = log_2n1_fact[l-1] - log((4*l-1)*(4*l+1));
      log_2n1_fact[l] = log_2n1_fact[l-1] - log(16*l*l-1);
      pow2n[l] = 0.25*pow2n[l-1];
    }
    if (TESTING) {
      printf("L:  ln(1/(2l+1)!!)  1/(2l+1)!!  1/2^(l)\n");
      for (l=0; l<=(Lmax+4); ++l)
      printf("%2d %20.12le %20.12le  %20.12le\n", 2*l,
              log_2n1_fact[l], exp(log_2n1_fact[l]), pow2n[l]);
    }
    // try evaluating using just the series expansion and the
    // envelope function:
    for (l=0; l<=Lmax; ++l){
      n=2*l;
      // Currently just the FIRST term in each series expansion...
      double pre = 2./M_PI*exp((n+1)*log(xmax)+log_2n1_fact[l]);
      // leading term...
      fval_series[l] = pre*int_value(l, pow2n);
      // next term...
      fval_series[l] -= pre*xmax*xmax*int_value(l+1, pow2n)/(2*(2*n+3));
    
      // we use x^(n+2), but SAME 1/(2*n+1)!! prefactor:
      f2val_series[l] = pre*xmax*xmax*int_value(l+1, pow2n);
      f2val_series[l] -= pre*xmax*xmax*xmax*xmax
        *int_value(l+2, pow2n)/(2*(2*n+3));
    }  
  */

  double x;
  for (x=dx; x<=(xmax + 0.5*dx); x+=dx) {
    if (! CYLINDRICAL) {
      Fn_int(Lmax, x, fval, TESTING, 1024);
      Fn_x2_int(Lmax, x, f2val, TESTING, 1024);
    }
    else {
      Fn_2d_int(Lmax, x, fval, TESTING, 1024);
      Fn_2d_x2_int(Lmax, x, f2val, TESTING, 1024);
    }
    
    printf("All integrals evaluated with fenv(x/xmax), from x=0.. %.5le\n",
	   x);
    //    printf("L:    int(2/Pi*j_l(x))   int(2/Pi*x^2*j_l(x))          ln comp  ln comp\n");
    printf("L:    int(2/Pi*j_l(x))   int(2/Pi*x^2*j_l(x))\n");
    for (l=0; l<=Lmax; ++l) {
      printf("%2d %20.12le %20.12le\n", 2*l, fval[l], f2val[l]);
      //      printf("%2d %20.12le %20.12le --series:", 
      //	     2*l, fval_series[l], f2val_series[l]);
      //      printf(" %8.4lf %8.4lf\n",
      //	     ((2*l+1)*log(xmax)+log_2n1_fact[l]),
      //	     ((2*l+3)*log(xmax)+log_2n1_fact[l]));
    }
  }
  if (! CYLINDRICAL) {
    printf("Envelope integrals: %18.12le %18.12le\n", 
	   kpoint_envelope_int(), kpoint_x2_envelope_int());
  }
  else {
    printf("Envelope integrals: %18.12le %18.12le\n", 
	   kpoint_envelope_2d_int()+log(xmax), kpoint_x2_envelope_2d_int());
  }

  if (TESTING && CYLINDRICAL) {
    printf("## Hypergeometric series test:\n");
    for (x=0; x<=100; x += 4.) {
      printf("## %.5lf %.8le\n", x, hypergeom_2d_eval(x));
    }
  }

  return 0;
}
