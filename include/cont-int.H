#ifndef __CONT_INT__
#define __CONT_INT__

// Our routines to numerically evaluate the integrals:
//
// 2/Pi * int( j_2l (x) * exp( -(x*b)^2 ), x=0..a) = F(a,b)
//  and
// 2/Pi * int( k_2l (x) * (1 - 3(x/a)^2 + 2(x/a)^3), x=0..a) = Fn(a)
// 
// numerically.  j_2l is the 2l-th order spherical bessel function.
//
// both a and b will always be finite, though both can be large, since
// a = r k_m, and b = Rc / r where r is the distance from the origin.
// We should only use this *if* r is non-zero.  If r *is* zero, there's
// a whole different integration to do, which is much easier.

// We proceed to evaluate these integrals by finding all of the roots
// of j_2l(x) = 0 that are less than a.  These are exactly the
// roots of J_(2l+1/2)(x), and are determined with gsl_sf_bessel_zero_Jnu.
// Between two roots, we use the gsl_integration_qag to integrate;
// it's a simple adaptive integration procedure.

// This is used by semicont; it is part of the continuous piece.

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>   // Bessel functions
#include <gsl/gsl_integration.h> // Numerical integration
#include <math.h>
#include "io.H"

// ========================= ENVELOPE FUNCTION =========================

// We can use this outside of our integration...
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

// ===================== THREE-DIMENSIONAL INTEGRALS =================

// **** Integration type switches:
const int CONT_INT_UNITY = 0;
const int CONT_INT_SQUARED = 1;
// **** Integration functions:
double kpoint_envelope_int (int INT_TYPE);

inline double kpoint_envelope_int () 
  { return kpoint_envelope_int(CONT_INT_UNITY); }
inline double kpoint_x2_envelope_int () 
  { return kpoint_envelope_int(CONT_INT_SQUARED); }

// our default exponential threshold; basically means values below 1e-18
// will just return 0.
const double base_threshold=-32.;
  
void Fn_int (int Lmax, double aparam, double *fval,
	     int TESTING, int MEMORY, int INT_TYPE, double EXP_THRESH);

inline void Fn_int (int Lmax, double aparam, double *fval,
		    int TESTING, int MEMORY) 
  { Fn_int(Lmax, aparam, fval, TESTING, MEMORY, CONT_INT_UNITY,
	   base_threshold); }
inline void Fn_x2_int (int Lmax, double aparam, double *fval,
		    int TESTING, int MEMORY) 
  { Fn_int(Lmax, aparam, fval, TESTING, MEMORY, CONT_INT_SQUARED,
	   base_threshold); }
  

// **** Code to do integrations...

// Integral of kpoint_envelope for x=0..1
double __kpoint_integrand_function__ (double x, void * params) 
{
  return kpoint_envelope(x);
}
// Integral of x^2*kpoint_envelope for x=0..1
double __kpoint_x2_integrand_function__ (double x, void * params) 
{
  return x*x*kpoint_envelope(x);
}

double kpoint_envelope_int (int INT_TYPE) 
{
  // Integration pieces
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1024);
  double result, error;

  // Make our function:
  gsl_function fint;

  if (INT_TYPE == CONT_INT_UNITY)
    fint.function = &__kpoint_integrand_function__;
  else
    fint.function = &__kpoint_x2_integrand_function__;
  fint.params = NULL; 

  gsl_integration_qag(&fint, 0, 1, 0, 1e-9, 1024, 
		      GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free(w);
  
  return result;
}


typedef struct 
{
  int n;
  double xm1;
} __Fn_param_type__;

double __Fn_integrand_function__ (double x, void * params) 
{
  int n = ((__Fn_param_type__ *)params)->n;
  double y = (x*((__Fn_param_type__ *)params)->xm1);
  
  return gsl_sf_bessel_jl(n, x)*kpoint_envelope(y);
}

double __Fn_x2_integrand_function__ (double x, void * params) 
{
  int n = ((__Fn_param_type__ *)params)->n;
  double y = (x*((__Fn_param_type__ *)params)->xm1);
  
  return x*x*gsl_sf_bessel_jl(n, x)*kpoint_envelope(y);
}


// NEWLY MODIFIED!
// We attempt (and I think succeed) at taking out those cases which
// would suffer from unrecoverable underflow.  We simply return 0 for
// those elements which we predict to be below a given exponential
// threshold; I recommend -32 (see base_threshold above)
int Lmax_calc=0;
const int Lmax_UPPER=100; // we'll never actually need this many...
// log_2n1_fact[n] = ln (1/(2*(2n)+1)!!)
double log_2n1_fact[Lmax_UPPER+5]={0.}; // our initial value

void init_log_2n1_fact(int Lmax) 
{
  if (Lmax > Lmax_UPPER) {
    fprintf(stderr, "Lmax=%d larger than Lmax_UPPER=%d\n", Lmax, Lmax_UPPER);
    fprintf(stderr, "please choose lower Lmax or recompile.\n");
    exit(-1);
  }
  if (Lmax<=Lmax_calc) return ;
  // fill in the values that have not been previously filled...
  for (int l=(Lmax_calc+1); l<=Lmax; ++l)
    log_2n1_fact[l] = log_2n1_fact[l-1] - log(16*l*l-1);
  Lmax_calc=Lmax;
}

void Fn_int (int Lmax, double aparam, double *fval,
	     int TESTING, int MEMORY, int INT_TYPE, double EXP_THRESH)
{
  int l, s;
  
  if (TESTING) {
    printf("## a = %.5le\n", aparam);
    printf("## Lmax = %d\n", Lmax);
  }
  
  if ((aparam < 0) || (Lmax < 0)) {
    fprintf(stderr, "Questionable a or Lmax values.\n");
    return;
  }

  double log_inc; // just to make comparison with unity or x^2  
  init_log_2n1_fact(Lmax); // just make sure we have everything we need.

  // Integration pieces
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (MEMORY);
  double result, error;

  // Make our function:
  gsl_function fint;
  __Fn_param_type__ fn_param;

  if (INT_TYPE == CONT_INT_UNITY) {
    fint.function = &__Fn_integrand_function__;
    log_inc = 1.;
  }
  else {
    fint.function = &__Fn_x2_integrand_function__;
    log_inc = 3.;
  }
  fint.params = &fn_param; // Now, if we change fn_param, it goes to fint.

  // 1. Determine which l's are above our EXP_THRESH value:
  int Lthresh, above=1;
  for (Lthresh=0; (Lthresh<=Lmax) && (above); ++Lthresh) {
    above = ( ((2*Lthresh+log_inc)*log(aparam)) > 
	      (EXP_THRESH-log_2n1_fact[Lthresh]) );
  }
  --Lthresh; // we always overshoot by one...
  
  double x0, x1; // Our integration limits
  const double Fscale = 2./M_PI;
  double fvalp, fvaln; // Our positive and negative values
  
  // Set our parameters:
  fn_param.xm1 = 1./aparam; 
  for (l=0; l<=Lthresh; ++l) {
    if (TESTING) {
      printf("## ++ integrating 2l = %d\n", 2*l);
    }
    // Set our parameters:
    fn_param.n = 2*l;
    // Initialization:
    x0 = 0;
    s = 1; // which root are we on...
    fvalp = 0;
    fvaln = 0;
    do {
      // Get our end limit...
      x1 = gsl_sf_bessel_zero_Jnu(2*l+0.5, s);
      if ( (x1+1e-2) > aparam) x1 = aparam; // Never go beyond the end
      // Integrate!
      // a, b, epsabs, epsrel, worksize, key, workspace, result, abserr
      //  --key determines which Gauss-Kronrod rule to use.
      //      gsl_integration_qag(&fint, x0, x1, 1e-18, 1e-9, MEMORY, 
      gsl_integration_qag(&fint, x0, x1, 0, 1e-9, MEMORY, 
      			  GSL_INTEG_GAUSS61, w, &result, &error);
      if (TESTING) {
	printf("## integrate x = %.5le .. %.5le = %.5le +- %.5le\n", 
	       x0, x1, result, error);
      }
      // add it to the tally.      
      if (result > 0) fvalp += result;
      else            fvaln += result;
      // Next limits:
      x0 = x1;
      ++s;      
    } while (x1 < aparam);
    fval[l] = Fscale*(fvalp+fvaln);
    if (TESTING) {
      printf("## final: %.5le\n", fval[l]);
    }
  }

  for (l=(Lthresh+1); l<=Lmax; ++l) {
    fval[l] = 0.;
    if (TESTING) {
      printf("## ++ integrating 2l = %d -- below threshold\n", 2*l);
      printf("## final: %.5le\n", fval[l]);
    }
  }
  // Garbage collection:
  gsl_integration_workspace_free(w);
}


// ====================== TWO-DIMENSIONAL INTEGRALS ==================
double kpoint_envelope_2d_int (int INT_TYPE);

inline double kpoint_envelope_2d_int () 
{ return kpoint_envelope_2d_int(CONT_INT_UNITY); }
inline double kpoint_x2_envelope_2d_int () 
{ return kpoint_envelope_2d_int(CONT_INT_SQUARED); }


void Fn_2d_int (int Nmax, double aparam, double *fval,
		int TESTING, int MEMORY, int INT_TYPE);

inline void Fn_2d_int(int Nmax, double aparam, double *fval,
		      int TESTING, int MEMORY) 
{ Fn_2d_int(Nmax, aparam, fval, TESTING, MEMORY, CONT_INT_UNITY);}

inline void Fn_2d_x2_int(int Nmax, double aparam, double *fval,
			 int TESTING, int MEMORY) 
{ Fn_2d_int(Nmax, aparam, fval, TESTING, MEMORY, CONT_INT_SQUARED);}



// Integral of kpoint_envelope for x=0..1
double __kpoint_2d_integrand_function__ (double x, void * params) 
{
  return kpoint_envelope(x)/x;
}
// Integral of x^2*kpoint_envelope for x=0..1
double __kpoint_2d_x2_integrand_function__ (double x, void * params) 
{
  return x*kpoint_envelope(x);
}

double kpoint_envelope_2d_int (int INT_TYPE) 
{
  // Integration pieces
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1024);
  double result, error;

  // Make our function:
  gsl_function fint;
  double x0;

  if (INT_TYPE == CONT_INT_UNITY) {
    fint.function = &__kpoint_2d_integrand_function__;
    x0 = kpoint_envelope_alpha;
  } else {
    fint.function = &__kpoint_2d_x2_integrand_function__;
    x0 = 0.;
  }
  fint.params = NULL; 

  gsl_integration_qag(&fint, x0, 1, 0, 1e-9, 1024, 
                      GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free(w);
  
  if (INT_TYPE == CONT_INT_UNITY)
    // Note: the user needs to add in log(kmax):
    result += log(kpoint_envelope_alpha) + M_EULER - M_LN2;
  return result;
}



// returns Jn(x) / x, but is careful to avoid the problem with
// limit x->0
inline double bessel_Jn_divx (const int &n, const double &x) 
{
  if (x < 1e-5) {
    // return 1/2 * (x/2)^(n-1)
    return 0.5*gsl_pow_int(0.5*x, n-1);
  }
  // else...
  return gsl_sf_bessel_Jn(n, x)/x;
}


double __Fn_2d_integrand_function__ (double x, void * params) 
{
  int n = ((__Fn_param_type__ *)params)->n;
  double y = (x*((__Fn_param_type__ *)params)->xm1);
  
  return bessel_Jn_divx(n, x)*kpoint_envelope(y);
}

double __Fn_2d_x2_integrand_function__ (double x, void * params) 
{
  int n = ((__Fn_param_type__ *)params)->n;
  double y = (x*((__Fn_param_type__ *)params)->xm1);
  
  return x*gsl_sf_bessel_Jn(n, x)*kpoint_envelope(y);
}

// evaluate our hypergeometric function:
//             n  n+1
// 1 inf   (-1)  a
// - SUM --------------
// 2 n=0 (n+1) (n+1)!^2
//
// ... where a = 1/4 * x^2.
//
// This is done by summing over every N terms, where N is even
// and N is approx x.  We group the terms as polynomials in a.
// We do this to (hopefully) avoid numerical errors for large x values

// The current upper limit for accuracy of the series seems to be x=35,
// which has an inaccuracy of 1e-5.  By x=40, the error is 4e-3, and by 
// x=43 the error is so large that the curve is actually no longer growing
// with x.

// To alleviate *that* error, we only use the series expansion up to
// x=30.  After that, we do numerical integration of the function
// 1/x(1-J0(x)) (which is the derivative of the hypergeom. function).
// We do that in Pi-sized chunks to take into account the oscillatory
// nature of J0(x), even though there are no roots of the derivative.

const double MAX_HYPER = 30;
const double INT_STEP_HYPER = M_PI;

double __hypergeometric_integrand_function__ (double x, void *param) 
{ return (1-gsl_sf_bessel_J0(x))/x; }

double hypergeom_2d_eval (double x) 
{
  int i, n, N, nN;
  double sum, term, baseterm;
  double prod;
  double x0;
  // we *should* use a factor of 2 here, but to extend the
  // applicability, it helps to use larger values here
  
  if (x > MAX_HYPER) x0 = MAX_HYPER;
  else               x0 = x;
  double a = x0*x0*0.25;
  if (x0 < 2)
    N = 2*(int)(0.5*sqrt(a)); // makes N even
  else
    N = 2*(int)(0.5*log(x0)*sqrt(a)); // produce larger N values...
  if (N == 0) N=2;
  sum = 0;
  n=0;
  baseterm = 1;
  do {
    nN = n*N;
    prod = 1./(double)(nN + N);
    for (i=N-1; i>0; --i)
      prod = 1/(double)(nN+i) - prod*a/gsl_pow_2( (double)(nN+i+1));

    term = baseterm*prod;
    sum += term;
    // now build up "base term" for next time through:
    for (i=2; i<=(N+1); ++i)
      baseterm *= a/gsl_pow_2( (double)(nN+i));
    ++n;
  } while ( fabs(term) > fabs(sum*1e-12) );
  sum *= a*0.5;
  // Now, see if we have to do numerical integration for the rest:
  if (x > x0) {
    // Integration pieces
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1024);
    double result, error, x1;
    // Make our function:
    gsl_function fint;
    fint.function = &__hypergeometric_integrand_function__;
    fint.params = NULL; 

    // work our way through, step by step:
    x1 = x0;
    do {
      x0 = x1;
      x1 += INT_STEP_HYPER;
      if (x1 > x) x1 = x;
      gsl_integration_qag(&fint, x0, x1, 0, 1e-9, 1024, 
			  GSL_INTEG_GAUSS61, w, &result, &error);
      sum += result;
    } while (x1 < x);
    // garbage collection:
    gsl_integration_workspace_free(w);
  }

  return sum;
}

void Fn_2d_int (int Nmax, double aparam, double *fval,
		int TESTING, int MEMORY, int INT_TYPE)
{
  int n, s;
  
  if (TESTING) {
    printf("## a = %.5le\n", aparam);
    printf("## Nmax = %d\n", Nmax);
  }
  
  if ((aparam < 0) || (Nmax < 0)) {
    fprintf(stderr, "Questionable a or Lmax values.\n");
    return;
  }

  // Integration pieces
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (MEMORY);
  double result, error;

  // Make our function:
  int nstart;
  gsl_function fint;
  __Fn_param_type__ fn_param;

  if (INT_TYPE == CONT_INT_UNITY) {
    fint.function = &__Fn_2d_integrand_function__;
    nstart=1; // we DON'T integrate the n=0 case like everything else
  }
  else {
    fint.function = &__Fn_2d_x2_integrand_function__;
    nstart=0;
  }
  fint.params = &fn_param; // Now, if we change fn_param, it goes to fint.

  double x0, x1; // Our integration limits
  const double Fscale = 1.;
  double fvalp, fvaln; // Our positive and negative values
  
  // Set our parameters:
  fn_param.xm1 = 1./aparam;
  if (nstart == 1) {
    if (TESTING) {
      printf("## ++ integrating 2n = 0  (special case)\n");
    }
    // we need to handle the n=0 case specially.
    // 0. constant factors:
    fval[0] = log(kpoint_envelope_alpha*aparam) - M_LN2 + M_EULER;
    // NOTE: the user must subtract log(R) from this to be correct.

    // 1. numerical integration:
    // Set our parameters:
    fn_param.n = 0;
    // Initialization:
    x0 = kpoint_envelope_alpha*aparam;
    // which root are we on...
    for (s=1; gsl_sf_bessel_zero_J0(s) < x0; ++s) ;
    fvalp = 0;
    fvaln = 0;
    do {
      // Get our end limit...
      x1 = gsl_sf_bessel_zero_J0(s);
      if ( (x1+1e-2) > aparam) x1 = aparam; // Never go beyond the end
      // Integrate!
      // a, b, epsabs, epsrel, worksize, key, workspace, result, abserr
      //  --key determines which Gauss-Kronrod rule to use.
      //      gsl_integration_qag(&fint, x0, x1, 1e-18, 1e-9, MEMORY, 
      gsl_integration_qag(&fint, x0, x1, 0, 1e-9, MEMORY, 
      			  GSL_INTEG_GAUSS61, w, &result, &error);
      if (TESTING) {
	printf("## integrate x = %.5le .. %.5le = %.5le +- %.5le\n", 
	       x0, x1, result, error);
      }
      // add it to the tally.      
      if (result > 0) fvalp += result;
      else            fvaln += result;
      // Next limits:
      x0 = x1;
      ++s;      
    } while (x1 < aparam);
    fval[0] += fvalp+fvaln;

    // 2. hypergeometric series:
    fval[0] += -hypergeom_2d_eval (kpoint_envelope_alpha*aparam);
    if (TESTING) {
      printf("## final value: %.5le\n", fval[0]);
    }
  }
 
  for (n=nstart; n<=Nmax; ++n) {
    if (TESTING) {
      printf("## ++ integrating 2n = %d\n", 2*n);
    }
    // Set our parameters:
    fn_param.n = 2*n;
    // Initialization:
    x0 = 0;
    s = 1; // which root are we on...
    fvalp = 0;
    fvaln = 0;
    do {
      // Get our end limit...
      x1 = gsl_sf_bessel_zero_Jnu(2*n, s);
      if ( (x1+1e-2) > aparam) x1 = aparam; // Never go beyond the end
      // Integrate!
      // a, b, epsabs, epsrel, worksize, key, workspace, result, abserr
      //  --key determines which Gauss-Kronrod rule to use.
      //      gsl_integration_qag(&fint, x0, x1, 1e-18, 1e-9, MEMORY, 
      gsl_integration_qag(&fint, x0, x1, 0, 1e-9, MEMORY, 
      			  GSL_INTEG_GAUSS61, w, &result, &error);
      if (TESTING) {
	printf("## integrate x = %.5le .. %.5le = %.5le +- %.5le\n", 
	       x0, x1, result, error);
      }
      // add it to the tally.      
      if (result > 0) fvalp += result;
      else            fvaln += result;
      // Next limits:
      x0 = x1;
      ++s;      
    } while (x1 < aparam);
    fval[n] = Fscale*(fvalp+fvaln);
    if (TESTING) {
      printf("## final: %.5le\n", fval[n]);
    }
  }

  // Garbage collection:
  gsl_integration_workspace_free(w);
}


#endif
