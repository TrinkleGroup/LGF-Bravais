// Little program to attempt to calculate the limit of series in the
// 2d EGF evaluation

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>   // Bessel functions

int main () 
{
  const int Nmax = 1<<30; // should be a multiple of 2...
  int n, N0;
  double x, j1, sum0, sum1;
  
  N0 = 1;
  sum0 = 0;
  sum1 = 0;
  for (n=1; n<=Nmax; ++n) {
    x = gsl_sf_bessel_zero_J0(n);
    j1 = gsl_sf_bessel_J1(x);
    sum1 += 2/(x*x*j1*j1);
    if (n == N0) {
      sum0 += sum1;
      sum1 = 0;
      N0 *= 2;
      printf("%d  %.15lf  %.15lf\n", n, log(x) - sum0,
	     log(gsl_sf_bessel_zero_J0(n+1)) - sum0);
    }
  }

  return 0;
}

    
