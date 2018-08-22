#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  int n, s;
  double ymax = 20*M_PI, y, j;
  const int Lmax=16;
  double xtest=60;

  for (y=(19*M_PI); y<=(20.001*M_PI); y += (0.01*M_PI))
    printf("%.5le %.8le\n", y, gsl_sf_bessel_jl(0,y));

  for (n=0; n<=Lmax; ++n)
    printf("Value of j_%d(%.5le) = %.8le\n", 2*n, xtest, gsl_sf_bessel_jl(2*n,xtest));

  for (n=0; n<=Lmax; ++n) {
    printf("Roots of j_%d(x) = 0\n", 2*n);
    s = 0;
    do {
      y = gsl_sf_bessel_zero_Jnu(2*n+0.5, s);
      j = gsl_sf_bessel_jl(2*n, y);
      printf("%3d: j_%d( %8.5le ) = %8.5le\n", s, 2*n, y, j);
      ++s;
    } while (y < ymax);
  }
  return 0;
}
