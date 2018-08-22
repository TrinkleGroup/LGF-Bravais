#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integrate.H"

int main (int argc, char **argv) {
  int N;
  int i, j;
  double dphi, dtheta;
  double phi, theta, r;

  if (argc <= 1) {
    fprintf(stderr, "Usage: %s <N>\n", argv[0]);
    fprintf(stderr, "  --outputs gaussian quadrature points and weights.\n");
    exit(1);
  }

  sscanf(argv[1], "%d", &N);
  if (N<=0) N = 1;

  double *x, *w;
  
  x = new double[N];
  w = new double[N];
  
  gauleg(-1.0, 1.0, x, w, N);

  // Output
  printf("%d  # Gaussian quadrature points and weight.\n", N);
  for (i=0; i<N; ++i) {
    printf("%17.15lf %22.15le\n", x[i], w[i]);
  }
  
  delete[] x;
  delete[] w;
  
  return 0;
}
