#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Outputs a triangular grid on the spherical octant */
int main (int argc, char **argv) {
  int N;
  int i, j, k;
  double x, y, z, dN;

  if (argc <= 1) {
    fprintf(stderr, "Usage: %s <N>\n", argv[0]);
    fprintf(stderr, "  --outputs a triangular grid projected onto the spherical octant.\n");
    exit(1);
  }

  sscanf(argv[1], "%d", &N);
  if (N<=0) N = 1;
  
  printf("%d  # Number of points\n", ((N+1)*(N+2))/2);
  dN = 1./(double)N;
  for (i=0; i<=N; ++i)
    for (j=0; j<=(N-i); ++j) {
      k = N - (i+j);
      /* Now, i+j+k = N, and all are between 0 and N. */
      x = sqrt ((double)i * dN);
      y = sqrt ((double)j * dN);
      z = sqrt ((double)k * dN);
      printf("%.12lf %.12lf %.12lf\n", x, y, z);
    }
  return 0;
}
