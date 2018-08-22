#include <stdio.h>
#include <stdlib.h>
#include "fourier.H"
// #include "cont-int.H"

int main () {
  int n[3];
  int N=4;

  for (n[0]=0; n[0]<=N; ++(n[0]))
  for (n[1]=0; n[1]<=N; ++(n[1]))
  for (n[2]=0; n[2]<=N; ++(n[2]))
    printf("%3d %3d %3d  gcd= %3d\n", n[0], n[1], n[2], gcd(n[0],n[1],n[2]));

  return 0;
}
