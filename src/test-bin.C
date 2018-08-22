#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sum-bin.H"

const double __SRAND_SCALE__ = (2./(double)( (1<<31) - 1));
inline double srand() 
{
  return __SRAND_SCALE__ * (double)(random()) - 1.;
}


int main (void)
{
  // setup our bin:
  //  bin_type *binp = init_bin(1.e-6, 1.e6, 10.);
  //  bin_type *binp = init_bin(1./1099511627776., 1099511627776., 16.);
  bin_type *binp = init_bin(1./1099511627776., 1., 2.);
  if (binp == NULL) exit(1);

  const int seed = 13;
  srandom(seed);
  
  const int MAXi = 10000000;
  int i;
  
  double sum = 0, x;
  for (i=0; i<MAXi; ++i) {
    //    x = cos(i) * srand() * exp(srand());
    x = cos(M_PI*0.099*i);
    sum += x;
    inc_bin(binp, x);
  }
  
  printf("bin contents:\n");
  for (i=0; i<binp->Nbins; ++i) {
    printf("%10.5le : %22.15le\n", binp->maxval[i], binp->val[i]);
  }

  double err, rerr, binsum, rbinsum;
  binsum = sum_bin(binp);
  rbinsum = recurse_sum_bin(binp);
  err = (sum - binsum)/binsum;
  rerr = (rbinsum - binsum)/rbinsum;
  printf("sum =  %.15le\n", sum);
  printf("bin =  %.15le\n", binsum);
  printf("rbin = %.15le\n", rbinsum);
  printf("err =        %.5le\n", err);
  printf("rerr =       %.5le\n", rerr);
  printf("1e-15*MAXi = %.5le\n", 1e-15*MAXi);

  
  free_bin(binp);

  return 0;
}
