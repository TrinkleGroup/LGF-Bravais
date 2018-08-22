#ifndef __SUM_BIN_H
#define __SUM_BIN_H

#include <stdio.h>
#include <stdlib.h>

const double MAX_DOUBLE = 1e300; // some upper limit ...

typedef struct 
{
  int Nbins;
  double *maxval, *val;
} bin_type;

// Takes a minbin and maxbin size, and a logarithmic step size:
// maxval[0] = 0.
// maxval[1] = minbin
// maxval[2] = minbin * step
// ...
// maxval[Nbins-2] = minbin * step^(Nbins-3) < maxbin
// maxval[Nbins-1] = MAX_DOUBLE
bin_type* init_bin (const double &minbin, const double &maxbin, 
		    const double &step) 
{
  if ( (minbin < 0) || (maxbin < 0) || (step <= 1) || (maxbin<minbin) )
    return NULL;

  int i;
  bin_type *binp = new bin_type;
  // round it to be on the safe side...
  binp->Nbins = (int)(3.5 + (log(maxbin)-log(minbin))/log(step));
  binp->maxval = new double[binp->Nbins];
  binp->val = new double[binp->Nbins];
  binp->maxval[0] = 0;       binp->val[0] = 0.;
  binp->maxval[1] = minbin;  binp->val[1] = 0.;
  for (i=2; i<binp->Nbins; ++i) {
    binp->maxval[i] = step * binp->maxval[i-1];
    binp->val[i] = 0.;
  }
  binp->maxval[binp->Nbins-1] = MAX_DOUBLE;
  
  return binp;
}

// Some default values... (never a good idea just to use defaults...)
inline bin_type* init_bin() 
{ return init_bin(1./1099511627776., 1099511627776., 1024); }


void free_bin (bin_type* &binp) 
{
  if (binp != NULL) {
    if (binp->maxval != NULL) delete[] binp->maxval;
    if (binp->val != NULL) delete[] binp->val;
    delete binp;
    binp = NULL;
  }
}

inline void clear_bin (const bin_type* binp) 
{
  for (int i=0; i<binp->Nbins; ++i) binp->val[i] = 0.;
}

inline double sum_bin (const bin_type* binp) 
{
  double result=0;
  for (int i=0; i<binp->Nbins; ++i) result += binp->val[i];
  return result;
}

// add x to bin; we do it with end recursion.
void inc_bin (const bin_type* binp, const double &x);
void inc_bin (const bin_type* binp, const double &x) 
{
  if (x==0) return;
  int i;
  double tx = fabs(x); // temporary value...
  for (i=0; (tx > binp->maxval[i]); ++i) ;
  tx = binp->val[i] + x;
  // see if stays in this bin...; else, zero out this bin and
  // recursively call inc_bin.
  if ( (fabs(tx) > binp->maxval[i]) || (fabs(tx) < binp->maxval[i-1]) ) {
    binp->val[i] = 0.;
    inc_bin(binp, tx);
    return;
  }
  else {
    // it belongs in this bin
    binp->val[i] = tx;
    return;
  }
}


// a recursive (destructive) method to sum up a bin--it folds down
// until there is but one element left.
double recurse_sum_bin (const bin_type* binp) 
{
  int i, j;
  double sum;
  for (i=0; (i<binp->Nbins) && (binp->val[i] == 0); ++i) ;
  for (j=(i+1); (j<binp->Nbins) && (binp->val[j] == 0); ++j) ;
  if (i == binp->Nbins) return 0; // on the off chance...
  if (j == binp->Nbins) return binp->val[i]; // only one value left.
  // else, sum up the i and j values, increment, and return:
  sum = binp->val[i] + binp->val[j];
  binp->val[i] = 0;
  binp->val[j] = 0;
  inc_bin(binp, sum);
  return recurse_sum_bin(binp);
}

#endif
