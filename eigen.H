#ifndef __EIGEN_H
#define __EIGEN_H
#include <math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_eigen.h>
#include "dcomp.H"

/*
  Program: eigen.H
  Author:  D. Trinkle
  Date:    Sept. 27, 2002
  Purpose: Define eigenvalues / vectors for 3x3 matrices.  We do this with
           gsl_poly.  See matrix.H for the simple (BLAS-like) manipulations.
*/

// even handles the sign correctly
inline double cube(double x) 
{
  if (x>0) return exp(1./3.*log(x));
  else return -exp(1./3.*log(-x));
}


#define __SWAP__(x,y) ({temp = x; x = y; y = temp;})
void sort3 (double a[3]) 
{
  double temp;
  if (a[0] > a[1]) __SWAP__(a[0], a[1]);
  if (a[1] > a[2]) __SWAP__(a[1], a[2]);
  if (a[0] > a[1]) __SWAP__(a[0], a[1]);
}
#undef __SWAP__

// EIGEN: Calculate the eigenvalues of a symmetric matrix.
//        We only do this for double matrices.
// We use the lapack-type routine for pulling this off.
// NOTE: We have to *copy* D over, because it gets modified in place
//       by the eigensolver routine.
void eigen_mat(double D[9], double* lambda) 
{
  int i;
  double Dcopy[9];
  for (i=0; i<9; ++i) Dcopy[i] = D[i];
  gsl_matrix_view m = gsl_matrix_view_array(Dcopy, 3, 3);
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(3);
  gsl_eigen_symm(&m.matrix, eval, w);
  gsl_eigen_symm_free(w);
  for (i=0; i<3; ++i) lambda[i] = gsl_vector_get(eval,i);
  gsl_vector_free(eval);
  sort3(lambda);
}


// EIGEN_VECT: Calculate the eigenvalues and vectors of a symmetric matrix.
//        We only do this for double matrices.
// We use the lapack-type routine for pulling this off.
// NOTE: We have to *copy* D over, because it gets modified in place
//       by the eigensolver routine, but then replace it with the vectors
// HORRIBLE ERROR previously; the nth column, dth component is n+3d NOT 3n+d
// Eigenvectors are the columns: D[n + 3*d] are the d components for value n
void eigen_vect(double D[9], double* lambda) 
{
  int i;
  double Dcopy[9];
  for (i=0; i<9; ++i) Dcopy[i] = D[i];
  gsl_matrix_view m = gsl_matrix_view_array(Dcopy, 3, 3);
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
  gsl_matrix_view eig = gsl_matrix_view_array(D, 3, 3);
  gsl_eigen_symmv(&m.matrix, eval, &eig.matrix, w);
  gsl_eigen_symmv_sort(eval, &eig.matrix, GSL_EIGEN_SORT_VAL_ASC);
  gsl_eigen_symmv_free(w);
  for (i=0; i<3; ++i) lambda[i] = gsl_vector_get(eval,i);
  gsl_vector_free(eval);
  //  sort3(lambda);
}


// EIGEN: Calculate the eigenvalues of a symmetric matrix.
//        We only do this for double matrices.
// Calculates the eigenvalues by solving directly the third order
// characteristic equation.  It assumes that D is symmetric.
void eigen_root(double D[9], double* lambda)
{
  double q, p, r;     // coeffecients of char. eqn
  
  // This follows the method of solution outlined in CRC 16th ed., p. 93
  // Our equation is written:
  // lambda^3 + p*lambda^2 + q*lambda + r = 0

  p = -(D[0]+D[4]+D[8]);  // -Tr D

  q = D[0]*D[4] + D[0]*D[8] + D[4]*D[8]
    - D[1]*D[1] - D[2]*D[2] - D[5]*D[5];

  r = D[0]*D[5]*D[5] + D[4]*D[2]*D[2] + D[8]*D[1]*D[1]
    - D[0]*D[4]*D[8] - 2*D[1]*D[2]*D[5];  // -det D

  gsl_poly_solve_cubic(p, q, r, lambda, lambda+1, lambda+2);
}

inline void eigen(double D[9], double* lambda) 
{
  // here we choose our prefered eigen solver:
  eigen_mat(D, lambda);
}



// Determine the corresponding eigenvector for a given eigenvalue.
// We assume:
//   (a) D is symmetric
//   (b) lambda is an eigenvalue
// The results of this procedure are not defined for
// lambda outside the eigenvalue set.
//
// Additionally, we assume that our eigenvalue lambda is unique;
// I'm not really sure what this procedure will do for a value
// of lambda which is degenerate.

void eigenvect(double D[9], double lambda, double vect[3]) 
{
  double a, c;
  double sqrta, sqrtc;
  int signb, signab;
  double D22minusl;
  double magn;
  
  D22minusl = D[8] - lambda;
  // If D22minusl = 0, we do a different solution:
  if (dcomp(D22minusl, 0.)) {
    a = D[2]*D[5];
    vect[0] = a - D[5]*D[5];
    vect[1] = a - D[2]*D[2];
    vect[2] = D[5]*(D[1] - D[0] + lambda) + D[2]*(D[1] - D[4] + lambda);
  }
  else {
    // Our normal solution
    a = D[2]*D[2] - (D[0]-lambda)*D22minusl;
    c = D[5]*D[5] - (D[4]-lambda)*D22minusl;
    // b = D[2]*D[5] - D[1]*(D[8] - lambda), so sign(b) =
    if ((D[2]*D[5]) > (D[1]*D22minusl))
      signb = 1;
    else
      signb = -1;
    // Now, we need to know sign(a) * sign(b)
    // Also, go ahead and calc. our sqrt's.
    if (dcomp(a, 0.)) {
      sqrta = 0.;
      if (c > 0) sqrtc = sqrt(c); else sqrtc = sqrt(-c);
      signab = 0;
    }
    else {
      if (a > 0.) {
	signab = signb;
	// This is for safety:
	sqrta = sqrt(a);
	if (dcomp(c, 0.)) sqrtc = 0.; else sqrtc = sqrt(c);
      }
      else {
	signab = -signb;
	// This is for safety:
	sqrta = sqrt(-a);
	if (dcomp(c, 0.)) sqrtc = 0.; else sqrtc = sqrt(-c);
      }
    }
    vect[0] = D22minusl * sqrtc;
    // Now the two sign dependent terms:
    if (signab == 1) {
      vect[1] = -sqrta*D22minusl;
      vect[2] = -(D[2]*sqrtc - D[5]*sqrta);
    }
    else {
      vect[1] = sqrta*D22minusl;
      vect[2] = -(D[2]*sqrtc + D[5]*sqrta);
    }
  }
  
  // Now, normalize it:
  magn = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
  if (dcomp(magn, 0.)) {
    vect[2] = 1.;
  }
  else {
    magn = 1./magn;
    vect[0] *= magn;
    vect[1] *= magn;
    vect[2] *= magn;
  }
}
#endif
