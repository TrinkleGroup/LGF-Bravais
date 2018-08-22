#ifndef __LAMBDA_H
#define __LAMBDA_H

/*
  Header: lambda.H
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a dynamical matrix, and calculate:

                       1
	   [ab,cd] = - - SUM D   x x
                       2  x   ab  c d

                          1
	   [ab,cdef] = - -- SUM D   x x x x
                         24  x   ab  c d e f

                            1
	   [ab,cdefgh] = - --- SUM D   x x x x x x
                           720  x   ab  c d e f g h


	   From these, we construct the matrices:

	   lambda (k) = SUM [ab,cd]k k
	                 cd         c d

                 (2)
	   lambda  (k) =  SUM [ab,cd]k k k k 
	                 cdef         c d e f

                 (3)
	   lambda  (k) =  SUM [ab,cdgh]k k k k k k
	                 cdefgh         c d e f g h

	   where k is normalized.  We construct the matrices:


                      -1
	   [lambda(k)]

                      -1       (2)               -1
	   [lambda(k)]  [lambda  (k)] [lambda(k)]

                      -1       (2)               -1       (2)               -1
	   [lambda(k)]  [lambda  (k)] [lambda(k)]  [lambda  (k)] [lambda(k)]
                          -1       (3)               -1
  	     - [lambda(k)]  [lambda  (k)] [lambda(k)]

	   These are used to accurately compute the LGF in Fourier
	   space, especially for small k values.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "voigt.H"
#include "Dij.H"

//***************************** SUBROUTINES ****************************

// How to take k.lambda.k and k.k.lambda2.k.k:

inline void lambda_k (double kpt[3], double lambda[9][9], double res[9]) 
{
  for (int i=0; i<9; ++i)
    res[i] = innerprod(kpt, lambda[i], kpt);
}

inline void lambda2_k (double kpt[3], double lambda2[9][9][9], double res[9]) 
{
  double temp[9];
  for (int i=0; i<9; ++i) {
    for (int j=0; j<9; ++j)
      temp[j] = innerprod(kpt, lambda2[i][j], kpt);
    res[i] = innerprod(kpt, temp, kpt);
  }
}

inline void lambda3_k (double kpt[3], double lambda3[9][9][9][9], 
		       double res[9]) 
{
  double temp2[9][9], temp[9];
  for (int i=0; i<9; ++i) {
    for (int j=0; j<9; ++j) {
      for (int k=0; k<9; ++k) {
	temp2[j][k] = innerprod(kpt, lambda3[i][j][k], kpt);
      }
      temp[j] = innerprod(kpt, temp2[j], kpt);
    }
    res[i] = innerprod(kpt, temp, kpt);
  }
}

inline void lambda4_k (double kpt[3], double lambda4[9][9][9][9][9], 
		       double res[9]) 
{
  double temp3[9][9][9], temp2[9][9], temp[9];
  for (int i=0; i<9; ++i) {
    for (int j=0; j<9; ++j) {
      for (int k=0; k<9; ++k) {
	for (int l=0; l<9; ++l)
	  temp3[j][k][l] = innerprod(kpt, lambda4[i][j][k][l], kpt);
	temp2[j][k] = innerprod(kpt, temp3[j][k], kpt);
      }
      temp[j] = innerprod(kpt, temp2[j], kpt);
    }
    res[i] = innerprod(kpt, temp, kpt);
  }
}


void make_lambda(int Np, point_type* p, double lambda[9][9]) 
{
  // Use Voigt notation to generate the symmetrized elements:
  double lambda_i[9][Ni];
  int ab, cd, i;
  int c, d;
  int np;
  
  for (ab=0; ab<9; ++ab)
    for (i=0; i<Ni; ++i) lambda_i[ab][i] = 0.;
  
  for (np=0; np<Np; ++np) 
    for (cd=0; cd<Ni; ++cd) {
      c = i2rank2[cd][0];
      d = i2rank2[cd][1];
      for (ab=0; ab<9; ++ab)
	lambda_i[ab][cd] += (p[np].mat[ab])*(p[np].Rcart[c])*(p[np].Rcart[d]);
    }

  // scaling:
  for (ab=0; ab<9; ++ab)
    for (i=0; i<Ni; ++i) lambda_i[ab][i] *= -0.5;

  // Now convert to full matrices:
  for (cd=0; cd<9; ++cd)
    for (ab=0; ab<9; ++ab)
      lambda[ab][cd] = lambda_i[ab][m9toi6[cd]];
}

void make_lambda2(int Np, point_type* p, double lambda2[9][9][9]) 
{
  // Use Voigt notation to generate the symmetrized elements:
  double lambda2_ii[9][Nii];
  int ab, cd, ef, ii;
  int c, d, e, f;
  int np;
  
  for (ab=0; ab<9; ++ab) 
    for (ii=0; ii<Nii; ++ii) lambda2_ii[ab][ii] = 0.;
  
  for (np=0; np<Np; ++np) 
    for (ii=0; ii<Nii; ++ii) {
      c = ii2rank4[ii][0];
      d = ii2rank4[ii][1];
      e = ii2rank4[ii][2];
      f = ii2rank4[ii][3];
      for (ab=0; ab<9; ++ab)
	lambda2_ii[ab][ii] += (p[np].mat[ab])*(p[np].Rcart[c])*(p[np].Rcart[d])
	  *(p[np].Rcart[e])*(p[np].Rcart[f]);
    }

  // scaling:
  for (ab=0; ab<9; ++ab) 
    for (ii=0; ii<Nii; ++ii) lambda2_ii[ab][ii] *= -1./24.;
  
  // Now convert to full matrices:
  for (cd=0; cd<9; ++cd) 
    for (ef=0; ef<9; ++ef)
      for (ab=0; ab<9; ++ab)
	lambda2[ab][cd][ef] = lambda2_ii[ab][mmtoii[cd][ef]];
}



void make_lambda3(int Np, point_type* p, double lambda3[9][9][9][9]) 
{
  // We don't use Voigt notation here... it's too cumbersome (though
  // it would save us a lot of work...)
  int ab, cd, ef, gh;
  int a, b, c, d, e, f, g, h;
  int np;
  
  for (ab=0; ab<9; ++ab) 
    for (cd=0; cd<9; ++cd) 
      for (ef=0; ef<9; ++ef) 
	for (gh=0; gh<9; ++gh) 
	  lambda3[ab][cd][ef][gh] = 0.;
  
  for (np=0; np<Np; ++np) 
    for (ab=0; ab<9; ++ab) {
      a = ab/3; b = ab%3;
      for (cd=0; cd<9; ++cd) {
	c = cd/3; d = cd%3;
	for (ef=0; ef<9; ++ef) {
	  e = ef/3; f = ef%3;
	  for (gh=0; gh<9; ++gh) {
	    g = gh/3; h = gh%3;
	    lambda3[ab][cd][ef][gh] += (p[np].mat[ab])
	      *(p[np].Rcart[c])*(p[np].Rcart[d])
	      *(p[np].Rcart[e])*(p[np].Rcart[f])
	      *(p[np].Rcart[g])*(p[np].Rcart[h]);
	  }
	}
      }
    }

  for (ab=0; ab<9; ++ab) 
    for (cd=0; cd<9; ++cd) 
      for (ef=0; ef<9; ++ef) 
	for (gh=0; gh<9; ++gh) 
	  lambda3[ab][cd][ef][gh] *= -1./720.;
}


void make_lambda4(int Np, point_type* p, double lambda4[9][9][9][9][9]) 
{
  // We don't use Voigt notation here... it's too cumbersome (though
  // it would save us a lot of work...)
  int ab, cd, ef, gh, ij;
  int a, b, c, d, e, f, g, h, i, j;
  int np;
  
  for (ab=0; ab<9; ++ab) 
    for (cd=0; cd<9; ++cd) 
      for (ef=0; ef<9; ++ef) 
	for (gh=0; gh<9; ++gh) 
	  for (ij=0; ij<9; ++ij) 
	  lambda4[ab][cd][ef][gh][ij] = 0.;
  
  for (np=0; np<Np; ++np) 
    for (ab=0; ab<9; ++ab) {
      a = ab/3; b = ab%3;
      for (cd=0; cd<9; ++cd) {
	c = cd/3; d = cd%3;
	for (ef=0; ef<9; ++ef) {
	  e = ef/3; f = ef%3;
 	  for (gh=0; gh<9; ++gh) {
	    g = gh/3; h = gh%3;
	    for (ij=0; ij<9; ++ij) {
	      i = ij/3; j = ij%3;
	      lambda4[ab][cd][ef][gh][ij] += (p[np].mat[ab])
		*(p[np].Rcart[c])*(p[np].Rcart[d])
		*(p[np].Rcart[e])*(p[np].Rcart[f])
		*(p[np].Rcart[g])*(p[np].Rcart[h])
		*(p[np].Rcart[i])*(p[np].Rcart[j]);
	    }
	  }
	}
      }
    }

  for (ab=0; ab<9; ++ab) 
    for (cd=0; cd<9; ++cd) 
      for (ef=0; ef<9; ++ef) 
	for (gh=0; gh<9; ++gh) 
	  for (ij=0; ij<9; ++ij) 
	    lambda4[ab][cd][ef][gh][ij] *= -1./40320.;
}


inline void eval_lambda_inv(double lambda[9][9], double kpt[3],
			    double lk_inv[9]) 
{
  double lk[9];

  lambda_k(kpt, lambda, lk);
  careful_inverse(lk, lk_inv);
}

inline void eval_lambda2(double lambda[9][9], double lambda2[9][9][9],
			 double kpt[3], double res[9]) 
{
  double magn;
  magn = 1./sqrt(dot(kpt,kpt));
  double tkpt[3];
  for (int i=0; i<3; ++i) tkpt[i] = kpt[i] * magn;
  
  double lk_inv[9], temp[9];
  eval_lambda_inv(lambda, tkpt, lk_inv);

  lambda2_k(tkpt, lambda2, res);

  // multiply
  mult(lk_inv, res, temp);
  mult(temp, lk_inv, res);
}


// We normalize our kpt, then multiply by |kpt|^2 at the end:
inline void eval_lambda3(double lambda[9][9], double lambda2[9][9][9],
			 double lambda3[9][9][9][9], 
			 double kpt[3], double res[9]) 
{
  double magn;
  magn = 1./sqrt(dot(kpt,kpt));
  double tkpt[3];
  for (int i=0; i<3; ++i) tkpt[i] = kpt[i] * magn;
  
  double lk_inv[9], lk2[9], lk3[9], temp[9];
  eval_lambda_inv(lambda, tkpt, lk_inv);

  lambda2_k(tkpt, lambda2, lk2);
  lambda3_k(tkpt, lambda3, res);

  // multiply and sum:
  mult(lk_inv, lk2, temp); // temp = (L1)^-1 L2
  mult(temp, temp, lk2);   // lk2 =  (L1)^-1 L2 (L1)^-1 L2
  mult(lk2, lk_inv, temp); // temp = (L1)^-1 L2 (L1)^-1 L2 (L1)^-1

  mult(lk_inv, res, lk3);  // lk3 =  (L1)^-1 L3
  mult(lk3, lk_inv, res);  // res =  (L1)^-1 L3 (L1)^-1

  for (int i=0; i<9; ++i) res[i] = temp[i] - res[i];
  // scale:
  mult(res, dot(kpt,kpt), res);
}



// We normalize our kpt, then multiply by |kpt|^4 at the end:
inline void eval_lambda4(double lambda[9][9], double lambda2[9][9][9],
			 double lambda3[9][9][9][9], 
			 double lambda4[9][9][9][9][9],
			 double kpt[3], double res[9]) 
{
  double magn;
  magn = 1./sqrt(dot(kpt,kpt));
  double tkpt[3];
  for (int i=0; i<3; ++i) tkpt[i] = kpt[i] * magn;
  
  double lk_inv[9], lk2[9], lk3[9], lk4[9], temp1[9], temp2[9];
  eval_lambda_inv(lambda, tkpt, lk_inv);

  lambda2_k(tkpt, lambda2, lk2);
  lambda3_k(tkpt, lambda3, lk3);
  lambda4_k(tkpt, lambda4, lk4);

  double lk1_2[9], lk1_3[9], lk1_4[9];
  mult(lk_inv, lk2, lk1_2); // lk1_2 = (L1)^-1 L2
  mult(lk_inv, lk3, lk1_3); // lk1_3 = (L1)^-1 L3
  mult(lk_inv, lk4, lk1_4); // lk1_4 = (L1)^-1 L4
  
  mult(lk1_4, lk_inv, res); // res = (L1)^-1 L4 (L1)^-1

  mult(lk1_2, lk1_3, temp1); // temp1 = (L1)^-1 (L2) (L1)^-1 (L3) 
  mult(lk1_3, lk1_2, temp2); // temp2 = (L1)^-1 (L3) (L1)^-1 (L2) 
  double lk1_23_1[9], lk1_32_1[9];
  mult(temp1, lk_inv, lk1_23_1); // lk1_23_1= (L1)^-1 (L2) (L1)^-1 (L3) (L1)^-1
  mult(temp2, lk_inv, lk1_32_1); // lk1_32_1= (L1)^-1 (L3) (L1)^-1 (L2) (L1)^-1
  
  mult(lk1_2, lk1_2, temp1); // temp1=  (L1)^-1 (L2) (L1)^-1 (L2) 
  mult(lk1_2, temp1, temp2); // temp2=  (L1)^-1 (L2) (L1)^-1 (L2)  (L1)^-1 (L2)
  mult(temp2, lk_inv, temp1); // temp1 = 1 2 1 2 1 2 1

  // sum:
  for (int i=0; i<9; ++i) res[i] += lk1_23_1[i]+lk1_32_1[i]+temp1[i];
  // scale:
  mult(res, dot(kpt,kpt)*dot(kpt,kpt), res);
}




void print_lambda(double lambda[9][9], const char *pre) 
{
  for (int a=0; a<3; ++a)
    for (int b=0; b<3; ++b) {
      int ab = index(a,b);
      for (int c=0; c<3; ++c)
	for (int d=0; d<3; ++d) {
	  int cd = index(c,d);
	  if (! zero(lambda[ab][cd]) )
	    printf("%s[%c%c,%c%c] = %.12le\n", pre,
		   coord[a], coord[b], coord[c], coord[d], lambda[ab][cd]);
	}
    }
}

void print_lambda2(double lambda2[9][9][9], const char *pre) 
{
  for (int a=0; a<3; ++a)
    for (int b=0; b<3; ++b) {
      int ab = index(a,b);
      for (int c=0; c<3; ++c)
	for (int d=0; d<3; ++d) {
	  int cd = index(c,d);
	  for (int e=0; e<3; ++e)
	    for (int f=0; f<3; ++f) {
	      int ef = index(e,f);
	      if (! zero(lambda2[ab][cd][ef]) )
		printf("%s[%c%c,%c%c%c%c] = %.12le\n", pre,
		       coord[a],coord[b],coord[c],coord[d],coord[e],coord[f],
		       lambda2[ab][cd][ef]);
	    }
	}
    }
}

void print_lambda3(double lambda3[9][9][9][9], const char *pre) 
{
  for (int a=0; a<3; ++a)
    for (int b=0; b<3; ++b) {
      int ab = index(a,b);
      for (int c=0; c<3; ++c)
	for (int d=0; d<3; ++d) {
	  int cd = index(c,d);
	  for (int e=0; e<3; ++e)
	    for (int f=0; f<3; ++f) {
	      int ef = index(e,f);
	      for (int g=0; g<3; ++g)
		for (int h=0; h<3; ++h) {
		  int gh = index(g,h);
		  if (! zero(lambda3[ab][cd][ef][gh]) )
		    printf("%s[%c%c,%c%c%c%c%c%c] = %.12le\n", pre,
			   coord[a],coord[b],
			   coord[c],coord[d],
			   coord[e],coord[f],
			   coord[g],coord[h],
			   lambda3[ab][cd][ef][gh]);
		}
	    }
	}
    }
}

void print_lambda4(double lambda4[9][9][9][9][9], const char *pre) 
{
  for (int a=0; a<3; ++a)
    for (int b=0; b<3; ++b) {
      int ab = index(a,b);
      for (int c=0; c<3; ++c)
	for (int d=0; d<3; ++d) {
	  int cd = index(c,d);
	  for (int e=0; e<3; ++e)
	    for (int f=0; f<3; ++f) {
	      int ef = index(e,f);
	      for (int g=0; g<3; ++g)
		for (int h=0; h<3; ++h) {
		  int gh = index(g,h);
		  for (int i=0; i<3; ++i)
		    for (int j=0; j<3; ++j) {
		      int ij = index(i,j);
		      double val = lambda4[ab][cd][ef][gh][ij];
		      if (! zero(val) )
			printf("%s[%c%c,%c%c%c%c%c%c%c%c] = %.12le\n", pre,
			       coord[a],coord[b],
			       coord[c],coord[d],
			       coord[e],coord[f],
			       coord[g],coord[h],
			       coord[i],coord[j],
			       val);
		    }
		}
	    }
	}
    }
}



#endif
