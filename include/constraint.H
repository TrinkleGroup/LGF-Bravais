#ifndef __CONSTRAINT_H
#define __CONSTRAINT_H

/*
  Program: constraint.H
  Author:  D. Trinkle
  Date:    June 9, 2004
  Purpose: General algorithm to solve a constrained least squares problem,
           now with (or without) roughness criterion.

           We have Nfree free variables, Nconst constraints, and Neval
	   evaluation points.  We specify three (two) matrices and a
	   weight vector:

	   alpha[Neval, Nfree]: eval[i] = alpha[i,j]*f[j]
	   C[Nconst, Nfree]: constraint -> C[i,j]*f[j] = c[j]
	   R[Nfree, Nfree]: roughness -> r = f*R*f (optional)
	   weight[Neval]: relative weight of each evaluation point.
	   weight_R:      relative weight of roughness

	   We return two matrices: FC[Nfree, Nconst] and FE[Nfree, Neval]
	   so that the solution for constraint values c[j] and evaluation
	   values e[j] is f[i] = FC[i,j]c[j] + FE[i,j]*e[j]

  Algo.:   Testing... We're gonna try our hands at GSL's interface to
           BLAS; at some point, we'll probably have to bite the bullet
	   and really *do* BLAS and LAPACK.  Might also need to rewrite
	   the code in a more efficient / accurate manner as needed.

  Output:  Whatever I need to output.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcomp.H"
#include <gsl/gsl_blas.h>   // BLAS baby!
#include <gsl/gsl_eigen.h>  // eigensolver (replaced by LAPACK in the future)
// #include <gsl/gsl_linalg.h> // LU decomp solver

//****************************** STRUCTURES ****************************

// Return codes for constrain_fit:
const int CONST_EXTRAC = 1; // extraneous constraints
const int CONST_UNDER = 2;  // underdetermined problem


//***************************** SUBROUTINES ****************************

// C[Nconst*Nfree], alpha[Neval*Nfree], weight[Neval]
// R[Nfree*Nfree], Rweight
// FC[Nfree*Nconst], FE[Nfree*Neval] -- should already be allocated.
int constrain_fit (int Nfree, int Nconst, int Neval,
		   double* Cmat, double* alpha, double* weight,
		   double* Rmat, double Rweight,
		   double* FC, double* FE, int TESTING);

inline int constrain_fit (int Nfree, int Nconst, int Neval,
			  double* Cmat, double* alpha, double* weight,
			  double* FC, double* FE, int TESTING) 
{
  return constrain_fit(Nfree, Nconst, Neval, Cmat, alpha, weight, NULL, 0.,
		       FC, FE, TESTING);
}


int constrain_fit (int Nfree, int Nconst, int Neval,
		   double* Cmat, double* alpha, double* weight,
		   double* Rmat, double Rweight,
		   double* FC, double* FE, int TESTING) 
{
  int i, j; // General counting variables.

  // ***************************** ANALYSIS **************************
  gsl_matrix_view alphaview = gsl_matrix_view_array(alpha, Neval, Nfree);
  gsl_matrix_view Cview = gsl_matrix_view_array(Cmat, Nconst, Nfree);
  if (TESTING) {
    printf("## Nfree= %d  Nconst= %d  Neval= %d\n", Nfree, Nconst, Neval);
    printf("## alpha (%d x %d):\n", Neval, Nfree);
    for (i=0; i<Neval; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %6.2lf", alpha[i*Nfree+j]);
      printf(" |  w%d = %6.2lf\n", i, weight[i]);
    }
    printf("## Cmat (%d x %d):\n", Nconst, Nfree);
    for (i=0; i<Nconst; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %6.2lf", Cmat[i*Nfree+j]);
      printf(" |\n");
    }
  }
  // 1. Rotate to the constraint space
  // 1.a. Square Cmat
  double* Csq = new double[Nfree*Nfree];
  gsl_matrix_view Csqview = gsl_matrix_view_array(Csq, Nfree, Nfree);
  
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., &Cview.matrix, &Cview.matrix,
    		 0, &Csqview.matrix);
  if (TESTING) {
    printf("## Csq: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %6.2lf", Csq[i*Nfree+j]);
      printf(" |\n");
    }
  }
  // 1.b. Diagonalize
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(Nfree);
  double* theta = new double[Nfree*Nfree]; // rotation matrix
  double* lambda = new double[Nfree];
  gsl_matrix_view thetaview = gsl_matrix_view_array(theta, Nfree, Nfree);
  gsl_vector_view lambdaview = gsl_vector_view_array(lambda, Nfree);

  gsl_eigen_symmv(&Csqview.matrix, &lambdaview.vector,
		  &thetaview.matrix, w);
  if (TESTING) {
    printf("## theta: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", theta[i*Nfree+j]);
      printf(" |\n");
    }
    printf("## lambda: (%d)\n##", Nfree);
    for (i=0; i<Nfree; ++i) printf(" %12.5le", lambda[i]);
    printf("\n");

  }
  // 1.c. Construct theta_proj -- zero out any the vectors with non-zero
  // eigenvalues.
  double* theta_proj = new double[Nfree*Nfree];
  for (j=0; j<Nfree; ++j) {
    if (zero(lambda[j])) 
      for (i=0; i<Nfree; ++i) theta_proj[i*Nfree+j]=theta[i*Nfree+j];
    else
      for (i=0; i<Nfree; ++i) theta_proj[i*Nfree+j]=0;
  }
  if (TESTING) {
    printf("## theta_proj: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", theta_proj[i*Nfree+j]);
      printf(" |\n");
    }
  }
  
  // 2. Determine the chi^2 minimum
  // 2.a. Rotated E matrix:
  double* Emat = new double[Neval*Nfree];
  gsl_matrix_view Eview = gsl_matrix_view_array(Emat, Neval, Nfree);
  gsl_matrix_view thetapview = gsl_matrix_view_array(theta_proj, Nfree, Nfree);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.,
		 &alphaview.matrix, &thetapview.matrix,
    		 0, &Eview.matrix);
  if (TESTING) {
    printf("## Emat (rot): (%d x %d)\n", Neval, Nfree);
    for (i=0; i<Neval; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", Emat[i*Nfree+j]);
      printf(" |\n");
    }
  }
  
  // 2.b. Square Emat
  double* Esq = new double[Nfree*Nfree];
  gsl_matrix_view Esqview = gsl_matrix_view_array(Esq, Nfree, Nfree);

  // We need to scale Emat by the weights, but we only want to do it
  // ONCE, so Etemp holds the original Emat before scaling, but
  // can be erased once our construction is finished
  double* Etemp = new double[Neval*Nfree]; // how to handle weights...
  gsl_matrix_view Etempview = gsl_matrix_view_array(Etemp, Neval, Nfree);
  for (i=0; i<(Neval*Nfree); ++i) Etemp[i] = Emat[i];
  // Now, scale:
  for (i=0; i<Neval; ++i)
    for (j=0; j<Nfree; ++j)
      Emat[i*Nfree+j] *= weight[i];

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., 
		 &Etempview.matrix, &Eview.matrix,
    		 0, &Esqview.matrix);
  if (TESTING) {
    printf("## Esq: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", Esq[i*Nfree+j]);
      printf(" |\n");
    }
  }
  // GC:
  delete[] Etemp;

  // 2.b.... Add in roughness to Esq:
  gsl_matrix_view Rview = gsl_matrix_view_array(Rmat, Nfree, Nfree);
  if ( (Rmat != NULL) && (!zero(Rweight)) ) {
    double* Rtemp = new double[Nfree*Nfree];
    gsl_matrix_view Rtempview = gsl_matrix_view_array(Rtemp, Nfree, Nfree);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., 
		   &Rview.matrix, &thetapview.matrix,
		   0, &Rtempview.matrix);
    // and add to Esq: Rweight scales it, and the 1. makes Esq be
    // incremented:
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, Rweight, 
		   &thetapview.matrix, &Rtempview.matrix, 
		   1., &Esqview.matrix);
    if (TESTING) {
      printf("## Esq with roughness: (%d x %d)\n", Nfree, Nfree);
      for (i=0; i<Nfree; ++i) {
	printf("## |");
	for (j=0; j<Nfree; ++j)
	  printf(" %12.5le", Esq[i*Nfree+j]);
	printf(" |\n");
      }
    }
    // GC:
    delete[] Rtemp;
  }
  
  
  // 2.c. Diagonalize
  double* Pi = new double[Nfree*Nfree]; // rotation matrix
  double* Lambda = new double[Nfree];
  gsl_matrix_view Piview = gsl_matrix_view_array(Pi, Nfree, Nfree);
  gsl_vector_view Lambdaview = gsl_vector_view_array(Lambda, Nfree);

  gsl_eigen_symmv(&Esqview.matrix, &Lambdaview.vector,
		  &Piview.matrix, w);
  if (TESTING) {
    printf("## Pi: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", Pi[i*Nfree+j]);
      printf(" |\n");
    }
    printf("## Lambda: (%d)\n##", Nfree);
    for (i=0; i<Nfree; ++i) printf(" %12.5le", Lambda[i]);
    printf("\n");
  }

  // 3. Construct the solution
  // 3.a. make lambda_inv and Lambda_inv vectors
  double* lambda_inv = new double[Nfree];
  double* Lambda_inv = new double[Nfree];

  int n_c=0, n_f=0; // count how many degrees of freedom we have:
  for (i=0; i<Nfree; ++i) {
    if (zero(lambda[i])) {lambda_inv[i] = 0.; ++n_c;}
    else                 lambda_inv[i] = 1./lambda[i];
    if (zero(Lambda[i])) {Lambda_inv[i] = 0.; ++n_f;}
    else                 Lambda_inv[i] = 1./Lambda[i];
  }
  n_c = Nfree - n_c; // now, the number of constraints we had
  n_f = Nfree - n_f; // number of truly free variables we had
  
  // 3.b. make the constraint-solving piece:
  double* theta_inv = new double[Nfree*Nfree];
  double* C_inv = new double[Nfree*Nfree];
  gsl_matrix_view thetainv_view = gsl_matrix_view_array(theta_inv,Nfree,Nfree);
  gsl_matrix_view Cinv_view = gsl_matrix_view_array(C_inv, Nfree, Nfree);
  gsl_matrix_view FCview = gsl_matrix_view_array(FC, Nfree, Nconst);

  for (i=0; i<Nfree; ++i) 
    for (j=0; j<Nfree; ++j)
      theta_inv[i*Nfree+j] = theta[i*Nfree+j]*lambda_inv[j];
  // Cinv = theta lambdainv theta^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &thetainv_view.matrix, &thetaview.matrix,
    		 0, &Cinv_view.matrix);
  // FC = Cinv Cmat^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &Cinv_view.matrix, &Cview.matrix,
    		 0, &FCview.matrix);
  delete[] theta_inv;
  delete[] C_inv;
  if (TESTING) {
    printf("## FC--before EF: (%d x %d)\n", Nfree, Nconst);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nconst; ++j)
	printf(" %12.5le", FC[i*Nconst+j]);
      printf(" |\n");
    }
  }

  // 3.c. make the value-solving piece:
  double* Pi_inv = new double[Nfree*Nfree];
  double* E_inv = new double[Nfree*Nfree];
  double* tempmat = new double[Nfree*Nfree];
  double* Frot = new double[Nfree*Nfree];
  gsl_matrix_view Piinv_view = gsl_matrix_view_array(Pi_inv, Nfree, Nfree);
  gsl_matrix_view Einv_view = gsl_matrix_view_array(E_inv, Nfree, Nfree);
  gsl_matrix_view tempview = gsl_matrix_view_array(tempmat, Nfree, Nfree);
  gsl_matrix_view Frotview = gsl_matrix_view_array(Frot, Nfree, Nfree);
  gsl_matrix_view FEview = gsl_matrix_view_array(FE, Nfree, Neval);

  for (i=0; i<Nfree; ++i) 
    for (j=0; j<Nfree; ++j)
      Pi_inv[i*Nfree+j] = Pi[i*Nfree+j]*Lambda_inv[j];
  // Einv = Pi Lambdainv Pi^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &Piinv_view.matrix, &Piview.matrix,
    		 0, &Einv_view.matrix);
  // FE = thetap Einv Emat_rot^T
  // Frot = thetap Einv thetap^T
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.,
		 &thetapview.matrix, &Einv_view.matrix,
    		 0, &tempview.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &tempview.matrix, &Eview.matrix,
    		 0, &FEview.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &tempview.matrix, &thetapview.matrix, 
    		 0, &Frotview.matrix);
  if (TESTING) {
    printf("## FE: (%d x %d)\n", Nfree, Neval);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Neval; ++j)
	printf(" %12.5le", FE[i*Neval+j]);
      printf(" |\n");
    }
  }

  delete[] tempmat;
  delete[] Pi_inv;
  delete[] E_inv;

  // 3.d. Now, correct FC matrix using FE and roughness
  // FC = FC - (FE*alpha* Rweight*Frot*R)*FC
  double* FEalpha = new double[Nfree*Nfree];
  double* FCtemp = new double[Nfree*Nconst];
  gsl_matrix_view FEalphaview = gsl_matrix_view_array(FEalpha, Nfree, Nfree);
  gsl_matrix_view FCtempview = gsl_matrix_view_array(FCtemp, Nfree, Nconst);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.,
		 &FEview.matrix, &alphaview.matrix,
    		 0, &FEalphaview.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, Rweight,
		 &Frotview.matrix, &Rview.matrix,
    		 1, &FEalphaview.matrix);
  // add to FC (need temporary matrix to pull this off...
  for (i=0; i<(Nfree*Nconst); ++i) FCtemp[i] = FC[i];
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.,
		 &FEalphaview.matrix, &FCtempview.matrix,
		 1, &FCview.matrix);
  // GC:
  delete[] FCtemp;

  if (TESTING) {
    printf("## FC final: (%d x %d)\n", Nfree, Nconst);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nconst; ++j)
	printf(" %12.5le", FC[i*Nconst+j]);
      printf(" |\n");
    }
  }

  // Determine return value.
  int ret=0;
  if (n_c != Nconst) ret = ret | CONST_EXTRAC;
  if ( (n_c + n_f) != Nfree) ret = ret | CONST_UNDER;

  if (TESTING) {
    printf("## Nconst= %d  n_c= %d\n", Nconst, n_c);
    printf("## Nfree= %d   n_c+n_f= %d\n", Nfree, n_c+n_f);
    printf("## return= %d\n", ret);
    /*
    // THIS DOES NOT WORK AS IT IS WRITTEN.  I really needs an
    // orthogonalization after projecting out the constraint degrees
    // of freedom.
    if (ret & CONST_UNDER) {
      // try to output the unspecified directions
      double *P=new double[Nfree*Nfree];
      double *dotp = new double[Nfree*Nfree];
      // dotp[k,k'] = Pi(k).theta(k'), but only if 
      // Lambda[k]=0 and lambda[k'] != 0
      for (int k=0; k<Nfree; ++k)
	for (int kp=0; kp<Nfree; ++kp) {
	  double *td = &(dotp[k*Nfree+kp]);
	  *td = 0;
	  //	  if ( zero(Lambda[k]) && (! zero(lambda[kp])) )
	  //	    for (int l=0; l<Nfree; ++l) 
	  //	      *td += Pi[l*Nfree+k]*theta[l*Nfree+kp];
	}
      for (i=0; i<Nfree; ++i)
	for (j=0; j<Nfree; ++j) {
	  double* tP = &(P[i*Nfree+j]);
	  *tP = 0;
	  for (int k=0; k<Nfree; ++k)
	    if (zero(Lambda[k])) {
	      double vki=Pi[i*Nfree+k], vkj=Pi[j*Nfree+k];
	      for (int kp=0; kp<Nfree; ++kp)
		if ( ! zero(lambda[kp])) {
		  vki -= dotp[k*Nfree+kp]*theta[i*Nfree+kp];
		  vkj -= dotp[k*Nfree+kp]*theta[j*Nfree+kp];
		}
	      *tP += vki*vkj;
	    }
	}
      delete[] dotp;
      double* freevect = new double[Nfree*Nfree];
      double* free = new double[Nfree];
      gsl_matrix_view Pview = gsl_matrix_view_array(P, Nfree, Nfree);
      gsl_matrix_view freevectview = gsl_matrix_view_array(freevect, Nfree, Nfree);
      gsl_vector_view freeview = gsl_vector_view_array(free, Nfree);

      gsl_eigen_symmv(&Pview.matrix, &freeview.vector,
		      &freevectview.matrix, w);
      printf("## undetermined degrees of freedom:\n");
      for (i=0; i<Nfree; ++i) 
	if (! zero(free[i])) {
	  printf("## proj= %.5le |", free[i]);
	  for (j=0; j<Nfree; ++j) {
	    double val = freevect[i*Nfree+j];
	    if (! zero(val)) printf(" %+.5le f_%d", val, j);
	  }
	  printf("\n");
	}
      // GC
      delete[] P;
      delete[] freevect;
      delete[] free;
    }
    */
  }

  // ************************* GARBAGE COLLECTION ********************

  delete[] FEalpha;
  delete[] Esq;
  delete[] Pi;
  delete[] Lambda;

  delete[] Emat;
  delete[] theta_proj;
  
  delete[] theta;
  delete[] lambda;
  gsl_eigen_symmv_free(w);
  delete[] Csq;

  return ret;
}

inline void constrain_fit (int Nfree, int Nconst, int Neval,
			   double* Cmat, double* alpha, double* weight,
			   double* FC, double* FE) 
{
  constrain_fit(Nfree, Nconst, Neval, Cmat, alpha, weight, FC, FE, 0);
}

#endif
