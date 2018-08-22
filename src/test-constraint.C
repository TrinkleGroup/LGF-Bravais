/*
  Program: test-constraint.C
  Author:  D. Trinkle
  Date:    June 3, 2004
  Purpose: I'm using this to test out some of the pieces of the algorithm
           for constraining a tricubic spline.  It kinda looks like a little
           nightmare... hopefully it won't be so bad when it's all said
           and done.

  Param.:  as needed...

  Flags:   MEMORY:  not used
           VERBOSE: not used
           TESTING: usual screen diahrea

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
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include <gsl/gsl_blas.h>  // BLAS baby!
#include <gsl/gsl_eigen.h> // eigensolver (replaced by LAPACK in the future)
// #include "matrix.H"
// #include "cell.H"
// #include "pointgroup.H"  // Do we still use this?
// #include "Dij.H"
// #include "kpts.H"
// #include "shell.H"
// #include "tricubic.H"

//****************************** STRUCTURES ****************************

//***************************** SUBROUTINES ****************************

double func_eval (double x) 
{
  //  return 0.5*(1-cos(2*M_PI*x));
  return -(cos(2*M_PI*x) + 0.25*cos(4*M_PI*x) + 1./36.*cos(6*M_PI*x));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<N>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'w'};  // flag characters.

const char* ARGEXPL = 
"  N  number of points along a given direction for spline\n\
  -w weight the midpoint doubly";

const char* FILEEXPL =
"";

int main ( int argc, char **argv ) 
{
  int d, i, j, k; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
                            VERBOSE, TESTING, MEMORY, 
                            NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  int WEIGHT = flagon[0];

  // ****************************** INPUT ****************************
  int Ni;
  
  sscanf(args[0], "%d", &Ni);
  if (Ni<3) {
    fprintf(stderr, "As of yet, I don't know how to do this for Ni=2 or less...\n");
    exit(1);
  }


  // ***************************** ANALYSIS **************************
  // 0. Construct our function to fit, including the second derivative
  // at the origin
  double f0;
  double fpp0; // our second derivative at 0 SCALED by dx^2 (deriv. wrt "i")
  double* func = new double[Ni];
  double* weight = new double[Ni];
  // our usual function...
  f0 = func_eval(0); // subtracted off, just in case...  
  for (i=0; i<Ni; ++i)
    func[i] = func_eval(i / (double)Ni) - f0;
  // calculated numerically...
  double dx=0.0001;
  fpp0 = (func_eval(dx)+func_eval(-dx) - 2*f0)/(dx*dx*Ni*Ni);
  for (i=0; i<Ni; ++i) weight[i] = 1;
  if (WEIGHT && ((Ni%2) == 0) )
    // only if we have an even grid...
    weight[Ni/2] = 2; // weight the "edge" point more

  // 1. Compute our alpha and alpha' matrices and our a and a' vectors
  // These describe what degrees of freedom we have, and how our variables
  // determine the values and derivatives at each lattice site from 0 to N-1.
  // We have 2N-5 degrees of freedom, as well.
  int Nfree = 2*Ni-5;
  double* alpha = new double[Ni*Nfree];
  double* alphap = new double[Ni*Nfree];
  double* a = new double[Ni];
  double* ap = new double[Ni];

  // initialize to zero
  for (i=0; i<(Nfree*Ni); ++i) {
    alpha[i]=0; alphap[i]=0;
  }
  for (i=0; i<Ni; ++i) {
    a[i]=0; ap[i]=0;
  }
  // f[0] = 0  -- nothing to do; added for compatibility
  // f'[0] = 0
  
  // f[j+2] = fbar[2j] for j=0..N-4
  // f'[j+2] = fbar[2j+1] for j=0..N-4
  for (j=0; j<=(Ni-4); ++j) {
    alpha[(j+2)*Nfree + 2*j] = 1.;    // f[j+2] = fbar[2j]
    alphap[(j+2)*Nfree + 2*j+1] = 1.; // f[j+2] = fbar[2j+1]
  }
  // We're fitting a quartic: f(x) = 1/2 f''_0 i^2 + f(4)_0 i^4
  //  f[1] = 1/2 f''0 + fbar[Nfree-1]
  //  f'[1] = f''0 + 4fbar[Nfree-1]
  //  f[N-1] = 1/2 f''0 + fbar[Nfree-1]
  //  f'[N-1] = -f''0 - 4fbar[Nfree-1]
  alpha[1*Nfree+Nfree-1] = 1.;         a[1] = 0.5*fpp0;
  alphap[1*Nfree+Nfree-1] = 4.;        ap[1] = fpp0;
  alpha[(Ni-1)*Nfree+Nfree-1] = 1.;    a[Ni-1] = 0.5*fpp0;
  alphap[(Ni-1)*Nfree+Nfree-1] = -4.;  ap[Ni-1] = -fpp0;

  if (TESTING) {
    printf("## Ni= %d  Nfree= %d\n", Ni, Nfree);
    printf("## f''0 = %.5le\n", fpp0);
    for (i=0; i<Ni; ++i)
      printf("## f0_%d = %+.5le\n", i, func[i]);
    printf("##\n");
    for (i=0; i<Ni; ++i) {
      printf("## f_%d =  %+.5le", i, a[i]);
      for (j=0; j<Nfree; ++j)
	if (! zero(alpha[i*Nfree+j])) {
	  printf(" %+.5le", alpha[i*Nfree+j]);
	  if (j == (Nfree-1)) printf(" f4_0");
	  else 
	    if (j%2) printf(" f'_%d", j/2 + 2);
	    else     printf(" f_%d", j/2 + 2);
	}
      printf("\n");
      printf("## f'_%d = %+.5le", i, ap[i]);
      for (j=0; j<Nfree; ++j)
	if (! zero(alphap[i*Nfree+j])) {
	  printf(" %+.5le", alphap[i*Nfree+j]);
	  if (j == (Nfree-1)) printf(" f4_0");
	  else 
	    if (j%2) printf(" f'_%d", j/2 + 2);
	    else     printf(" f_%d", j/2 + 2);
	}
      printf("\n");
    }
  }

  // 2. Construct the H and Del matrices
  double* Hmat = new double[Ni*Ni];
  double* Del = new double[Ni*Ni];

  for (i=0; i<(Ni*Ni); ++i) {
    Hmat[i] = 0.;
    Del[i] = 0.;
  }
  // Hmat[i,i+-1]= 1, Hmat[i,i]= 4
  // Del[i,i+1]= 3, Del[i,i-1]= -3
  for (i=0; i<Ni; ++i) {
    Hmat[i*Ni+i] = 4;
    Hmat[i*Ni+(i+1)%Ni] = 1;
    Hmat[i*Ni+(i+Ni-1)%Ni] = 1;
    Del[i*Ni+(i+1)%Ni] = 3;
    Del[i*Ni+(i+Ni-1)%Ni] = -3;
  }
  
  if (TESTING) {
    printf("## H: (%d x %d)\n", Ni, Ni);
    for (i=0; i<Ni; ++i) {
      printf("## |");
      for (j=0; j<Ni; ++j)
	printf(" %4.1lf", Hmat[i*Ni+j]);
      printf(" |\n");
    }
    printf("## Del: (%d x %d)\n", Ni, Ni);
    for (i=0; i<Ni; ++i) {
      printf("## |");
      for (j=0; j<Ni; ++j)
	printf(" %4.1lf", Del[i*Ni+j]);
      printf(" |\n");
    }
  }

  // 3. Construct our constraint matrices and vectors.
  // We do this with a completely different algorithm now.
  // The real constraint is that the appropriate second derivatives
  // are continuous across the grid points.  We define two matrices and
  // vectors:
  // f(2)+ = beta+ fbar + b+
  // f(2)- = beta- fbar + b-
  // Where f(2)+- is the second derivative at point i using the i+-1
  // direction.  Then C(m) = beta+ - beta- and c' = b+ - b-
  // We don't bother constructing beta+- and b+- explicitly; rather,
  // we just construct Cmat and cvect
  // OLD: Cmat = Hmat*alpha' - Del*alpha
  // OLD: cvect = Hmat*a' - Del*a
  double* Cmat = new double[Ni*Nfree];
  double* cvect = new double[Ni];

  // initialize
  for (i=0; i<(Ni*Nfree); ++i) Cmat[i]=0;
  for (i=0; i<Ni; ++i) cvect[i]=0;

  // b+(0)=b-(0)=fpp0, so cvect[0]=0
  // b+: (only -1 has a value)
  cvect[Ni-1] = fpp0;
  // beta+_-1:
  Cmat[(Ni-1)*Nfree+Nfree-1] += 12;
  // b-: (only +1 has a value)
  cvect[1] = -fpp0;
  // beta-_1:
  Cmat[Nfree+Nfree-1] -= 12;
  
  // beta+: only i=0,+-1 must be treated differently:
  for (i=1; i<(Ni-1); ++i) {
    // use a,alpha,a',alpha' to get the pieces
    for (j=0; j<Nfree; ++j) 
      Cmat[i*Nfree+j] += -6*alpha[i*Nfree+j] +6*alpha[(i+1)*Nfree+j]
	-4*alphap[i*Nfree+j] -2*alphap[(i+1)*Nfree+j];
    cvect[i] += -6*a[i]+6*a[i+1] -4*ap[i] -2*ap[i+1];
  }
  // beta-: only i=0,+-1 must be treated differently:
  for (i=2; i<Ni; ++i) {
    // use a,alpha,a',alpha' to get the pieces
    for (j=0; j<Nfree; ++j) 
      Cmat[i*Nfree+j] -= -6*alpha[i*Nfree+j] +6*alpha[(i-1)*Nfree+j]
	+4*alphap[i*Nfree+j] +2*alphap[(i-1)*Nfree+j];
    cvect[i] -= -6*a[i]+6*a[i-1] +4*ap[i] +2*ap[i-1];
  }

  gsl_matrix_view alphaview = gsl_matrix_view_array(alpha, Ni, Nfree);
  gsl_matrix_view alphapview = gsl_matrix_view_array(alphap, Ni, Nfree);
  gsl_vector_view aview = gsl_vector_view_array(a, Ni);
  gsl_vector_view apview = gsl_vector_view_array(ap, Ni);
  gsl_matrix_view Cview = gsl_matrix_view_array(Cmat, Ni, Nfree);
  gsl_vector_view cview = gsl_vector_view_array(cvect, Ni);

  /*
    // OLD METHOD (wrong):
  gsl_matrix_view Hview = gsl_matrix_view_array(Hmat, Ni, Ni);
  gsl_matrix_view Delview = gsl_matrix_view_array(Del, Ni, Ni);

  // C = H alpha' ...
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 &Hview.matrix, &alphapview.matrix, 0.0, &Cview.matrix);
  // - Del alpha.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, 
		 &Delview.matrix, &alphaview.matrix, 1.0, &Cview.matrix);

  // c = H a' ...
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Hview.matrix, &apview.vector,
		 0, &cview.vector);
  // - Del a.
  gsl_blas_dgemv(CblasNoTrans, -1.0, &Delview.matrix, &aview.vector,
		 1, &cview.vector);
  */

  if (TESTING) {
    printf("## Cmat (%d x %d) + cvect (%d):\n", Ni, Nfree, Nfree);
    for (i=0; i<Ni; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %6.2lf", Cmat[i*Nfree+j]);
      printf(" |   | %6.2lf |\n", cvect[i]);
    }
  }

  // 4. Rotate to the constraint space
  // 4.a. Square Cmat
  double* Csq = new double[Nfree*Nfree];
  for (i=0; i<(Nfree*Nfree); ++i) Csq[i]=0;
  gsl_matrix_view Csqview = gsl_matrix_view_array(Csq, Nfree, Nfree);
  
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., &Cview.matrix, &Cview.matrix,
    		 0, &Csqview.matrix);
  // squares, and stores in the upper half: DOESN'T WORK...?
  //  gsl_blas_dsyrk(CblasUpper, CblasTrans, 1., &Cview.matrix, 
  //		 0., &Csqview.matrix);
  // transpose down to the lower half:
  //  for (i=0; i<Nfree; ++i)
  //    for (j=0; j<i; ++j)
  //      Csq[i*Nfree+j] = Csq[j*Nfree+i];

  if (TESTING) {
    printf("## Csq: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %6.2lf", Csq[i*Nfree+j]);
      printf(" |\n");
    }
  }
  // 4.b. Diagonalize
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

  // 4.c. Construct theta_proj -- zero out any the vectors with non-zero
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
  
  // 5. Determine the chi^2 minimum
  // 5.a. Rotated E matrix:
  double* Emat = new double[Ni*Nfree];
  for (i=0; i<(Ni*Nfree); ++i) Emat[i]=0;
  gsl_matrix_view Eview = gsl_matrix_view_array(Emat, Ni, Nfree);
  gsl_matrix_view thetapview = gsl_matrix_view_array(theta_proj, Nfree, Nfree);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.,
		 &alphaview.matrix, &thetapview.matrix,
    		 0, &Eview.matrix);
  if (TESTING) {
    printf("## Emat (rot): (%d x %d)\n", Ni, Nfree);
    for (i=0; i<Ni; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", Emat[i*Nfree+j]);
      printf(" |\n");
    }
  }
  
  // 5.b. Square Emat
  double* Esq = new double[Nfree*Nfree];
  double* Etemp = new double[Ni*Nfree]; // how to handle weights...
  gsl_matrix_view Etempview = gsl_matrix_view_array(Etemp, Ni, Nfree);
  // We need to scale Emat by the weights, but we only want to do it
  // ONCE, so Etemp holds the original Emat before scaling, but
  // can be erased once our construction is finished
  for (i=0; i<(Ni*Nfree); ++i) Etemp[i] = Emat[i];
  // Now, scale:
  for (i=0; i<Ni; ++i)
    for (j=0; j<Nfree; ++j)
      Emat[i*Nfree+j] *= weight[i];
  
  for (i=0; i<(Nfree*Nfree); ++i) Esq[i]=0;
  gsl_matrix_view Esqview = gsl_matrix_view_array(Esq, Nfree, Nfree);
  
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., 
		 &Etempview.matrix, &Eview.matrix,
    		 0, &Esqview.matrix);
  
  // GC:
  delete[] Etemp;
  
  if (TESTING) {
    printf("## Esq: (%d x %d)\n", Nfree, Nfree);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Nfree; ++j)
	printf(" %12.5le", Esq[i*Nfree+j]);
      printf(" |\n");
    }
  }
  // 5.c. Diagonalize
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

  // 6. Construct the solution
  // 6.a. make lambda_inv and Lambda_inv vectors
  double* lambda_inv = new double[Nfree];
  double* Lambda_inv = new double[Nfree];
  
  for (i=0; i<Nfree; ++i) {
    if (zero(lambda[i])) lambda_inv[i] = 0.;
    else                 lambda_inv[i] = 1./lambda[i];
    if (zero(Lambda[i])) Lambda_inv[i] = 0.;
    else                 Lambda_inv[i] = 1./Lambda[i];
  }
  
  // 6.b. make the constraint-solving piece:
  double* theta_inv = new double[Nfree*Nfree];
  double* C_inv = new double[Nfree*Nfree];
  double* FCmat = new double[Nfree*Ni];
  gsl_matrix_view thetainv_view = gsl_matrix_view_array(theta_inv,Nfree,Nfree);
  gsl_matrix_view Cinv_view = gsl_matrix_view_array(C_inv, Nfree, Nfree);
  gsl_matrix_view FCview = gsl_matrix_view_array(FCmat, Nfree, Ni);

  for (i=0; i<Nfree; ++i) 
    for (j=0; j<Nfree; ++j)
      theta_inv[i*Nfree+j] = theta[i*Nfree+j]*lambda_inv[j];
  // Cinv = theta lambdainv theta^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &thetainv_view.matrix, &thetaview.matrix,
    		 0, &Cinv_view.matrix);
  // FC = - Cinv Cmat^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.,
		 &Cinv_view.matrix, &Cview.matrix,
    		 0, &FCview.matrix);
  if (TESTING) {
    printf("## FC: (%d x %d)\n", Nfree, Ni);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Ni; ++j)
	printf(" %12.5le", FCmat[i*Ni+j]);
      printf(" |\n");
    }
  }
  delete[] theta_inv;
  delete[] C_inv;


  // 6.c. make the value-solving piece:
  double* Pi_inv = new double[Nfree*Nfree];
  double* E_inv = new double[Nfree*Nfree];
  double* tempmat = new double[Nfree*Ni];
  double* EFmat = new double[Nfree*Ni];
  gsl_matrix_view Piinv_view = gsl_matrix_view_array(Pi_inv, Nfree, Nfree);
  gsl_matrix_view Einv_view = gsl_matrix_view_array(E_inv, Nfree, Nfree);
  gsl_matrix_view tempview = gsl_matrix_view_array(tempmat, Nfree, Ni);
  gsl_matrix_view EFview = gsl_matrix_view_array(EFmat, Nfree, Ni);

  for (i=0; i<Nfree; ++i) 
    for (j=0; j<Nfree; ++j)
      Pi_inv[i*Nfree+j] = Pi[i*Nfree+j]*Lambda_inv[j];
  // Einv = Pi Lambdainv Pi^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &Piinv_view.matrix, &Piview.matrix,
    		 0, &Einv_view.matrix);
  // EF = theta Einv Emat_rot^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.,
		 &Einv_view.matrix, &Eview.matrix,
    		 0, &tempview.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.,
		 &thetaview.matrix, &tempview.matrix,
    		 0, &EFview.matrix);
  if (TESTING) {
    printf("## EF: (%d x %d)\n", Nfree, Ni);
    for (i=0; i<Nfree; ++i) {
      printf("## |");
      for (j=0; j<Ni; ++j)
	printf(" %12.5le", EFmat[i*Ni+j]);
      printf(" |\n");
    }
  }

  delete[] tempmat;
  delete[] Pi_inv;
  delete[] E_inv;

  // 7. Finally, solve our problem (namely, produce the coefficients)
  double* f_a = new double[Ni]; // f^0 - a
  for (i=0; i<Ni; ++i)
    f_a[i] = func[i] - a[i];
  double* fvec = new double[Nfree]; // our final coefficients!
  for (i=0; i<Nfree; ++i) fvec[i]=0;
  // cview : cvector already allocated
  gsl_vector_view f_aview = gsl_vector_view_array(f_a, Ni);
  gsl_vector_view fvecview = gsl_vector_view_array(fvec, Nfree);

  // fvec = FC*c + EF*(f-a - alpha*fc); fc = FC*c
  gsl_blas_dgemv(CblasNoTrans, 1, &FCview.matrix, &cview.vector,
		 0., &fvecview.vector);
  for (i=0; i<Ni; ++i) 
    for (j=0; j<Nfree; ++j) f_a[i] -= alpha[i*Nfree+j]*fvec[j];
  gsl_blas_dgemv(CblasNoTrans, 1, &EFview.matrix, &f_aview.vector,
		 1., &fvecview.vector);
  
  if (TESTING) {
    printf("## fvec final:\n");
    for (i=0; i<Nfree; ++i)
      printf("## fbar_%d = %+.5le\n", i, fvec[i]);
  }

  // 8. Convert to f and fp values for each point:
  double* f = new double[Ni];
  double* fp = new double[Ni];
  for (i=0; i<Ni; ++i) {
    f[i] = a[i];
    fp[i] = ap[i];
    for (j=0; j<Nfree; ++j) {
      f[i] += alpha[i*Nfree+j]*fvec[j];
      fp[i] += alphap[i*Nfree+j]*fvec[j];
    }
  }
  
  // ****************************** OUTPUT ***************************
  // turn back into values and derivatives at each point

  if (VERBOSE) {
    printf("# f4_0 = %.8le\n", fvec[Nfree-1]);
    for (i=0; i<Ni; ++i) {
      printf("# f_%d =  %.8le\n", i, f[i]);
      printf("# f'_%d = %.8le\n", i, fp[i]);
    }
  }

  double finter; // interpolation value
  for (double x=0; x<=1.0000001; x += 0.001) {
    double u = x*Ni;
    i = (int)u;
    u = u-i;
    if (i==(Ni-1)) u = 1-u;
    if ( (i==0) || (i==(Ni-1)) )
      // quartic case:
      finter = u*u*(0.5*fpp0 + u*u*fvec[Nfree-1]);
    else {
      // cubic case:
      // de Boor's formulae:
      double f2 = -3*f[i]-2*fp[i]+3*f[i+1]-fp[i+1];
      double f3 = 2*f[i]+fp[i]-2*f[i+1]+fp[i+1];
      finter = f[i]+u*(fp[i]+u*(f2+u*f3));
    }
    printf("%.5lf %.8le %.8le\n", x, func_eval(x)-f0, finter);
  }
  


  // ************************* GARBAGE COLLECTION ********************
  delete[] f_a;
  delete[] fvec;
  
  delete[] FCmat;
  delete[] EFmat;

  delete[] Esq;
  delete[] Pi;
  delete[] Lambda;

  delete[] Emat;
  delete[] theta_proj;
  
  delete[] theta;
  delete[] lambda;
  gsl_eigen_symmv_free(w);
  delete[] Csq;

  delete[] Cmat;
  delete[] cvect;
  
  delete[] Hmat;
  delete[] Del;

  delete[] a;
  delete[] ap;
  delete[] alpha;
  delete[] alphap;

  delete[] func;
  delete[] weight;
  
  return 0;
}
