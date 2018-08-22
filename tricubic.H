#ifndef __TRICUBIC_H
#define __TRICUBIC_H
#include "matrix.H"
#include "io.H"

/*
  Program: tricubic.H
  Author:  D. Trinkle
  Date:    June 4, 2004
  Purpose: Define and compute the expansion coeffecients for a tricubic
           spline interpolator.  In this case, we are restricted to fully
	   periodic functions in 3D, and a uniform mesh.

	   This work is primarily an extension of de Boor's ideas for
	   a bicubic spline interpolator:

	   Carl de Boor, Journal of Math. and Phys. _41_ (1962) p.212-218.

	   I'll spare you the details, but everything translates over pretty
	   easily.
*/

typedef struct 
{
  int Ni[3];
  double **C;
} tricubic_type;

const int TRI_DIM = 64; // 4^3


// We need
// 1. write out C (include a *matrix* of C's?)
// 2. read in C (include a *matrix* of C's?)
// 3. given coeff. and a indexing function, generate tricubic_type
// 4. alloc tricubic
// 5. free tricubic
// 6. given a function and an indexing function, generate tricubic_type
//    (include a set of functions?)
// 7. given a tricubic_type and a point in the unit cube [0,1]^3, return
//    the interpolated function value.

// indexing function
inline int tri_index(const int& i, const int& j, const int& k,
		     const int& N1, const int& N2) 
{ return (i*N1+j)*N2 + k;}

inline int tri_index(int n[3], int Ni[3]) 
{ return tri_index(n[0], n[1], n[2], Ni[1], Ni[2]);}

inline int tri_index(const int& i, const int& j, const int& k,
		     tricubic_type &t) 
{ return tri_index(i, j, k, t.Ni[1], t.Ni[2]);}

inline int tri_index(int n[3], tricubic_type &t) 
{ return tri_index(n[0], n[1], n[2], t.Ni[1], t.Ni[2]);}


void alloc_tricubic (const int& N0, const int& N1, const int& N2,
		     tricubic_type &t) 
{
  t.Ni[0] = N0;
  t.Ni[1] = N1;
  t.Ni[2] = N2;
  int N = N0*N1*N2;
  t.C = new double *[N];
  for (int i=0; i<N; ++i) t.C[i] = new double[TRI_DIM];
}


void free_tricubic (tricubic_type &t) 
{
  int N = t.Ni[0]*t.Ni[1]*t.Ni[2];
  for (int i=0; i<N; ++i) delete[] t.C[i];
  delete[] t.C;
  t.C = NULL;
  for (int i=0; i<3; ++i) t.Ni[i] = 0;
}


// File format is pretty simple:
// N1 N2 N3 [d] # d is optional dimension--set to 1 if not there
// C000[0][0][0][0] ... C333[0][0][0][0]
// C000[0][0][0][1] ... C333[0][0][0][1]
// ...
// C000[0][0][0][d-1] ... C333[0][0][0][d-1]
// C000[0][0][1][0] ... C333[0][0][1][0]
// ...
// C000[N1-1][N2-1][N3-1][d-1] ... C333[N1-1][N2-1][N3-1][d-1]
// Note: this requires *64* entries per line.  No small feat :)
int read_tricubic (FILE* infile, tricubic_type* &t, int &Ndim) 
{
  char dump[4096]; // note the large size needed for handling this input!
  int Ni[3], n[3];
  int d;
  
  // Get our initial line
  nextnoncomment(infile, dump, sizeof(dump));
  int pct = sscanf(dump, "%d %d %d %d", Ni, Ni+1, Ni+2, &Ndim);
  if (pct == 3) Ndim = 1;
  if ( (Ndim < 1) || (Ni[0] < 1) || (Ni[1] < 1) || (Ni[2] < 1) ) {
    fprintf(stderr, "Error reading first line of tricubic spline.  Found:\n");
    fprintf(stderr, "%s\n", dump);
    fprintf(stderr, "Expected <N1> <N2> <N3> [<dim>]\n");
    return 1;
  }
  // Allocate space
  t = new tricubic_type[Ndim];
  for (d=0; d<Ndim; ++d)
    alloc_tricubic(Ni[0], Ni[1], Ni[2], t[d]);
  
  // Begin the reading!
  int ind = 0;
  for (n[0]=0; n[0]<Ni[0]; ++(n[0]))
    for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
      for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	for (d=0; d<Ndim; ++d) {
	  double *Cp = t[d].C[ind];
	  nextnoncomment(infile, dump, sizeof(dump));
	  sscanf(dump, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		 Cp   , Cp+ 1, Cp+ 2, Cp+ 3, Cp+ 4, Cp+ 5, Cp+ 6, Cp+ 7, 
		 Cp+ 8, Cp+ 9, Cp+10, Cp+11, Cp+12, Cp+13, Cp+14, Cp+15, 
		 Cp+16, Cp+17, Cp+18, Cp+19, Cp+20, Cp+21, Cp+22, Cp+23, 
		 Cp+24, Cp+25, Cp+26, Cp+27, Cp+28, Cp+29, Cp+30, Cp+31, 
		 Cp+32, Cp+33, Cp+34, Cp+35, Cp+36, Cp+37, Cp+38, Cp+39, 
		 Cp+40, Cp+41, Cp+42, Cp+43, Cp+44, Cp+45, Cp+46, Cp+47, 
		 Cp+48, Cp+49, Cp+50, Cp+51, Cp+52, Cp+53, Cp+54, Cp+55, 
		 Cp+56, Cp+57, Cp+58, Cp+59, Cp+60, Cp+61, Cp+62, Cp+63);
	}
	++ind;
      }
  return 0;
}

void write_tricubic (FILE* outfile, tricubic_type* &t, int &Ndim) 
{
  int* Ni, n[3];
  int d;
  Ni = t->Ni; // read from the first tricubic...
  
  fprintf(outfile, "%d %d %d %d\n", Ni[0], Ni[1], Ni[2], Ndim);
  int ind = 0;
  for (n[0]=0; n[0]<Ni[0]; ++(n[0]))
    for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
      for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	for (d=0; d<Ndim; ++d) {
	  double *Cp = t[d].C[ind];
	  fprintf(outfile, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n",
		  Cp[ 0],Cp[ 1],Cp[ 2],Cp[ 3],Cp[ 4],Cp[ 5],Cp[ 6],Cp[ 7], 
		  Cp[ 8],Cp[ 9],Cp[10],Cp[11],Cp[12],Cp[13],Cp[14],Cp[15], 
		  Cp[16],Cp[17],Cp[18],Cp[19],Cp[20],Cp[21],Cp[22],Cp[23], 
		  Cp[24],Cp[25],Cp[26],Cp[27],Cp[28],Cp[29],Cp[30],Cp[31], 
		  Cp[32],Cp[33],Cp[34],Cp[35],Cp[36],Cp[37],Cp[38],Cp[39], 
		  Cp[40],Cp[41],Cp[42],Cp[43],Cp[44],Cp[45],Cp[46],Cp[47], 
		  Cp[48],Cp[49],Cp[50],Cp[51],Cp[52],Cp[53],Cp[54],Cp[55], 
		  Cp[56],Cp[57],Cp[58],Cp[59],Cp[60],Cp[61],Cp[62],Cp[63]);
	}
	++ind;
      }
}

// Copies coefficients in, using an index function:
// t[d].C[ijk][m=0..63] = coeff[index[i][j][k]][d][m]
// index[N0][N1][N2], and coeff[*][Ndim][64]
void copy_coeff(int N0, int N1, int N2, int Ndim, int*** index,
		double*** coeff, tricubic_type* &t) 
{
  int n[3];
  int d;
  int Ni[3] = {N0,N1,N2};
  
  // allocate
  t = new tricubic_type[Ndim];
  for (d=0; d<Ndim; ++d)
    alloc_tricubic(N0, N1, N2, t[d]);

  int ind = 0;
  for (n[0]=0; n[0]<Ni[0]; ++(n[0]))
    for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
      for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	double **coeff_pp = coeff[index[n[0]][n[1]][n[2]]];
	for (d=0; d<Ndim; ++d) {
	  double *coeff_p = coeff_pp[d];
	  double *Cp = t[d].C[ind];
	  for (int m=0; m<TRI_DIM; ++m) Cp[m] = coeff_p[m];
	}
	++ind;
      }
}

// pretty efficient with the cosines, so a table lookup doesn't save
// us any appreciable time, really.
void calc_HinvDel (int Ni, double* &HinvDel) 
{
  int i, j, k;
  double PI2_N = 2*M_PI / (double)Ni;
  // allocate:
  HinvDel = new double[Ni*Ni];
  
  // 1. calc HinvDel[0,k]:
  for (k=0; k<Ni; ++k) {
    HinvDel[k] = 0.;
    for (j=0; j<Ni; ++j)
      HinvDel[k] += ( cos(PI2_N*j*(k-1)) - cos(PI2_N*j*(k+1)) ) /
        ( 4 + 2*cos(PI2_N*j) );
    // scale final result: (note the CRUCIAL factor of 3 here!)
    HinvDel[k] *= 3./(double)Ni;
  }
  // 2. Permute to rest of matrix:
  for (i=1; i<Ni; ++i)
    for (k=0; k<Ni; ++k)
      HinvDel[i*Ni+k] = HinvDel[ (k-i+Ni) % Ni ]; // +Ni to make sure %Ni >= 0
}

// This is de Boor's translation matrix that takes the 1dimension derivative
// values into the expansion coefficients:
const double Atransmat[16] =
  {1, 0, 0, 0,
   0, 1, 0, 0,
  -3,-2, 3,-1,
   2, 1,-2, 1};

 

// This is the only real "analysis" piece of all our code--here
// we actually *generate* the full tricubic expansion coeffecients
// for an Ndimensional function evaluated on a grid.  index points
// to the grid values themselves:
// f(x,y,z)[d] = func[d] [ index[x][y][z] ]  <- note where d appears in func
void calc_tricubic(int N0, int N1, int N2, int Ndim, int*** index,
		   double** func, tricubic_type* &t) 
{
  int n[3], Ni[3] = {N0,N1,N2};
  int N = N0*N1*N2;
  int d, i, j, k, ind;
  
  // allocate
  t = new tricubic_type[Ndim];
  for (d=0; d<Ndim; ++d)
    alloc_tricubic(N0, N1, N2, t[d]);

  // Make HinvDel matrices for each grid direction:
  double** HinvDel = new double *[3];
  calc_HinvDel(N0, HinvDel[0]);
  if (N1 == N0) HinvDel[1] = HinvDel[0];
  else     calc_HinvDel(N2, HinvDel[1]);
  if (N2 == N0) HinvDel[2] = HinvDel[0];
  else {
    if (N2 == N1) HinvDel[2] = HinvDel[1];
    else     calc_HinvDel(N2, HinvDel[2]);
  }

  // We're actually going to treat each dimension one by one
  // allocate some temporary space:
  double * f = new double[N];
  double * fx = new double[N];
  double * fy = new double[N];
  double * fz = new double[N];
  double * fxy = new double[N];
  double * fyz = new double[N];
  double * fzx = new double[N];
  double * fxyz = new double[N];
  double K[4][4][4]; // holding place for the 64 cube corner values
  
  for (d=0; d<Ndim; ++d) {
    // 1. Copy values into f:
    for (ind=0, n[0]=0; n[0]<Ni[0]; ++(n[0]))
      for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
	for (n[2]=0; n[2]<Ni[2]; ++(n[2]), ++ind) 
	  f[ind] = func[d][index[n[0]][n[1]][n[2]]];

    // 2. first derivatives!
    // 2.a. f_x--run over y & z
    for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
      for (n[2]=0; n[2]<Ni[2]; ++(n[2])) 
	for (n[0]=0; n[0]<Ni[0]; ++(n[0])) {
	  ind=tri_index(n, Ni);
	  fx[ind]=0;
	  j = n[0]*Ni[0];
	  for (i=0; i<Ni[0]; ++i) 
	    fx[ind] += HinvDel[0][j+i]*f[tri_index(i,n[1],n[2],N1,N2)];
	}
    // 2.b. f_y--run over z & x
    for (n[2]=0; n[2]<Ni[2]; ++(n[2])) 
      for (n[0]=0; n[0]<Ni[0]; ++(n[0])) 
	for (n[1]=0; n[1]<Ni[1]; ++(n[1])) {
	  ind=tri_index(n, Ni);
	  fy[ind]=0;
	  j = n[1]*Ni[1];
	  for (i=0; i<Ni[1]; ++i) 
	    fy[ind] += HinvDel[1][j+i]*f[tri_index(n[0],i,n[2],N1,N2)];
	}
    // 2.c. f_z--run over x & y
    for (n[0]=0; n[0]<Ni[0]; ++(n[0])) 
      for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
	for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	  ind=tri_index(n, Ni);
	  fz[ind]=0;
	  j = n[2]*Ni[2];
	  for (i=0; i<Ni[2]; ++i) 
	    fz[ind] += HinvDel[2][j+i]*f[tri_index(n[0],n[1],i,N1,N2)];
	}

    // 3. second derivatives!
    // 3.a. f_xy--run over y & z
    for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
      for (n[2]=0; n[2]<Ni[2]; ++(n[2])) 
	for (n[0]=0; n[0]<Ni[0]; ++(n[0])) {
	  ind=tri_index(n, Ni);
	  fxy[ind]=0;
	  j = n[0]*Ni[0];
	  for (i=0; i<Ni[0]; ++i) 
	    fxy[ind] += HinvDel[0][j+i]*fy[tri_index(i,n[1],n[2],N1,N2)];
	}
    // 3.b. f_yz--run over z & x
    for (n[2]=0; n[2]<Ni[2]; ++(n[2])) 
      for (n[0]=0; n[0]<Ni[0]; ++(n[0])) 
	for (n[1]=0; n[1]<Ni[1]; ++(n[1])) {
	  ind=tri_index(n, Ni);
	  fyz[ind]=0;
	  j = n[1]*Ni[1];
	  for (i=0; i<Ni[1]; ++i) 
	    fyz[ind] += HinvDel[1][j+i]*fz[tri_index(n[0],i,n[2],N1,N2)];
	}
    // 3.c. f_zx--run over x & y
    for (n[0]=0; n[0]<Ni[0]; ++(n[0])) 
      for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
	for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	  ind=tri_index(n, Ni);
	  fzx[ind]=0;
	  j = n[2]*Ni[2];
	  for (i=0; i<Ni[2]; ++i) 
	    fzx[ind] += HinvDel[2][j+i]*fx[tri_index(n[0],n[1],i,N1,N2)];
	}
  
    // 4. third derivative!
    // 4.a. f_xyz--run over x & y
    for (n[0]=0; n[0]<Ni[0]; ++(n[0])) 
      for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
	for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	  ind=tri_index(n, Ni);
	  fxyz[ind]=0;
	  j = n[2]*Ni[2];
	  for (i=0; i<Ni[2]; ++i) 
	    fxyz[ind] += HinvDel[2][j+i]*fxy[tri_index(n[0],n[1],i,N1,N2)];
	}

    // 5. Convert our mixed derivatives into the expansion coefficients.
    int trans_ind[8];
    for (n[0]=0; n[0]<Ni[0]; ++(n[0])) 
      for (n[1]=0; n[1]<Ni[1]; ++(n[1]))
	for (n[2]=0; n[2]<Ni[2]; ++(n[2])) {
	  // 5.a. translation indices--point to the 8 corners of the cube.
	  trans_ind[0] = tri_index(n[0], n[1], n[2], N1, N2);
	  trans_ind[1] = tri_index(n[0], n[1], (n[2]+1)%N2, N1, N2);
	  trans_ind[2] = tri_index(n[0], (n[1]+1)%N1, n[2], N1, N2);
	  trans_ind[3] = tri_index(n[0], (n[1]+1)%N1, (n[2]+1)%N2, N1, N2);
	  trans_ind[4] = tri_index((n[0]+1)%N0, n[1], n[2], N1, N2);
	  trans_ind[5] = tri_index((n[0]+1)%N0, n[1], (n[2]+1)%N2, N1, N2);
	  trans_ind[6] = tri_index((n[0]+1)%N0, (n[1]+1)%N1, n[2], N1, N2);
	  trans_ind[7] = tri_index((n[0]+1)%N0, (n[1]+1)%N1, (n[2]+1)%N2, N1, N2);
	  
	  // 5.b. Make the "K" matrix K_abc
	  // i, i+2 -> translate by one in that direction (trans_ind)
	  // i+1 -> derivative in that direction
	  for (i=0; i<2; ++i)
	    for (j=0; j<2; ++j)
	      for (k=0; k<2; ++k) {
		K[2*i  ][2*j  ][2*k  ] = f   [trans_ind[4*i+2*j+k]];
		K[2*i  ][2*j  ][2*k+1] = fz  [trans_ind[4*i+2*j+k]];
		K[2*i  ][2*j+1][2*k  ] = fy  [trans_ind[4*i+2*j+k]];
		K[2*i  ][2*j+1][2*k+1] = fyz [trans_ind[4*i+2*j+k]];
		K[2*i+1][2*j  ][2*k  ] = fx  [trans_ind[4*i+2*j+k]];
		K[2*i+1][2*j  ][2*k+1] = fzx [trans_ind[4*i+2*j+k]];
		K[2*i+1][2*j+1][2*k  ] = fxy [trans_ind[4*i+2*j+k]];
		K[2*i+1][2*j+1][2*k+1] = fxyz[trans_ind[4*i+2*j+k]];
	      }
	  // 5.c. Matrix multiplications to get C_mno values
	  double* Cp = t[d].C[trans_ind[0]];
	  for (i=0; i<4; ++i)
	    for (j=0; j<4; ++j)
	      for (k=0; k<4; ++k) {
		ind = 16*i+4*j+k;
		Cp[ind] = 0;
		for (int a=0; a<4; ++a)
		  for (int b=0; b<4; ++b)
		    for (int c=0; c<4; ++c)
		      Cp[ind] += 
			Atransmat[4*i+a]*Atransmat[4*j+b]*Atransmat[4*k+c]
			* K[a][b][c];
	      }
	  // That's it!! (*whew*)
	}
  }

  // garbage collection:
  // be careful to delete unique ones only:
  delete[] HinvDel[0];
  if (N1 != N0) delete[] HinvDel[1];
  if ( (N2 != N0) && (N2 != N1) ) delete[] HinvDel[2];
  delete[] HinvDel;
  delete[] f;    delete[] fx;   delete[] fy;   delete[] fz;
  delete[] fxy;  delete[] fyz;  delete[] fzx;  delete[] fxyz;
}


// Want to generalize this too, to calculate a general derivative
// value, too:

// for the nth power, the mth derivative is deriv_coeff[m][n]*x^(n-m)
const double deriv_coeff[4][4] = {
  {1,1,1,1}, {0,1,2,3}, {0,0,2,6}, {0,0,0,6}};

// evaluates the function or it's mixed derivative (given by m[3])
double eval_tricubic(double x[3], int m[3], tricubic_type& t) 
{
  int n[3];
  double dx[3];
  int i;
  
  for (i=0; i<3; ++i) {
    // translate to [0,1):
    dx[i] = x[i] - (int)x[i];
    if (dx[i] < 0.) dx[i] += 1.;
    dx[i] *= t.Ni[i];
    n[i] = (int)(dx[i]);  // which point do we need?
    dx[i] = dx[i] - n[i]; // and what's the vector in it's cube?
  }
  double *C = t.C[tri_index(n,t)]; // coeff. list

  double dx_pow[3];
  double res = 0;
  int n0, n1;
  double pre0, pre1; // prefactors
  dx_pow[0] = 1.;
  for (n[0]=m[0]; n[0]<4; ++(n[0])) {
    n0 = 16*n[0];
    pre0 = deriv_coeff[m[0]][n[0]];
    dx_pow[1] = 1.;
    for (n[1]=m[1]; n[1]<4; ++(n[1])) {
      n1 = 4*n[1];
      pre1 = deriv_coeff[m[1]][n[1]];
      dx_pow[2] = 1.;
      for (n[2]=m[2]; n[2]<4; ++(n[2])) {
	res += C[n0+n1+n[2]] * dx_pow[0] * dx_pow[1] * dx_pow[2]
	  * pre0 * pre1 * deriv_coeff[m[2]][n[2]];
	dx_pow[2] *= dx[2];
      }
      dx_pow[1] *= dx[1];
    }
    dx_pow[0] *= dx[0];
  }

  // final scaling: N0^m0 * N1^n1 * N2^m2
  for (int d=0; d<3; ++d)
    for (i=0; i<m[d]; ++i)
      res *= t.Ni[d];
  
  return res;
}


// some simple "convience" wrappers:

inline double eval_tricubic(double x[3], tricubic_type& t) 
{
  int m[3] = {0,0,0};
  return eval_tricubic(x, m, t);
}

inline double eval_tricubic(double x, double y, double z, tricubic_type& t) 
{
  double xv[3] = {x,y,z};
  return eval_tricubic(xv, t);
}

inline double eval_tricubic(double x[3], int mx, int my, int mz,
			    tricubic_type& t) 
{
  int m[3] = {mx,my,mz};
  return eval_tricubic(x, m, t);
}

inline double eval_tricubic(double x, double y, double z,
			    int mx, int my, int mz,
			    tricubic_type& t) 
{
  double xv[3] = {x,y,z};
  int m[3] = {mx,my,mz};
  return eval_tricubic(xv, m, t);
}

      

#endif
