#ifndef __FOURIER_H
#define __FOURIER_H

// Code to do fourier transforms, including one that does "closed shell"
// FT for efficiency using voigt matrices -- actually, it's NOT efficient,
// sadly.  The only way to make these FT's faster is to use table-lookup.
// And we gots that going ON, boyee.

#include "io.H"
#include "matrix.H"
#include "dcomp.H"
#include "Dij.H"
#include "pointgroup.H"
#include "shell.H"
#include "voigt.H"
#include "sum-bin.H" // to allow FT with binning... not effecient...

// computes the FT of the .mat part of p
// TYPE tells in what form kpt is; either cartesian coord. or lattice.
const int KPT_CART = 0;
const int KPT_LATT = 1;

// Two constants only defined in this file (undef-ed at the end)
#define M_2PI 6.28318530717958647692
#define M_1_2PI 0.15915494309189533576

void fourier (int Np, point_type* p, double kpt[3], double pk[9], 
	      const int TYPE) 
{
  int n, d;
  double coskR;
  int LATT = TYPE & KPT_LATT;

  for (d=0; d<9; ++d) pk[d] = 0.;
  
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    if (LATT) coskR = cos(M_2PI*dot(tp->Runit, kpt));
    else      coskR = cos(dot(tp->Rcart, kpt));
    for (d=0; d<9; ++d) pk[d] += coskR*(tp->mat[d]);
  }
}

void inv_fourier (point_type* p, shell_type *tsh, double DR[9],
		  int Nop, int inv_index[MAXop], 
		  double gcart[MAXop][9],
		  double Nkpt, double** Dk, 
		  double** kpt, double* w, const int TYPE) 
{
  int n, d, k;
  double kgR, coskgR;
  int LATT = TYPE & KPT_LATT;

  // 1. zero out:
  for (d=0; d<9; ++d) DR[d] = 0;
  // 2. add it all up over k in IRZ and g[n] in our group
  for (k=0; k<Nkpt; ++k) {
    double *tkpt=kpt[k];
    double sumD[9], Drot[9], temp[9];
    for (d=0; d<9; ++d) sumD[d] = 0;
    for (n=0; n<Nop; ++n) {
      if (LATT) kgR = M_2PI*dot(p[tsh->gR0[n]].Runit, tkpt);
      else      kgR = dot(p[tsh->gR0[n]].Rcart, tkpt);
      coskgR = cos(kgR);
      mult(gcart[inv_index[n]], Dk[k], temp);
      mult(temp, gcart[n], Drot);
      for (d=0; d<9; ++d) sumD[d] += coskgR*Drot[d];
    }
    for (d=0; d<9; ++d) DR[d] += sumD[d]*w[k];
  }
  // now, scale out by 1/Nop:
  for (d=0; d<9; ++d) DR[d] *= 1./Nop;
}

// A much simpler inv_fourier--doesn't scale by weights at all
void inv_fourier (point_type &p, double Nkpt, double** pk, 
		  double** kpt, const int TYPE) 
{
  int d, k;
  double kR, coskR;
  int LATT = TYPE & KPT_LATT;

  // 1. zero out:
  for (d=0; d<9; ++d) p.mat[d] = 0;
  // 2. add it all up over k
  for (k=0; k<Nkpt; ++k) {
    double *tkpt=kpt[k];
    if (LATT) kR = M_2PI*dot(p.Runit, tkpt);
    else      kR = dot(p.Rcart, tkpt);
    coskR = cos(kR);
    for (d=0; d<9; ++d) p.mat[d] += coskR*pk[k][d];
  }
  mult(p.mat, 1./(double)Nkpt, p.mat); // divide it out.
}



// So-called "closed-shell" fourier transform.  We use the symmetry
// properties of lattice functions to do this more efficiently, plus
// we do table-lookup on cos(k.R) in each shell
inline double modPi (double x) 
{
  double y = fabs(x) - (int)(fabs(x*M_1_2PI)) * M_2PI;
  if (y > M_PI) return M_2PI - y;
  else          return y;
}

// a little easier since only one entry shows up in the table
// we return either where the value is, or where it *should* be.
int binary_find(double x, double list[MAXop+2], int nt) 
{
  int lo, hi, probe;
  lo = -1;
  hi = nt;
  // get the floor first
  while ( (hi-lo) > 1) {
    probe = (hi+lo)/2;
    if (list[probe] > (x-TOLER)) hi = probe;
    else                         lo = probe;
  }
  return hi;
}

// The actual code to do the closed-shell fourier transform.
// This assumes that our lattice functions:
// 1. is defined on closed shells (though it *technically* doesn't need that)
// 2. is composed of symmetric matrices
// 3. transforms as a second-rank tensor lattice function with the full
//    crystal point group; i.e.:
//      D(g*R) = g*D(R)*g^-1
//    where D is a 3x3 matrix, R is a lattice point, and g is a point
//    group operation on the crystal.
void fourier_cs (int Np, point_type* p, int Nsh, shell_type* sh,
		 int Nop, int inv_index[MAXop], double g_i_ginv[MAXop][Ni][Ni],
		 double kpt[3], double pk[9], const int TYPE) 
{
  int i, j, d;
  double kR[MAXop+2], coskR[MAXop+2];
  int nt, n;
  double coskRA[Ni][Ni];
  double kgR, coskgR, gig;
  int LATT = TYPE & KPT_LATT;
  double mat6[Ni], mult_mat6[Ni];
  
  for (d=0; d<9; ++d) pk[d] = 0.;
  
  for (int ns=0; ns<Nsh; ++ns) {
    shell_type *tsh = sh + ns;
    for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i) coskRA[j][i] = 0;
    for (j=0; j<Ni; ++j) mat6[j] = p[tsh->elem[0]].mat[i6tom9[j]];
    // init table lookup:
    nt=2;
    kR[0] = 0.;   coskR[0] =  1.;
    kR[1] = M_PI; coskR[1] = -1.;
    for (int npi=0; npi<tsh->Nelem; ++npi) {
      int np = tsh->elem[npi]; // npi is just the index; np is the point number
      point_type *tp = p+np;
      if (LATT) kgR = modPi(2.*M_PI*dot(tp->Runit, kpt));
      else      kgR = modPi(dot(tp->Rcart, kpt));
      i = binary_find(kgR, kR, nt);
      if (dcomp(kgR, kR[i])) coskgR = coskR[i];
      else {
	// compute it and insert:
	coskgR = cos(kgR);
	for (j=nt; j>=i; --j) {
	  kR[j+1] = kR[j]; coskR[j+1] = coskR[j];
	}
	kR[i] = kgR; coskR[i] = coskgR;
	++nt;
      }
      // Now to find which group operation we have:
      for (n=0; tsh->gR0[n] != np; ++n) ;
      // Now, actually add to our running total:
      for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i)
	// note carefully how we call g_i_ginv:
	if (! zero(gig=g_i_ginv[n][i][j]) )
	  coskRA[j][i] += coskgR*gig;
    }
    for (j=0; j<Ni; ++j) {
      mult_mat6[j] = 0.;
      for (i=0; i<Ni; ++i) mult_mat6[j] += coskRA[j][i]*mat6[i];
    }
    // add to our running FT:
    for (d=0; d<9; ++d) pk[d] += mult_mat6[m9toi6[d]];
  }
}


// The actual code to do the closed-shell INVERSE fourier transform.
// This assumes that our lattice functions:
// 1. is defined on closed shells--only used to get the point positions
// 2. is composed of symmetric matrices
// 3. transforms as a second-rank tensor lattice function with the full
//    crystal point group; i.e.:
//      D(g*k) = g*D(k)*g^-1
//    where D is a 3x3 matrix, k is a rlv, and g is a point group operation
//    on the crystal.
// In this case, we only use a single point in a given shell, but we
// pass the *entire* list of points, but just the shell we're interested in.
// The result for the first element of tsh goes into DR; we don't bother
// doing the symmetric rotation, or doing it in place for obvious reasons.

void inv_fourier_cs (point_type* p, shell_type* tsh, double DR[9],
		     int Nop, int inv_index[MAXop], 
		     double g_i_ginv[MAXop][Ni][Ni],
		     double Nkpt, double** Dk, 
		     double** kpt, double* w, const int TYPE) 
{
  int i, j, d;
  double kR[MAXop+2], coskR[MAXop+2];
  int nt;
  double coskRA[Ni][Ni];
  double kgR, coskgR, gig;
  int LATT = TYPE & KPT_LATT;
  double mat6[Ni], mult_mat6[Ni];
  double Nop_1 = 1./(double)Nop;

  // We could take this line out, since we're careful about how we
  // add everything; there's no scaling involved.
  for (d=0; d<9; ++d) DR[d] = 0.;
  
  for (int k=0; k<Nkpt; ++k) {
    double *tkpt = kpt[k];
    for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i) coskRA[j][i] = 0;
    for (j=0; j<Ni; ++j) mat6[j] = Dk[k][i6tom9[j]];
    // init table lookup:
    nt=2;
    kR[0] = 0.;   coskR[0] =  1.;
    kR[1] = M_PI; coskR[1] = -1.;
    for (int n=0; n<Nop; ++n) {
      point_type *tp = p+tsh->gR0[inv_index[n]]; // this is g^-1 R
      if (LATT) kgR = modPi(2.*M_PI*dot(tp->Runit, tkpt));
      else      kgR = modPi(dot(tp->Rcart, tkpt));
      i = binary_find(kgR, kR, nt);
      if (dcomp(kgR, kR[i])) coskgR = coskR[i];
      else {
	// compute it and insert:
	coskgR = cos(kgR);
	for (j=nt; j>=i; --j) {
	  kR[j+1] = kR[j]; coskR[j+1] = coskR[j];
	}
	kR[i] = kgR; coskR[i] = coskgR;
	++nt;
      }
      // Now, actually add to our running total:
      for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i)
	// note carefully how we call g_i_ginv:
	if (! zero(gig=g_i_ginv[n][i][j]) )
	  coskRA[j][i] += coskgR*gig;
    }
    for (j=0; j<Ni; ++j) {
      mult_mat6[j] = 0.;
      for (i=0; i<Ni; ++i) mult_mat6[j] += Nop_1*coskRA[j][i]*mat6[i];
    }
    // add to our running FT:
    for (d=0; d<9; ++d) DR[d] += w[k]*mult_mat6[m9toi6[d]];
  }
}


// **** Some introductory number theory algorithms:
// iterative gcd, using Euclid's algorithm:
int gcd (int a, int b) 
{ int r, x=abs(a), y=abs(b);
 while (y!=0) {r = (x%y);  x=y;  y=r;}
 return x;}

// gcd(a,b,c) is the gcd of gcd(a,b) and c...
int gcd (int a, int b, int c) 
{ return gcd(gcd(a,b),c);}

int lcm (int a, int b) 
{ int x=abs(a), y=abs(b);
 if (x==0) return y; if (y==0) return x;
 // else, we use Euclid's algorithm to calc the gcd:
 return (x*y)/gcd(x,y);}

// lcm(a,b,c) is the lcm of lcm(a,b) and c
int lcm (int a, int b, int c) 
{ int x=abs(a), y=abs(b), z=abs(c);
  if (x==0) return lcm(y,z);  if (y==0) return lcm(x,z);
  if (z==0) return lcm(x,y);
  return lcm(lcm(x,y),z);
}


// Code to do FT by table-lookup; really, we just make the code to
// setup the table, access it, and free it.
typedef struct 
{
  int Nt;
  double Nt_2Pi;
  double* coskR;
} coskR_table_type;


// Right now, we only handle the case of a uniform mesh in all directions;
// everything else is, frankly, too complicated (but this may change
// at a later date).
// Recently changes so that we no longer throw errors in cases where
// we have N[i] == 0...

// Given a set Ngrid value:
int init_coskR_table (int Ngrid, coskR_table_type &table) 
{
  table.Nt = Ngrid;
  table.Nt_2Pi = Ngrid * M_1_2PI;
  table.coskR = new double[Ngrid];
  // fill the entries:
  for (int i=0; i<table.Nt; ++i)
    table.coskR [i]= cos(i/table.Nt_2Pi);
  return 0;
}

// tries to determine Ngrid...
int init_coskR_table (int Nkpt, double** kpt, const int TYPE, 
		      double cart[9], coskR_table_type &table) 
{
  int i;
  int LATT = TYPE & KPT_LATT;
  int N[3] = {0,0,0}, betterN;
  double ka[3];
  double cart_trans[9];
  
  if (!LATT) transpose(cart, cart_trans);
  
  for (int k=0; k<Nkpt; ++k) {
    double *tkpt = kpt[k];
    if (LATT)
      for (i=0; i<3; ++i) ka[i] = fabs(tkpt[i]);
    else
      for (i=0; i<3; ++i) ka[i] = fabs(dot(tkpt, cart_trans+3*i)*M_1_2PI);
    for (i=0; i<3; ++i)
      if (! zero(ka[i]))
	if ( (betterN = lround(1./ka[i])) > N[i])
	  N[i] = betterN;
  }
  int Ngrid;
  Ngrid = lcm(N[0], N[1], N[2]);
  if (Ngrid==0) Ngrid = 1;
  
  // sanity checks:
  //  if ( (N[0] == 0) || (N[1] == 0) || (N[2] == 0) )
  //    return -1;
  //  if ( (N[0] != N[1]) || (N[1] != N[2]) || (N[2] != N[0]) )
  //    return -1;
  // Now, make some new table space:
  return init_coskR_table(Ngrid, table);
}

void free_coskR_table(coskR_table_type &table) 
{
  table.Nt = 0;
  delete[] table.coskR;
}

// does NO sanity checks to make sure we're bang on...
inline double mod2Pi (const double &x) 
{
  return fabs(x) - (int)(fabs(x*M_1_2PI)) * M_2PI;
}

inline int coskR_table_index(const double &x, const double &Nt_2Pi, 
			     const int &Nt) 
{
  if (x>=0) { return (int)(x*Nt_2Pi + 0.5) % Nt; }
  else      { return (int)(-x*Nt_2Pi + 0.5) % Nt; }
}

inline double coskR_table(const double &x, const coskR_table_type &table) 
{
  //  return table.coskR[lround(mod2Pi(x)*table.Nt_2Pi)];
  return table.coskR[coskR_table_index(x, table.Nt_2Pi, table.Nt)];
}


void fourier_table (int Np, point_type* p, double kpt[3], double pk[9], 
		    const int TYPE, const coskR_table_type &table) 
{
  int n, d;
  double coskR;
  int LATT = TYPE & KPT_LATT;

  for (d=0; d<9; ++d) pk[d] = 0.;
  
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    if (LATT) coskR = M_2PI*dot(tp->Runit, kpt);
    else      coskR = dot(tp->Rcart, kpt);
    coskR = coskR_table(coskR, table); // now we use table lookup...
    for (d=0; d<9; ++d) pk[d] += coskR*(tp->mat[d]);
  }
}

void inv_fourier_table (point_type* p, shell_type *tsh, double DR[9],
			int Nop, int inv_index[MAXop], 
			double gcart[MAXop][9],
			double Nkpt, double** Dk, 
			double** kpt, double* w, const int TYPE,
			const coskR_table_type &table) 
{
  int n, d, k;
  double kgR, coskgR;
  int LATT = TYPE & KPT_LATT;

  // 1. zero out:
  for (d=0; d<9; ++d) DR[d] = 0;
  // 2. add it all up over k in IRZ and g[n] in our group
  for (k=0; k<Nkpt; ++k) {
    double *tkpt=kpt[k];
    double sumD[9], Drot[9], temp[9];
    for (d=0; d<9; ++d) sumD[d] = 0;
    for (n=0; n<Nop; ++n) {
      if (LATT) kgR = M_2PI*dot(p[tsh->gR0[n]].Runit, tkpt);
      else      kgR = dot(p[tsh->gR0[n]].Rcart, tkpt);
      coskgR = coskR_table(kgR, table);
      mult(gcart[inv_index[n]], Dk[k], temp);
      mult(temp, gcart[n], Drot);
      for (d=0; d<9; ++d) sumD[d] += coskgR*Drot[d];
    }
    for (d=0; d<9; ++d) DR[d] += sumD[d]*w[k];
  }
  // now, scale out by 1/Nop:
  for (d=0; d<9; ++d) DR[d] *= 1./Nop;
}

// A much simpler inv_fourier--doesn't scale by weights at all
void inv_fourier_table (point_type &p, double Nkpt, double** pk, 
			double** kpt, const int TYPE,
			const coskR_table_type &table) 
{
  int d, k;
  double kR, coskR;
  int LATT = TYPE & KPT_LATT;

  // 1. zero out:
  for (d=0; d<9; ++d) p.mat[d] = 0;
  // 2. add it all up over k
  for (k=0; k<Nkpt; ++k) {
    double *tkpt=kpt[k];
    if (LATT) kR = M_2PI*dot(p.Runit, tkpt);
    else      kR = dot(p.Rcart, tkpt);
    coskR = coskR_table(kR, table);
    for (d=0; d<9; ++d) p.mat[d] += coskR*pk[k][d];
  }
  mult(p.mat, 1./(double)Nkpt, p.mat); // divide it out.
}


// The actual code to do the closed-shell fourier transform: WITH TABLES!
// This assumes that our lattice functions:
// 1. is defined on closed shells (though it *technically* doesn't need that)
// 2. is composed of symmetric matrices
// 3. transforms as a second-rank tensor lattice function with the full
//    crystal point group; i.e.:
//      D(g*R) = g*D(R)*g^-1
//    where D is a 3x3 matrix, R is a lattice point, and g is a point
//    group operation on the crystal.
void fourier_cs_table (int Np, point_type* p, int Nsh, shell_type* sh,
		       int Nop, int inv_index[MAXop], 
		       double g_i_ginv[MAXop][Ni][Ni],
		       double kpt[3], double pk[9], const int TYPE,
		       const coskR_table_type &table) 
{
  int i, j, d;
  int n;
  double coskRA[Ni][Ni];
  double kgR, coskgR, gig;
  int LATT = TYPE & KPT_LATT;
  double mat6[Ni], mult_mat6[Ni];
  
  for (d=0; d<9; ++d) pk[d] = 0.;
  
  for (int ns=0; ns<Nsh; ++ns) {
    shell_type *tsh = sh + ns;
    for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i) coskRA[j][i] = 0;
    for (j=0; j<Ni; ++j) mat6[j] = p[tsh->elem[0]].mat[i6tom9[j]];
    for (int npi=0; npi<tsh->Nelem; ++npi) {
      int np = tsh->elem[npi]; // npi is just the index; np is the point number
      point_type *tp = p+np;
      if (LATT) kgR = 2.*M_PI*dot(tp->Runit, kpt);
      else      kgR = dot(tp->Rcart, kpt);
      coskgR = coskR_table(kgR, table);
      // Now to find which group operation we have:
      for (n=0; tsh->gR0[n] != np; ++n) ;
      // Now, actually add to our running total:
      for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i)
	// note carefully how we call g_i_ginv:
	if (! zero(gig=g_i_ginv[n][i][j]) )
	  coskRA[j][i] += coskgR*gig;
    }
    for (j=0; j<Ni; ++j) {
      mult_mat6[j] = 0.;
      for (i=0; i<Ni; ++i) mult_mat6[j] += coskRA[j][i]*mat6[i];
    }
    // add to our running FT:
    for (d=0; d<9; ++d) pk[d] += mult_mat6[m9toi6[d]];
  }
}


// The actual code to do the closed-shell INVERSE fourier transform:
// WITH TABLES!
// This assumes that our lattice functions:
// 1. is defined on closed shells--only used to get the point positions
// 2. is composed of symmetric matrices
// 3. transforms as a second-rank tensor lattice function with the full
//    crystal point group; i.e.:
//      D(g*k) = g*D(k)*g^-1
//    where D is a 3x3 matrix, k is a rlv, and g is a point group operation
//    on the crystal.
// In this case, we only use a single point in a given shell, but we
// pass the *entire* list of points, but just the shell we're interested in.
// The result for the first element of tsh goes into DR; we don't bother
// doing the symmetric rotation, or doing it in place for obvious reasons.

void inv_fourier_cs_table (point_type* p, shell_type* tsh, double DR[9],
			   int Nop, int inv_index[MAXop], 
			   double g_i_ginv[MAXop][Ni][Ni],
			   double Nkpt, double** Dk, 
			   double** kpt, double* w, const int TYPE,
			   const coskR_table_type &table) 
{
  int i, j, d;
  double coskRA[Ni][Ni];
  double kgR, coskgR, gig;
  int LATT = TYPE & KPT_LATT;
  double mat6[Ni], mult_mat6[Ni];
  double Nop_1 = 1./(double)Nop;

  // We could take this line out, since we're careful about how we
  // add everything; there's no scaling involved.
  for (d=0; d<9; ++d) DR[d] = 0.;
  
  for (int k=0; k<Nkpt; ++k) {
    double *tkpt = kpt[k];
    for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i) coskRA[j][i] = 0;
    for (j=0; j<Ni; ++j) mat6[j] = Dk[k][i6tom9[j]];
    for (int n=0; n<Nop; ++n) {
      point_type *tp = p+tsh->gR0[inv_index[n]]; // this is g^-1 R
      if (LATT) kgR = 2.*M_PI*dot(tp->Runit, tkpt);
      else      kgR = dot(tp->Rcart, tkpt);
      coskgR = coskR_table(kgR, table);
      // Now, actually add to our running total:
      for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i)
	// note carefully how we call g_i_ginv:
	if (! zero(gig=g_i_ginv[n][i][j]) )
	  coskRA[j][i] += coskgR*gig;
    }
    for (j=0; j<Ni; ++j) {
      mult_mat6[j] = 0.;
      for (i=0; i<Ni; ++i) mult_mat6[j] += Nop_1*coskRA[j][i]*mat6[i];
    }
    // add to our running FT:
    for (d=0; d<9; ++d) DR[d] += w[k]*mult_mat6[m9toi6[d]];
  }
}



// ************************** BINNING ROUTINES ***********************
// Test: try doing the sums with binning instead...

typedef struct 
{
  bin_type* bin[9];
  int alloc;
} mat_bin_type;

mat_bin_type* init_mat_bin (const double &minbin, const double &maxbin, 
                    const double &step) 
{
  if ( (minbin < 0) || (maxbin < 0) || (step <= 1) || (maxbin<minbin) )
    return NULL;
  mat_bin_type *mbinp = new mat_bin_type;
  for (int d=0; d<9; ++d)
    mbinp->bin[d] = init_bin(minbin, maxbin, step);
  mbinp->alloc=1;
  return mbinp;
}

inline mat_bin_type* init_mat_bin() 
{ return init_mat_bin(1./1099511627776., 1099511627776., 1024); }


void free_mat_bin (mat_bin_type* &mbinp) 
{
  if (mbinp != NULL) {
    for (int d=0; d<9; ++d) free_bin(mbinp->bin[d]);
    delete mbinp;
    mbinp = NULL;
  }
}

inline void clear_bin (const mat_bin_type* mbinp) 
{
  if (mbinp->alloc)
    for (int d=0; d<9; ++d) clear_bin(mbinp->bin[d]);
}

inline void sum_bin (const mat_bin_type* mbinp, double res[9]) 
{
  if (mbinp->alloc) {
    for (int d=0; d<9; ++d) {
      res[d] = 0;
      bin_type *binp = mbinp->bin[d];
      for (int i=0; i<binp->Nbins; ++i) res[d] += binp->val[i];
    }
  }
}

// add x[] to mbin
inline void inc_bin (const mat_bin_type* mbinp, double x[9]) 
{ for (int d=0; d<9; ++d) inc_bin(mbinp->bin[d], x[d]); }

// allows us to put a scale factor on...
inline void inc_bin (const mat_bin_type* mbinp, double x[9], double scale) 
{ for (int d=0; d<9; ++d) inc_bin(mbinp->bin[d], scale*x[d]); }

inline void inc_bin (const mat_bin_type* mbinp, double scale, double x[9]) 
{ inc_bin(mbinp, x, scale); }



void fourier_table_bin (int Np, point_type* p, double kpt[3], double pk[9], 
			const int TYPE, const coskR_table_type &table,
			const mat_bin_type *mbinp) 
{
  int n;
  double coskR;
  int LATT = TYPE & KPT_LATT;

  if (mbinp == NULL) return;
  if (! mbinp->alloc) return;
  clear_bin(mbinp); // initialize
  
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    if (LATT) coskR = M_2PI*dot(tp->Runit, kpt);
    else      coskR = dot(tp->Rcart, kpt);
    coskR = coskR_table(coskR, table); // now we use table lookup...
    inc_bin(mbinp, tp->mat, coskR);
  }

  sum_bin(mbinp, pk);
}

void inv_fourier_table_bin (point_type* p, shell_type *tsh, double DR[9],
			    int Nop, int inv_index[MAXop], 
			    double gcart[MAXop][9],
			    double Nkpt, double** Dk, 
			    double** kpt, double* w, const int TYPE,
			    const coskR_table_type &table,
			    const mat_bin_type *mbinp) 
{
  int n, d, k;
  double kgR, coskgR;
  int LATT = TYPE & KPT_LATT;

  if (mbinp == NULL) return;
  if (! mbinp->alloc) return;
  clear_bin(mbinp); // initialize
  
  // 1. add it all up over k in IRZ and g[n] in our group
  for (k=0; k<Nkpt; ++k) {
    double *tkpt=kpt[k];
    double sumD[9], Drot[9], temp[9];
    for (d=0; d<9; ++d) sumD[d] = 0;
    for (n=0; n<Nop; ++n) {
      if (LATT) kgR = M_2PI*dot(p[tsh->gR0[n]].Runit, tkpt);
      else      kgR = dot(p[tsh->gR0[n]].Rcart, tkpt);
      coskgR = coskR_table(kgR, table);
      mult(gcart[inv_index[n]], Dk[k], temp);
      mult(temp, gcart[n], Drot);
      for (d=0; d<9; ++d) sumD[d] += coskgR*Drot[d];
    }
    inc_bin(mbinp, sumD, w[k]);
  }
  // now, scale out by 1/Nop:
  sum_bin(mbinp, DR);
  for (d=0; d<9; ++d) DR[d] *= 1./Nop;
}




void fourier_cs_table_bin (int Np, point_type* p, int Nsh, shell_type* sh,
			   int Nop, int inv_index[MAXop], 
			   double g_i_ginv[MAXop][Ni][Ni],
			   double kpt[3], double pk[9], const int TYPE,
			   const coskR_table_type &table) 
{
  int i, j, d;
  int n;
  double coskRA[Ni][Ni];
  double kgR, coskgR, gig;
  int LATT = TYPE & KPT_LATT;
  double mat6[Ni], mult_mat6[Ni];
  
  bin_type* sum[Ni];
  for (i=0; i<Ni; ++i) sum[i] = init_bin(); // use the default values...
  
  for (int ns=0; ns<Nsh; ++ns) {
    shell_type *tsh = sh + ns;
    for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i) coskRA[j][i] = 0;
    for (j=0; j<Ni; ++j) mat6[j] = p[tsh->elem[0]].mat[i6tom9[j]];
    for (int npi=0; npi<tsh->Nelem; ++npi) {
      int np = tsh->elem[npi]; // npi is just the index; np is the point number
      point_type *tp = p+np;
      if (LATT) kgR = 2.*M_PI*dot(tp->Runit, kpt);
      else      kgR = dot(tp->Rcart, kpt);
      coskgR = coskR_table(kgR, table);
      // Now to find which group operation we have:
      for (n=0; tsh->gR0[n] != np; ++n) ;
      // Now, actually add to our running total:
      for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i)
	// note carefully how we call g_i_ginv:
	if (! zero(gig=g_i_ginv[n][i][j]) )
	  coskRA[j][i] += coskgR*gig;
    }
    for (j=0; j<Ni; ++j) {
      mult_mat6[j] = 0.;
      for (i=0; i<Ni; ++i) mult_mat6[j] += coskRA[j][i]*mat6[i];
    }
    // add to our running FT:
    for (i=0; i<Ni; ++i) inc_bin(sum[i], mult_mat6[i]);
  }

  // put back into normal matrix form:
  for (d=0; d<9; ++d) pk[d] = sum_bin(sum[m9toi6[d]]);

  // garbage collection
  for (i=0; i<Ni; ++i) free_bin(sum[i]);
}


void inv_fourier_cs_table_bin (point_type* p, shell_type* tsh, double DR[9],
			       int Nop, int inv_index[MAXop], 
			       double g_i_ginv[MAXop][Ni][Ni],
			       double Nkpt, double** Dk, 
			       double** kpt, double* w, const int TYPE,
			       const coskR_table_type &table) 
{
  int i, j, d;
  double coskRA[Ni][Ni];
  double kgR, coskgR, gig;
  int LATT = TYPE & KPT_LATT;
  double mat6[Ni], mult_mat6[Ni];
  double Nop_1 = 1./(double)Nop;

  bin_type* sum[Ni];
  for (i=0; i<Ni; ++i) sum[i] = init_bin(); // use the default values...
  
  for (int k=0; k<Nkpt; ++k) {
    double *tkpt = kpt[k];
    for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i) coskRA[j][i] = 0;
    for (j=0; j<Ni; ++j) mat6[j] = Dk[k][i6tom9[j]];
    for (int n=0; n<Nop; ++n) {
      point_type *tp = p+tsh->gR0[inv_index[n]]; // this is g^-1 R
      if (LATT) kgR = 2.*M_PI*dot(tp->Runit, tkpt);
      else      kgR = dot(tp->Rcart, tkpt);
      coskgR = coskR_table(kgR, table);
      // Now, actually add to our running total:
      for (j=0; j<Ni; ++j) for (i=0; i<Ni; ++i)
	// note carefully how we call g_i_ginv:
	if (! zero(gig=g_i_ginv[n][i][j]) )
	  coskRA[j][i] += coskgR*gig;
    }
    for (j=0; j<Ni; ++j) {
      mult_mat6[j] = 0.;
      for (i=0; i<Ni; ++i) mult_mat6[j] += Nop_1*coskRA[j][i]*mat6[i];
    }
    for (i=0; i<Ni; ++i) inc_bin(sum[i], w[k]*mult_mat6[i]);
  }

  // put back into normal matrix form:
  for (d=0; d<9; ++d) DR[d] = sum_bin(sum[m9toi6[d]]);

  // garbage collection
  for (i=0; i<Ni; ++i) free_bin(sum[i]);
}


// ************************** Small k expansions *********************
// We're interested in the FT for small k... so here they is.
// We actually DON'T calc. the FT; instead of expanding with cos(k.R),
// we expand with 1-1/2(k.R)^2-cos(k.R).  We do this with table-lookup
// as well, but perhaps most importantly is that we also use a Taylor
// series expansion for small k.R values so that we don't die with
// accuracy.

// we use these for doing the Taylor series; they're 1/4! 1/6! 1/8! and
// 1/10! out to 20 digits (thanks to bc)
#define __QUARTIC_4  4.16666666666666666667e-2
#define __QUARTIC_6  1.38888888888888888889e-3
#define __QUARTIC_8  2.48015873015873015873e-5
#define __QUARTIC_10 2.75573192239858906525e-7
#define __QUARTIC_THRESHOLD 0.5

inline double cos_quartic_expand (const double &x2) 
{ return x2*x2*(-__QUARTIC_4
		+ x2*(__QUARTIC_6
		      + x2*(-__QUARTIC_8
			   + x2*__QUARTIC_10)));
}

inline double cos_quartic (const double &x) 
{
  double x2 = x*x;
  if (x<__QUARTIC_THRESHOLD) {
    // Taylor expansion:
    return cos_quartic_expand(x2);
  }
  // else:
  return 1-0.5*x2-cos(x);
}

inline double cos_quartic_table (const double &x, 
				 const coskR_table_type &table) 
{
  double x2 = x*x;
  if (x<__QUARTIC_THRESHOLD) {
    // Taylor expansion:
    return cos_quartic_expand(x2);
  }
  // else:
  return 1-0.5*x2-coskR_table(x, table);
}

#undef __QUARTIC_4
#undef __QUARTIC_6
#undef __QUARTIC_8
#undef __QUARTIC_10
#undef __QUARTIC_THRESHOLD


void fourier_smallk (int Np, point_type* p, double kpt[3], double pk[9], 
		     const int TYPE) 
{
  int n, d;
  double kR;
  int LATT = TYPE & KPT_LATT;

  for (d=0; d<9; ++d) pk[d] = 0.;
  
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    if (LATT) kR = M_2PI*dot(tp->Runit, kpt);
    else      kR = dot(tp->Rcart, kpt);
    kR = cos_quartic(kR);
    for (d=0; d<9; ++d) pk[d] += kR*(tp->mat[d]);
  }
}

void fourier_smallk_table (int Np, point_type* p, double kpt[3], double pk[9], 
			   const int TYPE, const coskR_table_type &table) 
{
  int n, d;
  double kR;
  int LATT = TYPE & KPT_LATT;

  for (d=0; d<9; ++d) pk[d] = 0.;
  
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    if (LATT) kR = M_2PI*dot(tp->Runit, kpt);
    else      kR = dot(tp->Rcart, kpt);
    kR = cos_quartic_table(kR, table);
    for (d=0; d<9; ++d) pk[d] += kR*(tp->mat[d]);
  }
}



#undef M_2PI
#undef M_1_2PI

#endif
