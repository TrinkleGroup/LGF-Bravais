#ifndef __SUPERCELL__
#define __SUPERCELL__

/* A few routines to deal with creating supercells.

  free_super (Natoms, atomlist)
  --frees space allocated for our supercell list (fill_super or 
    make_super_shifts)

  fill_super (cart, n, Natoms, atomlist)
  --makes a Wigner-Seitz supercell, and sorts the list by distance
    from origin

  free_index (maxn, index)
  --frees space allocated for our index

  make_index (Np, points, maxn, index)
  --takes a set of points in direct coordinates, and makes an "index"
    of them.  That is, index[maxn+i][maxn+j][maxn+k] = -1 if the point
    isn't in points; else, it's the index value in points[].

  get_index(maxn, index, point)
  --returns the index value for a point; does range-checking

  check_point(maxn, index, point)
  --returns 0 or -1 if point is not or is in index

  make_super_shifts (maxn, index, n, Nshift, shiftlist)
  --calculates all of the supercell commensurate shifts needed to sum
    up over points provided we're *only* using points inside the supercell
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcomp.H"
#include "matrix.H"
#include "unit.H"
#include "Dij.H"
#include "eigen.H"
#include "pointgroup.H"

void free_super (int Natoms, int** &atomlist) 
{
  int n;
  if (atomlist != NULL) {
    for (n=0; n<Natoms; ++n)
      delete[] atomlist[n];
    delete[] atomlist;
    atomlist = NULL;
  }
}

// We make two versions of fill_super, and call one from the other:

void fill_super (double cart[9], int supern[9], 
		 int &Natoms, int** &atomlist);

void fill_super (double cart[9], int supern[9],
		 int &Natoms, point_type* &p) 
{
  int n, d;
  int** atomlist;
  
  fill_super(cart, supern, Natoms, atomlist);
  p = new point_type[Natoms];
  for (n=0; n<Natoms; ++n) {
    point_type *tp = p+n;
    for (d=0; d<3; ++d) tp->Runit[d] = atomlist[n][d];
    mult_vect(cart, tp->Runit, tp->Rcart);
    tp->Rmagn = sqrt(dot(tp->Rcart, tp->Rcart));
    //    for (d=0; d<9; ++d) tp->lat[d] = 0.;
  }
  free_super(Natoms, atomlist);
}

void fill_super (double cart[9], int supern[9], 
		 int &Natoms, int** &atomlist) 
{
  int i, na, n;
  int m[3];
  int mmax;
  int invn[9];
  double cartsq[9];
  unit_point_type *work;  // comes from unit.H
  
  square(cart, cartsq);

  Natoms = inverse(supern, invn); // Natoms = det(supern)
  if (Natoms < 0) Natoms = -Natoms;

  work = new unit_point_type[Natoms];

  // First, construct our Natoms list
  for (mmax=0, i=0; i<9; ++i)
    if (mmax < abs(supern[i])) mmax = abs(supern[i]);
  mmax *= 2; // safety...
  na = 0;
  for (m[0]=-mmax; m[0] <= mmax; ++(m[0])) {
    for (m[1]=-mmax; m[1] <= mmax; ++(m[1])) {
      for (m[2]=-mmax; m[2] <= mmax; ++(m[2])) {
	int temp[3];
	mult_vect(invn, m, temp);
	if ( (temp[0] >= 0) && (temp[0] < Natoms) &&
	     (temp[1] >= 0) && (temp[1] < Natoms) &&
	     (temp[2] >= 0) && (temp[2] < Natoms) ) {
	  // We've got a new atom in our supercell!
	  for (i=0; i<3; ++i) work[na].R[i] = m[i];
	  ++na;
	}
	if (na == Natoms) break;
      }
      if (na == Natoms) break;
    }
    if (na == Natoms) break;
  }
  if (na != Natoms) {
    fprintf(stderr, "HORRIBLE ERROR in fill_super!  Bailing out.\n");
    exit(-1);
  }
  // Now fold down into the Wigner-Seitz cell
  // First construct the shifts:
  const int Nshift = 27;
  int alpha[3], nshift[27][3];
  for (alpha[0]=-1; alpha[0] <= 1; ++(alpha[0]))
    for (alpha[1]=-1; alpha[1] <= 1; ++(alpha[1]))
      for (alpha[2]=-1; alpha[2] <= 1; ++(alpha[2]))
	mult_vect(supern, alpha, nshift[13+alpha[0]+3*(alpha[1]+3*alpha[2])]);

  // We do this by adding supercell vectors to each:
  int trym[3], m0[3];
  double bestm2, magn2;
  for (na=0; na<Natoms; ++na) {
    int* mp = work[na].R;
    for (i=0; i<3; ++i) m0[i] = mp[i];
    bestm2 = magnsq(cartsq, m0);
    for (n=0; n<Nshift; ++n) {
      for (i=0; i<3; ++i) trym[i] = m0[i] + nshift[n][i];
      magn2 = magnsq(cartsq, trym);
      if (magn2 < bestm2) {
	bestm2 = magn2;
	for (i=0; i<3; ++i) mp[i] = trym[i];
      }
    }
    work[na].Rmagn = magnsq(cartsq, work[na].R);
  }

  // Final step: sort by Rmagn
  qsort(work, Natoms, sizeof(unit_point_type), __COMPARE__UNIT__);
  
  // Copy to atomlist:
  atomlist = new int*[Natoms];  
  for (na=0; na<Natoms; ++na) {
    atomlist[na] = new int[3];
    for (i=0; i<3; ++i) atomlist[na][i] = work[na].R[i];
  }
  delete[] work;
}


// frees space allocated for our index
void free_index (int maxn, int*** &index) 
{
  int i, j;
  if (index != NULL) {
    for (i=0; i<=(2*maxn); ++i) {
      for (j=0; j<=(2*maxn); ++j) 
	delete[] index[i][j];
      delete[] index[i];
    }
    delete[] index;
    index = NULL;
  }
}

    
//takes a set of points in direct coordinates, and makes an "index"
// of them.  That is, index[maxn+i][maxn+j][maxn+k] = -1 if the point
// isn't in points; else, it's the index value in points[].
// Two flavors: the first determines maxn; other assumes it.
void make_index (int Np, int** points, int &maxn, int*** &index);
void make_index (int maxn, int*** &index);

void make_index (int Np, point_type *p, int &maxn, int*** &index) 
{
  int n;
  int **atomlist;

  atomlist = new int*[Np];
  for (n=0; n<Np; ++n) {
    atomlist[n] = p[n].Runit; // very sneaky...
    //    atomlist[n] = new int[3];
    //    for (d=0; d<3; ++d)
    //      atomlist[n][d] = p[n].Runit[d];
  }
  make_index(Np, atomlist, maxn, index);
  // Garbage collection...
  delete[] atomlist;
}


void make_index (int Np, int** points, int &maxn, int*** &index) 
{
  int i, j, k;

  // determine maxn first:
  for (i=0, maxn=0; i<Np; ++i) {
    if (maxn < abs(points[i][0])) maxn = abs(points[i][0]);
    if (maxn < abs(points[i][1])) maxn = abs(points[i][1]);
    if (maxn < abs(points[i][2])) maxn = abs(points[i][2]);
  }
  
  // allocate space and initialize
  index = new int**[2*maxn+1];
  for (i=0; i<=(2*maxn); ++i) {
    index[i] = new int*[2*maxn+1];
    for (j=0; j<=(2*maxn); ++j) {
      index[i][j] = new int[2*maxn+1];
      // Initial value:
      for (k=0; k<=(2*maxn); ++k)
	index[i][j][k] = -1;
    }
  }
  
  // Now, fill up index:
  for (i=0; i<Np; ++i)
    index[maxn+points[i][0]][maxn+points[i][1]][maxn+points[i][2]] = i;
}

void make_index (int maxn, int*** &index) 
{
  // allocate space and initialize
  index = new int**[2*maxn+1];
  for (int i=0; i<=(2*maxn); ++i) {
    index[i] = new int*[2*maxn+1];
    for (int j=0; j<=(2*maxn); ++j) {
      index[i][j] = new int[2*maxn+1];
      // Initial value:
      for (int k=0; k<=(2*maxn); ++k)
	index[i][j][k] = -1;
    }
  }
}


// Also does some range-checking
inline int get_index (int maxn, int*** index, int point[3]) 
{
  if ( (abs(point[0]) > maxn) ||
       (abs(point[1]) > maxn) ||
       (abs(point[2]) > maxn) )
    return -1;
  return index[maxn+point[0]][maxn+point[1]][maxn+point[2]];
}


inline void set_index (int maxn, int*** index, int point[3], int val) 
{
  if ( (abs(point[0]) > maxn) ||
       (abs(point[1]) > maxn) ||
       (abs(point[2]) > maxn) )
    return;
  index[maxn+point[0]][maxn+point[1]][maxn+point[2]] = val;
}

inline int check_point (int maxn, int*** index, int point[3]) 
{
  return (get_index(maxn, index, point) != -1);
}


// Symmetric version of the indexing functions:
typedef struct 
{
  int index, op;
} symm_index_type;


// Also does some range-checking
inline symm_index_type *get_symm_index (const int maxn, 
					symm_index_type*** index,
					const int point[3]) 
{
  if ( (abs(point[0]) > maxn) ||
       (abs(point[1]) > maxn) ||
       (abs(point[2]) > maxn) )
    return NULL;
  return &(index[maxn+point[0]][maxn+point[1]][maxn+point[2]]);
}

inline int check_symm_point (const int maxn, symm_index_type*** index, 
			     const int point[3]) 
{return (get_symm_index(maxn, index, point)->index != -1);}


void make_symm_index(int Np, point_type *GL, int Nop, int gunit[MAXop][9],
		     int &maxn, symm_index_type*** &index) 
{
  int i, j, k, n;
  int gR[3];

  // determine maxn first:
  // now fixed to get maxn correctly (had trouble with symmetry producing
  // index values *larger* than maxn)--now we run through the symmetric
  // options as well.
  for (i=0, maxn=0; i<Np; ++i) {
    int *tp = GL[i].Runit;
    for (n=0; n<Nop; ++n) {
      mult_vect(gunit[n], tp, gR);
      for (int d=0; d<3; ++d)
	if (maxn < abs(gR[d])) maxn = abs(gR[d]);
    }
  }
  
  // allocate space and initialize
  index = new symm_index_type **[2*maxn+1];
  for (i=0; i<=(2*maxn); ++i) {
    index[i] = new symm_index_type *[2*maxn+1];
    for (j=0; j<=(2*maxn); ++j) {
      index[i][j] = new symm_index_type[2*maxn+1];
      // Initial value:
      for (k=0; k<=(2*maxn); ++k) 
	index[i][j][k].index = -1;
    }
  }
  
  // Now, fill up index:
  for (i=0; i<Np; ++i) {
    int *tp = GL[i].Runit;
    for (n=0; n<Nop; ++n) {
      mult_vect(gunit[n], tp, gR);
      symm_index_type *ti = get_symm_index(maxn, index, gR);
      if (ti->index == -1) {
	ti->index = i;
	ti->op = n;
      }
    }
  }
}
  
// frees space allocated for our index
void free_symm_index (int maxn, symm_index_type*** &index) 
{
  int i, j;
  if (index != NULL) {
    for (i=0; i<=(2*maxn); ++i) {
      for (j=0; j<=(2*maxn); ++j) 
	delete[] index[i][j];
      delete[] index[i];
    }
    delete[] index;
    index = NULL;
  }
}




// calculates all of the supercell commensurate shifts needed to sum
// up over points provided we're *only* using points inside the supercell
// NOTE: these points are in NO ORDER
void make_super_shifts (int maxn, int*** index, int supern[9],
			int &Nshift, int** &shiftlist) 
{
  int i, j;
  
  const int Nbase=27;
  int baseshift[Nbase][3];
  int alpha[3];
  
  // Calculate all of our possible lattice shifts:
  for (alpha[0]=-1; alpha[0]<=1; ++(alpha[0]))
    for (alpha[1]=-1; alpha[1]<=1; ++(alpha[1]))
      for (alpha[2]=-1; alpha[2]<=1; ++(alpha[2]))
	mult_vect(supern, alpha, 
		  baseshift[13+alpha[0]+3*(alpha[1]+3*alpha[2])]);
  // Now, calc. mmax by looking at the eigenvalues of [n]^T[n]:
  double nsq[9], lambda[3];
  square(supern, nsq); // also takes care of converting int to double
  eigen(nsq, lambda); // Note: lambda[0] is the smallest value
  int mmax;
  mmax = lround(maxn / sqrt(lambda[0])) + 2;

  // Next, let's try to make 'em:
  int m[3];
  int MAXshift; // (2*mmax+1)^3
  MAXshift = 1+mmax*(6+mmax*(12+mmax*8));
  shiftlist = new int*[MAXshift];
  Nshift=0;
  for (m[0]=-mmax; m[0]<=mmax; ++(m[0]))
    for (m[1]=-mmax; m[1]<=mmax; ++(m[1]))
      for (m[2]=-mmax; m[2]<=mmax; ++(m[2])) {
	int nmvect[3], nmshift[3];
	int found;
	mult_vect(supern, m, nmvect);
	for (found=0, i=0; (i<Nbase) && (!found); ++i) {
	  for (j=0; j<3; ++j) nmshift[j]=nmvect[j]+baseshift[i][j];
	  found = check_point(maxn, index, nmshift);
	}
	if (found) {
	  // Add our point!
	  shiftlist[Nshift] = new int[3];
	  for (j=0; j<3; ++j)
	    shiftlist[Nshift][j] = nmvect[j];
	  ++Nshift;
	}
      }
  if (Nshift > MAXshift) {
    fprintf(stderr, "Something went horribly wrong in make_super_shifts.\n");
  }
}


#endif
