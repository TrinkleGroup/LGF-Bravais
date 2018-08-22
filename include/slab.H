#ifndef __SLAB_H
#define __SLAB_H

/*
  Program: slab.H
  Author:  D. Trinkle
  Date:    October 6, 2004
  Purpose: Construct a cylindrical slab.  Pretty icky stuff, but it's here.
           The algorithm is stolen from the Dislocation code, and modified
	   to work with our terrible data structures.
*/

#include <stdio.h>
#include <math.h>
#include "dcomp.H"
#include "matrix.H"
#include "Dij.H"
#include "dislocation.H"
#include "shell.H" // for sortin'

// Makes a cylindrical slab, where:
// t is the vertical axis (will be the z axis)
// m is the cut axis (will be the x axis)
// n is the mutual perp (will be the y axis)
// Rcut is the radius in the xy plane
// c is the cart. coord. of the center of the dislocation
// |m|=|n| = 1, while |t| is the thickness of the slab
// [a] is the cartesian coordinates of the lattice
// u is the matrix of atoms, in unit cell coord.
// Natoms is the number of atoms in the unit cell
//   We output:
// Nslab: number of atoms in the slab
// xyz: xyz positions of atoms in the slab (xyz[n][1,2,3])
int construct_slab (double t[3], double m[3], double n[3], double c[3],
		    double Rcut,
		    double a[9], double** u, int Natoms,
		    int& Nslab, double** &xyz);

// a "point_type"-friendly version of construct_slab:
int construct_slab (double cart[9], double Rcut,
		    int t_unit[3], int m_unit[3], int n_unit[3],
		    int& Nslab, point_type* &p) 
{
  double** xyz;
  double** u, c[3] = {0,0,0};
  double tvect[3], mvect[3], nvect[3];
  int ERROR=0;
  
  u = new double*[1];
  u[0] = new double[3];
  for (int i=0; i<3; ++i) u[0][i] = 0.;

  // get our coordinates:
  make_dislcoord(cart, tvect, mvect, nvect, t_unit, m_unit, n_unit,
		 DISLCOORD_UNIT + DISLCOORD_NOCHANGE);
  // we need to "redo" tvect, because make_dislcoord scales it to
  // magnitude 1:
  mult_vect(cart, t_unit, tvect);
  ERROR = construct_slab(tvect, mvect, nvect, c, Rcut, cart, u, 1,
			 Nslab, xyz);
  // garbage collection:
  delete[] u[0];  delete[] u;

  if (ERROR) return ERROR;
  // Now, put into our data structure, and rotate back to Cart. coord:
  p = new point_type[Nslab];

  // rotation matrix:
  double theta[9], tmagn = 1./sqrt(dot(tvect, tvect));
  theta[0] = mvect[0];  theta[1] = nvect[0];  theta[2] = tvect[0]*tmagn;
  theta[3] = mvect[1];  theta[4] = nvect[1];  theta[5] = tvect[1]*tmagn;
  theta[6] = mvect[2];  theta[7] = nvect[2];  theta[8] = tvect[2]*tmagn;
  double cart_inv[9], dump[3];
  careful_inverse(cart, cart_inv);
  for (int n=0; n<Nslab; ++n) {
    point_type *tp = p + n;
    double* r = xyz[n];
    tp->Rmagn = sqrt(r[0]*r[0] + r[1]*r[1]); // only keep dist. in plane
    mult_vect(theta, r, tp->Rcart);
    mult_vect(cart_inv, tp->Rcart, dump);
    for (int i=0; i<3; ++i) tp->Runit[i] = lround(dump[i]);
    // garbage collection...
    delete[] xyz[n];
  }
  delete[] xyz;
  // finally, sort our points by magnitude:
  sort_points(p, Nslab);
  return ERROR;
}


int construct_slab (double t[3], double m[3], double n[3], double c[3],
		    double Rcut,
		    double a[9], double** u, int Natoms,
		    int& Nslab, double** &xyz) 
{
  int ERROR = 0;
  int i, j, k;
  // S = (m,n,t): transformation matrix
  // Sa = S^-1*[a]
  // aS = [a]^-1*S
  double S[9], Sa[9], aS[9];
  double ainv[9], Sinv[9];
  double deter;
  double Rcut2;
  double tmagn;

  Rcut2 = Rcut*Rcut;

  // Construct the rotation matrix
  S[0] = m[0];  S[1] = n[0];  S[2] = t[0];
  S[3] = m[1];  S[4] = n[1];  S[5] = t[1];
  S[6] = m[2];  S[7] = n[2];  S[8] = t[2];
  tmagn = det(S);

  deter = 1./inverse(S, Sinv);
  for (i=0; i<9; ++i) Sinv[i] *= deter;
  mult(Sinv, a, Sa);

  deter = 1./inverse(a, ainv);
  for (i=0; i<9; ++i) ainv[i] *= deter;
  mult(ainv, S, aS); // aS = [a]^-1 * S
  // We use this to make our maximums for i0, i1, i2:
  int imax[3];
  for (i=0; i<3; ++i) {
    imax[i] = (int)(Rcut*(fabs(aS[3*i]) + fabs(aS[3*i+1])) + fabs(aS[3*i+2])
		    + 1.999);
  }
  
  int Napprox;
  int N;

  // Guess how many atoms we'll end up with at the end of the day...
  Napprox = (int)(2*Natoms * M_PI * Rcut*Rcut * tmagn / det(a));
  if (Napprox < 128) Napprox = 128;
  double** s;
  s = new double*[3];
  for (i=0; i<3; ++i) s[i] = new double[Napprox];

  // Calculate the atom positions for the unit cell, with the shift
  // built in:
  double** s_atom;
  double cu[3]; // Center, in unit coord.

  mult_vect(ainv, c, cu);
  for (i=0; i<3; ++i) cu[i] = insidecell(cu[i]);

  s_atom = new double*[Natoms];
  for (j=0; j<Natoms; ++j) {
    s_atom[j] = new double[3];
    double ushift[3];
    // Start inside the unit cell:
    for (i=0; i<3; ++i) ushift[i] = insidecell(u[j][i] - cu[i]);
    mult_vect(Sa, ushift, s_atom[j]);
    // Make sure we're still inside the cylinder in the z direction:
    s_atom[j][2] = insidecell(s_atom[j][2]);
  }

  // Now, try to put all the atoms in.
  int cell[3];
  double scell[3], stry[3];
  int found; // Collision detection
  N = 0;
  for (cell[0]=-imax[0]; cell[0]<=imax[0]; ++cell[0])
    for (cell[1]=-imax[1]; cell[1]<=imax[1]; ++cell[1])
      for (cell[2]=-imax[2]; cell[2]<=imax[2]; ++cell[2]) {
	mult_vect(Sa, cell, scell);
	for (j=0; j<Natoms; ++j) {
	  for (i=0; i<3; ++i) stry[i] = scell[i] + s_atom[j][i];
	  if ( (stry[0]*stry[0] + stry[1]*stry[1]) < Rcut2 ) {
	    // Attempt to check for "collisions"
	    found = dcomp(stry[2], 0.) || dcomp(stry[2], 1.);
	    if (found) {
	      // Do the *real* collision detection
	      found = 0;
	      for (k=0; (k<N) && (!found); ++k)
		found = dcomp(stry[0], s[0][k]) && dcomp(stry[1], s[1][k]);
	      stry[2] = insidecell(stry[2]);
	    }
	    else
	      // Else, only add if we're between 0 and 1 along z axis:
	      found = (stry[2] < 0.) || ( stry[2] > 1.);
	    // We only add if *not* found.
	    if (!found) {
	      if (N >= Napprox)
		{++N; break;}
	      for (i=0; i<3; ++i)
		s[i][N] = stry[i];
	      ++N;
	    }
	  }
	}
      }
  if (N>Napprox) {
    fprintf(stderr, "Error constructing cylinder... something went horribly awry\n");
    fprintf(stderr, "Found %d atoms, when estimated %d atoms.\n", N, Napprox);
    fprintf(stderr, "Continuing, but atoms will be missing.  Change construct_slab() routine.\n");
    ERROR = -1;
    N = Napprox;
  }
  // Now, convert into xyz coordinates:
  Nslab = N;
  xyz = new double *[N];
  for (i=0; i<N; ++i) {
    xyz[i] = new double[3];
    xyz[i][0] = s[0][i];
    xyz[i][1] = s[1][i];
    xyz[i][2] = s[2][i] * tmagn;
  }
  // Garbage collection
  for (j=0; j<Natoms; ++j)
    delete[] s_atom[j];
  delete[] s_atom;
  
  for (i=0; i<3; ++i) 
    delete[] s[i];
  delete[] s;

  return ERROR;
}

#endif
