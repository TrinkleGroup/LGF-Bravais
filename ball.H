#ifndef __BALL_H
#define __BALL_H

/*
  Program: ball.H
  Author:  D. Trinkle
  Date:    Feb 13, 2004
  Purpose: Provide code to generate a spherical ball of atoms around
           given a lattice.  We also handle the smooth smearing of the
	   sums.

	   This head file contains the following subroutines:

	   refine_Rcut (cart, Rcut)
	     --takes in our initial value for Rcut and a lattice given
	       by cart, and tries to produce an equivalent value for
	       Rcut that lies midway between the two shells at that
	       distance.  In this way, the new Rcut doesn't "hit"
	       any atoms.

	  smear(Rcut, a, R)
	    --returns a smearing factor for R, given Rcut (outer
	      limit) and a width a.  We use a cubic smearing.  For 
	      R > Rcut, smear is 0; R < Rcut-a, smear is 1, and the
	      intermediate values smoothly interpolate betwixt.

	  atom_sphere(cart, Rcut, Natoms, Rlist)
	    --takes in a lattice defined by cart, and a cutoff distance
	      of Rcut, and makes a sphere of atoms that are closer
	      than that distance.  This list is sorted from largest
	      distance to smallest, and returned as an array:
	      Rlist[n][0..2]: normalized direction
	      Rlist[n][3]:    magnitude
	     
 	  atom_sphere_unit(cart, Rcut, Natoms, Runit)
	    --same as atom_sphere, except the returned list are
	      unit cell coordinates (i.e., integers).  This list is
	      sorted from smallest distance to largest and returned
	      as an array (no information about magnitude, however).
	      Runit[n][0..2]: unit cell coordinates
*/

#include "matrix.H"
#include "unit.H"
#include "eigen.H"

/*============================== FUNCTIONS ===========================*/

//================================ smear ===============================
// Our smearing function:
// smear = 0 for R>=Rcut
// smear = 1 for R<=Rcut-a
// smear = 1/2 for R=Rcut-a/2
// ... and we smoothly interpolate in between using a cubic polynomial
double smear (double Rcut, double a, double R) 
{
  if (R>=Rcut) return 0;
  if (R<=(Rcut-a)) return 1;
  double u=(2*(R-Rcut)/a+1); // value between -1 and 1.
  return 0.5 - u*(0.75 - 0.25*u*u);
}

//============================= refine_Rcut ============================
// Takes in an initial guess for Rcut, and places it so that it's not
// "on top of" an atom position.
double refine_Rcut (double cart[9], double Rcut) 
{
  int n[3];
  double cartsq[9];
  int nmax;
  double rlo2, rmid2, rhi2;
  double Rcut2, Rcut2_1;
  double magn;
  
  square(cart, cartsq);
  double lambda[3];
  eigen(cartsq, lambda);
  nmax = lround(Rcut/sqrt(lambda[0]))+1; // lambda[0] is the smallest eigenv.

  Rcut2 = Rcut*Rcut; Rcut2_1 = 1./Rcut2;
  rlo2 = 0.; // Need to bracket our cutoff
  rhi2 = Rcut2*100.0;
  rmid2 = 0; // Assume no mid is found...
  for (n[0]=-nmax; n[0]<=nmax; ++(n[0])) 
    for (n[1]=-nmax; n[1]<=nmax; ++(n[1])) 
      for (n[2]=-nmax; n[2]<=nmax; ++(n[2])) {
	magn = magnsq(cartsq, n);
	// First, treat the special case that we have Rcut which
	// falls on an atom:
	if (dcomp(magn*Rcut2_1, 1.)) rmid2 = magn;
	else {
	  // Next, search for highs and lows.
	  if ((magn > Rcut2) && (magn < rhi2)) rhi2 = magn;
	  if ((magn < Rcut2) && (magn > rlo2)) rlo2 = magn;
	}
      }
  // Now, Rcut2 should be right between rhi2 and rlo2; so return the
  // *average* of the square root of the two
  // Check if we found an atom at distance Rcut; if we did, overestimate:
  if (rmid2 == 0)
    return 0.5*(sqrt(rlo2) + sqrt(rhi2));
  else
    return 0.5*(sqrt(rmid2) + sqrt(rhi2));
}

// Comparison routine used for sorting in atom_sphere.
// NOTE: it's "backwards"--since we want the list sorted in descending
// order
int __COMPARE__RVECT__ (const void* a, const void* b) 
{
  double amag, bmag;
  amag = ((double **)a)[0][3]; // I know--it looks weird :)
  bmag = ((double **)b)[0][3];
  if (dcomp(amag,bmag)) return 0;
  if (amag<bmag) return 1;
  return -1;
}

// Note: creates an atom at the origin!
int atom_sphere(double cart[9], double Rcut, 
		int &Natoms, double** &Rlist) 
{
  int n[3], i;
  double cartsq[9];
  int nmax;
  double Rcut2;
  double magn;

  int ERROR = 0;
  
  square(cart, cartsq);
  double lambda[3], rvec[3];
  eigen(cartsq, lambda);
  nmax = lround(Rcut/sqrt(lambda[0]))+1; // lambda[0] is the smallest eigenv.
  int nc; // number of atom
  
  Natoms = 8*nmax*nmax*nmax; // Our upper limit...
  Rlist = new double*[Natoms];
  
  Rcut2 = Rcut*Rcut;
  nc = 0;
  for (n[0]=-nmax; n[0]<=nmax; ++(n[0])) {
    for (n[1]=-nmax; n[1]<=nmax; ++(n[1])) {
      for (n[2]=-nmax; n[2]<=nmax; ++(n[2])) {
	magn = magnsq(cartsq, n);
	if (magn<Rcut2) {
	  // Add our atom!
	  Rlist[nc] = new double[4];
	  mult_vect(cart, n, rvec);
	  // Scale out magn:
	  magn = sqrt(magn);
	  if (nc == Natoms) {
	    fprintf(stderr, "Error in atom_sphere; too many atoms??\n");
	    ERROR = -1;
	    break;
	  }
	  if (dcomp(magn, 0.)) {
	    Rlist[nc][0] = 0;
	    Rlist[nc][1] = 0;
	    Rlist[nc][2] = 0;
	    Rlist[nc][3] = 0;
	  }
	  else {
	    Rlist[nc][3] = magn;
	    for (i=0; i<3; ++i) Rlist[nc][i] = rvec[i]/magn;
	  }
	  ++nc;
	}
      }
      if (ERROR) break;
    }
    if (ERROR) break;
  }
  // Now, NULL out the "extra" atoms
  for (i=nc; i<Natoms; ++i) Rlist[i] = NULL;
  Natoms = nc;

  // Final step: sort by Rlist[n][4]
  qsort((void *)Rlist, Natoms, sizeof(double *), __COMPARE__RVECT__);
  return ERROR;
}


// Note: creates an atom at the origin!
int atom_sphere_unit(double cart[9], double Rcut, 
		     int &Natoms, int** &Runit) 
{
  int n[3], i;
  double cartsq[9];
  int nmax;
  double Rcut2;
  double magn;
  unit_point_type *work;  // comes from unit.H

  int ERROR = 0;
  
  square(cart, cartsq);
  double lambda[3];
  eigen(cartsq, lambda);
  nmax = lround(Rcut/sqrt(lambda[0]))+1; // lambda[0] is the smallest eigenv.
  int nc; // number of atom
  
  Natoms = 8*nmax*nmax*nmax; // Our upper limit...
  work = new unit_point_type[Natoms];
  
  Rcut2 = Rcut*Rcut;
  nc = 0;
  for (n[0]=-nmax; n[0]<=nmax; ++(n[0])) {
    for (n[1]=-nmax; n[1]<=nmax; ++(n[1])) {
      for (n[2]=-nmax; n[2]<=nmax; ++(n[2])) {
	magn = magnsq(cartsq, n);
	if (magn<Rcut2) {
	  // Add our atom!
	  for (i=0; i<3; ++i) work[nc].R[i] = n[i];
	  work[nc].Rmagn = magn; // Can keep squared magnitude
	  ++nc;
	}
      }
      if (ERROR) break;
    }
    if (ERROR) break;
  }
  // Reset our atom count:
  Natoms = nc;

  // Final step: sort by Rlist[n][4]
  qsort(work, Natoms, sizeof(unit_point_type), __COMPARE__UNIT__);

  // Copy to Runit:
  Runit = new int*[Natoms];
  for (nc=0; nc<Natoms; ++nc) {
    Runit[nc] = new int[3];
    for (i=0; i<3; ++i) Runit[nc][i] = work[nc].R[i];
  }
    
  delete[] work;
  return ERROR;
}


#endif
