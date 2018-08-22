#ifndef __PLANE_H
#define __PLANE_H

/*
  Program: plane.H
  Author:  D. Trinkle
  Date:    2012 June 20
  Purpose: Code for the planar coordinate system.

*/

#include "dcomp.H"
#include "matrix.H"
#include "io.H"
#include "dislocation.H"

//****************************** SUBROUTINES ****************************

// read the triad of planar coordinates, in unit cell coord.
inline int read_plane (FILE* infile, int t1_unit[3], int t2_unit[3], int n_unit[3]) {
  return read_disl(infile, t1_unit, t2_unit, n_unit);
}

int make_planecoord (double cart[9], 
		     double t1v[3], double t2v[3], double nv[3],
		     int t1u[3], int t2u[3], int nu[3], int TYPE);


// Makes planar coordinate system (normalizes, orthogonalizes, etc.)
// TYPE determines whether it read in unit cell or cart coords.  It
// overwrites all vectors accordingly.
// :: modified from make_dislcoord
int make_planecoord (double cart[9], 
		     double t1v[3], double t2v[3], double nv[3],
		     int t1u[3], int t2u[3], int nu[3], int TYPE) 
{
  int d, ERROR=0;
  int MAKECART = ( (TYPE & DISLCOORD_UNIT) != 0);
  int KEEPUNIT = ( (TYPE & DISLCOORD_NOCHANGE) != 0);
  double v_magn;
  
  if (MAKECART) {
    mult_vect(cart, t1u, t1v);
    mult_vect(cart, t2u, t2v);
    mult_vect(cart, nu, nv);
  }
  // Normalize:
  v_magn = sqrt(dot(t1v, t1v));
  for (d=0; d<3; ++d) t1v[d] /= v_magn;
  v_magn = sqrt(dot(t2v, t2v));
  for (d=0; d<3; ++d) t2v[d] /= v_magn;
  if (zero(fabs(dot(t1v,t2v))-1)) {
    fprintf(stderr, "t1 and t2 are parallel...?  t1_unit=(%d,%d,%d)  t2_unit=(%d,%d,%d)\n",
            t1u[0],t1u[1],t1u[2], t2u[0],t2u[1],t2u[2]);
    return(1); // get out now.
  }
  // make n as the cross product of t1 and t2:
  crossprod(t1v, t2v, nv);
  
  // finally, put in the new unit cell coord:
  if (!KEEPUNIT) {
    // We do not make_unit on t1 and t2--because their lengths MATTER.
    //    make_unit(cart, t1v, t1u);
    //    make_unit(cart, t2v, t2u);
    make_unit(cart, nv, nu);
  }
  return ERROR;
}

#endif
