#ifndef __POINTGROUP_H
#define __POINTGROUP_H

/*
  Header file for point group generation.
*/

const int MAXop = 48;


// Returns the number of operations; returns 0 if there was an error.
int gen_pointgroup (double cart[9], double gcart[MAXop][9], 
		    int gunit[MAXop][9], int inv_index[MAXop]) 
{
  int Nop;
  int d, i, n;
  
  int maxn = 1, MORE;
  int gtry[9];
  double gtrycart[9];
  const double ident[9] = {1,0,0, 0,1,0, 0,0,1};

  for (i=0; i<9; ++i) gtry[i] = -maxn;
  MORE = 1;
  Nop = 0;
  double cartinv[9], cellvol;
  cellvol = inverse(cart, cartinv);
  mult(cartinv, 1./cellvol, cartinv);
  while (MORE) {
    if (abs(det(gtry)) == 1) {
      double temp[9], gsq[9];
      // cartesian component...
      mult(cart, gtry, temp);
      mult(temp, cartinv, gtrycart);
      // square it...
      square(gtrycart, gsq);
      // Compare to the identity matrix:
      if (equal(gsq, ident)) {
        // We got a live one!
        for (d=0; d<9; ++d) {
          gunit[Nop][d] = gtry[d];
          gcart[Nop][d] = gtrycart[d];
          inv_index[Nop] = -1;
        }
        ++Nop;
      }
    }
    // increment counters
    for (i=8, ++(gtry[i]); (i>=0) && (gtry[i] > maxn); --i) {
      gtry[i] = -maxn;
      if (i>0) ++(gtry[i-1]);
    }
    MORE = ( (i>=0) && (Nop < MAXop) );
  }
  // Now, calc all the inverses:
  for (n=0; n<Nop; ++n) {
    int ginv[9], gsign;
    if (inv_index[n] < 0) {
      gsign = inverse(gunit[n], ginv);
      mult(ginv, gsign, ginv);
      for (i=n; !(equal(ginv, gunit[i])); ++i) ;
      inv_index[n] = i;
      inv_index[i] = n;
    }
  }
  
  return Nop;
}

/* analyzes a point-group operation in terms of the Euler angles;
   returns 0 if a proper rotation, or 1 if an improper rotation.
   Return a -1 if an error (as in, doesn't look like a group op). */
const double inv[9]={-1,0,0, 0,-1,0, 0,0,-1};
int euler_rot(double grot[9],
	      double &alpha, double &beta, double &gamma) 
{
  double g[9], detg;
  int improp;

  detg = det(grot); // has to be +-1 
  if (! dcomp(fabs(detg), 1) ) return -1;

  if (detg < 0) improp = 1;
  else          improp = 0;
  if (improp) mult(inv, grot, g);
  else for (int d=0; d<9; ++d) g[d] = grot[d];
  
  beta = acos(g[8]); // acos (g_zz)
  double sinbeta = sqrt(1-g[8]*g[8]);
  if (! zero(sinbeta) ) {
    // atan2(sin(), cos()):
    alpha = atan2(g[7], g[6]);  // sin(a)=g_zy/sin(b), cos(a)=g_zx/sin(b)
    gamma = atan2(g[5], -g[2]); // sin(g)=g_yz/sin(b), cos(g)=-g_xz/sin(b)
  }
  else {
    // treat the special case: (rotation around z)
    alpha = atan2(-g[3], g[4]); // sin(a)=-g_yx, cos(a)=g_yy
    gamma = 0;                  // sin(g)=0,     cos(g)=1
  }
  return improp;
}


#endif
