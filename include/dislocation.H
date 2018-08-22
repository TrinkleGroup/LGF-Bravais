#ifndef __DISLOCATION_H
#define __DISLOCATION_H

/*
  Program: dislocation.H
  Author:  D. Trinkle
  Date:    Sept. 9, 2004
  Purpose: Code for the dislocation coordinate system, and handling the
           fourier series stuff.

*/

#include "dcomp.H"
#include "matrix.H"
#include "io.H"

//****************************** SUBROUTINES ****************************

// read the triad of dislocation coordinates, in unit cell coord.
int read_disl (FILE* infile, int t_unit[3], int m_unit[3], int n_unit[3]);

// Free our 2degf
void free_2degf(int N, double** &gncos, double** &gnsin);

// Read in our 2degf from a file (needs only even components)
int read_2degf (FILE* infile, int& N, double** &gncos, double** &gnsin);

// Takes in a vector and finds the smallest unit cell vector pointing
// in that direction
void make_unit (double cart[9], double v[3], int u[3]);

// Makes dislocation coordinate system (normalizes, orthogonalizes, etc.)
// TYPE determines whether it read in unit cell or cart coords.  It
// overwrites the cartesian vectors accordingly, and will change the unit
// cell coord if DISLCOORD_CHANGE is set.
const int DISLCOORD_CART = 0;
const int DISLCOORD_UNIT = 1;
const int DISLCOORD_CHANGE = 0;
const int DISLCOORD_NOCHANGE = 2;

int make_dislcoord (double cart[9], 
		    double tv[3], double mv[3], double nv[3],
		    int tu[3], int mu[3], int nu[3], int TYPE);


int read_disl (FILE* infile, int t_unit[3], int m_unit[3], int n_unit[3]) 
{
  int ERROR = 0;
  char dump[512];
  
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", t_unit, t_unit+1, t_unit+2);
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", m_unit, m_unit+1, m_unit+2);
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", n_unit, n_unit+1, n_unit+2);

  if ( (t_unit[0]==0) && (t_unit[1]==0) && (t_unit[2]==0) ) {
    fprintf(stderr, "t vector is zero\n");
    ERROR = ERROR_BADFILE;
  }
  if ( (m_unit[0]==0) && (m_unit[1]==0) && (m_unit[2]==0) ) {
    fprintf(stderr, "m vector is zero\n");
    ERROR = ERROR_BADFILE;
  }
  if ( (n_unit[0]==0) && (n_unit[1]==0) && (n_unit[2]==0) ) {
    fprintf(stderr, "n vector is zero\n");
    ERROR = ERROR_BADFILE;
  }
  return ERROR;
}

void free_2degf(int N, double** &gncos, double** &gnsin) 
{
  for (int i=0; i<=N; ++i) {
    if (gncos[i] != NULL) delete[] gncos[i];
    if (gnsin[i] != NULL) delete[] gnsin[i];
  }
  delete[] gncos;
  delete[] gnsin;
}

// Now with safe allocations:
const int ALLOC_LIST = 0;
const int ALLOC_ALL = 1;

int alloc_2degf(int N, double** &gncos, double** &gnsin, int TYPE) 
{
  int DOALL = (TYPE == ALLOC_ALL);
  if ( (gncos = new double*[N+1]) == NULL ) return ERROR_MEMORY;
  if ( (gnsin = new double*[N+1]) == NULL ) return ERROR_MEMORY;
  if (! DOALL)
    for (int n=0; n<=N; ++n) {
      gncos[n] = NULL;
      gnsin[n] = NULL;
    }
  else
    for (int n=0; n<=N; ++n) {
      if ( (gncos[n] = new double[9]) == NULL ) return ERROR_MEMORY;
      if ( (gnsin[n] = new double[9]) == NULL ) return ERROR_MEMORY;
    }
  return 0;
}

// default value:
inline int alloc_2degf(int N, double** &gncos, double** &gnsin) 
{ return alloc_2degf(N, gncos, gnsin, ALLOC_ALL); }


// Read in our egf from a file:
int read_2degf (FILE* infile, int& N, double** &gncos, double** &gnsin,
		char *firstline);

inline int read_2degf (FILE* infile, int& N, double** &gncos, double** &gnsin) 
{ 
  char *null = NULL;
  return read_2degf(infile, N, gncos, gnsin, null);
}


int read_2degf (FILE* infile, int& N, double** &gncos, double** &gnsin,
		char *firstline) 
{
  int ERROR = 0;
  int i, j, d;
  char dump[512];
  
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d", &N);

  if (firstline != NULL)
    strcpy(firstline, dump);
  
  N /= 2; // should be even anyway...

  ERROR = alloc_2degf(N, gncos, gnsin, ALLOC_LIST); // don't allocate entries
  if (ERROR != 0) {
    fprintf(stderr, "Error allocating in read_2degf.\n");
    return(ERROR);
  }
  while ( (!feof(infile)) && (!ERROR) ) {
    double gin_cos[9], gin_sin[9];
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &i,
           gin_cos, gin_cos+1, gin_cos+2,
           gin_cos+3, gin_cos+4, gin_cos+5,
           gin_cos+6, gin_cos+7, gin_cos+8);
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &j,
           gin_sin, gin_sin+1, gin_sin+2,
           gin_sin+3, gin_sin+4, gin_sin+5,
           gin_sin+6, gin_sin+7, gin_sin+8);
    if ( (i%2)!= 0 ) {
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Should only have even FT components.\n");
    }
    if ( i != j ) {
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Read a line where real and imag components do not have correct index?\n");
    }
    if (i > (2*N)) {
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Index (%d) too large (>%d)\n", i, 2*N);
    }
    if (!ERROR) {
      i /= 2;
      // If i < 0, the gin_sin components have the right sign;
      // else, we negate.
      if (i < 0)
        i = -i;
      else
        for (d=0; d<9; ++d) gin_sin[d] = -gin_sin[d];
      if (gncos[i] == NULL) {
        // allocate and initialize
        gncos[i] = new double[9]; for (d=0; d<9; ++d) gncos[i][d]=0;
        gnsin[i] = new double[9]; for (d=0; d<9; ++d) gnsin[i][d]=0;
      }
      for (d=0; d<9; ++d) {
        gncos[i][d] += gin_cos[d];
        gnsin[i][d] += gin_sin[d];
      }
    }
  }
  
  if (!ERROR) {
    // find correct N value:
    for ( ; (N>0) && (gncos[N] == NULL); --N) ;
    for (i=0; i<=N; ++i)
      if (gncos[i] == NULL) {
        // allocate and initialize
        gncos[i] = new double[9]; for (d=0; d<9; ++d) gncos[i][d]=0;
        gnsin[i] = new double[9]; for (d=0; d<9; ++d) gnsin[i][d]=0;
      }
  }  
  else {
    // garbage collect whatever we may have allocated
    free_2degf(N, gncos, gnsin);
  }
  return ERROR;
}



inline void eval_cossin_expansion(int Nmax, double theta,
				  double **gncos, double **gnsin,
				  double g[9]) 
{
  int d;
  for (d=0; d<9; ++d) g[d] = 0;
  for (int n=0; n<=Nmax; ++n) {
    double cost = cos(2*n*theta);
    double sint = sin(2*n*theta);
    for (d=0; d<9; ++d) g[d] += gncos[n][d]*cost + gnsin[n][d]*sint;
  }
}


// Makes dislocation coordinate system (normalizes, orthogonalizes, etc.)
// TYPE determines whether it read in unit cell or cart coords.  It
// overwrites all vectors accordingly.
int make_dislcoord (double cart[9], 
		    double tv[3], double mv[3], double nv[3],
		    int tu[3], int mu[3], int nu[3], int TYPE) 
{
  int d, ERROR=0;
  int MAKECART = ( (TYPE & DISLCOORD_UNIT) != 0);
  int KEEPUNIT = ( (TYPE & DISLCOORD_NOCHANGE) != 0);
  double v_magn;
  
  if (MAKECART) {
    mult_vect(cart, tu, tv);
    mult_vect(cart, mu, mv);
    mult_vect(cart, nu, nv);
  }
  // Normalize:
  v_magn = sqrt(dot(tv, tv));
  for (d=0; d<3; ++d) tv[d] /= v_magn;
  // make sure m is normal to t:
  v_magn = dot(mv, tv);
  for (d=0; d<3; ++d) mv[d] -= v_magn * tv[d];
  // normalize
  v_magn = sqrt(dot(mv, mv));
  if (zero(v_magn)) {
    fprintf(stderr, "t and m are parallel...?  t_unit=(%d,%d,%d)  m_unit=(%d,%d,%d)\n",
            tu[0],tu[1],tu[2], mu[0],mu[1],mu[2]);
    return(1); // get out now.
  }
  for (d=0; d<3; ++d) mv[d] /= v_magn;
  // make n as the cross product of t and m:
  crossprod(tv, mv, nv);
  
  // finally, put in the new unit cell coord:
  if (!KEEPUNIT) {
    // We no longer run make_unit on t--because it's length MATTERS.
    //    make_unit(cart, tv, tu);
    make_unit(cart, mv, mu);
    make_unit(cart, nv, nu);
  }
  return ERROR;
}



inline double min_magn (double a[3]) 
{
  double min=1e30;
  for (int i=0; i<3; ++i) 
    if (! zero(a[i]))
      if (fabs(a[i]) < min) min = fabs(a[i]);
  return min;
}

const double DIV_TOL = 1e-3;

double divisor (double a[3]) 
{
  int m, i, fail;
  double am, am_i;

  m=0; fail=1;
  while (fail) {
    ++m;
    fail=0;
    for (i=0; i<3; ++i) {
      am = a[i] * m;
      am_i = round(am);
      if ( fabs(am-am_i) > DIV_TOL) fail = 1;
    }
  }
  return (double)m;
}

void make_unit (double cart[9], double v[3], int u[3]) 
{
  double cart_inv[9], scale;
  double uinit[3];
  int i;
  
  careful_inverse(cart, cart_inv);
  mult_vect(cart_inv, v, uinit);
  
  scale = min_magn(uinit);
  for (i=0; i<3; ++i) uinit[i] *= 1./scale;
  scale = divisor(uinit);
  for (i=0; i<3; ++i) u[i] = lround(uinit[i]*scale);
}


#endif
