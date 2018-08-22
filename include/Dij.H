#ifndef __DIJ_H
#define __DIJ_H

#include <string.h>
#include "io.H"
#include "matrix.H"

typedef struct 
{
  int Runit[3];
  double Rcart[3];
  double Rmagn;
  double mat[9];
} point_type;

void copy_point (const point_type &src, point_type &dest) 
{
  for (int d=0; d<3; ++d) dest.Runit[d] = src.Runit[d];
  for (int d=0; d<3; ++d) dest.Rcart[d] = src.Rcart[d];
  dest.Rmagn = src.Rmagn;
  for (int d=0; d<9; ++d) dest.mat[d] = src.mat[d];
}
  

void write_Dij (FILE* outfile, int Np, point_type* p, char *firstline) 
{
  int n;
  if (firstline == NULL) printf("%d 9 # number of points\n", Np);
  else                   printf("%s\n", firstline);
  for (n=0; n<Np; ++n) {
    point_type *tp = p+n;
    fprintf(outfile, "%3d %3d %3d %.12le %.12le %.12le %.12le %.12le %.12le %.12le %.12le %.12le\n",
	   tp->Runit[0], tp->Runit[1], tp->Runit[2], 
	   tp->mat[0], tp->mat[1], tp->mat[2],
	   tp->mat[3], tp->mat[4], tp->mat[5],
	   tp->mat[6], tp->mat[7], tp->mat[8]);
  }
}

inline void write_Dij (FILE* outfile, int Np, point_type* p) 
{
  write_Dij(outfile, Np, p, NULL);
}


// determines whether we read in the matrix values or not.
const int READ_DIJ_NOMAT = 0;
const int READ_DIJ_MAT = 1;

// general one:
int read_Dij (FILE* infile, double cart[9], int& Np, point_type* &p,
	      int READ_TYPE, char *firstline);

// no first line:
inline int read_Dij (FILE* infile, double cart[9], int& Np, point_type* &p,
	      int READ_TYPE) 
{ char *null = NULL;
  return read_Dij(infile, cart, Np, p, READ_TYPE, null);}

// always read the matrix:
inline int read_Dij (FILE* infile, double cart[9], int& Np, point_type* &p,
		     char *firstline) 
{ return read_Dij(infile, cart, Np, p, READ_DIJ_MAT, firstline);}  

// always read the matrix, no first line:
inline int read_Dij (FILE* infile, double cart[9], int& Np, point_type* &p) 
{ char *null = NULL;
  return read_Dij(infile, cart, Np, p, READ_DIJ_MAT, null);}  


int read_Dij (FILE* infile, double cart[9], int& Np, point_type* &p,
	      int READ_TYPE, char *firstline) 
{
  int n;
  int ERROR = 0;
  char dump[512];
  double cart2[9];
  int DIJ = READ_TYPE & READ_DIJ_MAT;

  square(cart, cart2);

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Np);
  
  if (firstline != NULL)
    strcpy(firstline, dump);

  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in read_Dij\n", Np);
    ERROR = ERROR_BADFILE;
  }
  else {
    //    if (p != NULL) delete[] p;
    p = new point_type[Np];
    if (p == NULL) ERROR = ERROR_MEMORY;
  }
  
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    point_type *tp = p + n; // this point
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    if (DIJ)
      sscanf(dump, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	     &(tp->Runit[0]), &(tp->Runit[1]), &(tp->Runit[2]),
	     &(tp->mat[0]), &(tp->mat[1]), &(tp->mat[2]),
	     &(tp->mat[3]), &(tp->mat[4]), &(tp->mat[5]),
	     &(tp->mat[6]), &(tp->mat[7]), &(tp->mat[8]));
    else
      sscanf(dump, "%d %d %d", 
	     &(tp->Runit[0]), &(tp->Runit[1]), &(tp->Runit[2]));
    // convert to cartesian coordinates:
    mult_vect(cart, tp->Runit, tp->Rcart);
    tp->Rmagn = sqrt(magnsq(cart2, tp->Runit));
  }
  if (n!=Np) {
    fprintf(stderr, "Didn't read enough points in read_Dij\n");
    ERROR = ERROR_BADFILE;
  }
  return ERROR;
}


#endif
