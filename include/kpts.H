#ifndef __KPTS_H
#define __KPTS_H
#include "io.H"
#include "matrix.H"

// We include *both* versions of the very popular "fold_down", 
// and a routine for reading kpt meshes, both with and without weights.
double fold_down(double kpt_orig[3], double cart_b[9], double kpt[3]);
inline void fold_down(double kpt[3], double cart_b[9]) 
{
  fold_down (kpt, cart_b, kpt);
}

void free_kpts(int Nk, double** &kpt)
{ 
  int i;
  if (kpt != NULL) {
    for (i=0; i<Nk; ++i) delete[] kpt[i];
    delete[] kpt;
  }
}

// MAGN vs. NOMAGN: decides whether kpt[n][3] is the magnitude or not
// CART vs. LATT: determines how we return the kpts read in
// WEIGHTS: whether to expect weights or not
const int READ_KPTS_NOMAGN = 0;
const int READ_KPTS_MAGN = 1;
const int READ_KPTS_CART = 0;
const int READ_KPTS_LATT = 2;
const int READ_KPTS_NOWEIGHTS = 0;
const int READ_KPTS_WEIGHTS = 4;

int read_kpts (FILE* infile, double cart_rlv[9], int &Nkpt, 
	       double** &kpt, double* &w, const int READTYPE) 
{
  int i, n;
  int MAGN, WEIGHTS, RETURN_LATT;
  double k0, tempvec[3]; 
  double cell_vol, cart_rlv_inv[9];
  char type_ident;
  int LATT; // how they're specified in the *file*

  cell_vol = inverse(cart_rlv, cart_rlv_inv);
  for (i=0; i<9; ++i) cart_rlv_inv[i] *= 1./cell_vol;
  
  char dump[512];
  int ERROR = 0;

  // pull out *exactly* what we're supposed to do...
  MAGN = READTYPE & READ_KPTS_MAGN;
  RETURN_LATT = READTYPE & READ_KPTS_LATT;
  WEIGHTS = READTYPE & READ_KPTS_WEIGHTS;

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %lf %c", &Nkpt, &k0, &type_ident);
  
  if (Nkpt < 1) {
    fprintf(stderr, "No points listed in read_kpts\n");
    ERROR = ERROR_BADFILE;
  }
  if (k0 <= 0.) {
    fprintf(stderr, "Bad scale factor in read_kpts\n");
    ERROR = ERROR_BADFILE;
  }
  type_ident = tolower(type_ident);
  if ( (type_ident != 'l') && (type_ident != 'c') ) {
    fprintf(stderr, "Invalid type character (%c) in read_kpts; should be C or L\n", 
            type_ident);
    ERROR = ERROR_BADFILE;
  }
  else {
    LATT = (type_ident == 'l'); // convert from rlv to cartesian
  }
  if (ERROR) {
    myclose(infile);
    exit(ERROR);
  }
  else {
    // Allocate!
    //    if (kpt != NULL) delete[] kpt;
    kpt = new double *[Nkpt];
    if (kpt == NULL) ERROR = ERROR_MEMORY;
    if (!ERROR) {
      if (MAGN)
	for (i=0; i<Nkpt; ++i) kpt[i] = new double[4];
      else 
	for (i=0; i<Nkpt; ++i) kpt[i] = new double[3];
      if (kpt[Nkpt-1] == NULL) ERROR = ERROR_MEMORY;
      if (WEIGHTS) {
	//      if (w != NULL) delete[] w;
	w = new double[Nkpt];
	if (w == NULL) ERROR = ERROR_MEMORY;
      }
    }
  }
  
  // Now, do the readin'!
  for (n=0; (!ERROR) && (n<Nkpt) && (!feof(infile)); ++n) {
    double* kpt_p = kpt[n];
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;

    if (WEIGHTS)
      sscanf(dump, "%lf %lf %lf %lf", &(kpt_p[0]), &(kpt_p[1]), &(kpt_p[2]),
	     &(w[n]));
    else
      sscanf(dump, "%lf %lf %lf", &(kpt_p[0]), &(kpt_p[1]), &(kpt_p[2]));
    for (i=0; i<3; ++i) kpt_p[i] *= k0; // scale!
    // Convert to lattice coord
    if (! LATT) {
      mult_vect(cart_rlv_inv, kpt_p, tempvec);
    }
    // Fold into BZ:
    if (MAGN)
      kpt_p[3] = fold_down(kpt_p, cart_rlv, kpt_p);
    else
      fold_down(kpt_p, cart_rlv);
    if (!RETURN_LATT) {
      // Translate back into cartesian coordinates
      mult_vect(cart_rlv, kpt_p, tempvec);
      for (i=0; i<3; ++i) kpt_p[i] = tempvec[i];
    }
  }
  // Make sure we read enough points...
  if (n != Nkpt) {
    fprintf(stderr, "Not enough points in read_kpt\n");
    ERROR = ERROR_BADFILE;
  }
  return ERROR;
}
  
  
// not in place, and returns the magnitude...
double fold_down(double kpt_orig[3], double cart_b[9], double kpt[3]) 
{
  double bb[9];
  int a[3], a0[3], a1[3];
  double q[3], k0[3];
  double kmagn, qmagn;
  int i;
  int better;

  square(cart_b, bb);
  
  for (i=0; i<3; ++i) {
    kpt[i] = kpt_orig[i] - (int)(kpt_orig[i]);
    if (dcomp(kpt[i],-1.0)) kpt[i] += 1.0;
    if (dcomp(kpt[i],1.0)) kpt[i] -= 1.0;
    // Now, guess our possible shifts
    a0[i] = -1; a1[i] = 1;
    if (! dcomp(kpt[i], 0.0)) {
      if (kpt[i] > 0.0) a1[i] = 0;
      if (kpt[i] < 0.0) a0[i] = 0;
    }
    k0[i] = kpt[i];
  }
  // Magnitudes:
  kmagn = magnsq(bb, k0); // k0 -- our best guess.
  for (a[0]=a0[0]; a[0]<=a1[0]; ++(a[0]))
    for (a[1]=a0[1]; a[1]<=a1[1]; ++(a[1]))
      for (a[2]=a0[2]; a[2]<=a1[2]; ++(a[2])) {
        // Loop through possible shifts:
        for (i=0; i<3; ++i) q[i] = kpt[i] + a[i];
        qmagn = magnsq(bb, q);
        better = (qmagn < kmagn); // Better k-point?
        if (dcomp(qmagn, kmagn)) {
          // Complicated tests...
          better = (q[0] > k0[0]);
          if (dcomp(q[0], k0[0])) {
            better = (q[1] > k0[1]);
            if (dcomp(q[1], k0[1]))
              better = (q[2] > k0[2]);
          }
        }
        if (better) {
          kmagn = qmagn;
          for (i=0; i<3; ++i) k0[i] = q[i];
        }
      }
  // Now we've got our "best" possible k-point:
  for (i=0; i<3; ++i) kpt[i] = k0[i];
  return sqrt(kmagn);
}

#endif
