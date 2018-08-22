#ifndef __SPHERE_HARM__
#define __SPHERE_HARM__

/*
  Header: sphere-harm.H --routines to calculate values for a given spherical
          harmonic expansion
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include "dcomp.H"
#include "io.H" // just to get error codes

//============================ allocation routines =====================

void free_Ylm (int Lmax, double*** &RYlm, double*** &IYlm) 
{
  int l, m;
  if (RYlm != NULL) {
    for (l=0; l<=Lmax; ++l) 
      if (RYlm[l] != NULL) {
	for (m=0; m<=(4*l); ++m) 
	  if (RYlm[l][m] != NULL) delete[] RYlm[l][m];
	delete[] RYlm[l];
      }
    delete[] RYlm;
    RYlm = NULL;
  }
  if (IYlm != NULL) {
    for (l=0; l<=Lmax; ++l) 
      if (IYlm[l] != NULL) { 
	for (m=0; m<=(4*l); ++m) 
	  if (IYlm[l][m] != NULL) delete[] IYlm[l][m];
	delete[] IYlm[l];
      }
    delete[] IYlm;
    IYlm = NULL;
  }
}

// Now with ERROR checking!
int alloc_Ylm (int Lmax, int Ndim, double*** &RYlm, double*** &IYlm) 
{
  int ERROR = 0;
  int l, m, d;

  RYlm = new double **[Lmax+1];
  IYlm = new double **[Lmax+1];
  ERROR = (IYlm == NULL);
  for (l=0; (l<=Lmax) && (!ERROR); ++l) {
    RYlm[l] = new double *[4*l+1];
    IYlm[l] = new double *[4*l+1];
    ERROR = (IYlm[l] == NULL);
    if (ERROR) continue;
    for (m=0; (m<=(4*l)) && (!ERROR); ++m) {
      RYlm[l][m] = new double[Ndim];
      IYlm[l][m] = new double[Ndim];
      ERROR = (IYlm[l][m] == NULL);
      if (ERROR) continue;
      for (d=0; d<Ndim; ++d) {
	// Initialize
	RYlm[l][m][d] = 0.;
	IYlm[l][m][d] = 0.;
      }
    }
  }
  // if we threw an error, just undo all allocation:
  if (ERROR) {
    free_Ylm(Lmax, RYlm, IYlm);
    ERROR = ERROR_MEMORY;
  }
  return ERROR;
}



//============================= I/O routines ===========================

int read_Ylm_even (FILE* infile, int &Lmax, int &Ndim,
		   double*** &RYlm, double*** &IYlm) 
{
  int ERROR = 0;
  int d, i, j, k, l, m, n;
  char dump[512];

    fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d", &Lmax, &Ndim);
  if ( (Lmax < 0) || (Ndim < 1) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value\n", 
            Lmax, Ndim);
    ERROR = ERROR_BADFILE;
  }
  else {
    // Turn Lmax into Lmax/2:
    Lmax = Lmax/2;
    // Allocate: (handles deallocating on error)
    ERROR = alloc_Ylm(Lmax, Ndim, RYlm, IYlm);
  }

  int no_more;
  int l_index, m_index, l_max_0;
  l_max_0 = 0; // Find the *true* Lmax value...
  for (no_more=0, n=0; (!ERROR) && (!no_more) && (!feof(infile)); ++n) {
    // Read the real elements, and ignore # comments:
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    sscanf(dump, "%d %d", &l, &m);
    if ( l < 0 ) {
      no_more = -1; // Our end of data marker...
      break;
    }
    // Check against:
    // a. l > 2*Lmax
    // b. |m| > l
    // c. l odd
    if ( (l > 2*Lmax) || (abs(m) > l) || ( (l%2) == 1 ) ) {
      fprintf(stderr, "Entry %d has bad lm components = %d %d\n", n+1, l, m);
      ERROR = ERROR_BADFILE;
      break;
    }
    if (l_max_0 < (l/2)) l_max_0 = (l/2);
    m_index = l+m;
    l_index = l/2; // l/2 :)
    
    // Code to pull off each elem one by one...
    char *startp, *endp;
    startp = dump; // Start at the beginning...
    // Need to "prime the pump" by going past first two entries.
    for (k=0; k<2; ++k) {
      strtol(startp, &endp, 10); // use l because we've got two integer elems.
      startp = endp;
    }
    for (d=0; d<Ndim; ++d) {
      RYlm[l_index][m_index][d] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough entries?
    if (d != Ndim) { 
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Didn't read enough entries on real entry %d.\n", n+1);
    }
    
    // Now, repeat with imag piece:
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d", &i, &j);
    if ( (i != l) || (j != m) ) {
      fprintf(stderr,
              "Entry %d's imag. lm's don't match real lm's = %d %d vs. %d %d\n",
              n+1, l, m, i, j);
      ERROR = ERROR_BADFILE;
      break;
    }
    
    // Code to pull off each elem one by one...
    startp = dump; // Start at the beginning...
    // Need to "prime the pump" by going past first two entries.
    for (k=0; k<2; ++k) {
      strtol(startp, &endp, 10); // use l because we've got two integer elems.
      startp = endp;
    }
    for (d=0; d<Ndim; ++d) {
      IYlm[l_index][m_index][d] = strtod(startp, &endp);
      // check for conversion:
      if (startp == endp) break;
      startp = endp;
    }
    // TEST: Did we read enough entries?
    if (d != Ndim) { 
      ERROR = ERROR_BADFILE;
      fprintf(stderr, "Didn't read enough entries on imag entry %d.\n", n+1);
    }
  }

  if (ERROR) return ERROR; // time to bail out yet?

  // Now, remove the "extra" l elements, if l_max_0 < Lmax.
  if (l_max_0 < Lmax) {
    // Free up the space:
    for (l=(l_max_0+1); l < Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) { // Since l is half it's value...
        delete[] RYlm[l][m];
        delete[] IYlm[l][m];
      }
      delete[] RYlm[l];
      delete[] IYlm[l];
    }
    Lmax = l_max_0; // Correct Lmax value
  }

  return ERROR;
}


void verbose_outYlm (int Lmax, int Ndim, double*** &RYlm, double*** &IYlm,
		     const char* leadin) 
{
  int l, m, d;
  
  // Print out *all* the elements we read in...
  printf("%sYlm components:\n", leadin);
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(4*l); ++m) {
      // Real part
      printf("%s%3d %3d R ", leadin, 2*l, m-2*l);
      for (d=0; d<Ndim; ++d)
	printf(" %.12le", RYlm[l][m][d]);
      printf("\n");
      // Imag part
      printf("%s%3d %3d I ", leadin, 2*l, m-2*l);
      for (d=0; d<Ndim; ++d)
	printf(" %.12le", IYlm[l][m][d]);
      printf("\n");
    }
  }

  // print out information about *size* of each component.
  printf("%s Magnitude of terms:\n", leadin);
  double zmagn;
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(4*l); ++m) {
      zmagn = 0.;
      for (d=0; d<Ndim; ++d)
	zmagn += RYlm[l][m][d]*RYlm[l][m][d] + 
	  IYlm[l][m][d]*IYlm[l][m][d];
      zmagn = sqrt(zmagn /(double)Ndim);
      printf("%s %3d %3d %.8le\n", leadin, 2*l, m-2*l, zmagn);
    }
  }
}



//========================== Evaluation routines =======================

// Our workspace
typedef struct 
{
  double **RsphP, **IsphP;
  double *result_temp;
  int Lmax;
} SH_workspace_type;


SH_workspace_type * allocate_SH_workspace (int Lmax) 
{
  if (Lmax < 0) return NULL;

  int l;
  SH_workspace_type * w = new SH_workspace_type;

  w->Lmax = Lmax;
  w->RsphP = new double *[Lmax+1];
  w->IsphP = new double *[Lmax+1];
  for (l=0; l<=Lmax; ++l) {
    w->RsphP[l] = new double[2*l+1]; // Only calc for m=0..2*l
    w->IsphP[l] = new double[2*l+1];
  }
  w->result_temp = new double[2*Lmax+1]; // It calcs from m=0..2*Lmax

  return w;
}

void free_SH_workspace (SH_workspace_type* &w) 
{
  int l;
  
  if (w != NULL) {
    delete[] w->result_temp;
    for (l=0; l<=w->Lmax; ++l) {
      delete[] w->RsphP[l];
      delete[] w->IsphP[l];
    }
    delete[] w->RsphP;
    delete[] w->IsphP;
    delete w;
    w = NULL;
  }
}

inline void calc_angle(double lmn[3], double &cost, double &phi) 
{
  cost = lmn[2];
  if (cost > 1) cost = 1.;
  if (cost < -1) cost = -1.;
  phi = atan2(lmn[1], lmn[0]);
}

// Evaluates Ylm expansion at x; Lmax is the maximum element in the
// Ylm expansion, and Ndim is the dimension of the vector
// Note: *just* outputs the real part
// RYlm[l/2 = 0..Lmax][m+2l = 0..4l][d=0..Ndim-1]
int eval_Ylm_expansion_R (int Lmax, double x[3], int Ndim, 
			  double ***RYlm, double ***IYlm,
			  double *Rval, SH_workspace_type *w) 
{
  int i, l, m, d;
  int m_index;
  // Name our workspace
  double *result_temp = w->result_temp;
  double **RsphP = w->RsphP, **IsphP = w->IsphP;
 
  // Check that we've got adequate workspace
  if (Lmax > w->Lmax) return -1;

  // Check for the case of x = 0
  double xmagn, xv[3];
  xmagn = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  if (dcomp(xmagn, 0)) {
    // Only Y(0,0) matters...
    const double Y00 = M_2_SQRTPI * 0.25;  // Y00 = 1/sqrt(4Pi)
    for (d=0; d<Ndim; ++d)
      Rval[d] = Y00 * RYlm[0][0][d];
    return 0;
  }
  else {
    xmagn = 1./xmagn;
    for (i=0; i<3; ++i) xv[i] = x[i]*xmagn;
  }
  
  double phi, cost; // cos(theta) for kpt
  calc_angle(xv, cost, phi);

  // First, calc the spherical harmonic values --
  for (m=0; m<=2*Lmax; ++m) {
    double cosm, sinm;
    // exp(I*m*phi):
    sinm = sin(m*phi);
    cosm = cos(m*phi);
    // Plm(cost):
    gsl_sf_legendre_sphPlm_array(2*Lmax, m, cost, result_temp);
    // Now, put in the appropriate spots:
    for (l=2*Lmax; l>=m; l -= 2) {
      // result_temp goes from m to Lmax
      RsphP[l/2][m] = cosm*result_temp[l-m];
      IsphP[l/2][m] = sinm*result_temp[l-m];
    }
  }
  // Now, do our sum:
  int mplus, mneg;
  for (d=0; d<Ndim; ++d) Rval[d] = 0;
  for (l=0; l<=Lmax; ++l) {
    // m = 0:
    // note: xYlm[l][m_index], m_index = 2*l+m, so for m=0, we use 2*l.
    m_index = 2*l;
    for (d=0; d<Ndim; ++d) 
      Rval[d] += RYlm[l][m_index][d]*RsphP[l][0];
    // Next, rest of m terms:
    for (m=1; m<=2*l; ++m) {
      mplus = m+2*l;
      mneg = -m+2*l;
      for (d=0; d<Ndim; ++d)
	Rval[d] += (RYlm[l][mplus][d] + RYlm[l][mneg][d])*RsphP[l][m]
	  - (IYlm[l][mplus][d] - IYlm[l][mneg][d])*IsphP[l][m];
    }
  }
  return 0;
}
#endif
