/*
  Program: phonon-DOs.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a kpoint set (only to get the weights) and phonons
           evaluated the kpoint set and compute the phonon DOS.

  Param.:  <kpt> <phonon>
	   kpt:    list of kpoints to calculate our phonons (WITH weights)
	   phonon: phonon frequencies in THz for each kpt

	   ==== kpts ====
	   N k0 <type>     # Number of points, scale factor, C/L for units
	   k1.1 k1.2 k1.3 w1
	   ...
	   kN.1 kN.2 kN.3 wN
	   ==== kpts ====

	   ==== phonon ====
	   N  # Number of points; must match Nkpts
	   * w1.1 w1.2 w1.3  # ignored column, followed by phonon freq.
	   ...
	   * wN.1 wN.2 wN.3
	   ==== phonon ====
	   

  Flags:   MEMORY:  number of bins used--default = 1/3 Nkpt
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and calculate.

  Output:  We output the points in unit coord.  (Switch for lattice coord?)
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"
#include "kpts.H"
#include "eigen.H"

// We're given our frequencies in THz; we want to convert them to a 
// temperature.  That means multiplying nu by h/kB = 4.79921566e-11 K/s
// Converting that to THz gives:
const double hdivkB = 47.9921566;

// determines sigma = (prefact)*dx
// erf(1/(2 pre)) determines how much of the centered delta function lies
// completely in a single bin.
// erf(1/2) = 0.52 --> pre = 1
// erf(1)   = 0.84 --> pre = 0.5
const double sigma_prefactor = 1;

//***************************** SUBROUTINES ****************************

void add_DOS (double* DOS, const int &Nbins, 
	      const double &dx, const double &sigma, 
	      const double &nu, const double &w);

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<kpt> <phonon>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'k'}; // flag characters.

const char* ARGEXPL = 
"  kpts:   list of kpoints to calculate our phonons (must have weights)\n\
  phonon: set of phonon frequencies calculated using kpts\n\
  -k: convert from linear frequencies nu in THz to temperature in K\n\
  MEMORY value determines number of bins used; default = 1/30 Nkpt";

const char* FILEEXPL =
"==== kpts ====\n\
N k0 <type>     # Number of points, scale factor, C/L for units\n\
k1.1 k1.2 k1.3 w1\n\
...\n\
kN.1 kN.2 kN.3 wN\n\
==== kpts ====\n\
\n\
==== phonon ====\n\
N  # Number of points; must match Nkpts\n\
* w1.1 w1.2 w1.3  # ignored column, followed by phonon freq.\n\
...\n\
* wN.1 wN.2 wN.3\n\
==== phonon ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // number of bins

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (i=0; i<NFLAGS; ++i) flagon[i] = 0;
  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // Flags
  int TEMP = flagon[0]; // do we want to convert to kelvin? (add eV later)

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *kpts_name, *phonon_name;
  
  // Let's pull off the args:
  kpts_name = args[0];
  phonon_name = args[1];

  //++ ==== kpts ====
  int Nkpt;
  double** kpt, *w;
  
  infile = myopenr(kpts_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", kpts_name);
    exit(ERROR_NOFILE);
  }
  // the last entry specifies what we read and return:
  double cart[9] = {1,0,0,0,1,0,0,0,1};
  ERROR = read_kpts(infile, cart, Nkpt, kpt, w, 
                    READ_KPTS_WEIGHTS);
  myclose(infile);
  //-- ==== kpts ====

  if (TESTING) {
    printf("## %d kpoints.\n", Nkpt);
    printf("##  i    k_1      k_2      k_3\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## %2d %8.5lf %8.5lf %8.5lf  %8.5lf\n",
             k+1, kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
    }
  }
  // garbage collect kpt... we didn't need it anyway.
  free_kpts(Nkpt, kpt);

  if (ERROR) exit(ERROR);
  
  //++ ==== phonon ====
  // double nu[Nkpt][3];
  double** nu = new double*[Nkpt];
  for (int n=0; n<Nkpt; ++n) nu[n] = new double[3];

  int imag_freq = 0; // check for imaginary frequencies...
  infile = myopenr(phonon_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", phonon_name);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &j);
  if (j!=Nkpt) {
    ERROR = -1;
    fprintf(stderr, "File %s doesn't match the kpt set?\n", phonon_name);
  }
  else 
    for (k=0; (k<Nkpt) && (!feof(infile)); ++k) {
      do {
	fgets(dump, sizeof(dump), infile);
      } while ((!feof(infile)) && (dump[0] == '#'));
      if (feof(infile)) break;
      
      // %*d : does conversion, suppresses assignment
      sscanf(dump, "%*d %lf %lf %lf", 
	     &(nu[k][0]), &(nu[k][1]), &(nu[k][2]));
      sort3(nu[k]);
      for (d=0; d<3; ++d)
	if (nu[k][d] < -TOLER) imag_freq = -1;
    }
  myclose(infile);
  //-- ==== phonon ====

  if (TESTING) {
    printf("## %d kpoints.\n", Nkpt);
    if (imag_freq)
      printf("## IMAGINARY FREQUENCIES--difficult to do DOS.\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## %2d", k+1);
      printf(" %8.5lf %8.5lf %8.5lf",
	     nu[k][0], nu[k][1], nu[k][2]);
      printf("\n");
    }
  }
  if (ERROR) exit(ERROR);


  // Scale by temperature:
  if (TEMP) {
    for (k=0; k<Nkpt; ++k)
      for (d=0; d<3; ++d)
	nu[k][d] *= hdivkB;
  }
  
  // ***************************** ANALYSIS **************************
  // Determine number of bins, and step size
  int Nbins;
  
  if (MEMORY != 0) Nbins = MEMORY;
  else             Nbins = (Nkpt/30);
  if (Nbins < 10) {
    fprintf(stderr, "Warning: number of bins (%d) is very low (<10)\n", Nbins);
    fprintf(stderr, "Consider rerunning with more kpoints, or setting MEMORY\n");
    Nbins = 10;
  }

  double max_nu, dx, sigma;
  for (max_nu=0, k=0; k<Nkpt; ++k)
    for (d=0; d<3; ++d)
      if (nu[k][d] > max_nu) max_nu = nu[k][d];
  
  dx = max_nu/(Nbins-5);
  sigma = sigma_prefactor * dx;

  double DOS[Nbins+1];
  for (i=0; i<=Nbins; ++i) DOS[i] = 0.;

  // add up the DOS:
  for (k=0; k<Nkpt; ++k) 
    for (d=0; d<3; ++d) 
      add_DOS(DOS, Nbins, dx, sigma, nu[k][d], w[k]);

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.
  // theta's relative errors:
  for (i=0; i<=Nbins; ++i) {
    printf("%13.8lf %.8le\n", (i+0.5)*dx, DOS[i]); // print the CENTER
    //    printf("%13.8lf %.8le\n", i*dx, DOS[i]);
    //    printf("%13.8lf %.8le\n", (i+1)*dx, DOS[i]);
  }
  
  // ************************* GARBAGE COLLECTION ********************
  for (int n=0; n<Nkpt; ++n) nu[n] = new double[3];
  delete[] nu;
  delete[] w;
  return 0;
}

// x1x0 = x1 - x0
inline double delta (const double &x1x0, const double &sigma,
		     const double &dx) 
{
  return 0.5*(erf( (x1x0+dx)/sigma ) - erf( (x1x0)/sigma ) );
}


void add_DOS (double* DOS, const int &Nbins, 
	      const double &dx, const double &sigma, 
	      const double &nu, const double &w) 
{
  int i, ilo, ihi;
  
  ilo = (int) ((nu-5*sigma)/dx);
  ihi = (int) ((nu+5*sigma)/dx);
  if (ilo < 0) ilo = 0;
  if (ihi > Nbins) ihi = Nbins;
  // quick sanity check on ranges:
  if ( (ilo > Nbins) || (ihi < 0) ) return;
  
  for (i=ilo; i<=ihi; ++i) {
    double x1 = i*dx;
    DOS[i] += w*delta(x1-nu, sigma, dx);
  }
}

