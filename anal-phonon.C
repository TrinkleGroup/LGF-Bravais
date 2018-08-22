/*
  Program: anal-phonon.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a kpoint set (only to get the weights) and two different
           phonons evaluated on each and compute:

	   1. High temperature Debye temps Theta_0 and Theta_2 for each
	   2. RMS relative error in the phonons (assuming the first is
	      the "correct" one)

  Param.:  <kpt> <phonon1> <phonon2>
	   kpt:     list of kpoints to calculate our phonons (WITH weights)
	   phononi: phonon frequencies in THz for each kpt

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
	   

  Flags:   MEMORY:  not used
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

//***************************** SUBROUTINES ****************************

inline double square (double x) 
{
  return x*x;
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<kpt> <phonon1> <phonon2>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'o'}; // flag characters.

const char* ARGEXPL = 
"  kpts:    list of kpoints to calculate our phonons (must have weights)\n\
  phonon1: set of phonon frequencies calculated using kpts (the correct ones)\n\
  phonon2: set of phonon frequencies calculated using kpts (calc. rel. error)\n\
  -o: use old RMS error definition";

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
  int MEMORY = 0;   // not used

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
  int BADRMS = flagon[0]; // how do we want to calc. the rms deviation?

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  const int Nphonon = 2;
  char *kpts_name, *phonon_name[Nphonon];
  
  // Let's pull off the args:
  kpts_name = args[0];
  for (i=0; i<Nphonon; ++i) 
    phonon_name[i] = args[i+1];

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
             k+1, kpt[k][0], kpt[k][1], kpt[k][2]);
    }
  }
  // garbage collect kpt... we didn't need it anyway.
  free_kpts(Nkpt, kpt);

  if (ERROR) exit(ERROR);
  
  //++ ==== phonon ====
  double nu[Nphonon][Nkpt][3];
  int imag_freq = 0; // check for imaginary frequencies...
  for (i=0; i<Nphonon; ++i) {
    infile = myopenr(phonon_name[i]);
    if (infile == NULL) {
      fprintf(stderr, "Couldn't open %s for reading.\n", phonon_name[i]);
      exit(ERROR_NOFILE);
    }
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d", &j);
    if (j!=Nkpt) {
      ERROR = -1;
      fprintf(stderr, "File %s doesn't match the kpt set?\n", phonon_name[i]);
    }
    else 
      for (k=0; (k<Nkpt) && (!feof(infile)); ++k) {
	do {
	  fgets(dump, sizeof(dump), infile);
	} while ((!feof(infile)) && (dump[0] == '#'));
	if (feof(infile)) break;

	// %*d : does conversion, suppresses assignment
	sscanf(dump, "%*d %lf %lf %lf", 
	       &(nu[i][k][0]), &(nu[i][k][1]), &(nu[i][k][2]));
	sort3(nu[i][k]);
	for (d=0; d<3; ++d)
	  if (nu[i][k][d] < -TOLER) imag_freq = -1;
      }
    myclose(infile);
  }
  //-- ==== phonon ====

  if (TESTING) {
    printf("## %d kpoints.\n", Nkpt);
    if (imag_freq)
      printf("## IMAGINARY FREQUENCIES--computation of Debye temp. will be questionable.\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## %2d", k+1);
      for (i=0; i<Nphonon; ++i)
	printf(" %8.5lf %8.5lf %8.5lf",
	       nu[i][k][0], nu[i][k][1], nu[i][k][2]);
      printf("\n");
    }
  }
  if (ERROR) exit(ERROR);
  
  // ***************************** ANALYSIS **************************
  // Debye temps, and rlv error:
  double theta0[Nphonon] = {0.,0.}, theta2[Nphonon] = {0.,0.};
  double relerr = 0.;
  
  for (k=0; k<Nkpt; ++k) {
    for (d=0; d<3; ++d) {
      // Relative error first:
      if (BADRMS) {
	if (! zero(nu[0][k][d]) )
	  relerr += w[k] * square(1 - nu[1][k][d]/nu[0][k][d]);
      }
      else
	relerr += w[k] * square(nu[1][k][d] - nu[0][k][d]);
      // theta's:
      for (i=0; i<Nphonon; ++i) {
	// theta0:
	if (! zero(nu[i][k][d]))
	  theta0[i] += w[k]*log(fabs(nu[i][k][d]));
	// theta2:
	if (nu[i][k][d] > 0.)
	  theta2[i] += w[k]*square(nu[i][k][d]);
	else
	  theta2[i] -= w[k]*square(nu[i][k][d]);
      }
    }
  }

  // now, scale 'em:
  relerr = sqrt(relerr);
  if (! BADRMS) relerr *= 1./sqrt(theta2[0]); // scale by squared phonon
  for (i=0; i<Nphonon; ++i) {
    theta0[i] = hdivkB*exp(theta0[i]);
    theta2[i] = (5./3.)*hdivkB*sqrt(theta2[i]);
  }

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.
  // theta's relative errors:
  printf("%.8lf", fabs(1-theta0[1]/theta0[0]));
  for (i=0; i<Nphonon; ++i)
    printf(" %.8le", theta0[i]);
  printf(" # relative error, theta0");
  if (imag_freq) printf(": imaginary freq!\n");
  else           printf("\n");
    
  printf("%.8lf", fabs(1-theta2[1]/theta2[0]));
  for (i=0; i<Nphonon; ++i)
    printf(" %.8le", theta2[i]);
  printf(" # relative error, theta2\n");
  
  printf("%.8lf # relative RMS error in phonons\n", relerr);

  // ************************* GARBAGE COLLECTION ********************
  return 0;
}
