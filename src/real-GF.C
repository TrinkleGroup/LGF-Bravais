/*
  Program: real-GF.C
  Author:  D. Trinkle
  Date:    February 9, 2004
  Purpose: Given the elastic greens function (calculated with egf-FT)
           we evaluate the greens function in real space at a series
	   of *lattice* sites.  This is useful for... well, I dunno just
	   yet... but gives us something more to play around with.

  Param.:  <cell> <Gk(lm)> <points>
           cell:   cell file describing our lattice
           Gk(lm): spherical harmonic components of elastic GF in k-space
	   points: set of lattice points to eval GF

           ==== cell ====
           a0                            # Scale factor for unit cell
           a1.x a1.y a1.z                # Cartesian coord of unit cell
           a2.x a2.y a2.z
           a3.x a3.y a3.z
           crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
           Natoms                        # Number of atoms in first unit cell
           u1.1 u1.2 u1.3                # Atom locations, in direct coord.
           ...
           uN.1 uN.2 uN.3
           ==== cell ====

           ==== Gk(lm) ====
           Lmax d       # Maximum L value, dimensionality of vector
           l1 m1 R(...) # l,m pair -- real components
           l1 m1 I(...) # l,m pair -- imag components
           l2 m2 R(...) 
           l2 m2 I(...)
           ...
           # Read until EOF or l=-1.
           # NOTE: only EVEN l values should appear!!
           ==== Gk(lm) ====

	   ==== points ====
	   N               # Number of points
	   n1.1 n1.2 n1.3  # unit coord. of lattice points
	   ...
	   nN.1 nN.2 nN.3
	   ==== points ====

  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in the Gk(lm) and our points, and calc the GF at each
           point using the spherical harmonic expansion.

  Output:  I figure for now, I'll just output G_ij for each lattice
	   point, and that should suffice.
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
#include "cell.H"  // Do we still use this?
#include "sphere-harm.H" // Spherical harmonic evaluation
#include "cont-int.H" // included so that we can do the origin


//***************************** SUBROUTINES ****************************

double max_k_BZ (double cart_b[9]) 
{
  int n[3];
  const int nmax = 3; // How far out should we search for nn?
  double bb[9];
  double kmin, kmagn;
  
  square(cart_b, bb);
  kmin = bb[0]; // simple starting guess.

  for (n[0] = -nmax; n[0] <= nmax; ++(n[0]))
    for (n[1] = -nmax; n[1] <= nmax; ++(n[1]))
      for (n[2] = -nmax; n[2] <= nmax; ++(n[2])) {
        if ( (n[0] == 0) && (n[1] == 0) && (n[2] == 0))
          continue;
        kmagn = magnsq(bb, n);
        if (kmagn < kmin) kmin = kmagn;
      }
  return 0.5*sqrt(kmin);
}


// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <Gk(lm)> <points>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  Gk(lm): spherical harmonic components of elastic GF in k-space\n\
  points: set of lattice points to eval GF";

const char* FILEEXPL =
"==== cell ====\n\
a0                            # Scale factor for unit cell\n\
a1.x a1.y a1.z                # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.\n\
Natoms                        # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3                # Atom locations, in direct coord.\n\
...\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== Gk(lm) ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...)\n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
# NOTE: only EVEN l values should appear!!\n\
==== Gk(lm) ====\n\
\n\
==== points ====\n\
N               # Number of points\n\
n1.1 n1.2 n1.3  # unit coord. of lattice points\n\
...\n\
nN.1 nN.2 nN.3\n\
==== points ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *points_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  points_name = args[2];

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input -- we read it, but don't use it.
  int Natoms = NO_ATOMS;
  double** u_atoms = NULL;
  double atomic_mass;

  //++ ==== cell ====
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms);
  // Read in the atomic mass:
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%lf", &atomic_mass);
  myclose(infile);
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }
  if (Natoms != 1) {
    fprintf(stderr, "Sorry.  Currently we can only do single atom cells.\n");
    exit(1);
  }
  if (TESTING) {
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    printf("# atomic mass = %.3lf\n", atomic_mass);
  }
  //-- ==== cell ====

  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cart_rlv[9];
  double cell_vol;
  cell_vol = inverse(cart, cart_rlv);        // invert
  self_transpose(cart_rlv);                  // transpose in place
  mult(cart_rlv, 2*M_PI/cell_vol, cart_rlv); // scale

  //++ ==== G(lm) ====
  // NOTE: this is different than other versions, since we're *forced*
  // to only have EVEN l values.  So [l] is really 2*l.
  int Lmax;  // Maximum l value / 2
  int Ndim;  // Dimensionality -- must be 9
  double ***RYlm, ***IYlm; // Separate our real and imag. components

  infile = myopenr(Ylm_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Ylm_name);
    exit(ERROR_NOFILE);
  }

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d", &Lmax, &Ndim);
  if ( (Lmax < 0) || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
	    Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    // Turn Lmax into Lmax/2:
    Lmax = Lmax/2;
    // Start allocating...
    RYlm = new double **[Lmax+1];
    IYlm = new double **[Lmax+1];
    for (l=0; l<=Lmax; ++l) {
      RYlm[l] = new double *[4*l+1];
      IYlm[l] = new double *[4*l+1];
      for (m=0; m<=(4*l); ++m) {
	RYlm[l][m] = new double[Ndim];
	IYlm[l][m] = new double[Ndim];
	for (d=0; d<Ndim; ++d) {
	  // Initialize
	  RYlm[l][m][d] = 0.;
	  IYlm[l][m][d] = 0.;
	}
      }
    }
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
  myclose(infile);
  //-- ==== G(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

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
    if (TESTING) {
      printf("## Reset Lmax from %d to %d.\n", Lmax, l_max_0);
    }
    Lmax = l_max_0; // Correct Lmax value
  }

  if (TESTING) {
    // Print out *all* the elements we read in...
    printf("## Ylm components:\n");
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) {
	// Real part
	printf("## %3d %3d R ", 2*l, m-2*l);
	print_mat(RYlm[l][m]); printf("\n");
	// Imag part
	printf("## %3d %3d I ", 2*l, m-2*l);
	print_mat(IYlm[l][m]); printf("\n");
      }
    }
  }

  // VERBOSE: print out information about *size* of each component.
  if (TESTING) {
    printf("## Magnitude of terms:\n");
    double zmagn;
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=(4*l); ++m) {
	zmagn = 0.;
	for (d=0; d<Ndim; ++d)
	  zmagn += RYlm[l][m][d]*RYlm[l][m][d] + 
	    IYlm[l][m][d]*IYlm[l][m][d];
	zmagn = sqrt(zmagn /(double)Ndim);
	printf("## %3d %3d %.8le\n", 2*l, m-2*l, zmagn);
      }
    }
  }

  //++ ==== points ====
  int Np;
  int **Runit;
  
  infile = myopenr(points_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", points_name);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Np);
  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in %s\n", Np, points_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    Runit = new int *[Np];
  }
  
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    Runit[n] = new int[3];
    sscanf(dump, "%d %d %d", &(Runit[n][0]), &(Runit[n][1]), &(Runit[n][2]));
  }
  myclose(infile);
  if (n!=Np) {
    fprintf(stderr, "Didn't read enough points in %s\n", points_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== points ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  double Rvect[3]; // Cartesian coord. of vector
  double Rmagn;

  double gc[Np][9]; // Our GF at each point.  
  
  // Spherical harmonics
  SH_workspace_type * w = allocate_SH_workspace(Lmax);

  // Now eval. the inverse FT of our g~_lm -> xYlm_R
  double ***RYlm_R, ***IYlm_R;
  double Fl_val, g_scale = 1./(4*M_PI); // All that'll be left is 1./Rmagn
  RYlm_R = new double **[Lmax+1];
  IYlm_R = new double **[Lmax+1];
  Fl_val = 1.; // This will be prod_i=(1..l/2) (- (2i-1)/2i)
  for (l=0; l<=Lmax; ++l) {
    RYlm_R[l] = new double *[4*l+1];
    IYlm_R[l] = new double *[4*l+1];
    for (m=0; m<=(4*l); ++m) {
      RYlm_R[l][m] = new double[Ndim];
      IYlm_R[l][m] = new double[Ndim];
      for (d=0; d<Ndim; ++d) {
	// Initialize
	RYlm_R[l][m][d] = RYlm[l][m][d] * g_scale * Fl_val;
	IYlm_R[l][m][d] = IYlm[l][m][d] * g_scale * Fl_val;
      }
    }
    // Update Fl_val:
    Fl_val *= - (double)(2*l+1) / (double)(2*l+2);
  }

  // Now, loop over our points:
  for (n=0; n<Np; ++n) {
    if ( (Runit[n][0] == 0) && (Runit[n][1] == 0) && (Runit[n][2] == 0) ) {
      // Do the zero point using integration... VERY approximate,
      // but better than setting it to 0, I've found.  Taken from gf-disc.C
      
      double kmax = max_k_BZ(cart_rlv);
      double g0_scale = M_2_SQRTPI/(8*M_PI*M_PI); // strange (exact) number
      g0_scale *= kmax * kpoint_envelope_int();  // and our maximum k value...
      for (d=0; d<Ndim; ++d)
        gc[n][d] = g0_scale * RYlm[0][0][d];
    }
    else {
      mult_vect(cart, Runit[n], Rvect);
      Rmagn = 1./sqrt(dot(Rvect, Rvect));
      // spherical harmonic:
      eval_Ylm_expansion_R (Lmax, Rvect, 9, RYlm_R, IYlm_R, gc[n], w);
      for (d=0; d<9; ++d) gc[n][d] *= Rmagn;
    }
  }
	
  // garbage collection:
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(4*l); ++m) {
      delete[] RYlm_R[l][m];
      delete[] IYlm_R[l][m];
    }
    delete[] RYlm_R[l];
    delete[] IYlm_R[l];
  }
  delete[] RYlm_R;
  delete[] IYlm_R;
  

  // ****************************** OUTPUT ***************************
  // Human readable (sorta) first:
  
  printf("%d 9 # Number of lattice points\n", Np);
  for (n=0; n<Np; ++n) {
    printf("%3d %3d %3d", Runit[n][0], Runit[n][1], Runit[n][2]);
    for (d=0; d<9; ++d) printf(" %.12le", gc[n][d]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  free_SH_workspace(w);

  // Grid of points
  for (n=0; n<Np; ++n)
    delete[] Runit[n];
  delete[] Runit;

  // spherical harmonic components
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(4*l); ++m) {
      delete[] RYlm[l][m];
      delete[] IYlm[l][m];
    }
    delete[] RYlm[l];
    delete[] IYlm[l];
  }
  delete[] RYlm;
  delete[] IYlm;

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
