/*
  Program: phonon.C
  Author:  D. Trinkle
  Date:    February 12, 2004
  Purpose: Given a spherical harmonic expansion for our elastic Greens
           function, and a cell file, we go off and calculate the
	   phonons (!) for different k-points.

	   Sounds like a dream come true, eh?

	   Later we'll have to add in the lattice GF correction terms
	   to our calc... but that day is not today.

  Param.:  <cell> <G(lm)> <kpts> <Rcut>
           cell:  cell file describing our lattice
           G(lm): file of spherical harmonic components of elastic GF
           kpts:  grid of k-points to evaluate phonons
	   Rcut:  cutoff for transition from discrete to continuum (in Ang)

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

	   ==== G(lm) ====
	   Lmax d       # Maximum L value, dimensionality of vector
	   l1 m1 R(...) # l,m pair -- real components
	   l1 m1 I(...) # l,m pair -- imag components
	   l2 m2 R(...) 
	   l2 m2 I(...)
	   ...
	   # Read until EOF or l=-1.
	   # NOTE: only EVEN l values should appear!!
	   ==== G(lm) ====

	   ==== kpts ====
	   N k0 <type>     # Number of points, scale factor, C/L for units
	   k1.1 k1.2 k1.3
	   ...
	   kN.1 kN.2 kN.3
	   ==== kpts ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Given our lm spherical harmonic components, we calc our
           function at each grid point given.

	   If a (l,m) component pair is *not* given, it is assumed to
	   be zero.  However--we must *always* give both real and imag.
	   parts, even if the imag part is all zero.

  Output:  Just the phonons at each point -- scaling to get meV or THz???
           ==== output ====
           N      # Number of grid points, dimensionality
	   k1.1 k1.2 k1.3  w1.1 w1.2 w1.3  # k-points in Cart, frequencies
	   ...
	   kN.1 kN.2 kN.3  w1.1 w1.2 w1.3  # k-points in Cart, frequencies
           ==== output ====
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "ball.H"  // How to make a sphere of atoms, and smear is out
#include "eigen.H"

// If we are fed our GF in units of A^3/eV, then after we diagonalize
// invert and square root, we multiply by 1/(2Pi)*(ec*Na/10)^(1/2)
// to get the linear frequency nu in THz, once we divide by the
// square root of the mass in amu.  Whew.
const double THz_scale = 15.6333;


//****************************** SUBROUTINES ****************************

inline void calc_angle(double lmn[3], double &cost, double &phi) 
{
  cost = lmn[2];
  if (cost > 1) cost = 1.;
  if (cost < -1) cost = -1.;
  phi = atan2(lmn[1], lmn[0]);
}

// Fold into the BZ
void fold_down(double kpt[3], double cart_b[9]);


// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<cell> <G(lm)> <kpts> <Rcut>";

const int NFLAGS = 6;
const char USERFLAGLIST[NFLAGS] = {'k', 'n', 'd', 'c', 'z', 's'}; 
// flag characters.

const char* ARGEXPL = 
"  cell:  cell file describing our lattice\n\
  G(lm): file of spherical harmonic components of elastic GF\n\
  kpts:  grid of k-points to evaluate phonons\n\
  Rcut:  cutoff for transition from discrete to continuum (in Ang)\n\
  ** output flags: (default: -k)\n\
  -k: output the cartesian kpoint values\n\
  -n: output the kpoint number\n\
  ** fourier transform flags:\n\
  -d: only discrete part for GF\n\
  -c: only continuum part, cutoff by Rcut (set Rcut=0 and -z for full continuum)\n\
  -z: no calculation of GF at R=0; note--this flag shouldn't be set with -d\n\
  -s: set smearing parameter to 0";

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
==== G(lm) ====\n\
Lmax d       # Maximum L value, dimensionality of vector\n\
l1 m1 R(...) # l,m pair -- real components\n\
l1 m1 I(...) # l,m pair -- imag components\n\
l2 m2 R(...) \n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
# NOTE: only EVEN l values should appear!!\n\
==== G(lm) ====\n\
\n\
==== kpts ====\n\
N k0 <type>     # Number of points, scale factor, C/L for units\n\
k1.1 k1.2 k1.3\n\
...\n\
kN.1 kN.2 kN.3\n\
==== kpts ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // Not used
  int KPOINT_CART, KPOINT_NUM; // Describe how we output the kpoint info
  int FT_DISCRETE, FT_CONT;  // Describe which pieces of FT to compute 
  int FT_ZERO;               // Calculate R=0 (needed by discrete piece)
  int SMEARING;              // Use smearing or not?

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      // Note: we don't print out the elastic const. info... mainly
      // because we just ignore all that stuff anyway.
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // Flags
  KPOINT_CART = flagon[0];
  KPOINT_NUM = flagon[1];
  if ( (! KPOINT_CART) && (! KPOINT_NUM) )
    // default, if no flags set.
    KPOINT_CART = 1;

  FT_DISCRETE = flagon[2];
  FT_CONT = flagon[3];
  FT_ZERO = !(flagon[4]);
  if ( (! FT_DISCRETE) && (! FT_CONT) ) {
    // default, if no flags set.
    FT_DISCRETE = 1;
    FT_CONT = 1;
  }

  SMEARING = !(flagon[5]);

  if ( (! FT_ZERO) && FT_DISCRETE && VERBOSE ) {
    printf("##!! Warning: being asked to compute the discrete part WITHOUT R=0 term.\n");
    printf("##!! Results are going to be practically meaningless.\n");
  }
  
  if (TESTING) {
    printf("## Calculating FT using");
    if (FT_CONT) printf(" continuous");
    if (FT_DISCRETE) printf(" discrete");
    printf(" pieces");
    if (FT_ZERO) printf(" including R=0 term.\n");
    else printf(".\n");
    if (! SMEARING) printf("## No smearing.\n");
    printf("## Will output k-point");
    if (KPOINT_NUM) printf(" number");
    if (KPOINT_CART) printf(" position");
    printf("\n");
  }

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *kpts_name;
  double Rcut;
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  kpts_name = args[2];
  sscanf(args[3], "%lf", &Rcut);
  
  if (Rcut < 0.) {
    fprintf(stderr, "Rcut less than 0?\n");
    exit(1);
  }
  if ( (Rcut == 0) && (! FT_DISCRETE) && VERBOSE ) {
    printf("##!! Warning: we're doing the purely elastic, i.e., Rcut=0 version.\n");
    printf("##!! Thus our values near the BZ edges will be terrible.\n");
  }

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

  //++ ==== kpts ====

  int Nkpt;   // Total number of grid points
  double k0;  // Our scale
  char type_ident;  // Character specifying the kpoint type
  int LATT;   // = 1 if we've got RLV, = 0 if cartesian
  
  double** kpt; // Our points -- note: has dimensionality of 4
  // [0..2] is the unit vector, and [3] is the magnitude.
  
  infile = myopenr(kpts_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", kpts_name);
    exit(ERROR_NOFILE);
  }

  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %lf %c", &Nkpt, &k0, &type_ident);

  if (Nkpt < 1) {
    fprintf(stderr, "No points listed in %s\n", kpts_name);
    ERROR = ERROR_BADFILE;
  }
  if (k0 <= 0.) {
    fprintf(stderr, "Bad scale factor in %s\n", kpts_name);
    ERROR = ERROR_BADFILE;
  }
  type_ident = tolower(type_ident);
  if ( (type_ident != 'l') && (type_ident != 'c') ) {
    fprintf(stderr, "Invalid type character (%c) in %s; should be C or L\n", 
	    type_ident, kpts_name);
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
    kpt = new double *[Nkpt];
    for (i=0; i<Nkpt; ++i) kpt[i] = new double[4]; // [0..2] = khat, [3] = k
  }

  // Now, do the readin'!
  for (n=0; (!ERROR) && (n<Nkpt) && (!feof(infile)); ++n) {
    double* kpt_p;
    kpt_p = kpt[n];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", &(kpt_p[0]), &(kpt_p[1]), &(kpt_p[2]));
    for (i=0; i<3; ++i) kpt_p[i] *= k0; // scale!
    // Convert to lattice coord
    double tempvec[3];
    if (! LATT) {
      mult_vect(cart, kpt_p, tempvec);
      for (i=0; i<3; ++i) kpt_p[i] = tempvec[i] / (2*M_PI);
    }
    // Fold into BZ:
    fold_down(kpt_p, cart_rlv);
    // Convert to Cartesian coord:
    mult_vect(cart_rlv, kpt_p, tempvec);
    for (i=0; i<3; ++i) kpt_p[i] = tempvec[i];
    // Magnitude
    kpt_p[3] = sqrt(kpt_p[0]*kpt_p[0] + kpt_p[1]*kpt_p[1] + kpt_p[2]*kpt_p[2]);
    // Scale out: (if k = 0, then set the normalized part to 0)
    if (dcomp(kpt_p[3], 0.)) 
      for (i=0; i<3; ++i) kpt_p[i] = 0.;
    else 
      for (i=0; i<3; ++i) kpt_p[i] *= 1./kpt_p[3];
  }
  // Make sure we read enough points...
  if (n != Nkpt) {
    fprintf(stderr, "Not enough points in %s\n", kpts_name);
    ERROR = ERROR_BADFILE;
  }
  myclose(infile);
  //-- ==== kpts ====

  if (TESTING) {
    printf("## %d kpoints.\n", Nkpt);
    printf("##  i    k_x      k_y      k_z       |k|\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## %2d %8.5lf %8.5lf %8.5lf  %8.5lf\n",
	     k+1, kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
    }
  }
  

  // ***************************** ANALYSIS **************************

  // Allocate space for our GF FT: G_ab(k_vec):
  double **Gabk;
  Gabk = new double*[Nkpt];
  for (k=0; k<Nkpt; ++k) {
    Gabk[k] = new double[Ndim];
    for (d=0; d<Ndim; ++d) Gabk[k][d] = 0.; // initial value.
  }
  
  // To calculate our FT, we have to do three separate calculations
  // and put the pieces together.
  // 1. Calculate G_E with cutoff corrections
  // 2. Calculate G(R) inside the cutoff sphere, and FT
  // 3. Calculate the lattice GF correction inside the sphere, and FT
  
  // For now, we're not going to bother doing 3.

  if (FT_CONT) {
    // *************************** FT_CONT *************************

    // Go through and calculate our spherical harmonic elements at each
    // grid point.
    // Note: RsphP[2l][-m] = RsphP[2l][m], and IsphP[2l][-m] = -IsphP[2l][m]
    //       so we only calc for m=0..2*l
    double ***RsphP, ***IsphP;
    // Allocate first...
    RsphP = new double **[Lmax];
    IsphP = new double **[Lmax];
    for (l=0; l<=Lmax; ++l) {
      RsphP[l] = new double *[2*l+1]; // Only calc for m=0..l
      IsphP[l] = new double *[2*l+1];
      for (m=0; m<=(2*l); ++m) {
	RsphP[l][m] = new double[Nkpt];
	IsphP[l][m] = new double[Nkpt];
      }
    }
    
    // Now fill in the elements: (use gsl_sf_legendre_sphPlm_array)
    double* result_temp;
    result_temp = new double[2*Lmax+1];
    for (m=0; m<=2*Lmax; ++m) {
      for (n=0; n<Nkpt; ++n) {
	double phi, cosm, sinm;
	double cost; // cos(theta) for kpt
	calc_angle(kpt[n], cost, phi);
	// exp(I*m*phi):
	sinm = sin(m*phi);
	cosm = cos(m*phi);
	// Plm(cost):
	gsl_sf_legendre_sphPlm_array(2*Lmax, m, cost, result_temp);
	// Now, put in the appropriate spots:
	for (l=2*Lmax; l>=m; l -= 2) {
	  // result_temp goes from m to Lmax
	  RsphP[l/2][m][n] = cosm*result_temp[l-m];
	  IsphP[l/2][m][n] = sinm*result_temp[l-m];
	}
      }
    }
    delete[] result_temp; // Garbage collection...
    
    
    // To do (1), we need to compute (-1)^l F_l(kRc) by an expansion technique.
    // The key to this is to evaluate all of the coefficients f^(n)_l
    // up to some max n, defined by kmax*Rc and our error tolerance.
    // We include (-1)^l in f^(n)_l, so we have just a polynomial in kRc:
    //   (-1)^l F_l(kRc) = sum_n=0^nmax f^(n)_l * (kRc)^(2n)
    
    int nmax;
    double kmax, toler = 1e-5; // Change this toler as desired
    
    for (kmax=0., n=0; n<Nkpt; ++n)
      if (kmax < kpt[n][3]) kmax = kpt[n][3];
    { double ratio = (kmax*Rcut)*(kmax*Rcut)/toler;
    int n1, n2; // Our two "guesses" for what a reasonable nmax is...
    n1 = lround(exp(log(ratio/16.)/3));
    n2 = lround(0.5*exp(log(ratio)/4));
    if (n1 < n2) nmax = n2;
    else nmax = n1;
    if (nmax <= (Lmax+1)) nmax = Lmax+2; // We always do at least this.
    }
    double **fln;
    double fdenom = 0.5; // fdenom = 1/( (2l+2)! * (4l+1)!! ); init l=0
    fln = new double *[Lmax+1];
    if (TESTING) {
      printf("## f^(n)_l; nmax = %d  toler = %.3le\n", nmax, toler);
      printf("## Rcut = %.5lf\n", Rcut);
    }
    for (l=0; l<=Lmax; ++l) {
      // Allocate
      fln[l] = new double[nmax];
      // Now, eval:
      fln[l][0] = 1.;
      for (i=1; i<=l; ++i) fln[l][0] *= (double)i / ((double)i-0.5);
      if ((l%2)==1) fln[l][0] *= -1; // (-1)^l
      // That was easy... let's do the derivatives:
      for (n=1; n<=l; ++n) fln[l][n] = 0.;
      // only non-zero for l<n:
      fln[l][l+1] = (2*l+1)*fdenom;
      if ( (l%2) == 1) fln[l][l+1] *= -1;
      for (n=l+2; n<nmax; ++n) {
	// Use our recursion relation:
	fln[l][n] = -fln[l][n-1]/(2*n*(2*n-3)*(2*n-2*l-2)*(2*n+2*l-2));
      }
      if (TESTING) {
	printf("## f^(n)_%d =", l);
	for (n=0; n<nmax; ++n) {
	  if (fln[l][n] != 0) printf(" %.4le", fln[l][n]);
	  else printf(" 0");
	  printf("(%d)", n);
	}
	printf("\n");
      }
      // Update our factorial counters (note: use l+1 in formula)
      fdenom *= 1./((2*l+4)*(2*l+3)*(4*l+5)*(4*l+3));
    }
    
    if (TESTING) {
      // Show what values for F_l we get for each kpoint:
      printf("## F_l(kpt):\n");
      for (k=0; k<Nkpt; ++k) {
	double x;
	x = kpt[k][3]*Rcut;
	x *= x;  // x = (kRc)^2
	printf("## |k|= %8.5lf ;", kpt[k][3]);
	for (l=0; l<=Lmax; ++l) {
	  double Fval;
	  Fval = fln[l][nmax-1];
	  for (n=(nmax-2); n>=0; --n) {
	    Fval = fln[l][n] + x*Fval;
	  }
	  printf(" F_%d %.3le", l, Fval);
	}
	printf("\n");
      }
    }
    
    if (TESTING) {
      printf("## **** Continuum calculation ****\n");
    }
    // Now, we can evaluate the elastic - center correction part of G_ab(k).
    double gab_scale;
    gab_scale = 4*M_PI / det(cart); // = 4Pi/cell_volume
    for (k=0; k<Nkpt; ++k) {
      double x;
      x = kpt[k][3]*Rcut;
      x *= x;  // x = (kRc)^2
      double Gk[9];
      for (d=0; d<Ndim; ++d) Gk[d] = 0;
      for (l=0; l<=Lmax; ++l) {
	// Eval Fl(kRc):
	double Fval;
	Fval = fln[l][nmax-1];
	for (n=(nmax-2); n>=0; --n) {
	  Fval = fln[l][n] + x*Fval;
	}
	// Now we just need to eval g_ab(khat):
	// This is the l contribution:
	int mplus, mneg;
	double Rval[Ndim]; // will hold sum_m g_ab(2l,m)*Y_(2l,m)(khat)
	for (d=0; d<Ndim; ++d) {
	  Rval[d] = 0.;
	  // m = 0:
	  // note: xYlm[l][m_index], m_index = 2*l+m, so for m=0, we use 2*l.
	  m_index = 2*l;
	  Rval[d] += RYlm[l][m_index][d]*RsphP[l][0][k];
	  // Next, rest of m terms:
	  for (m=1; m<=2*l; ++m) {
	    mplus = m+2*l;
	    mneg = -m+2*l;
	    Rval[d] += (RYlm[l][mplus][d] + RYlm[l][mneg][d])*RsphP[l][m][k]
	      - (IYlm[l][mplus][d] - IYlm[l][mneg][d])*IsphP[l][m][k];
	  }
	  // Multiply by Fval and gab_scale and add to our Gabk:
	  Gk[d] += Fval*gab_scale*Rval[d];
	}
      }
      // Set our final value:
      for (d=0; d<Ndim; ++d) Gabk[k][d] = Gk[d];
      if (TESTING) {
	printf("## Gk(%8.5lf %8.5lf %8.5lf %8.5lf) = ",
	       kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
	print_mat(Gk);
	printf(" / %8.5lf\n", kpt[k][3]*kpt[k][3]);
      }
    }
    
    // **** Garbage collection
    for (l=0; l<=Lmax; ++l) delete[] fln[l];
    delete[] fln;
    
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=2*l; ++m) {
	delete[] RsphP[l][m];
	delete[] IsphP[l][m];
      }
      delete[] RsphP[l];
      delete[] IsphP[l];
    }
    delete[] RsphP;
    delete[] IsphP;
  }


  if (FT_DISCRETE) {
    // ************************* FT_DISCRETE ***********************
    // This is the part where we sum up the discrete part of the
    // FT inside a sphere cutoff by Rcut

    // For now, just test our refine_Rcut method, and evaluate asmear
    // We set asmear according to the volume per atom
    // and get a new Rcut by refining, including half of our smear
    // -- in theory, our old Rcut will have a value of 0.5 from
    //    smear(); not exactly true, because refine_Rcut will move
    //    it slightly.
    double asmear;
    asmear = cube(det(cart)*0.75/M_PI/(double)Natoms);
    if (!SMEARING) asmear = 0;
    if (TESTING) {
      printf("## smear parameter a = %.5lf\n", asmear);
      printf("## original Rcut = %.5lf\n", Rcut);
    }
    Rcut = refine_Rcut(cart, Rcut+asmear*0.5);

    double **Rlist;
    int Natoms;

    // Make a ball of atoms...
    ERROR = atom_sphere(cart, Rcut, Natoms, Rlist);

    if (TESTING) {
      printf("## new Rcut (added in asmear/2) = %.5lf\n", Rcut);
      // Print out our atoms:
      printf("## Natoms = %d\n", Natoms);
      printf("##  i    R_x      R_y      R_z        |R|      smear(R)\n");
      for (n=0; n<Natoms; ++n) {
	printf("## %2d  %8.5lf %8.5lf %8.5lf  %10.5lf  %.5lf\n", n+1,
	       Rlist[n][0], Rlist[n][1], Rlist[n][2], Rlist[n][3],
	       smear(Rcut, asmear, Rlist[n][3]));
      }
    }
    // NOTE: the point at Natoms-1 is 0,0,0; currently, we're
    // treating it "specially" -- i.e., ignoring it.

    // Now we pull a stunt to evaluate the Ylm at each Rhat, like
    // we did for k's.  Should look familiar at least.
    double ***RsphP, ***IsphP;
    // Allocate first...
    RsphP = new double **[Lmax];
    IsphP = new double **[Lmax];
    for (l=0; l<=Lmax; ++l) {
      RsphP[l] = new double *[2*l+1]; // Only calc for m=0..l
      IsphP[l] = new double *[2*l+1];
      for (m=0; m<=(2*l); ++m) {
	RsphP[l][m] = new double[Natoms-1];
	IsphP[l][m] = new double[Natoms-1];
      }
    }
    
    // Now fill in the elements: (use gsl_sf_legendre_sphPlm_array)
    double* result_temp;
    result_temp = new double[2*Lmax+1];
    for (m=0; m<=2*Lmax; ++m) {
      for (n=0; n<(Natoms-1); ++n) {
	double phi, cosm, sinm;
	double cost; // cos(theta) for kpt
	calc_angle(Rlist[n], cost, phi);
	// exp(I*m*phi):
	sinm = sin(m*phi);
	cosm = cos(m*phi);
	// Plm(cost):
	gsl_sf_legendre_sphPlm_array(2*Lmax, m, cost, result_temp);
	// Now, put in the appropriate spots:
	for (l=2*Lmax; l>=m; l -= 2) {
	  // result_temp goes from m to Lmax
	  RsphP[l/2][m][n] = cosm*result_temp[l-m];
	  IsphP[l/2][m][n] = sinm*result_temp[l-m];
	}
      }
    }
    delete[] result_temp; // Garbage collection...

    if (TESTING) {
      printf("## **** Discrete calculation ****\n");
    }
    // Now the g_ab(Rhat)*cos(kR(khat.Rhat))*smear(R) / R for each:
    for (n=0; n<(Natoms-1); ++n) {
      // calc g_ab(Rhat):
      double Rval[Ndim]; // will hold sum_lm g_ab(2l,m)*Y_(2l,m)(Rhat)
      double GR[Ndim];
      for (d=0; d<Ndim; ++d) Rval[d] = 0.;
      for (d=0; d<Ndim; ++d) GR[d] = 0;
      for (l=0; l<=Lmax; ++l) {
	// This is the l contribution:
	int mplus, mneg;
	for (d=0; d<Ndim; ++d) {
	  // m = 0:
	  // note: xYlm[l][m_index], m_index = 2*l+m, so for m=0, we use 2*l.
	  m_index = 2*l;
	  Rval[d] += RYlm[l][m_index][d]*RsphP[l][0][n];
	  // Next, rest of m terms:
	  for (m=1; m<=2*l; ++m) {
	    mplus = m+2*l;
	    mneg = -m+2*l;
	    Rval[d] += (RYlm[l][mplus][d] + RYlm[l][mneg][d])*RsphP[l][m][n]
	      - (IYlm[l][mplus][d] - IYlm[l][mneg][d])*IsphP[l][m][n];
	  }
	}
	// include the smear(R)/R contribution:
	double Rscale, gtemp[Ndim];
	Rscale = smear(Rcut, asmear, Rlist[n][3])/Rlist[n][3];
	mult(Rval, Rscale, gtemp);
	for (d=0; d<Ndim; ++d) GR[d] += gtemp[d];
	// loop over k-points, and include the k^2 cos(k.R) term:
	for (k=0; k<Nkpt; ++k) {
	  double k2_cos_kR; // k^2 * cos(k.R)
	  k2_cos_kR = cos(kpt[k][3]*Rlist[n][3]*dot(kpt[k], Rlist[n]))
	    * kpt[k][3] * kpt[k][3];
	  for (d=0; d<Ndim; ++d) Gabk[k][d] += k2_cos_kR*gtemp[d];
        }
      }
      if (TESTING) {
	printf("## GR(%8.5lf %8.5lf %8.5lf %8.5lf) = ",
	       Rlist[n][0], Rlist[n][1], Rlist[n][2], Rlist[n][3]);
	print_mat(GR);
	printf("\n");
      }
    }
    
    // **** Garbage collection
    for (l=0; l<=Lmax; ++l) {
      for (m=0; m<=2*l; ++m) {
	delete[] RsphP[l][m];
	delete[] IsphP[l][m];
      }
      delete[] RsphP[l];
      delete[] IsphP[l];
    }
    delete[] RsphP;
    delete[] IsphP;

    for (n=0; n<Natoms; ++n)
      delete[] Rlist[n];
    delete[] Rlist;
  }


  if (FT_ZERO) {
    // *************************** FT_ZERO *************************
    // We treat the Greens function value at 0 in a special way.
    double G0[Ndim];
    const double PRE_FACTOR = 0.2228139187;  // 6^(1/3) / Pi^(11/6)
    double a_cell;
    a_cell = cube(det(cart)); // (Cell volume)^(1/3)
    // G_ab(0) = pre * g_ab(0,0) / a -- very simple.
    for (d=0; d<Ndim; ++d) 
      G0[d] = (PRE_FACTOR/a_cell) * RYlm[0][0][d];
    if (TESTING) {
      printf("## **** R=0 calculation ****\n");
    }
    if (TESTING) {
      printf("## GF at R=0 from FT:\n");
      printf("## G(R=0) = ");
      print_mat(G0); printf("\n");
    }
    // NOTE: We include a factor of k^2 in G_ab(k), because
    // we divide out k^2 in the continuum part.
    for (k=0; k<Nkpt; ++k) {
      for (d=0; d<Ndim; ++d) Gabk[k][d] += kpt[k][3]*kpt[k][3]*G0[d];
    }
  }
  
  if (TESTING) {
    printf("## **** Final Green's function FT ****\n");
    for (k=0; k<Nkpt; ++k) {
      printf("## G~(%8.5lf %8.5lf %8.5lf %8.5lf) = ",
	     kpt[k][0], kpt[k][1], kpt[k][2], kpt[k][3]);
      print_mat(Gabk[k]);
      printf(" / %8.5lf\n", kpt[k][3]*kpt[k][3]);
    }
  }

  // *********************** FREQUENCY CALCULATION *******************
  // Compute the frequencies:
  double omega[Nkpt][3], omega_scale;
  omega_scale = THz_scale/sqrt(atomic_mass);
  if (TESTING) {
    printf("## **** Frequency calculation ****\n");
  }
  for (k=0; k<Nkpt; ++k) {
    double lambda[3];
    eigen(Gabk[k], lambda); // Get the eigenvalues.
    if (TESTING) {
      printf("## lambda_%d = %.5lf %.5lf %.5lf\n", k+1, lambda[0],lambda[1],lambda[2]);
    }
    for (d=0; d<3; ++d) {
      if (dcomp(lambda[d], 0.)) {
	if (! dcomp(kpt[k][3], 0.)) {
	  fprintf(stderr, "At kpt %d %.3lf %.3lf %.3lf, got zero inverse frequency without zero k-point magnitude\n",
		  k+1, kpt[k][0]*kpt[k][3], kpt[k][1]*kpt[k][3], 
		  kpt[k][2]*kpt[k][3]);
	}
	omega[k][d] = 0.;
	continue;
      }
      if (lambda[d] < 0) 
	omega[k][d] = -omega_scale*kpt[k][3]/sqrt(-lambda[d]);
      else
	omega[k][d] =  omega_scale*kpt[k][3]/sqrt(lambda[d]);
    }
    // Now, we'll order them by size:
    sort3(omega[k]);
  }

  // ****************************** OUTPUT ***************************
  // Now, we output EVERYTHING.

  printf("%d # Nkpt\n", Nkpt);
  for (k=0; k<Nkpt; ++k) {
    if (KPOINT_NUM)
      printf("%3d ", k);
    if (KPOINT_CART)
      printf("%15.12lf %15.12lf %15.12lf ", 
	     kpt[k][0]*kpt[k][3], kpt[k][1]*kpt[k][3], kpt[k][2]*kpt[k][3]);
    // Output our "frequencies":
    printf(" %15.12lf %15.12lf %15.12lf\n", 
	   omega[k][0], omega[k][1], omega[k][2]);
  }

  // ************************* GARBAGE COLLECTION ********************

  // Greens function
  for (n=0; n<Nkpt; ++n) delete[] Gabk[n];
  delete[] Gabk;

  // Grid points
  for (i=0; i<Nkpt; ++i) delete[] kpt[i];
  delete[] kpt;

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


// Fold into the BZ
void fold_down(double kpt[3], double cart_b[9]) 
{
  double bb[9];
  int a[3], a0[3], a1[3];
  double q[3], k0[3];
  double kmagn, qmagn;
  int i;
  int better;

  square(cart_b, bb);
  
  for (i=0; i<3; ++i) {
    kpt[i] = kpt[i] - (int)(kpt[i]);
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
	    if (dcomp(q[2], k0[2]))
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
}
