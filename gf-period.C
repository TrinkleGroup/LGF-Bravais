/*
  Program: gf-period.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. elastic GF in k-space spherical harmonic components
	   3. discrete correction evaluated at a set of points
	   4. a periodic supercell

	   and we compute
	   1. periodic GF evaluated for a test force applied at the
	   origin (with constant compensating forces all around)

  Param.:  <cell> <Gk(lm)> <Gd(R)> <supercell>
           cell:    cell file describing our lattice
           Gk(lm):  spherical harmonic components of elastic GF in k-space
	   Gd(R):   discrete correction to GF evaluated at series of lattice sites
	   supercell: supercell geometry

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

	   ==== Gd(R) ====
	   N kmax                      # number of points, kmax
	   n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG
	   ...
	   nN.1 nN.2 nN.3  Gxx .. Gzz
	   ==== Gd(R) ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Currently just testing the gridding of k-points

  Output:  Going to output a series of lattice points inside the supercell, 
           and the periodic GF evaluated at those points.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> // GSL error handling
#include <gsl/gsl_sf_legendre.h> // Spherical legendre polynomials
#include <gsl/gsl_sf_bessel.h>   // Spherical bessel functions
#include <gsl/gsl_integration.h> // Needed for doing integration
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "cont-int.H" // For doing the integration over the continuous piece
#include "sphere-harm.H" // Spherical harmonic evaluation
#include "supercell.H" // supercell 

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
double max_k_BZ (double cart_b[9]);

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
const char* ARGLIST = "<cell> <Gk(lm)> <Gd(R)> <supercell>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Gk(lm):  spherical harmonic components of elastic GF in k-space\n\
  Gd(R):   discrete correction to GF evaluated at series of lattice sites\n\
  supercell: our supercell geometry";

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
l2 m2 R(...) \n\
l2 m2 I(...)\n\
...\n\
# Read until EOF or l=-1.\n\
# NOTE: only EVEN l values should appear!!\n\
==== Gk(lm) ====\n\
\n\
==== Gd(R) ====\n\
N kmax                      # number of points, kmax\n\
n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GF\n\
...\n\
nN.1 nN.2 nN.3  Gxx .. Gzz\n\
==== Gd(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *Gd_name, *super_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  Gd_name = args[2];
  super_name = args[3];
  
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

  //++ ==== Gk(lm) ====
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
  
  ERROR = read_Ylm_even(infile, Lmax, Ndim, RYlm, IYlm);
  myclose(infile);

  if ( ERROR || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
	    Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== Gk(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

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

  //++ ==== Gd ====
  int Np;
  double kmax;
  int **Runit;
  double **Gd;
  
  infile = myopenr(Gd_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Gd_name);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %lf", &Np, &kmax);
  if (Np < 1) {
    fprintf(stderr, "Bad Np (%d) value in %s\n", Np, Gd_name);
    ERROR = ERROR_BADFILE;
  }
  else {
    Runit = new int *[Np];
    Gd = new double *[Np];
  }
  
  for (n=0; (!ERROR) && (!feof(infile)) && (n<Np); ++n) {
    do {
      fgets(dump, sizeof(dump), infile);
    } while ((!feof(infile)) && (dump[0] == '#'));
    if (feof(infile)) break;
    // Parse it.
    Runit[n] = new int[3];
    Gd[n] = new double[Ndim];
    sscanf(dump, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	   &(Runit[n][0]), &(Runit[n][1]), &(Runit[n][2]),
	   &(Gd[n][0]), &(Gd[n][1]), &(Gd[n][2]),
	   &(Gd[n][3]), &(Gd[n][4]), &(Gd[n][5]),
	   &(Gd[n][6]), &(Gd[n][7]), &(Gd[n][8]));
  }
  myclose(infile);
  if (n!=Np) {
    fprintf(stderr, "Didn't read enough points in %s\n", Gd_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== Gd ====

  if (TESTING) {
    printf("## Discrete correction:\n");
    for (np=0; np<Np; ++np) {
      printf("## %3d %3d %3d  ", Runit[np][0],  Runit[np][1],  Runit[np][2]);
      print_mat(Gd[np]);
      printf("\n");
    }
  }

  //++ ==== supercell ====
  int super_n[Ndim];
  infile = myopenr(super_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", super_name);
    exit(ERROR_NOFILE);
  }
  // For now, just read the matrix (later on, read the cartesian coord?)
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d %d %d %d %d %d %d",
	 &(super_n[0]), &(super_n[1]), &(super_n[2]), 
	 &(super_n[3]), &(super_n[4]), &(super_n[5]), 
	 &(super_n[6]), &(super_n[7]), &(super_n[8]));
  myclose(infile);
  
  //-- ==== supercell ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }


  // ***************************** ANALYSIS **************************
  // Make our supercell
  int **atomlist;
  int Nsuper;
  fill_super(cart, super_n, Nsuper, atomlist);
  if (TESTING) {
    printf("## Supercell geometry (%d atoms):\n", Nsuper);
    for (n=0; n<Nsuper; ++n) {
      printf("## %3d : %3d %3d %3d\n", n+1, 
	     atomlist[n][0], atomlist[n][1], atomlist[n][2]);
    }
  }

  // Now, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[Ndim], cartB[Ndim];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

  // Make our set of gamma equivalent kpoints:
  int** klist;
  int Nkpt;
  int transn[Ndim]; // transpose of super_n
  transpose(super_n, transn);
  int inv_transn[Ndim];
  Nkpt = inverse(transn, inv_transn);
  mult(cartb, inv_transn, cartB);
  mult(cartB, 1./Nkpt, cartB);

  // Now use fill_super to make the points
  fill_super(cartB, transn, Nkpt, klist);
  // Finally, write in terms of rlv and cartesian coordinates:
  // take out gamma:
  double klist_unit[Nkpt-1][3], klist_cart[Nkpt-1][4]; // [3] = magn
  for (n=0; n<(Nkpt-1); ++n) {
    int temp[3];
    mult_vect(inv_transn, klist[n+1], temp);
    for (i=0; i<3; ++i) klist_unit[n][i] = temp[i]*(1./Nkpt);
    // Now, make cartesian:
    mult_vect(cartb, klist_unit[n], klist_cart[n]);
    // magnitude:
    klist_cart[n][3] = sqrt(dot(klist_cart[n], klist_cart[n]));
  }
  // Garbage collection...
  free_super(Nkpt, klist);
  --Nkpt;

  if (TESTING) {
    printf("## %d kpoints equivalent to gamma:\n", Nkpt);
    for (n=0; n<Nkpt; ++n) {
      printf("## %3d : %8.5lf %8.5lf %8.5lf | %8.5lf %8.5lf %8.5lf | %8.5lf\n", 
	     n+1,
	     klist_unit[n][0], klist_unit[n][1], klist_unit[n][2], 
	     klist_cart[n][0], klist_cart[n][1], klist_cart[n][2],
	     klist_cart[n][3]);
    }
  }

  // Make our list of supercell translations we need to sum up Gd:
  int*** index, max_index; // Our index of points
  int** shiftlist, Nshift; // the shifts we need

  // make our index
  make_index(Np, Runit, max_index, index);
  // make our list of shifts
  make_super_shifts(max_index, index, super_n, Nshift, shiftlist);
  if (TESTING) {
    printf("## Indexing array: max_index = %d\n", max_index);
    int testp[3] = {0,0,0}, ind;
    ind = get_index(max_index, index, testp);
    printf("## 000 index = %d  Runit[] = %3d %3d %3d\n",
	   ind, Runit[ind][0], Runit[ind][1], Runit[ind][2]);
    printf("## checking index...\n");
    for (np=0; np<Np; ++np)
      if (get_index(max_index, index, Runit[np]) != np) {
	printf("## FAILURE for point %d\n", np);
      }
    printf("## ...check complete.\n");
    printf("## Nshifts = %d\n", Nshift);
    for (n=0; n<Nshift; ++n)
      printf("## %3d : %3d %3d %3d\n", n+1,
	     shiftlist[n][0], shiftlist[n][1], shiftlist[n][2]);
  }
  
  // NOW!! The real analysis can begin.
  double** Gperiod; // for each supercell point...
  Gperiod = new double *[Nsuper];
  for (np=0; np<Nsuper; ++np) Gperiod[np] = new double[Ndim];

  // First, evaluate Gsc at each k_gamma point:
  double** Gsc_kgamma;
  Gsc_kgamma = new double *[Nkpt];
  for (n=0; n<Nkpt; ++n) Gsc_kgamma[n] = new double[Ndim];

  int Nkpt_inc; // how many kpooints are less than kmax
  // Workspace for spherical harmonic evaluation
  SH_workspace_type * w = allocate_SH_workspace(Lmax);
  Nkpt_inc = Nkpt;
  for (k=0; k<Nkpt_inc; ++k) {
    // Our klist is sorted by magnitude, so once we find one point
    // above kmax, all points are above kmax.
    if (klist_cart[k][3] >= kmax)
      Nkpt_inc = k; // our way of shortcutting out...
    else {
      eval_Ylm_expansion_R (Lmax, klist_cart[k], Ndim, RYlm, IYlm, 
			    Gsc_kgamma[k], w);
      // Multiply by our envelope function, divide by a^3*k^2:
      double g_scale = kpoint_envelope(klist_cart[k][3] / kmax) / 
	(cell_vol * klist_cart[k][3] * klist_cart[k][3]);
      for (d=0; d<Ndim; ++d) Gsc_kgamma[k][d] *= g_scale;
    }
  }
  // Now evaluate the 1/N sum G_d(R) term:
  double G_dsum[Ndim];
  for (d=0; d<Ndim; ++d) G_dsum[d] = 0;
  for (np=0; np<Np; ++np)
    for (d=0; d<Ndim; ++d) G_dsum[d] += Gd[np][d];
  for (d=0; d<Ndim; ++d) G_dsum[d] *= 1./Nsuper;

  if (TESTING) {
    printf("## Gsc_kgamma: (%d points)\n", Nkpt_inc);
    for (k=0; k<Nkpt_inc; ++k) {
      printf("## %d %8.5lf %8.5lf %8.5lf ", k+1, 
	     klist_unit[k][0], klist_unit[k][1], klist_unit[k][2]);
      print_mat(Gsc_kgamma[k]);
      printf("\n");
    }
    printf("## G_dsum: ");
    print_mat(G_dsum);
    printf("\n");
  }
  
  
  // ==== Finally, put it all together in one package:
  if (TESTING) {
    printf("## evaluation:\n");
  }
  for (n=0; n<Nsuper; ++n) {
    double *gperiod = Gperiod[n];
    int *atom = atomlist[n];
    // zero out Gperiod:
    for (d=0; d<Ndim; ++d) gperiod[d] = 0.;
    // kgamma sum:
    for (k=0; k<Nkpt_inc; ++k) {
      // Do dot product in direct coord:
      double cos_kR = cos(2.*M_PI*dot(klist_unit[k], atom));
      for (d=0; d<Ndim; ++d) gperiod[d] += cos_kR*Gsc_kgamma[k][d];
    }
    // Scale:
    for (d=0; d<Ndim; ++d) gperiod[d] *= 1./Nsuper;

    if (TESTING) {
      printf("## Runit[ %3d ] = %3d %3d %3d\n", n+1, atom[0],atom[1],atom[2]);
      printf("## Gsc: ");
      print_mat(gperiod);
      printf("\n");
    }
    
    // S:[A] sum:
    for (i=0; i<Nshift; ++i) {
      int ind, point[3];
      for (k=0; k<3; ++k) point[k] = atom[k] + shiftlist[i][k];
      ind = get_index(max_index, index, point);
      if (ind != -1) {
	// Add contribution
	for (d=0; d<Ndim; ++d) gperiod[d] += Gd[ind][d];

	if (TESTING) {
	  printf("##  + Gd( %3d %3d %3d ) =\n", point[0],point[1],point[2]);
	  printf("##      ");
	  print_mat(Gd[ind]);
	  printf("\n");
	}
      }
    }
    
    // Finally, subtract 1/N sum G_d(R) term:
    for (d=0; d<Ndim; ++d) gperiod[d] -= G_dsum[d];

    if (TESTING) {
      printf("## TOTAL: ");
      print_mat(gperiod);
      printf("\n\n");
    }
  }

  // ****************************** OUTPUT ***************************
  // A slightly different header file than usual:
  for (d=0; d<Ndim; ++d)
    printf("%3d ", super_n[d]);
  printf("# supercell definition\n");
  printf("%d # Number of points\n", Nsuper);
  for (n=0; n<Nsuper; ++n) {
    printf("%3d %3d %3d", atomlist[n][0], atomlist[n][1], atomlist[n][2]);
    for (d=0; d<Ndim; ++d)
      printf(" %.8le", Gperiod[n][d]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Nkpt; ++n) delete[] Gsc_kgamma[n];
  delete[] Gsc_kgamma;

  for (np=0; np<Nsuper; ++np) delete[] Gperiod[np];
  delete[] Gperiod;

  free_index(max_index, index);
  free_super(Nshift, shiftlist);
  
  free_SH_workspace(w);
  free_super(Nsuper, atomlist);

  for (np=0; np<Np; ++np) {
    delete[] Runit[np];
    delete[] Gd[np];
  }
  delete[] Runit;
  delete[] Gd;
  
  // spherical harmonic components
  free_Ylm(Lmax, RYlm, IYlm);

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
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
