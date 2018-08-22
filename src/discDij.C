/*
  Program: discDij.C
  Author:  D. Trinkle
  Date:    February 18, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. elastic GF in k-space spherical harmonic components
	   3. second order correction to GF in k-space Ylm
	   4. dynamical matrix evaluated at a set of points
	   5. a symmetrized kpoint mesh for integration
	   6. a set of points to evaluate the discretization correction

	   and we compute
	   1. the cutoff kmax -- maximum k value in our sphere
	   2. the discrete correction for the GF at each of the points

	   Changed: for small k points, we calculate the deviation from
	   D(k)^-1 very carefully.  We have to do this, else
	   the error in Gk(lm) and G0(lm) makes us want to kill ourselves.
	   Seriously: it will induce homi-/sui-cidal feelings when you
	   realize how bad it is.  Grep krealsmall to see when we switch
	   over.

	   It's worth mentioning that the small k expansion is EXACT
	   for ALL k values PROVIDED we are BELOW where the envelope
	   function kicks in... it's just not that efficient, so it's
	   worth our while to only use it when we absolutely have to.

  Param.:  <cell> <Gk(lm)> <G0(lm)> <Dij(R)> <kpts> <latt>
           cell:    cell file describing our lattice
           Gk(lm):  spherical harmonic components of elastic GF in k-space
           G0(lm):  spherical harmonic components of discontiunity correction
	   Dij(R):  dynamical matrix evaluated at a series of points
	   kpts:    symmetrized kpt mesh (with weights)
	   latt:    list of lattice points to eval Gdisc

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

           ==== Dij(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij
           ...
           nN.1 nN.2 nN.3  Dxx .. Dzz
           ==== Dij(R) ====

           ==== kpts ====
           N k0 <type>       # Number of points, scale factor, C/L for units
           k1.1 k1.2 k1.3 w1 # kpoint and weight
           ...
           kN.1 kN.2 kN.3 wN
           ==== kpts ====

           ==== latt ====
           N               # number of points
           n1.1 n1.2 n1.3  # unit cell coord
           ...
           nN.1 nN.2 nN.3
           ==== latt ====
	   

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   Currently just testing the gridding of k-points

  Output:  The discretization correction of the lattice GF.
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
#include "pointgroup.H"
#include "Dij.H"
#include "kpts.H"
#include "shell.H"
#include "eigen.H"
#include "fourier.H"
#include "lambda.H" // for evaluation (on the fly) of lambda^-1 and lambda^(2)

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
		   double lambda[9][9], double lambda2[9][9][9],
		   const int TYPE, const coskR_table_type* table);

inline void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
			  double lambda[9][9], double lambda2[9][9][9],
			  const int TYPE, const coskR_table_type &table) 
{ return smallk_expand(Np, p, kpt, pk, lambda, lambda2, TYPE, &table);}
  
inline void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
			  double lambda[9][9], double lambda2[9][9][9],
			  const int TYPE) 
{ return smallk_expand(Np, p, kpt, pk, lambda, lambda2, TYPE, NULL);}



// Determine the maximum magnitude k sphere we can inscribe in
// our BZ
double max_k_BZ (double cart_b[9]);

// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  //  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
  printf("xx= %15.8le yy= %15.8le zz= %15.8le  xy= %15.8le yz= %15.8le zx= %15.8le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 6;
const char* ARGLIST = "<cell> <Gk(lm)> <G0(lm)> <Dij(R)> <kpts> <latt>";

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'l','s','b'}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Gk(lm):  spherical harmonic components of elastic GF in k-space\n\
  G0(lm):  spherical harmonic components of discontiunity correction\n\
  Dij(R):  dynamical matrix evaluated at a series of points\n\
  kpts:    symmetrized kpt mesh (with weights)\n\
  latt:    list of lattice points to eval Gdisc\n\
  -l: implement a lookup table for more efficient FT--highly recommended\n\
  -s: enforce symmetrization of GF at each k-point\n\
  -b: do FT's with binning for more accurate (?) summation (turns -l on too)";

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
==== Dij(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Dxx .. Dzz  # unit cell coord, Dij\n\
...\n\
nN.1 nN.2 nN.3  Dxx .. Dzz\n\
==== Dij(R) ====\n\
\n\
==== kpts ====\n\
N k0 <type>     # Number of points, scale factor, C/L for units\n\
k1.1 k1.2 k1.3\n\
...\n\
kN.1 kN.2 kN.3\n\
==== kpts ====\n\
\n\
==== latt ====\n\
N               # number of points\n\
n1.1 n1.2 n1.3  # unit cell coord\n\
...\n\
nN.1 nN.2 nN.3\n\
==== latt ====";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter
  int LOOKUP = 0;   // use a lookup table for FT--more efficient
  int SYMM = 0;     // enforce symmetrization at each kpoint
  int BINSUM = 0;   // use bin summation technique with FT

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

  // flags
  LOOKUP = flagon[0];
  SYMM = flagon[1];
  BINSUM = flagon[2];
  if (BINSUM) LOOKUP = 1; // need lookup ...

  if (TESTING) {
    if (LOOKUP) printf("## Using lookup tables for FT\n");
    if (SYMM) printf("## Enforcing symmetrization at each kpoint\n");
    if (BINSUM) printf("## Using bin summation (requires lookup)\n");
  }
  
  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *Ylm_disc_name, 
    *Dij_name, *kpt_name, *latt_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  Ylm_disc_name = args[2];
  Dij_name = args[3];
  kpt_name = args[4];
  latt_name = args[5];

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

  ERROR = read_Ylm_even(infile, Lmax, Ndim, RYlm, IYlm);
  myclose(infile);

  if ( ERROR || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
            Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== G(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm, IYlm, "##");
  }


  //++ ==== G0(lm) ====
  // NOTE: this is different than other versions, since we're *forced*
  // to only have EVEN l values.  So [l] is really 2*l.
  int Lmax_0;  // Maximum l value / 2
  double ***RYlm_0, ***IYlm_0; // Separate our real and imag. components

  infile = myopenr(Ylm_disc_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Ylm_disc_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Ylm_even(infile, Lmax_0, Ndim, RYlm_0, IYlm_0);
  myclose(infile);

  if ( ERROR || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
            Lmax_0, Ndim, Ylm_disc_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== G0(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax_0, Ndim, RYlm_0, IYlm_0, "##");
  }


  //++ ==== Dij(R) ====
  int Np;
  point_type* Dij;
  
  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Dij(infile, cart, Np, Dij, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== Dij(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", Dij_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (n=0; n<Np; ++n) {
      printf("## %3d %3d %3d ", 
	     Dij[n].Runit[0], Dij[n].Runit[1], Dij[n].Runit[2]);
      print_mat(Dij[n].mat);
      printf("\n");
    }
  }


  //++ ==== kpts ====
  int Nkpt;
  double** kpt, *w;
  
  infile = myopenr(kpt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", kpt_name);
    exit(ERROR_NOFILE);
  }
  // the last entry specifies what we read and return:
  ERROR = read_kpts(infile, cart_rlv, Nkpt, kpt, w, 
		    READ_KPTS_MAGN | READ_KPTS_CART | READ_KPTS_WEIGHTS);
  myclose(infile);
  //-- ==== kpts ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", kpt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## %d kpoints\n", Nkpt);
    for (i=0; i<Nkpt; ++i)
      printf("## %8.5lf %8.5lf %8.5lf %.5le\n", 
	     kpt[i][0], kpt[i][1], kpt[i][2], w[i]);
  }


  //++ ==== latt ====
  // Finally, the set of lattice points for which we want to calculate
  // our discretization correction
  point_type *Gd; // our discretization correction; held in .mat
  int Nlatt; // number of lattice points

  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Nlatt, Gd, READ_DIJ_NOMAT);
  myclose(infile);
  //-- ==== latt ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", latt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## lattice points to calc Gd:\n");
    for (n=0; n<Nlatt; ++n)
      printf("## %3d %3d %3d\n", 
	     Gd[n].Runit[0], Gd[n].Runit[1], Gd[n].Runit[2]);
  }

  // ***************************** ANALYSIS **************************
  // lookup table for FT:
  coskR_table_type coskR;
  mat_bin_type *mbinp;
  if (LOOKUP) {
    ERROR = init_coskR_table(Nkpt, kpt, KPT_CART, cart, coskR);
    if (ERROR) {
      fprintf(stderr, "Problem creating lookup table... perhaps kpt mesh isn't uniform?\n"
	      );
      exit(ERROR);
    }
  }
  if (BINSUM) {
    mbinp = init_mat_bin(1./1099511627776., 1099511627776., 16);
    ERROR = (mbinp == NULL);
    if (!ERROR) ERROR = (! mbinp->alloc);
    if (ERROR) {
      fprintf(stderr, "Problem creating bins... this shouldn\'t happen.\n");
      exit(ERROR);
    }
  }
      
    

  // Let's get our point group information:
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  // Symmetrization information:
  int* kpt_type; // specifies which type each kpoint is
  int Nkpt_type;
  double** Sk;   // the symmetrization matrix for each type
  
  if (SYMM) {
    // 0. allocate
    kpt_type = new int[Nkpt];
    // 1. make the sets for each k-point
    double kdiff[3], kint[3];
    // Harold Stokes says there are 98 subgroups for the cubic point group.
    const int MAX_type = 98;
    int** gk_equal = new int*[MAX_type];
    int* gk_count = new int[MAX_type];
    int* t_gk = new int[Nop]; // a temporary holding place...
    int found;
    Nkpt_type=0;
    for (k=0; k<Nkpt; ++k) {
      double *tkpt = kpt[k];
      int count=0;
      for (n=0; n<Nop; ++n) {
	mult_vect(gcart[n], tkpt, kdiff);
	for (d=0; d<3; ++d) kdiff[d] -= tkpt[d];
	for (d=0; d<3; ++d) kint[d] = round(kdiff[d]);
	t_gk[n] = equal_vect(kdiff, kint);
	if (t_gk[n]) ++count;
      }
      // See if we have a new type or not...
      found=0;
      for (i=0; (i<Nkpt_type) && (!found); ++i)
	if (count == gk_count[i]) {
	  found = 1;
	  for (n=0; (n<Nop) && found; ++n)
	    found = (gk_equal[i][n] == t_gk[n]);
	}
      // if we've found it, we're done...
      if (found)
	kpt_type[k] = (i-1); // have to decrement, because we shoot past
      else {
	if (Nkpt_type == MAX_type) {
	  fprintf(stderr, "ERROR!  Too many subgroups... Harold Stokes owes you a beer.\n");
	  ERROR = -1;
	  break;
	}
	// we need to add it:
	kpt_type[k] = Nkpt_type;
	gk_equal[Nkpt_type] = t_gk;
	gk_count[Nkpt_type] = count;
	// "increment":
	t_gk = new int[Nop];
	++Nkpt_type;
      }
    }
    if (! found) delete[] t_gk; // garbage coll. (in case we made a new t_gk)

    // 2. make the g_i_ginv matrices:
    double g_i_ginv[MAXop][Ni][Ni];
    for (n=0; n<Nop; ++n) 
      for (i=0; i<Ni; ++i)
	gen_g_i_ginv (i, gcart[n], gcart[inv_index[n]], g_i_ginv[n][i]);

    // 3. make the symmetrization matrices:
    Sk = new double*[Nkpt_type];
    double Skt[Ni*Ni];
    for (int t=0; t<Nkpt_type; ++t) {
      Sk[t] = new double[9*9];
      for (d=0; d<(Ni*Ni); ++d) Skt[d] = 0.;
      for (n=0; n<Nop; ++n) 
	if (gk_equal[t][n]) {
	  for (i=0; i<Ni; ++i)
	    for (j=0; j<Ni; ++j)
	      Skt[i*Ni+j] += g_i_ginv[n][j][i];
	}
      // scale:
      for (d=0; d<(Ni*Ni); ++d) Skt[d] *= 1./gk_count[t];
      // expand out to Sk[t]:
      int e;
      for (d=0; d<9; ++d) {
	for (e=0; e<9; ++e) Sk[t][d*9+e] = 0.;
	for (j=0; j<Ni; ++j) {
	  Sk[t][d*9+i6tom9[j]] += 0.5*Skt[m9toi6[d]*Ni+j];
	  Sk[t][d*9+trans[i6tom9[j]]] += 0.5*Skt[m9toi6[d]*Ni+j];
	}
      }
    }
    
    if (TESTING) {
      printf("## Symmetrization matrices:\n");
      printf("## Nkpt_type = %d\n", Nkpt_type);
      for (int t=0; t<Nkpt_type; ++t) {
	printf("## type: %d  Nop: %d\n", t, gk_count[t]);
	printf("##");
	for (n=0; n<Nop; ++n) if (gk_equal[t][n]) printf(" g%d", n);
	printf("\n");
	for (i=0; i<9; ++i) {
	  printf("## %d%d =", i/3, i%3);
	  for (j=0; j<9; ++j) 
	    if (! zero(Sk[t][i*9+j]))
	      printf(" %+.5lf %d%d", Sk[t][i*9+j], j/3, j%3);
	  printf("\n");
	}
	printf("##\n");
      }
    }
    
    // garbage collection
    for (i=0; i<Nkpt_type; ++i)
      delete[] gk_equal[i];
    delete[] gk_equal;
    delete[] gk_count;
  }
  if (ERROR) exit(ERROR);

  // Now, we sort Gd and shell it; we use this to reduce the number of
  // Gd's we have to calculate.
  shell_type* sh;
  int Nsh;
  
  sort_points(Gd, Nlatt);
  ERROR = gen_shell_list_open(Gd, Nlatt, gunit, Nop, sh, Nsh);

  // Now we set kmax; so that when we do the integration using a
  // sphere inscribed inside our BZ, the result is exact.
  double kmax;
  
  kmax = max_k_BZ(cart_rlv);
  if (TESTING) {
    printf("## ++ Analytic spherical piece integration.\n");
    printf("## kmax = %.5lf ; volume fraction = %.5lf\n", 
	   kmax, 4*M_PI*kmax*kmax*kmax/(3.*det(cart_rlv)));
  }

  // Test to see if we've got enough k-points for our set of points:
  double mink, maxR;
  maxR = Gd[Nlatt-1].Rmagn;
  for (k=0, mink=kmax; k<Nkpt; ++k)
    if ( (mink > kpt[k][3]) && (! zero(kpt[k][3])) )
      mink = kpt[k][3];
  if (TESTING) {
    printf("## maxR = %.8le  mink = %.8le\n", maxR, mink);
    printf("## pi/maxR*mink = %.5lf\n", M_PI/(maxR*mink));
  }
  // Now, our criterion is that mink*maxR must be less than Pi/6:
  if ( (mink*maxR) > (M_PI/5.99) )
    fprintf(stderr, "WARNING: maxR = %.4le  mink = %.4le  pi/maxR*mink = %.5lf ; FT might not be accurate.\n", 
	    maxR, mink, M_PI/(maxR*mink));

  // value for the dc:
  double kmax_0;
  //  kmax_0 = 24*mink;
  kmax_0 = kmax*0.25;
  if (kmax_0 > kmax) kmax_0 = kmax;
  if (TESTING) {
    printf("## ++ Analytic spherical piece integration: discontinuity corr.\n");
    printf("## kmax_0 = %.5lf ; volume fraction = %.5lf\n", 
	   kmax_0, 4*M_PI*kmax_0*kmax_0*kmax_0/(3.*det(cart_rlv)));
  }  

  double ksmall, krealsmall;
  //  ksmall = kmax*0.5;     // the limit below which we don't use Ylm
  ksmall = kmax;     // the limit below which we don't use Ylm
  //  krealsmall = 0.25/maxR; // the limit below which we use expansion of k^2
  krealsmall = kmax*0.1; // the limit below which we really get "clever"
  if (krealsmall > (kpoint_envelope_alpha*kmax))
    krealsmall = (kpoint_envelope_alpha*kmax);
  if (krealsmall > (kpoint_envelope_alpha*kmax_0))
    krealsmall = (kpoint_envelope_alpha*kmax_0);
  
  // We no longer bother with lambda3 and lambda4... because we're
  // ridiculously careful about everything now.  We've got this calculation
  // in a bike-helmet with knee, elbow, and wrist pads.  Not to mention
  // the mylar vest.  It's that careful.
  double lambda[9][9], lambda2[9][9][9];
  //  double lambda3[9][9][9][9], lambda4[9][9][9][9][9];
  make_lambda(Np, Dij, lambda);
  make_lambda2(Np, Dij, lambda2);
  //  make_lambda3(Np, Dij, lambda3);
  //  make_lambda4(Np, Dij, lambda4);
  if (TESTING) {
    printf("## ksmall =     %.5lf :below which we don\'t use Ylm.\n", ksmall);
    printf("## krealsmall = %.5lf :below which we very carefully evaluated Gdk.\n", 
	   krealsmall);
    printf("## lambda: [ab,cd]: --non-zero only\n");
    print_lambda(lambda, "## ");
    printf("## lambda(2): [ab,cdef]: --non-zero only\n");
    print_lambda2(lambda2, "## ");
    //    printf("## lambda(3): [ab,cdefgh]: --non-zero only\n");
    //    print_lambda3(lambda3, "## ");
    //    printf("## lambda(4): [ab,cdefghij]: --non-zero only\n");
    //    print_lambda4(lambda4, "## ");
  }

  // Some workspace:
  // Spherical harmonics
  SH_workspace_type * work;
  // we'll use the same workspace for both, so pick the largest Lmax:
  if (Lmax >= Lmax_0) work = allocate_SH_workspace(Lmax);
  else                work = allocate_SH_workspace(Lmax_0);

  int np; // Loop over our points:

  //  double Gdk[Nkpt][9]; // the FT of our discrete correction
  double** Gdk; // the FT of our discrete correction
  Gdk = new double *[Nkpt];
  for (n=0; n<Nkpt; ++n)
    Gdk[n] = new double[9];

  if (TESTING) {
    printf("## Discretization correction FT: Gdk:\n");
  }
  for (k=0; k<Nkpt; ++k) {
    double* tkpt = kpt[k]; // this kpoint
    if (TESTING) {
      printf("## kpt: %8.5lf %8.5lf %8.5lf", tkpt[0], tkpt[1], tkpt[2]);
      if (SYMM) printf("  type: %d\n", kpt_type[k]);
      else printf("\n");
    }
    
    if (dcomp(tkpt[3], 0.)) {
      // Treat the origin as a special point:
      for (d=0; d<9; ++d) Gdk[k][d] = 0;
    }
    else {
      // **** WARNING!! ****
      // This code is now streamlined... so we only calculate those things
      // that we need.  If you look at older versions of the code, there
      // are versions where each and every piece is calculated, and only
      // the "correct" pieces are used.  So if you want to do some verbose
      // testing, you might want to pull that back out.

      if (tkpt[3] < krealsmall) {
	// NOTE: this is *exact* for all k where fenv=1... it's just that
	// it's really tuned to work for *small* k values...
	if (LOOKUP) smallk_expand(Np, Dij, tkpt, Gdk[k], lambda, lambda2,
				  KPT_CART, coskR);
	else        smallk_expand(Np, Dij, tkpt, Gdk[k], lambda, lambda2,
				  KPT_CART);
      }
      else {
	// 1. Fourier transform Dij...
	double Dk[9];
	if (LOOKUP) {
	  if (BINSUM) fourier_table_bin(Np, Dij, tkpt, Dk, KPT_CART,
					coskR, mbinp);
	  else        fourier_table(Np, Dij, tkpt, Dk, KPT_CART, coskR);
	}
	else        fourier(Np, Dij, tkpt, Dk, KPT_CART);
	double Dkinv[9], err;
	// Use our (more stable?) symmetric inverse:
	careful_inverse(Dk, Dkinv, err);
	if (err > 1e-13) {
	  fprintf(stderr, "Somehow managed to get an appreciable error in Dkinv: err= %.5le\n", err);
	  fprintf(stderr, "this most likely indicated that Dk:\n");
	  for (d=0; d<9; ++d) fprintf(stderr, " %.15le", Dk[d]);
	  fprintf(stderr, "\nis singular for kpt= %.12lf %.12lf %.12lf.\n",
		  tkpt[0], tkpt[1], tkpt[2]);
	  fprintf(stderr, "det(Dk) = %.15le", det(Dk));
	}
	// 2. Calculate Gsc(k) and Gdisc(k)
 	double fenv = kpoint_envelope(tkpt[3] / kmax);
 	double fenv_0 = kpoint_envelope(tkpt[3] / kmax_0);
	double Gsck[9], Gdisc[9];
	if (tkpt[3] > ksmall) {
	  // Ylm expansions for each:
	  double g_scale =  1./(cell_vol * tkpt[3] * tkpt[3]);
	  eval_Ylm_expansion_R(Lmax, tkpt, Ndim, RYlm, IYlm, Gsck, work);
	  mult(Gsck, g_scale*fenv, Gsck);

	  g_scale =  1./ cell_vol;
	  eval_Ylm_expansion_R(Lmax_0, tkpt, Ndim, RYlm_0, IYlm_0, Gdisc,work);
	  mult(Gdisc, g_scale*fenv_0, Gdisc);
	}
	else {
	  eval_lambda_inv(lambda, tkpt, Gsck);
	  mult(Gsck, fenv, Gsck);
	  eval_lambda2(lambda, lambda2, tkpt, Gdisc);
	  mult(Gdisc, fenv_0, Gdisc);
	}
	if (TESTING) {
	  printf("##   Dk:    "); print_mat(Dk); printf("\n");
	  printf("##   Dk^-1: "); print_mat(Dkinv); printf(" err = %.5le\n",err);
	  printf("##   Gsck:  "); print_mat(Gsck); printf("\n");
	  printf("##   Gdisc: "); print_mat(Gdisc); printf("\n");
	}

	// 3. Subtract:
	for (d=0; d<9; ++d)
	  Gdk[k][d] = Dkinv[d] - Gsck[d] - Gdisc[d];
      }
    }
    if (TESTING) {
      printf("##   Gdk:   "); print_mat(Gdk[k]); printf("\n");
    }
    if (SYMM) {
      // need to apply our symmetrization op:
      double Gi[9];
      int t=kpt_type[k];
      for (i=0; i<9; ++i) {
	Gi[i]=0.;
	for (j=0; j<9; ++j) 
	  Gi[i] += Sk[t][i*9+j]*Gdk[k][j];
      }
      for (d=0; d<9; ++d) Gdk[k][d] = Gi[d];
      if (TESTING) {
	printf("##   Gdk.s: ");	print_mat(Gdk[k]); printf("\n");
      }
    }
  }
  
  // Now all that remains is to *invert* the FT.  This takes advantage
  // of the symmetry that must have been built into our kpoint mesh
  // and the symmetry in Gd; we run through our shells one by one
  double ***RYlm_R, ***IYlm_R; // for handling the discontinuity
  double Fl_val[Lmax_0+1];
  double gint[9];
  ERROR = alloc_Ylm(Lmax_0, Ndim, RYlm_R, IYlm_R);
  if (ERROR) {
    fprintf(stderr, "Error allocating our Ylm expansion...\n");
    exit(ERROR);
  }
   
  for (int nsh=0; nsh<Nsh; ++nsh) {
    shell_type *tsh = sh + nsh;
    
    // first point in the shell:
    point_type *tGd = Gd + tsh->elem[0];
    if (TESTING) {
      printf("## point: %3d %3d %3d\n", 
	     tGd->Runit[0], tGd->Runit[1], tGd->Runit[2]);}

    // 1. inverse FT the discrete piece:
    if (LOOKUP) {
      if (BINSUM)
	inv_fourier_table_bin(Gd, tsh, tGd->mat, Nop, inv_index, gcart, 
			      Nkpt, Gdk, kpt, w, KPT_CART, coskR, mbinp);
      else 
	inv_fourier_table(Gd, tsh, tGd->mat, Nop, inv_index, gcart, 
			  Nkpt, Gdk, kpt, w, KPT_CART, coskR);
    }
    else
      inv_fourier(Gd, tsh, tGd->mat, Nop, inv_index, gcart, 
		  Nkpt, Gdk, kpt, w, KPT_CART);

    // 2. inverse FT the discontinuity correction
    // 2.a. calc the FT: 
    // we've pulled this algorithm from gf-disc;
    if (tGd->Rmagn > 0) {
      Fn_x2_int(Lmax_0, kmax_0*tGd->Rmagn, Fl_val, TESTING, 1024);
      // Fn_x2_int(Lmax_0, kmax_0*tGd->Rmagn, Fl_val, 0, 1024);
      if (TESTING) {
	printf("## x^2 int-- Fl_val:\n");
        for (l=0; l<=Lmax_0; ++l) {
          printf("## %d %.12le\n", l*2, Fl_val[l]);
        }
      }
      double g_scale = 1./(4*M_PI*(tGd->Rmagn*tGd->Rmagn*tGd->Rmagn)), l_scale;
      int l_sign;
      for (l=0, l_sign = 1; l<=Lmax_0; ++l, l_sign = -l_sign) {
        l_scale = g_scale * l_sign * Fl_val[l];
        for (m=0; m<=(4*l); ++m) {
          for (d=0; d<Ndim; ++d) {
            // Initialize
            RYlm_R[l][m][d] = RYlm_0[l][m][d] * l_scale;
            IYlm_R[l][m][d] = IYlm_0[l][m][d] * l_scale;
          }
        }
      }
      // Get the continuum piece:
      eval_Ylm_expansion_R (Lmax_0, tGd->Rcart, Ndim, RYlm_R, IYlm_R, gint,
			    work);
    }
    else {
      // For R = 0, the continuum piece is really simple:
      double g0_scale = M_2_SQRTPI/(8*M_PI*M_PI); // strange (exact) number
      g0_scale *= kmax_0*kmax_0*kmax_0*
	kpoint_x2_envelope_int();  // and our maximum k value...
      for (d=0; d<Ndim; ++d)
        gint[d] = g0_scale * RYlm_0[0][0][d];
    }
    if (TESTING) {
      printf("## point: %3d %3d %3d\n", 
	     tGd->Runit[0], tGd->Runit[1], tGd->Runit[2]);
      printf("## discrete FT: "); print_mat(tGd->mat); printf("\n");
      printf("## discont. GF: "); print_mat(gint); printf("\n");
    }
    // 2.b. add correction to the .mat:
    for (d=0; d<9; ++d) tGd->mat[d] += gint[d];
    if (TESTING) {
      printf("## final. GF:   "); print_mat(tGd->mat); printf("\n");
    }

    // 3. rotate to the rest of the points in the shell:
    for (j=1; j<tsh->Nelem; ++j) {
      np = tsh->elem[j];
      // find the group operation that takes R0 to our R:
      for (n=0; tsh->gR0[n] != np; ++n) ;
      double temp[9];
      mult(gcart[n], tGd->mat, temp);
      mult(temp, gcart[inv_index[n]], Gd[np].mat);
    }
  }

  // Finally, estimate the discretization length l_d:
  // ******* WARNING ********
  // we'll have to modify this later if we EVER do a non-cubic system
  // --there you'd want to take the ratio of the eigenvalues; in a
  // cubic system, this is trivial.
  double l_d = 0.;
  double gdisc_eigen = kpoint_x2_envelope_int() * RYlm_0[0][0][0];
  double gelas_eigen = kpoint_envelope_int() * RYlm[0][0][0];

  if (! dcomp(gelas_eigen, 0))
    l_d = sqrt(fabs(gdisc_eigen/gelas_eigen));

  // ****************************** OUTPUT ***************************

  // Now, we output EVERYTHING.
  if (VERBOSE) {
    // Human readable:
    for (np=0; np<Nlatt; ++np) {
      printf("# |R| = %8.5lf  Rhat = %8.5lf %8.5lf %8.5lf", 
	     Gd[np].Rmagn, Gd[np].Rcart[0], Gd[np].Rcart[1], Gd[np].Rcart[2]);
      printf("# "); print_mat(Gd[np].mat); printf("\n");
    }
  }

  printf("%d %.12le  # Number of points, kmax; ld = %.5le\n", Nlatt, kmax,l_d);
  for (np=0; np<Nlatt; ++np) {
    point_type *tp = Gd+np;
    printf("%3d %3d %3d", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    for (d=0; d<Ndim; ++d)
      printf(" %.12le", tp->mat[d]);
    printf("\n");
  }


  // ************************* GARBAGE COLLECTION ********************
  if (LOOKUP) free_coskR_table(coskR);
  if (BINSUM) free_mat_bin(mbinp);
  if (SYMM) {
    delete[] kpt_type;
    for (int t=0; t<Nkpt_type; ++t) delete[] Sk[t];
    delete[] Sk;
  }

  free_Ylm(Lmax, RYlm, IYlm);
  free_Ylm(Lmax_0, RYlm_0, IYlm_0);
  free_Ylm(Lmax_0, RYlm_R, IYlm_R);
  
  free_SH_workspace(work);

  free_kpts(Nkpt, kpt);

  for (n=0; n<Nkpt; ++n)
    delete[] Gdk[n];
  delete[] Gdk;

  free_shell(sh, Nsh);

  delete[] Dij;
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


// We use this to calculate Gsck for small kpoints.  We return
// (D(k))^-1 - (lambda(k))^-1 - (lambda(k))^-1(lambda2(k))(lambda(k))^-1
// calculated for small k.
// NOTE: We DON'T use an envelope function on either of the lambda
// matrices, so this can ONLY be used in kpoint regimes where the
// envelope functions would return 1--which is generally a good idea
// anyway...
//
// The algorithm works like this:
// We calculate X(k) which is sum_R (1-(k.R)^2-cos(k.R))D(R) by using
// fourier_smallk().  Then we calculate (lambda(k))^-1 and 
// lambda(k)^-1 lambda2(k) lambda(k)^-1 using eval_lambda_inv and
// eval_lambda2.  We then determine the eigenvalues of X and lambda_inv;
// define b = max(x_i)max(l_i); we make nseries = ln tol / ln b (tol = 1e-8?)
// This value determines if we use a series expansion, or if we directly
// do the inversion.  If nseries < NMAX, we expand out (1-B)^-1 as 
// B+B^2+B^3+... where B = lambda(k)^-1 X(k).  We subtract the lambda(2)
// piece from B*lambda^-1 before adding it, and left multiply everything
// by lambda^-1.
// If nseries > NMAX, we write A=1-B, and calculate (1-A^-1)lambda^-1.
// From that, we subtract the lambda(2) piece and return it.
//
// Whew.  Hopefully, it will all work.

const double SMALLK_TOL = 1e-8;
const double LOG_SMALLK_TOL = log(SMALLK_TOL);
const int SMALLK_NMAX = 10;

void smallk_expand(int Np, point_type* p, double kpt[3], double pk[9], 
		   double lambda[9][9], double lambda2[9][9][9],
		   const int TYPE, const coskR_table_type *table) 
{
  int i, d;
  double Xk[9], lambdak_inv[9], lambda2k[9];
  double B[9];
  int nseries;

  // 1. Calculate the matrices we need to work with:
  // call according to whether we have a lookup table or not...
  if (table == NULL) fourier_smallk(Np, p, kpt, Xk, TYPE);
  else               fourier_smallk_table(Np, p, kpt, Xk, TYPE, *table);
  eval_lambda_inv(lambda, kpt, lambdak_inv);
  eval_lambda2(lambda, lambda2, kpt, lambda2k);
  
  double xi[3], li[3], b;
  eigen(Xk, xi);
  eigen(lambdak_inv, li);
  // the eigenvalues are sorted, so the LARGEST ones are the last;
  // however, that doesn't consider the signs... so xk could be negative.
  // but the li eigenvalues are ALL positive--because we have a
  // stable lattice.
  if (fabs(xi[0]) > fabs(xi[2])) b = fabs(xi[0]);
  else                           b = fabs(xi[2]);
  b *= li[2];
  // Now, we take some logs:
  if (b != 0.0) nseries = lround(LOG_SMALLK_TOL / log(b));
  else nseries = 1;

  // Our B matrix:
  mult(lambdak_inv, Xk, B);

  double A[9], Ainv[9], scaleA;
  if (nseries > SMALLK_NMAX) {
    // do the inversion directly:
    for (d=0; d<9; ++d) A[d] = ident[d] - B[d];  // A = 1-B
    careful_inverse(A, Ainv);
    for (d=0; d<9; ++d) A[d] = ident[d]-Ainv[d]; // A = 1-(1-B)^1
    mult(A, lambdak_inv, pk);                    // pk = (1-(1-B)^1)lambda^-1
    // finally, subtract off the lambda2k piece:
    for (d=0; d<9; ++d) pk[d] -= lambda2k[d];
    // our work is done here
  }
  else {
    // else... do the inversion with a power series in B
    double tempA[9];
    for (d=0; d<9; ++d) tempA[d] = B[d]; // tempA = B
    for (int n=0; n<nseries; ++n) {
      for (d=0; d<9; ++d) A[d] = ident[d]+tempA[d]; // A = 1+tempA
      mult(B, A, tempA);                            // tempA = B(1+tempA)
    }
    // lastly, we'll multiply by B and lambda^-1... 
    mult(B, tempA, A);          // A = B^2(1+B(1+B(1+...)))
    mult(A, lambdak_inv, tempA); // tempA = B^2(1+B(1+B(1+...)))lambda^-1
    // ... and add B lambda^-1 - lambda2 piece:
    mult(B, lambdak_inv, A);     // A = B lambda^-1
    for (d=0; d<9; ++d) pk[d] = (A[d]-lambda2k[d]) + tempA[d];
  }
}

