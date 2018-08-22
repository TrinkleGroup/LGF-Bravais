/*
  Program: gf-polish.C
  Author:  D. Trinkle
  Date:    May 16, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. Ylm expansion for EGF
	   3. elastic GF discreziation correction
	   4. lattice GF discreziation correction (same set of points)
	   5. the dynamical matrix 

	   Through a very complicated algorithm, we attempt to iteratively
	   improve the GF.  There are two main ways to do this:
	   
	   1. Iterate using the self-consistent equation
	   2. Iterate using the quadratic polishing equation (2G-GDG)

  Param.:  <cell> <Gk(lm)> <GE.d(R)> <GL.d(R)> <Dij>
           cell:    cell file describing our lattice
           Gk(lm):  spherical harmonic components of elastic GF in k-space
	   GE.d(R): discrete correction to EGF
	   GL.d(R): discrete correction to LGF
	   Dij:     dynamical matrix

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

  Algo.:   

  Output:  
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
#include "Dij.H"
#include "cell.H"  // Do we still use this?
#include "cont-int.H" // For doing the integration over the continuous piece
#include "sphere-harm.H" // Spherical harmonic evaluation
#include "supercell.H" // supercell 
#include "ball.H"
#include "pointgroup.H"
#include "shell.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

inline double mat_error(const double a[9], const double b[9]) 
{
  double err=0;
  for (int d=0; d<9; ++d) err += fabs(a[d]-b[d]);
  return err;
}

inline double mat_zero(const double a[9]) 
{
  double err=0;
  for (int d=0; d<9; ++d) err += fabs(a[d]);
  return err;
}

inline void error_point(const int a[3], const char *loc) 
{
  fprintf(stderr, "In %s couldn\'t find a match for %3d %3d %3d\n", loc, a);
}

// Print out our 3x3 symmetric matrix elements
inline void print_mat (double a[9]) 
{
  printf("xx= %11.4le yy= %11.4le zz= %11.4le  xy= %11.4le yz= %11.4le zx= %11.4le",
	 a[0], a[4], a[8],
	 0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}


/*================================= main ==================================*/


// This is the default error tolerance before we kick out:
const double errorTOL = 1e-7;


// Arguments first, then flags, then explanation.
const int NUMARGS = 5;
const char* ARGLIST = "<cell> <Gk(lm)> <GE.d(R)> <GL.d(R)> <Dij>";


const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'q'}; 
// flag characters.

const char* ARGEXPL = 
"  cell:    cell file describing our lattice\n\
  Gk(lm):  spherical harmonic components of elastic GF in k-space\n\
  GE.d(R): discrete correction to EGF\n\
  GL.d(R): discrete correction to LGF\n\
  Dij:     dynamical matrix\n\
  MEMORY:  maximum number of iterations before we bail out\n\
  -q       use quadratic polishing (as opposed to self-consistent polishing)";

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
  int MEMORY = 1024;// Our default maximum number of iterations.

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
  int QUAD = flagon[0];
  if (VERBOSE) {
    if (QUAD) printf("# Quadratic polishing method.\n");
    else      printf("# Self-consistent equation polishing method.\n");
  }

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *Ylm_name, *GEd_name, *GLd_name, *Dij_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  Ylm_name = args[1];
  GEd_name = args[2];
  GLd_name = args[3];
  Dij_name = args[4];
  
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

  if (ERROR) {
    fprintf(stderr, "Something wrong with %s\n", Ylm_name);
    exit(ERROR); // time to bail out yet?
  }
  

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm, IYlm, "## ");
  }

  //++ ==== GEd ====
  int NEp;
  double kmax;
  point_type *GEd;
  
  infile = myopenr(GEd_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", GEd_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, NEp, GEd, READ_DIJ_MAT, dump);
  myclose(infile);
  sscanf(dump, "%*d %lf", &kmax); // %* to skip the first entry
  //-- ==== GEd ====
  if (ERROR) {
    fprintf(stderr, "Something wrong with %s\n", GEd_name);
    exit(ERROR); // time to bail out yet?
  }


  //++ ==== GLd ====
  int NLp;
  point_type *GLd;
  
  infile = myopenr(GLd_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", GLd_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, NLp, GLd, READ_DIJ_MAT, dump);
  myclose(infile);
  ERROR = (ERROR) || (NLp != NEp);
  if (ERROR) {
    fprintf(stderr, "Np= %d from %s but %d from %s\n", NEp,
	    GEd_name, NLp, GLd_name);
    fprintf(stderr, "files must match.\n");
  }
  if (!ERROR) {
    double kmaxL;
    sscanf(dump, "%*d %lf", &kmaxL); // %* to skip the first entry
    ERROR = (! dcomp(kmax, kmaxL));
    if (ERROR) {
      fprintf(stderr, "kmax = %.8le from %s but %.8le from %s\n", kmax,
	      GEd_name, kmaxL, GLd_name);
      fprintf(stderr, "files must match.\n");
    }
  }
  //-- ==== GLd ====
  if (ERROR) {
    fprintf(stderr, "Something wrong with %s\n", GLd_name);
    exit(ERROR); // time to bail out yet?
  }


  //++ ==== Dij ====
  int NDp;
  point_type *Dij;
  
  infile = myopenr(Dij_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Dij_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, NDp, Dij, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== Dij ====
  if (ERROR) {
    fprintf(stderr, "Something wrong with %s\n", Dij_name);
    exit(ERROR); // time to bail out yet?
  }


  // ***************************** ANALYSIS **************************
  double gcart[MAXop][9];
  int gunit[MAXop][9], inv_index[MAXop], Nop;
  // point group operations:
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  double RcutD, RcutGF;
  
  // 1. Sort and determine cutoffs:
  if (VERBOSE) {
    printf("# Sorting points, getting cutoffs...\n");
  }
  sort_points(Dij, NDp); RcutD = Dij[NDp-1].Rmagn;
  sort_points(GEd, NEp);
  sort_points(GLd, NLp); RcutGF = GLd[NLp-1].Rmagn;
  if (VERBOSE) {
    printf("# RcutD= %.5le  RcutGF= %.5le\n", RcutD, RcutGF);
  }
  
  
  // 2. Generate a ball out to our max. cutoff, and shell
  // 2.a. make the sphere...
  double Rcut; 
  int Ncut;
  int **Runit;
  // note: different Rcut's for different algorithms...
  Rcut = refine_Rcut(cart, RcutD+RcutGF);
  if (QUAD) {
    // use the larger cutoff for quadratic polishing:
    Rcut = refine_Rcut(cart, 2*RcutGF+RcutD);
    if ( (RcutGF+2*RcutD) > Rcut) 
      Rcut = refine_Rcut(cart, RcutGF+2*RcutD);
  }
  
  if (VERBOSE) {
    printf("# Making ball out to max cutoff %.5le\n", Rcut);
  }
  atom_sphere_unit(cart, Rcut, Ncut, Runit);
  // 2.b. shell it into representative points:
  point_type *GL;
  int Np;
  gen_representative_points(cart, Ncut, Runit, Nop, gunit, Np, GL);
  // 2.c. make an index into the representative point list
  symm_index_type*** index; // Our index of points
  int max_index;
  make_symm_index(Np, GL, Nop, gunit, max_index, index);
  if (VERBOSE) {
    printf("# %d representative points out to cutoff\n", Np);
  }
  if (TESTING) {
    for (np=0; np<Np; ++np)
      printf("## %3d : %3d %3d %3d | %10.5lf %10.5lf %10.5lf  %9.6lf\n",
	     np+1, GL[np].Runit[0], GL[np].Runit[1], GL[np].Runit[2],
	     GL[np].Rcart[0], GL[np].Rcart[1], GL[np].Rcart[2], 
	     GL[np].Rmagn);
  }
  

  // 3. fill in the "representative" versions of D and GL.
  if (VERBOSE) {
    printf("# Filling up D:\n");
  }
  point_type *Drep;
  int NDrep;
  // 3.a. find out how many points we need for D:
  for (NDrep=0; (GL[NDrep].Rmagn < (RcutD+TOLER)); ++NDrep) ;
  Drep = new point_type[NDrep];
  int jinit=0;
  for (i=0; i<NDrep; ++i) {
    // always make sure jinit points at the first member of the 
    // shell
    // CHANGED by D. Trinkle to *initialize* jinit=0
    //    for (; (Dij[jinit].Rmagn < GL[i].Rmagn); ++jinit) ;
    for (jinit=0; (Dij[jinit].Rmagn < (GL[i].Rmagn - TOLER)); ++jinit) ;
    
    for (d=0; d<3; ++d) Drep[i].Runit[d] = GL[i].Runit[d];
    for (d=0; d<3; ++d) Drep[i].Rcart[d] = GL[i].Rcart[d];
    Drep[i].Rmagn = GL[i].Rmagn;
    int found=0;
    // Changed: we may have cases where points in Drep exist
    // but aren't in Dij (since Drep is a sphere).  This requires
    // us to *pad out* entries in Drep accordingly.
    // So, we might *not* find our entry... if so, we make the
    // matrix elements zero.
    // NB: a smarter way to do this would be to just pull representative
    // elements from Dij itself by shelling it... but there might be
    // reasons in the code below to not do it this way.
    for (j=jinit; (! found) && (j<NDp); ++j) {
      if (dcomp(Dij[j].Rmagn, Drep[i].Rmagn)) 
	found = equal_vect(Dij[j].Runit, Drep[i].Runit);
    }
    --j;
    // Now j is the index for the corresponding point:
    if (found)
      for (d=0; d<9; ++d) Drep[i].mat[d] = Dij[j].mat[d];
    else
      // pad entries
      for (d=0; d<9; ++d) Drep[i].mat[d] = 0;
  }
  if (TESTING) {
    printf("## representative elements for Dij:\n");
    printf("## %d points\n", NDrep);
    for (np=0; np<NDrep; ++np) {
      printf("## %3d %3d %3d ", Drep[np].Runit[0], Drep[np].Runit[1], 
	     Drep[np].Runit[2]);
      print_mat(Drep[np].mat);
      printf("\n");
    }
  }

  // 3.b. work through the points for GE...
  if (VERBOSE) {
    printf("# Filling up GE:\n");
  }
  double ***RYlm_R, ***IYlm_R;
  double Fl_val, g_scale = 1./(4*M_PI); // All that'll be left is 1./Rmagn
  SH_workspace_type *w = allocate_SH_workspace(Lmax);
  alloc_Ylm(Lmax, 9, RYlm_R, IYlm_R);
  Fl_val = 1.; // This will be prod_i=(1..l/2) (- (2i-1)/2i)
  for (l=0; l<=Lmax; ++l) {
    for (m=0; m<=(4*l); ++m) 
      for (d=0; d<Ndim; ++d) {
        // Initialize
        RYlm_R[l][m][d] = RYlm[l][m][d] * g_scale * Fl_val;
        IYlm_R[l][m][d] = IYlm[l][m][d] * g_scale * Fl_val;
      }
    // Update Fl_val:
    Fl_val *= - (double)(2*l+1) / (double)(2*l+2);
  }

  // stolen from real-GF:
  double Gsc0[9];
  for (np=0; np<Np; ++np) {
    int *Runit = GL[np].Runit;
    double *Rcart = GL[np].Rcart;
    double Rmagn = GL[np].Rmagn;
    double *mat = GL[np].mat;
    
    if ( zero(Rmagn) ) {
      // Do the zero point using integration... VERY approximate,
      // but better than setting it to 0, I've found.  Taken from gf-disc.C
      // kmax we already know...
      double g0_scale = M_2_SQRTPI/(8*M_PI*M_PI); // strange (exact) number
      g0_scale *= kmax * kpoint_envelope_int();  // and our maximum k value...
      for (d=0; d<Ndim; ++d)
        mat[d] = g0_scale * RYlm[0][0][d];
      for (d=0; d<Ndim; ++d) Gsc0[d] = mat[d]; // store this for later.
    }
    else {
      // spherical harmonic:
      eval_Ylm_expansion_R (Lmax, Rcart, 9, RYlm_R, IYlm_R, mat, w);
      mult(mat, 1./Rmagn, mat);
    }
  }
  // garbage collection
  free_Ylm(Lmax, RYlm_R, IYlm_R);
  free_SH_workspace(w);
  
  // 3.c. ... then work through GL
  if (VERBOSE) {
    printf("# Filling up GL:\n");
  }
  int NGrep, NDG;
  point_type *Gsc, *DG, *GLnew;
  for (NGrep=0; (GL[NGrep].Rmagn < (RcutGF+TOLER)); ++NGrep) ;
  Gsc = new point_type[NGrep]; // also put these elements in...
  GLnew = new point_type[NGrep];
  if (QUAD) {
    // for quad, we need more DG points than we do with self-consist.
    for (NDG=NGrep; (GL[NDG].Rmagn < (RcutGF+RcutD+TOLER)); ++NDG) ;
  }
  else NDG = NGrep;
  DG = new point_type[NDG]; // also put these elements in...
  if (VERBOSE) {
    printf("# GL uses %d points.\n", NGrep);
    printf("# DG uses %d points.\n", NDG);
  }
  jinit=0;
  for (i=0; i<NGrep; ++i) {
    // always make sure jinit points at the first member of the 
    // shell
    // CHANGED by D. Trinkle to *initialize* jinit=0
    //    for (; (GLd[jinit].Rmagn < GL[i].Rmagn); ++jinit) ;
    for (jinit=0; (GLd[jinit].Rmagn < (GL[i].Rmagn-TOLER)); ++jinit) ;

    for (d=0; d<3; ++d) DG[i].Runit[d] = GL[i].Runit[d];
    for (d=0; d<3; ++d) DG[i].Rcart[d] = GL[i].Rcart[d];
    DG[i].Rmagn = GL[i].Rmagn;

    for (d=0; d<3; ++d) GLnew[i].Runit[d] = GL[i].Runit[d];
    for (d=0; d<3; ++d) GLnew[i].Rcart[d] = GL[i].Rcart[d];
    GLnew[i].Rmagn = GL[i].Rmagn;
    
    int found=0;
    for (j=jinit; (! found); ++j) {
      if (dcomp(GLd[j].Rmagn, GL[i].Rmagn)) 
	found = equal_vect(GLd[j].Runit, GL[i].Runit);
    }
    --j;
    // Now j is the index for the corresponding point:
    if (! equal_vect(GLd[j].Runit, GEd[j].Runit)) {
      fprintf(stderr, "ERROR while trying to reconstruct GL using GEd.\n");
      fprintf(stderr, "GLd point %3d %3d %3d doesn\'t match GEd %3d %3d %3d\n",
	      GLd[j].Runit[0], GLd[j].Runit[1], GLd[j].Runit[2],
	      GEd[j].Runit[0], GEd[j].Runit[1], GEd[j].Runit[2]);
      fprintf(stderr, "Probably means that the two sets aren't over the same set of points.\n");
      ERROR=-1;
    }
    for (d=0; d<9; ++d) {
      Gsc[i].mat[d] = GL[i].mat[d]-GEd[j].mat[d];  // semicontinuum at i
      GL[i].mat[d] += GLd[j].mat[d]-GEd[j].mat[d]; // full LGF at i
    }
  }
  // if QUAD, we have to fill in a few more points:
  for (i=NGrep; i<NDG; ++i) {
    for (d=0; d<3; ++d) DG[i].Runit[d] = GL[i].Runit[d];
    for (d=0; d<3; ++d) DG[i].Rcart[d] = GL[i].Rcart[d];
    DG[i].Rmagn = GL[i].Rmagn;
  }
  if (TESTING) {
    printf("## representative elements for GL:\n");
    printf("## %d points\n", Np);
    for (np=0; np<Np; ++np) {
      point_type *tp = GL+np;
      printf("## %3d %3d %3d ", tp->Runit[0], tp->Runit[1], 
	     tp->Runit[2]);
      if (np<NGrep) printf("L: ");
      else          printf("E: ");
      print_mat(tp->mat);
      printf("\n");
    }
  }
  // finally, fix the origin:
  for (d=0; d<9; ++d) Gsc[0].mat[d] = Gsc0[d];


  // 4. Now, it's time to get our iteration *ON*.
  // 4.a. First step is to fix the value at zero, before we bother
  //      with anything else, because it's the one most likely to be
  //      screwed up (because of the way GEd is defined)
  // we iterate over the points in Dij, but *skip* the origin:
  int elem[Nop][3], Nelem; // we use this to make sure we don't double count
  double G0[9] = {0,0,0, 0,0,0, 0,0,0}; // initialized, too.
  double Drot[9], Grot[9], DG0[9];
  if (VERBOSE) {
    printf("# Correcting R=0 value.\n");
  }
  for (np=1; np<NDrep; ++np) {
    Nelem=0;
    for (n=0; n<Nop; ++n) {
      int gx[3];
      mult_vect(gunit[n], Drep[np].Runit, gx);
      for (j=0; (j<Nelem) && (! equal_vect(gx, elem[j])); ++j) ;
      if (j == Nelem) {
	for (d=0; d<3; ++d) elem[Nelem][d] = gx[d];
	++Nelem;
	// we haven't added this point's contribution yet...
	symm_index_type *ti = get_symm_index(max_index, index, gx);
	if (ti == NULL) {
	  error_point(gx, "G(0) refinement");
	  ERROR = -1;
	}
	m = ti->op;
	// rotate D to be correct:
	mult(gcart[m], Drep[ti->index].mat, gcart[inv_index[m]], Drot);
	// rotate G to be correct: (easy, since we're using the same point)
	mult(gcart[m], GL[ti->index].mat, gcart[inv_index[m]], Grot);
	// multiply and add
	mult(Drot, Grot, DG0);
	for (d=0; d<9; ++d) G0[d] += DG0[d];
	if (TESTING) {
	  printf("## point %3d %3d %3d\n", gx[0], gx[1], gx[2]);
	  printf("##   Dij: "); print_mat(Drot); printf("\n");
	  printf("##   Gij: "); print_mat(Grot); printf("\n");
	  printf("##   DG0: "); print_mat(DG0); printf("\n");
	}
      }
    }
  }
  // negate, and add unity:
  for (d=0; d<9; ++d) G0[d] = ident[d] - G0[d];
  double newG[9], invD0[9];
  careful_inverse(Drep[0].mat, invD0);
  mult(invD0, G0, newG);
  if (VERBOSE) {
    printf("# Old GL(0): "); print_mat(GL[0].mat); printf("\n");
    printf("# New GL(0): "); print_mat(newG); printf("\n");
    printf("# error: %.5le\n", mat_error(GL[0].mat, newG));
  }
  mult(newG, 1, GL[0].mat);

  // 4.b. Let the iteration begin!
  int MAXiter = MEMORY;
  int iter, iter_print;
  double err = 2*errorTOL;
  iter_print = 10; // one day, we'll want to change this
  if (TESTING) {
    printf("## **** Iterations beginning.\n");
  }
  for (iter=0; (iter<MAXiter) && (err > errorTOL); ++iter) {
    err = 0.;
    // iterate over R; go out as far as NDG requires:
    for (int nGp=0; nGp<NDG; ++nGp) {
      int *Runit = DG[nGp].Runit;
      double *mat = DG[nGp].mat;
      for (d=0; d<9; ++d) newG[d] = 0.; // zero out
      
      if (TESTING) {
	printf("## ** iterating GL( %3d %3d %3d )\n", 
	       Runit[0], Runit[1], Runit[2]);
      }
      
      // Let's do our thing.
      // iterate over x:
      for (np=1; np<NDrep; ++np) {
	Nelem=0;
	for (n=0; n<Nop; ++n) {
	  int gx[3], Rx[3];
	  mult_vect(gunit[n], Drep[np].Runit, gx);
	  for (j=0; (j<Nelem) && (! equal_vect(gx, elem[j])); ++j) ;
	  if (j == Nelem) {
	    for (d=0; d<3; ++d) elem[Nelem][d] = gx[d];
	    ++Nelem;
	    // we haven't added this point's contribution yet...
	    symm_index_type *ti = get_symm_index(max_index, index, gx);
	    if (ti == NULL) {
	      error_point(gx, "DG calc");
	      ERROR = -1;
	    }
	    m = ti->op;
	    // rotate D to be correct:
	    mult(gcart[m], Drep[ti->index].mat, gcart[inv_index[m]], Drot);
	    for (d=0; d<3; ++d) Rx[d] = Runit[d]+gx[d];
	    ti = get_symm_index(max_index, index, Rx);
	    m = ti->op;
	    // rotate G to be correct:
	    mult(gcart[m], GL[ti->index].mat, gcart[inv_index[m]], Grot);
	    // multiply and add
	    mult(Drot, Grot, DG0);
	    for (d=0; d<9; ++d) newG[d] += DG0[d];
	    if (TESTING) {
	      printf("## point D( %3d %3d %3d ) GL( %3d %3d %3d )\n", 
		     gx[0], gx[1], gx[2], Rx[0], Rx[1], Rx[2]);
	      printf("##   Dij: "); print_mat(Drot); printf("\n");
	      printf("##   Gij: "); print_mat(Grot); printf("\n");
	      printf("##   DG0: "); print_mat(DG0); printf("\n");
	    }
	  }
	}
      }
      if (QUAD) {
	// fix DG--stored in "mat"--so that it's sum_y D(y) G(x+y) for 
	// each x:
	mult(Drep[0].mat, GL[nGp].mat, mat);   // mat = D(0)G(R)
	for (d=0; d<9; ++d) mat[d] += newG[d]; // mat = sum_y D(y) G(x+y)
      }

      // Now, negate and (possibly) shift by unity, then multiply:
      for (d=0; d<9; ++d) newG[d] = -newG[d];
      if (nGp == 0)
	for (d=0; d<9; ++d) newG[d] += ident[d];
      // Calculate the true error:
      double DGfull[9];
      mult(Drep[0].mat, GL[nGp].mat, DGfull); // D(0)G(R)
      for (d=0; d<9; ++d) DGfull[d] = newG[d] - DGfull[d];
      // DGfull = 1delta(R) - sum_x D(x)G(R+x), so it SHOULD be zero:
      double thiserr = mat_zero(DGfull); // accumulate error
      err += thiserr;

      if (!QUAD) {
	// only bother with this if we're not doing quadratic polishing
	mult(invD0, newG, GLnew[nGp].mat); // this is our new value for SC
	if (TESTING) {
	  printf("## R: %3d %3d %3d\n", Runit[0], Runit[1], Runit[2]);
	  printf("## Old GL(R): "); print_mat(GL[nGp].mat); printf("\n");
	  printf("## New GL(R): "); print_mat(GLnew[nGp].mat); printf("\n");
	  printf("## error: %.5le\n", thiserr);
	  printf("##\n");
	}
      }
    }
    // scale error by number of points (the fair thing to do)
    err *= 1./NGrep;
    
    if (QUAD) {
      // We actually have MORE work to do for quadratic convergence...
      // In reality, GLnew(R) = 2*G(R) - sum_x G(R+x)DG(x), so we just
      // do a similar type of sum as before, but NOW over G and DG.
      // Also it runs R out to RGF_c (i.e., NGrep) and x runs out to
      // RGF_c+RD_c (i.e., NDG):
      
      if (TESTING) {
	printf("## **** quadratic step\n");
      }
      for (int nGp=0; nGp<NGrep; ++nGp) {
	int *Runit = GLnew[nGp].Runit;
	double *mat = GLnew[nGp].mat;
	for (d=0; d<9; ++d) newG[d] = 0.; // zero out
	
	if (TESTING) {
	  printf("## ** iterating GL( %3d %3d %3d )\n", 
		 Runit[0], Runit[1], Runit[2]);
	}
	
	// Let's do our thing.
	// iterate over ALL x--now runs over RGF_c+RD_c
	for (np=0; np<NDG; ++np) {
	  Nelem=0;
	  for (n=0; n<Nop; ++n) {
	    int gx[3], Rx[3];
	    mult_vect(gunit[n], DG[np].Runit, gx);
	    for (j=0; (j<Nelem) && (! equal_vect(gx, elem[j])); ++j) ;
	    if (j == Nelem) {
	      for (d=0; d<3; ++d) elem[Nelem][d] = gx[d];
	      ++Nelem;
	      // we haven't added this point's contribution yet...
	      symm_index_type *ti = get_symm_index(max_index, index, gx);
	      if (ti == NULL) {
		error_point(gx, "quadratic step");
		ERROR = -1;
	      }
	      m = ti->op;
	      // rotate DG to be correct:
	      mult(gcart[m], DG[ti->index].mat, gcart[inv_index[m]], Drot);
	      for (d=0; d<3; ++d) Rx[d] = Runit[d]+gx[d];
	      ti = get_symm_index(max_index, index, Rx);
	      m = ti->op;
	      // rotate G to be correct:
	      mult(gcart[m], GL[ti->index].mat, gcart[inv_index[m]], Grot);
	      // multiply and add
	      mult(Grot, Drot, DG0); // now GDG
	      for (d=0; d<9; ++d) newG[d] += DG0[d];
	      if (TESTING) {
		printf("## point DG( %3d %3d %3d ) GL( %3d %3d %3d )\n", 
		       gx[0], gx[1], gx[2], Rx[0], Rx[1], Rx[2]);
		printf("##   Gij:  "); print_mat(Grot); printf("\n");
		printf("##   DGij: "); print_mat(Drot); printf("\n");
		printf("##   GDG0: "); print_mat(DG0); printf("\n");
	      }
	    }
	  }
	}
	// Now, negate and add 2G: (mat=GLnew[nGp].mat)
	for (d=0; d<9; ++d) mat[d] = 2*GL[nGp].mat[d] - newG[d];
	if (TESTING) {
	  printf("## R: %3d %3d %3d\n", Runit[0], Runit[1], Runit[2]);
	  printf("## Old GL(R): "); print_mat(GL[nGp].mat); printf("\n");
	  printf("## New GL(R): "); print_mat(GLnew[nGp].mat); printf("\n");
	}
      }
    }
    

    // Update:
    for (int nGp=0; nGp<NGrep; ++nGp)
      for (d=0; d<9; ++d) GL[nGp].mat[d] = GLnew[nGp].mat[d];
    if (VERBOSE) {
      printf("# Iteration %d total accumulated error= %.5le\n", iter+1, err);
    }
    // now so we get some kind of error reporting:
    if ( ((iter+1)%iter_print) == 0 ){
      fprintf(stderr, "Iteration %d  error= %.5le\n", iter+1, err);
    }
  }
  
  fprintf(stderr, "Final iteration %d  error= %.5le\n", iter+1, err);

  // ****************************** OUTPUT ***************************

  printf("%d %.15le  # number of points, kmax; %d iterations, error= %.5le\n",
	 NLp, kmax, iter+1, err);
  // Now run through the points and output:
  for (int nGp=0; nGp<NGrep; ++nGp) {
    int *Runit = GL[nGp].Runit;
    double mat[9]; // will be the discrete correction!
    for (d=0; d<9; ++d) mat[d] = GL[nGp].mat[d]-Gsc[nGp].mat[d];
    
    // Let's do our thing.
    // iterate over x:
    Nelem=0;
    for (n=0; n<Nop; ++n) {
      int gR[3];
      double Grot[9];
      mult_vect(gunit[n], Runit, gR);
      for (j=0; (j<Nelem) && (! equal_vect(gR, elem[j])); ++j) ;
      if (j == Nelem) {
	// add it to our list:
	for (d=0; d<3; ++d) elem[Nelem][d] = gR[d];
	++Nelem;
	// rotate the matrix
	mult(gcart[n], mat, gcart[inv_index[n]], Grot);
	// print it out:
	printf("%3d %3d %3d", gR[0], gR[1], gR[2]);
	for (d=0; d<9; ++d) printf(" %.15le", Grot[d]);
	printf("\n");
      }
    }
  }
  
  // ************************* GARBAGE COLLECTION ********************

  delete[] Gsc;
  delete[] DG;
  delete[] GLnew;
  free_symm_index(max_index, index);
  delete[] Drep;
  delete[] GL;
  
  delete[] GEd;
  delete[] GLd;
  delete[] Dij;
  
  // spherical harmonic components
  free_Ylm(Lmax, RYlm, IYlm);

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
