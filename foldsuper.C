/*
  Program: foldsuper.C
  Author:  D. Trinkle
  Date:    April 12, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. lattice function defined at a set of points
	   3. a periodic supercell

	   and we compute
	   1. periodic version of lattice function evaluated at points
	      inside of the supercell only
	   2. if symmetrization is turned *on*, we enforce the point
	      group symmetry correctly; this produces *more* points
	      in the folddown operation.

  Param.:  <cell> <latt(R))> <supercell>
           cell:    cell file describing our lattice
	   latt(R): matrix function evaluated at a set of points
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

	   ==== latt(R) ====
	   N kmax                      # number of points, kmax
	   n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij
	   ...
	   nN.1 nN.2 nN.3  Lxx .. Lzz
	   ==== latt(R) ====

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We just evaluate sum_S mat(S+R) for each R, where S runs over
           the supercell.

  Output:  Going to output a series of lattice points inside the supercell, 
           and the periodic GF evaluated at those points.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "cell.H"  // Do we still use this?
#include "supercell.H" // supercell 
#include "Dij.H"
#include "voigt.H"
#include "pointgroup.H"
#include "shell.H"
#include <gsl/gsl_linalg.h> // needed for inverting Aji in symmetrization

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************
int gen_symmrel(point_type* p, int super[9], int gunit[MAXop][9], int Nop,
		int** symmrel);

void solve_symm(point_type *p, double g_i_ginv[MAXop][Ni][Ni], 
		int** symmrel, int Nsymm);

inline void mult3 (const double a[9], const double b[3], const double c[3],
		   double result[9]) 
{
  double temp[9];
  mult(a, b, temp);
  mult(temp, c, result);
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
const char* ARGLIST = "<cell> <latt(R)> <supercell>";


const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'s'}; // symmetrization?
// flag characters.

const char* ARGEXPL = 
"  cell:      cell file describing our lattice\n\
  latt(R):   matrix function evaluated at a set of points\n\
  supercell: supercell geometry\n\
  -s:        enforce symmetrization (usually produces more points)";

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
==== latt(R) ====\n\
N kmax                      # number of points, kmax\n\
n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij\n\
...\n\
nN.1 nN.2 nN.3  Lxx .. lzz\n\
==== latt(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, l, m, n, np; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter
  int SYMMETRIZE = 0;

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
  SYMMETRIZE = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *latt_name, *super_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  latt_name = args[1];
  super_name = args[2];
  
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

  //++ ==== latt(R) ====
  int Np;
  point_type *latt;
  
  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, latt, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== latt(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", latt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (n=0; n<Np; ++n) {
      point_type *tp = latt + n;
      printf("## %3d %3d %3d ", 
             tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      print_mat(tp->mat);
      printf("\n");
    }
  }


  //++ ==== supercell ====
  int super_n[9];
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
  point_type *period;
  int Nsuper;
  fill_super(cart, super_n, Nsuper, period);
  if (TESTING) {
    printf("## Supercell geometry (%d atoms):\n", Nsuper);
    for (n=0; n<Nsuper; ++n) {
      printf("## %3d : %3d %3d %3d\n", n+1, 
	     period[n].Runit[0], period[n].Runit[1], period[n].Runit[2]);
    }
  }

  // Make our list of supercell translations we need to sum up Gd:
  int*** index, max_index; // Our index of points
  int** shiftlist, Nshift; // the shifts we need

  // make our index
  make_index(Np, latt, max_index, index);
  // make our list of shifts
  make_super_shifts(max_index, index, super_n, Nshift, shiftlist);
  if (TESTING) {
    printf("## Indexing array: max_index = %d\n", max_index);
    int testp[3] = {0,0,0}, ind;
    ind = get_index(max_index, index, testp);
    printf("## 000 index = %d  Runit[] = %3d %3d %3d\n",
	   ind, latt[ind].Runit[0], latt[ind].Runit[1], latt[ind].Runit[2]);
    printf("## checking index...\n");
    for (np=0; np<Np; ++np)
      if (get_index(max_index, index, latt[np].Runit) != np) {
	printf("## FAILURE for point %d\n", np);
      }
    printf("## ...check complete.\n");
    printf("## Nshifts = %d\n", Nshift);
    for (n=0; n<Nshift; ++n)
      printf("## %3d : %3d %3d %3d\n", n+1,
	     shiftlist[n][0], shiftlist[n][1], shiftlist[n][2]);
  }
  
  // NOW!! The real analysis can begin.
  for (np=0; np<Nsuper; ++np) {
    point_type *tp = period + np;
    // zero out Gperiod:
    for (d=0; d<9; ++d) tp->mat[d] = 0.;

    for (i=0; i<Nshift; ++i) {
      int ind, point[3];
      for (k=0; k<3; ++k) point[k] = tp->Runit[k] + shiftlist[i][k];
      ind = get_index(max_index, index, point);
      if (ind != -1) {
	// Add contribution
	for (d=0; d<9; ++d) tp->mat[d] += latt[ind].mat[d];
      }
      if (TESTING) {
	printf("##  mat( %3d %3d %3d ) =\n", 
	       tp->Runit[0], tp->Runit[1], tp->Runit[2]);
	printf("##      ");
	print_mat(tp->mat);
	printf("\n");
      }
    }
  }

  // This is the MESSY bit... designed to enforce point group symmetry
  // on our set of points.  This can actually be used "post mortem"
  // on a set of points that were folded down already *without*
  // symmetrization turned on.
  if (SYMMETRIZE) {
    // Point group analysis--only used if SYMMETRIZE is turned on.
    double gcart[MAXop][9];
    int gunit[MAXop][9];
    int inv_index[MAXop];
    int Nop;
    double g_i_ginv[MAXop][Ni][Ni];
    shell_type* sh = NULL;
    int Nsh;
    // 1. point group operations:
    Nop = gen_pointgroup(cart, gcart, gunit, inv_index);
    // 2. g [i] ginv, where [i] is a voigt matrix
    for (n=0; n<Nop; ++n) {
      for (i=0; i<Ni; ++i) {
	gen_g_i_ginv(i, gcart[n], gcart[inv_index[n]], g_i_ginv[n][i]);
      }
    }
    // 3. make open shell list
    sort_points(period, Nsuper);
    ERROR = gen_shell_list_open(period, Nsuper, gunit, Nop, sh, Nsh);

    if (TESTING) {
      printf("## g[i]g^-1:\n");
      for (n=0; n<Nop; ++n) {
	printf("## g(%d)[ ]g(%d)^-1 = (", n, n);
	for (i=0; i<Ni; ++i) {
	  for (j=0; j<Ni; ++j)
	    printf(" %.2le", g_i_ginv[n][i][j]);
	  printf(" |");
	}
	printf("\n");
      }
    }
    
    // Now, we need to estimate how many new points are needed:
    int Nnew = 0; // count up how many open shell points there are
    int opensh[Nsh];
    point_type *pnew;
    for (n=0; n<Nsh; ++n) {
      opensh[n] = 0;
      for (i=0; i<Nop; ++i)
	if (sh[n].gR0[i] == -1) {
	  ++Nnew;
	  opensh[n] = -1;
	}
    }
    pnew = new point_type[Nnew];
    
    // Now, let's make 'em one by one...
    int **symmrel = new int *[Nop]; // symm. related points
    int Nsymm;
    for (n=0; n<Nop; ++n) symmrel[n] = new int[4];
    
    Nnew = 0;
    for (n=0; n<Nsh; ++n) 
      if (opensh[n]) {
	shell_type *tsh = sh + n;
	for (np=0; np<tsh->Nelem; ++np) {
	  point_type *tp = period + (tsh->elem[np]);
	  Nsymm = gen_symmrel(tp, super_n, gunit, Nop, symmrel);
	  // solve the problem on the initial point in place:
	  solve_symm(tp, g_i_ginv, symmrel, Nsymm);
	  // Now, we need to add the "extra" points:
	  for (i=0; i<Nsymm; ++i) {
	    // don't bother adding our initial point...
	    if (! equal_vect(symmrel[i], tp->Runit)) {
	      // 1. position information:
	      for (d=0; d<3; ++d) pnew[Nnew].Runit[d] = symmrel[i][d];
	      mult_vect(cart, symmrel[i], pnew[Nnew].Rcart);
	      pnew[Nnew].Rmagn = sqrt(dot(pnew[Nnew].Rcart, pnew[Nnew].Rcart));
	      // 2. transform matrix:
	      mult3(gcart[symmrel[i][3]], tp->mat,
		    gcart[inv_index[symmrel[i][3]]], pnew[Nnew].mat);
	      ++Nnew;
	    }
	  }
	}
      }
    // NOW, we're slick... we merge the two lists:
    point_type* pcombine;
    pcombine = new point_type[Nsuper+Nnew];
    for (n=0; n<Nsuper; ++n) pcombine[n] = period[n];
    for (n=0; n<Nnew; ++n) pcombine[n+Nsuper] = pnew[n];
    // remove our old separate lists:
    delete[] period;
    delete[] pnew;
    period = pcombine;
    Nsuper += Nnew;
    sort_points(period, Nsuper);

    // garbage collection
    for (n=0; n<Nop; ++n) delete[] symmrel[n];
    delete[] symmrel;

    delete[] sh;
  }
  
  // ****************************** OUTPUT ***************************
  //  for (d=0; d<9; ++d)
  //    printf("%3d ", super_n[d]);
  //  printf("# supercell definition\n");
  printf("%d 9 # Number of points", Nsuper);
  if (SYMMETRIZE) printf(" symmetrized\n");
  else printf("\n");
  for (np=0; np<Nsuper; ++np) {
    point_type *tp = period + np;
    printf("%3d %3d %3d", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    for (d=0; d<9; ++d)
      printf(" %.12le", tp->mat[d]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  free_index(max_index, index);
  free_super(Nshift, shiftlist);

  delete[] period;
  delete[] latt;
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


// construct the list of points equivalent to p by *both* a point group
// operation as well as a supercell lattice vector... woohoo.
int gen_symmrel(point_type* p, int super[9], int gunit[MAXop][9], int Nop,
		int** symmrel) 
{
  int n, d, i;
  int Nsymm;
  int prot[3], pdiff[3], pdiff_inv[3];
  int super_inv[9], detsuper;
  detsuper = inverse(super, super_inv);
  
  Nsymm = 0;
  for (n=0; n<Nop; ++n) {
    mult_vect(gunit[n], p->Runit, prot);
    for (d=0; d<3; ++d) pdiff[d] = p->Runit[d] - prot[d];
    // Now, see if the difference is a supercell vector:
    mult_vect(super_inv, pdiff, pdiff_inv);
    if ( ((pdiff_inv[0] % detsuper) == 0) &&
	 ((pdiff_inv[1] % detsuper) == 0) &&
	 ((pdiff_inv[2] % detsuper) == 0) ) {
      // Now, we have to see if it's in our list...
      for (i=0; (i<Nsymm) && (! equal_vect(symmrel[i], prot)); ++i) ;
      if (i == Nsymm) {
	// New point!
	for (d=0; d<3; ++d) symmrel[Nsymm][d] = prot[d];
	symmrel[Nsymm][3] = n;
	++Nsymm;
      }
    }
  }
  return Nsymm;
}


// this constructs the linear problem to solve for a given set of
// symmetry related points, and solves it using the voigt construction
void solve_symm(point_type *p, double g_i_ginv[MAXop][Ni][Ni], 
		int** symmrel, int Nsymm) 
{
  int i, j, n;
  double* Amat = new double[Ni*Ni];
  double* bvect = new double[Ni];
  
  // construct our A_ij = sum_R g_i_ginv[g(R)][j][i]
  for (i=0; i<(Ni*Ni); ++i) Amat[i] = 0;
  for (n=0; n<Nsymm; ++n) {
    for (i=0; i<Ni; ++i)
      for (j=0; j<Ni; ++j) 
	Amat[i*Ni+j] += g_i_ginv[symmrel[n][3]][j][i];
  }
  for (i=0; i<Ni; ++i) bvect[i] = p->mat[i6tom9[i]];
  
  // Now, solve our equation.  We hav to use the SVD because
  // A is most likely singular due to symmetry... *sigh*
  gsl_matrix_view m = gsl_matrix_view_array(Amat, Ni, Ni);
  gsl_vector_view b = gsl_vector_view_array(bvect, Ni);
  gsl_vector *x = gsl_vector_alloc(Ni);
  // SVD needs:
  gsl_matrix *V = gsl_matrix_alloc(Ni, Ni);
  gsl_vector *S = gsl_vector_alloc(Ni);
  gsl_vector *w = gsl_vector_alloc(Ni);
  
  gsl_linalg_SV_decomp (&m.matrix, V, S, w);
  // Now, our matrix may be singular, but that's okay--it just
  // means that the symmetry is already built-in.  We use our
  // tolerance to make these values 0:
  for (i=0; i<Ni; ++i)
    if (dcomp( gsl_vector_get(S, i), 0.))
      gsl_vector_set(S, i, 0.);

  gsl_linalg_SV_solve (&m.matrix, V, S, &b.vector, x);

  // Now, put in the entries!
  for (i=0; i<9; ++i) 
    p->mat[i] = gsl_vector_get(x, m9toi6[i]);
  // garbage collection 
  gsl_vector_free(x);
  gsl_vector_free(S);
  gsl_matrix_free(V);
  delete[] Amat;
  delete[] bvect;
}
