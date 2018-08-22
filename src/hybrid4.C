/*
  Program: hybrid4.C
  Author:  D. Trinkle
  Date:    May 20, 2004
  
  NOTE:    This algorithm is not very good.  It fails because the use
           of symmetrization is simply incorrect--it attempts to fold out
	   the LGF, but I don't think that needs to be done.  It also
	   has a problem that because it uses a kpt mesh of k equivalent to
	   gamma (mod [n]), the IFT isn't very good for large R.  Though
	   we seem to avoid that, because I think we only use the IFT for
	   the set of points in the supercell.  I don't think that makes
	   a lot of sense... suffice it to say, I think while there were
	   some good ideas here, the final result does not produce good
	   answers.

  Purpose: This takes in:
           1. unit cell definition (including Cij--very important)
	   2. a periodic supercell
	   3. a folded dynamical matrix
	   4. the Ylm expansion of (L)^-1(L2)(L)^-1, where L2 is folded
	   5. the discretization correction for the elastic GF

	   and we attempt to compute the discretization correction for
	   the lattice GF at the set of points inside the supercell.

	   This is a so-called "hybrid method" where we attempt to combine
	   the correct aspects of the direct force method and the LGF
	   method.  The main difference now is that we attempt to subtract
	   off the discontinuity piece correctly... maybe later, we'll get
	   around to trying to fit that piece instead.

  Param.:  <cell> <supercell> <folded Dij> <Gdc.lm> <disc.elas>
           cell:       cell file describing our lattice
	   supercell:  supercell geometry
	   folded Dij: folded down dynamical matrix (i.e., direct force calc)
	   Gdc.lm:     Ylm expansion of discontuity correction (folded L2)
	   disc.elas:  discretization correction of EGF

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

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   This algorithm borrows from about 5 other programs that I've
           written in the course of trying to do this bloody calculation.
	   Basically we calculate GL.d at each k=0(mod[n]) in the BZ:

	   GL.d(k) = [D(k)]^-1 - GE(k)*fenv(k/km)

	   and

	   GL.d(0) = the spherical average of (lambda1)^-1(lambda2)(lambda1)^-1
	   
	   Then, we take GE.d(R) and FT to get GE.d(k); finally, we set
	   dG.d(k) = GL.d(k)-GE.d(k), and then inverse FT using the
	   k=0(mod[n]).  This gives dG.d(R); this is unfolded using
	   symmetry, and added to the original GE.d(R) to give our final
	   result.

	   Whew.

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
#include "sphere-harm.H"
#include "shell.H"
#include "fourier.H"
#include "kpts.H"
#include "cont-int.H"  // envelope function
#include "lambda.H"    // lambda2 matrix evaluation
#include "integrate.H" // used for doing the 4Pi average
#include <gsl/gsl_linalg.h> // needed for inverting Aji in symmetrization

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************
int gen_symmrel(point_type* p, int super[9], int gunit[MAXop][9], int Nop,
		int** symmrel);

void solve_symm(point_type *p, double g_i_ginv[MAXop][Ni][Ni], 
		int** symmrel, int Nsymm);


// Calculate (kk)
void a_mult_a (double a[3], double Cijkl[9][9], double aa[9]) 
{
  int i, j, k, l;
  for (i=0; i<3; ++i) {
    for (j=0; j<i; ++j)
      aa[index(i,j)] = aa[index(j,i)];
    for (   ; j<3; ++j) {
      aa[index(i,j)] = 0.;
      for (k=0; k<3; ++k)
        for (l=0; l<3; ++l)
          aa[index(i,j)] += a[k]*Cijkl[index(k,i)][index(j,l)]*a[l];
    }
  }
}


// Do a 4Pi average of (lambda1)^-1 lambda2 (lambda1)^-1
void avg_4pi(double Cijkl[9][9], double lambda2[9][9][9], double ave[9]);


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
const int NUMARGS = 5;
const char* ARGLIST = "<cell> <supercell> <folded Dij> <Gdc.lm> <disc.elas>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'k'}; // determine how we handle k=0
// flag characters.

const char* ARGEXPL = 
"  cell:       cell file describing our lattice\n\
  supercell:  supercell geometry\n\
  folded Dij: folded down dynamical matrix (i.e., direct force calc)\n\
  Gdc.lm:     Ylm expansion of discontuity correction (folded lambda2)\n\
  disc.elas:  discretization correction of EGF\n\
  -k          set the kmax value (ignore value in disc.elas)";

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
  int FIX_KMAX = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *super_name, *latt_name, *Ylm_name, *elas_name;
  
  // Let's pull off the args:
  cell_name = args[0];
  super_name = args[1];
  latt_name = args[2];
  Ylm_name = args[3];
  elas_name = args[4];
  
  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; 
  double Cijkl[9][9];
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
  //-- ==== cell ====
  if (TESTING) {
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    printf("# atomic mass = %.3lf\n", atomic_mass);
  }
  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);


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
    fprintf(stderr, "Error encountered with file %s\n", super_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Supercell definition:\n");
    for (i=0; i<3; ++i) {
      printf("## A%d =", i+1);
      for (j=0; j<3; ++j) printf(" %+3d a%d", super_n[j*3+i], j+1);
      printf("\n");
    }
    printf("## %d atoms\n", Natoms*det(super_n));
  }

  //++ ==== Dij(R) ====
  int Np;
  point_type *Dij;
  
  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, Dij, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== Dij(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", latt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (n=0; n<Np; ++n) {
      point_type *tp = Dij + n;
      printf("## %3d %3d %3d ", 
             tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      print_mat(tp->mat);
      printf("\n");
    }
  }

  //++ ==== G0(lm) ====
  // NOTE: this is different than other versions, since we're *forced*
  // to only have EVEN l values.  So [l] is really 2*l.
  int Lmax, Ndim;  // Maximum l value / 2
  double ***RYlm0, ***IYlm0; // Separate our real and imag. components

  infile = myopenr(Ylm_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", Ylm_name);
    exit(ERROR_NOFILE);
  }

  ERROR = read_Ylm_even(infile, Lmax, Ndim, RYlm0, IYlm0);
  myclose(infile);

  if ( ERROR || (Ndim != 9) ) {
    fprintf(stderr, "Bad Lmax (%d) or Ndim (%d) value in %s\n", 
            Lmax, Ndim, Ylm_name);
    ERROR = ERROR_BADFILE;
  }
  //-- ==== G0(lm) ====

  if (ERROR) exit(ERROR); // time to bail out yet?

  if (TESTING) {
    verbose_outYlm(Lmax, Ndim, RYlm0, IYlm0, "##");
  }


  //++ ==== disc.elas(R) ====
  int Nd;
  point_type *GEd;
  double kmax;
  
  infile = myopenr(elas_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", elas_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Nd, GEd, READ_DIJ_MAT, dump);
  myclose(infile);
  sscanf(dump, "%*d %lf", &kmax); // %* to skip the first entry
  //-- ==== disc.elas(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", elas_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Elastic discretization correction:\n");
    for (n=0; n<Nd; ++n) {
      point_type *tp = GEd + n;
      printf("## %3d %3d %3d ", 
             tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      print_mat(tp->mat);
      printf("\n");
    }
  }

  // Let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[9], cartB[9];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

  if (FIX_KMAX) {
    kmax = max_k_BZ(cartb);
  }

  // ***************************** ANALYSIS **************************
  // There are a lot of steps here, so we'll need to hold on tightly,
  // just to keep from getting thrown to the wolves.

  // 0. Generate point group info
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  double g_i_ginv[MAXop][Ni][Ni];
  // 0.a. point group operations:
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);
  // 0.b. g [i] ginv, where [i] is a voigt matrix
  for (n=0; n<Nop; ++n) {
    for (i=0; i<Ni; ++i) {
      gen_g_i_ginv(i, gcart[n], gcart[inv_index[n]], g_i_ginv[n][i]);
    }
  }
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

  // 1. Make our supercell
  point_type *dGR;
  int Nsuper;
  fill_super(cart, super_n, Nsuper, dGR);
  if (TESTING) {
    printf("## Supercell geometry (%d atoms):\n", Nsuper);
    for (n=0; n<Nsuper; ++n) {
      printf("## %3d : %3d %3d %3d\n", n+1, 
	     dGR[n].Runit[0], dGR[n].Runit[1], dGR[n].Runit[2]);
    }
  }

  // 2. Time to make k = 0 (mod[n]) and the weights

  // 2.b. Make our set of gamma equivalent kpoints:
  int** klist;
  int Nkpt;
  int transn[9]; // transpose of super_n
  transpose(super_n, transn);
  int inv_transn[9];
  Nkpt = inverse(transn, inv_transn);
  mult(cartb, inv_transn, cartB);
  mult(cartB, 1./Nkpt, cartB);

  // Now use fill_super to make the points
  fill_super(cartB, transn, Nkpt, klist);
  // Finally, write in terms of rlv and cartesian coordinates:
  // but *keep* gamma:
  double klist_unit[Nkpt][3], klist_cart[Nkpt][4]; // [3] = magn
  for (n=0; n<Nkpt; ++n) {
    int temp[3];
    mult_vect(inv_transn, klist[n], temp);
    for (i=0; i<3; ++i) klist_unit[n][i] = temp[i]*(1./Nkpt);
    // Now, make cartesian:
    mult_vect(cartb, klist_unit[n], klist_cart[n]);
    // magnitude:
    klist_cart[n][3] = sqrt(dot(klist_cart[n], klist_cart[n]));
  }
  // Garbage collection...
  free_super(Nkpt, klist);

  // 2.c. Calculate the weights for each kpoint
  double w[Nkpt];
  for (k=0; k<Nkpt; ++k) w[k] = 1./Nkpt;

  if (TESTING) {
    printf("## %d kpoints equivalent to gamma:\n", Nkpt);
    printf("## equal weights = %.8le\n", w[0]);
    for (n=0; n<Nkpt; ++n) {
      printf("## %3d : %8.5lf %8.5lf %8.5lf | %8.5lf %8.5lf %8.5lf | %8.5lf\n", 
             n+1,
             klist_unit[n][0], klist_unit[n][1], klist_unit[n][2], 
             klist_cart[n][0], klist_cart[n][1], klist_cart[n][2],
             klist_cart[n][3]);
    }
  }

  // 3. At this stage, we'll construct lambda2 again, so as to get
  //    more accurate numbers in fourier space.
  double lambda2[9][9][9];
  make_lambda2(Np, Dij, lambda2);
  
  // 4. If we've made it this far, we can construct the discretization
  //    correction for our LGF at all k = 0(mod[n]), and subtract off
  //    the elastic disc. correction
  double dGk[Nkpt][9];
  int kgamma; // index indicating which point is truly gamma

  for (k=0; k<Nkpt; ++k) {
    double *tkpt = klist_cart[k];
    if (TESTING) {
      printf("## kpt %3d %8.5lf %8.5lf %8.5lf | %8.5lf\n", 
             k+1,tkpt[0], tkpt[1], tkpt[2], tkpt[3]);
    }
    if (! zero(tkpt[3])) {
      double Dk[9], Dkinv[9], kCk[9], GEk[9], L2[9], Gdck[9];
      // a. D(k):
      fourier(Np, Dij, tkpt, Dk, KPT_CART);
      careful_inverse(Dk, Dkinv);
      // b. GE(k):
      a_mult_a(tkpt, Cijkl, kCk);
      careful_inverse(kCk, GEk);
      // c. Gdc(k):
      lambda2_k(tkpt, lambda2, L2);
      mult(GEk, L2, GEk, Gdck);
      // scale both GE and Gdc by the envelope function:
      double gscale = kpoint_envelope(tkpt[3]/kmax)/(cell_vol);
      mult(GEk, gscale, GEk);
      mult(Gdck, gscale, Gdck);
      for (d=0; d<9; ++d) dGk[k][d] = Dkinv[d] - GEk[d] - Gdck[d];
      if (TESTING) {
	printf("## Dk:       "); print_mat(Dk); printf("\n");
	printf("## Dk^-1:    "); print_mat(Dkinv); printf("\n");
	printf("## GEk:      "); print_mat(GEk); printf("\n");
	printf("## Gdck:     "); print_mat(Gdck); printf("\n");
      }
    }
    else {
      for (d=0; d<9; ++d) dGk[k][d] = 0.;
    }
    if (TESTING) {
      printf("## GL-GE-Gdc:"); print_mat(dGk[k]); printf("\n");
    }
    
    // GE.d(k): 
    double GEdk[9];
    fourier(Nd, GEd, tkpt, GEdk, KPT_CART);
    // subtract off GE.d:
    for (d=0; d<9; ++d) dGk[k][d] += -GEdk[d];
    if (TESTING) {
      printf("## GEdk:     "); print_mat(GEdk); printf("\n");
      printf("## dGk:      "); print_mat(dGk[k]); printf("\n");
    }
      
  }

  // 5. Inverse FT using the kpoints equivalent to gamma:
  if (TESTING) {
    printf("## inverse FT: (pre-symmetrization)\n");
  }
  for (n=0; n<Nsuper; ++n) {
    point_type *tp = dGR + n;
    for (d=0; d<9; ++d) tp->mat[d] = 0;
    for (k=0; k<Nkpt; ++k) {
      double *tkpt = klist_cart[k];
      double coskR = cos(dot(tp->Rcart, tkpt));
      for (d=0; d<9; ++d) tp->mat[d] += w[k]*coskR*dGk[k][d];
    }
    if (TESTING) {
      printf("## point: %3d %3d %3d\n", tp->Runit[0], tp->Runit[1], 
	     tp->Runit[2]);
      printf("## cart:  %8.5lf %8.5lf %8.5lf\n", tp->Rcart[0], 
	     tp->Rcart[1], tp->Rcart[2]);
      printf("## dGR: "); print_mat(tp->mat); printf("\n");
    }
  }

  // 6. This is the MESSY bit... designed to enforce point group symmetry
  // on our set of points.  This can actually be used "post mortem"
  // on a set of points that were folded down already *without*
  // symmetrization turned on. -- taken from fold-super:

  // 6.a. make open shell list
  shell_type* sh = NULL;
  int Nsh;
  sort_points(dGR, Nsuper);
  ERROR = gen_shell_list_open(dGR, Nsuper, gunit, Nop, sh, Nsh);

    
  // 6.b. Now, we need to estimate how many new points are needed:
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
    
  // 6.c. Now, let's make 'em one by one...
  int **symmrel = new int *[Nop]; // symm. related points
  int Nsymm;
  for (n=0; n<Nop; ++n) symmrel[n] = new int[4];
  
  Nnew = 0;
  for (n=0; n<Nsh; ++n) 
    if (opensh[n]) {
      shell_type *tsh = sh + n;
      for (np=0; np<tsh->Nelem; ++np) {
	point_type *tp = dGR + (tsh->elem[np]);
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
	    mult(gcart[symmrel[i][3]], tp->mat,
		 gcart[inv_index[symmrel[i][3]]], pnew[Nnew].mat);
	    ++Nnew;
	  }
	}
      }
    }
  // NOW, we're slick... we merge the two lists:
  point_type* pcombine;
  pcombine = new point_type[Nsuper+Nnew];
  for (n=0; n<Nsuper; ++n) pcombine[n] = dGR[n];
  for (n=0; n<Nnew; ++n) pcombine[n+Nsuper] = pnew[n];
  // remove our old separate lists:
  delete[] dGR;
  delete[] pnew;
  dGR = pcombine;
  Nsuper += Nnew;
  sort_points(dGR, Nsuper);
  
  // garbage collection
  for (n=0; n<Nop; ++n) delete[] symmrel[n];
  delete[] symmrel;
  
  delete[] sh;

  // 7. We combine the GEd with dGR.
  sort_points(GEd, Nd);
  Nnew = Nsuper + Nd; // try to add 'em up
  pnew = new point_type[Nnew];
  Nnew = 0;
  i=0; j=0;
  int combinetype;
  while ( (i<Nsuper) || (j<Nd) ) {
    combinetype = 0;
    if (i==Nsuper) combinetype = 1;
    if (j==Nd) combinetype = -1;
    if (! combinetype) 
      combinetype = __COMPARE_POINT_TYPE__((void*)(dGR+i), (void*)(GEd+j));
    switch (combinetype) {
    case -1: // copy in dGR:
      copy_point(dGR[i], pnew[Nnew]);
      ++i;
      break;
    case 1:  // copy in GEd:
      copy_point(GEd[j], pnew[Nnew]);
      ++j;
      break;
    case 0: 
      copy_point(dGR[i], pnew[Nnew]);
      for (d=0; d<9; ++d) pnew[Nnew].mat[d] += GEd[j].mat[d];
      ++i;
      ++j;
    }
    ++Nnew;
  }
  delete[] GEd;
  GEd = pnew;
  Nd = Nnew;


  // 8. We now have to put the inverse FT of the discontinuity correction
  //    into our points.
  //    Note: we have to do this *after* symmetrization, because we don't
  //    have an easy way to handle supercell points which might be "doubled"

  // allocate some workspace, and our transformed Ylm expansions:
  SH_workspace_type * work = allocate_SH_workspace(Lmax);
  double ***RYlm_R, ***IYlm_R; // for handling the discontinuity
  double Fl_val[Lmax+1];
  double gint[9];
  ERROR = alloc_Ylm(Lmax, Ndim, RYlm_R, IYlm_R);
  if (ERROR) {
    fprintf(stderr, "Error allocating our Ylm expansion...\n");
    exit(ERROR);
  }
  double lastRmagn = -1;
  for (n=0; n<Nd; ++n) {
    point_type *tGd = GEd + n;
    // inverse FT the discontinuity correction
    // we've pulled this algorithm from discDij;
    double gint[9];
    if (! zero(tGd->Rmagn)) {
      if (! dcomp(tGd->Rmagn, lastRmagn)) {
	// we need to redo this calc...
	Fn_x2_int(Lmax, kmax*tGd->Rmagn, Fl_val, TESTING, 1024);
	// Fn_x2_int(Lmax_0, kmax*tGd->Rmagn, Fl_val, 0, 1024);
	if (TESTING) {
	  printf("## x^2 int-- Fl_val:\n");
	  for (l=0; l<=Lmax; ++l) {
	    printf("## %d %.12le\n", l*2, Fl_val[l]);
	  }
	}
	double g_scale = 1./(4*M_PI*(tGd->Rmagn*tGd->Rmagn*tGd->Rmagn)), l_scale;
	int l_sign;
	for (l=0, l_sign = 1; l<=Lmax; ++l, l_sign = -l_sign) {
	  l_scale = g_scale * l_sign * Fl_val[l];
	  for (m=0; m<=(4*l); ++m) {
	    for (d=0; d<Ndim; ++d) {
	      // Initialize
	      RYlm_R[l][m][d] = RYlm0[l][m][d] * l_scale;
	      IYlm_R[l][m][d] = IYlm0[l][m][d] * l_scale;
	    }
	  }
	}
	lastRmagn = tGd->Rmagn;
      }
      // else, we calculated it last time, so keep using it.
      // Get the continuum piece:
      eval_Ylm_expansion_R (Lmax, tGd->Rcart, Ndim, RYlm_R, IYlm_R, gint,
                            work);
    }
    else {
      // For R = 0, the continuum piece is really simple:
      double g0_scale = M_2_SQRTPI/(8*M_PI*M_PI); // strange (exact) number
      g0_scale *= kmax*kmax*kmax*kpoint_x2_envelope_int();  // and our maximum k value...
      for (d=0; d<Ndim; ++d)
        gint[d] = g0_scale * RYlm0[0][0][d];
    }
    if (TESTING) {
      printf("## point: %3d %3d %3d\n", tGd->Runit[0], tGd->Runit[1], 
	     tGd->Runit[2]);
      printf("## gint: "); print_mat(gint); printf("\n");
      printf("## tGd:  "); print_mat(tGd->mat); printf("\n");
    }
    // add back to our real space piece:
    for (d=0; d<9; ++d) tGd->mat[d] += gint[d];
    if (TESTING) {
      printf("## sum:  "); print_mat(tGd->mat); printf("\n");
    }
  }
  // garbage collection:
  free_Ylm(Lmax, RYlm_R, IYlm_R);
  free_SH_workspace(work);
  


  // ****************************** OUTPUT ***************************
  //  for (d=0; d<9; ++d)
  //    printf("%3d ", super_n[d]);
  //  printf("# supercell definition\n");
  printf("%d %.12le # Number of points\n", Nd, kmax);
  for (np=0; np<Nd; ++np) {
    point_type *tp = GEd + np;
    printf("%3d %3d %3d", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    for (d=0; d<9; ++d)
      printf(" %.12le", tp->mat[d]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************

  free_Ylm(Lmax, RYlm0, IYlm0);

  delete[] GEd;
  delete[] dGR;
  delete[] Dij;
  
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


// Do a 4Pi average of (lambda1)^-1 lambda2 (lambda1)^-1
void avg_4pi(double Cijkl[9][9], double lambda2[9][9][9], double ave[9]) 
{
  int d;
  const int N=16;
  double x[N], w[N];
  double lambda1k[9], lambda1_inv[9], lambda2k[9], res[9];
  double kpt[3], costheta, sintheta;
  
  gauleg(-1.0, 1.0, x, w, N);
  
  for (d=0; d<9; ++d) ave[d] = 0.;
  for (int i=0; i<N; ++i) {
    costheta = x[i];
    kpt[2] = costheta; 
    sintheta = sqrt(1-costheta*costheta);
    double tave[9];
    for (d=0; d<9; ++d) tave[d] = 0.;
    for (int j=0; j<N; ++j) {
      kpt[0] = sintheta*cos(2*M_PI*j);
      kpt[1] = sintheta*sin(2*M_PI*j);
      a_mult_a(kpt, Cijkl, lambda1k);
      careful_inverse(lambda1k, lambda1_inv);
      lambda2_k(kpt, lambda2, lambda2k);
      mult(lambda1_inv, lambda2k, lambda1_inv, res);
      for (d=0; d<9; ++d) tave[d] += res[d];
    }
    mult(tave, 1./N, tave);
    for (d=0; d<9; ++d) ave[d] += tave[d] * w[i];
  }
}
      
