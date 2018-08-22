/*
  Program: fitsuper.C
  Author:  D. Trinkle
  Date:    June 9, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. the dynamical matrix defined at a set of points
	   3. a periodic supercell

	   and we compute
	   1. periodic version of lattice function evaluated at points
	      inside of the supercell only

	   This is fit to match the elastic constants in the unit
	   cell, and is symmetrized.

  Param.:  <cell> <Dij(R))> <supercell>
           cell:    cell file describing our lattice
	   Dij(R):  matrix function evaluated at a set of points
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

	   ==== Dij(R) ====
	   N kmax                      # number of points, kmax
	   n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij
	   ...
	   nN.1 nN.2 nN.3  Lxx .. Lzz
	   ==== Dij(R) ====

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
#include "fourier.H"
#include "kpts.H" // fold_down
#include "constraint.H"
#include <gsl/gsl_linalg.h> // needed for inverting Aji in symmetrization

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************
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
const char* ARGLIST = "<cell> <Dij(R)> <supercell>";


const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'f'}; // flag characters.

const char* ARGEXPL = 
"  cell:      cell file describing our lattice\n\
  Dij(R):    matrix function evaluated at a set of points\n\
  supercell: supercell geometry\n\
  -f         free fit (i.e., no constraints)";

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
==== Dij(R) ====\n\
N kmax                      # number of points, kmax\n\
n1.1 n1.2 n1.3  Lxx .. Lzz  # unit cell coord, Lij\n\
...\n\
nN.1 nN.2 nN.3  Lxx .. lzz\n\
==== Dij(R) ====\n";

int main ( int argc, char **argv ) 
{
  int d, i, j, k, m, n, np; // General counting variables.

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
  int FREE = flagon[0];

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
  double* Cmn_list; // elastic constant input
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
  if (TESTING) {
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    printf("# atomic mass = %.3lf\n", atomic_mass);
  }
  //-- ==== cell ====

  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);

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
  // 1.a. Make our supercell points
  point_type *period;
  int Nsuper;
  fill_super(cart, super_n, Nsuper, period);

  // 1.b. Symmetrize 'em and fill out the shell list.
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  double g_i_ginv[MAXop][Ni][Ni];
  shell_type* sh = NULL;
  int Nsh;
  // 1.b: point group operations:
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);
  // 1.b: g [i] ginv, where [i] is a voigt matrix
  for (n=0; n<Nop; ++n) {
    for (i=0; i<Ni; ++i) {
      gen_g_i_ginv(i, gcart[n], gcart[inv_index[n]], g_i_ginv[n][i]);
    }
  }
  // 1.c. make open shell list
  sort_points(period, Nsuper);
  ERROR = gen_shell_list_open(period, Nsuper, gunit, Nop, sh, Nsh);
  
  if (TESTING) {
    /*
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
    */
    printf("## Open shell: Nsh= %d  Np= %d\n", Nsh, Nsuper);
    for (n=0; n<Nsh; ++n) {
      printf("## shell %d :Nelem= %d\n", n, sh[n].Nelem);
      printf("##  ");
      for (i=0; i<sh[n].Nelem; ++i)
	printf(" (%3d%3d%3d)", 
	       period[sh[n].elem[i]].Runit[0], 
	       period[sh[n].elem[i]].Runit[1], 
	       period[sh[n].elem[i]].Runit[2]);
      printf("\n");
    }
  }
  // 1.d. close the shell list
  close_shell_list(period, Nsuper, cart, gunit, Nop, sh, Nsh);
  
  if (TESTING) {
    printf("## Closed shell: Nsh= %d  Np= %d\n", Nsh, Nsuper);
    for (n=0; n<Nsh; ++n) {
      printf("## shell %d :Nelem= %d\n", n, sh[n].Nelem);
      printf("##  ");
      for (i=0; i<sh[n].Nelem; ++i)
	printf(" (%3d%3d%3d)", 
	       period[sh[n].elem[i]].Runit[0], 
	       period[sh[n].elem[i]].Runit[1], 
	       period[sh[n].elem[i]].Runit[2]);
      printf("\n");
    }
    printf("## Supercell geometry (%d atoms):\n", Nsuper);
    for (n=0; n<Nsuper; ++n) {
      printf("## %3d : %3d %3d %3d  |  %9.5lf %9.5lf %9.5lf\n", n+1, 
	     period[n].Runit[0], period[n].Runit[1], period[n].Runit[2],
	     period[n].Rcart[0], period[n].Rcart[1], period[n].Rcart[2]);
    }
  }
  

  // 2. Time to make k = 0 (mod[n]) and the weights
  // 2.a. First, let's make our reciprocal lattice vectors:
  // [b] = 2Pi*([a]^-1)^T
  double cartb[9], cartB[9];
  double cell_vol;
  cell_vol = inverse(cart, cartb);     // invert
  self_transpose(cartb);               // transpose in place
  mult(cartb, 2*M_PI/cell_vol, cartb); // scale

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
  for (k=0; k<Nkpt; ++k) {
    int temp[3];
    mult_vect(inv_transn, klist[k], temp);
    for (i=0; i<3; ++i) klist_unit[k][i] = temp[i]*(1./Nkpt);
    fold_down(klist_unit[k], cartb); // folddown into the BZ
    // Now, make cartesian:
    mult_vect(cartb, klist_unit[k], klist_cart[k]);
    // magnitude:
    klist_cart[k][3] = sqrt(dot(klist_cart[k], klist_cart[k]));
  }
  // Garbage collection...
  free_super(Nkpt, klist);

  // 2.c. Calculate the weights for each kpoint
  // --we weight each point by 1, unless the inverse isn't there.
  double w[Nkpt];
  for (k=0; k<Nkpt; ++k) w[k]=0;
  for (k=0; k<Nkpt; ++k)
    if (w[k] == 0) {
      double kinv[3] = {-klist_unit[k][0],-klist_unit[k][1],-klist_unit[k][2]};
      fold_down(kinv, cartb);
      for (j=(k+1); j<Nkpt; ++j)
	if (w[j] == 0) 
	  if (equal_vect(klist_unit[j], kinv))
	    break;
      if (j != Nkpt) {
	w[k] = 1;
	w[j] = 1;
      }
      else w[k] = 2;
    }
  if (TESTING) {
    double scale=0;
    for (k=0; k<Nkpt; ++k) scale += w[k];
    scale *= 0.5;
    printf("## %d kpts: k = 0(mod [n])  1/2 sum w[k] = %.2lf\n", Nkpt, scale);
    for (k=0; k<Nkpt; ++k) {
      printf("## %8.5lf %8.5lf %8.5lf | %8.5lf %8.5lf %8.5lf  w= %1.0lf\n",
	     klist_unit[k][0], klist_unit[k][1], klist_unit[k][2],
	     klist_cart[k][0], klist_cart[k][1], klist_cart[k][2],
	     w[k]);
    }
  }


  // 3.a. Construct the constraint and alpha matrices for fitting
  // This isn't too terribly difficult, but it's harder than before,
  // because we do it over shells now.
  int Nconst = 7*Ni; // value at gamma, 6 second derivatives
  int Nfree = Nsh*Ni;
  int Neval = Nkpt*Ni;
  if (FREE) Nconst=Ni;
  double *Cmat = new double[Nconst*Nfree];
  double *alpha = new double[Neval*Nfree];
  double *weight = new double[Neval];

  for (k=0; k<Nkpt; ++k)
    for (j=0; j<Ni; ++j)
      weight[k*Ni+j] = w[k];
  
  // second derivatives: use voigt notation.
  if (!FREE) {
    for (int ab=0; ab<Ni; ++ab) {
      int a=i2rank2[ab][0], b=i2rank2[ab][1];
      for (n=0; n<Nsh; ++n) {
	shell_type *tsh = sh+n;
	for (int cd=0; cd<Ni; ++cd)
	  for (i=0; i<Ni; ++i)
	    Cmat[(ab*Ni+cd)*Nfree+(n*Ni+i)] = 0;
	for (int g=0; g<Nop; ++g)
	  for (int cd=0; cd<Ni; ++cd) 
	    for (i=0; i<Ni; ++i)
	      Cmat[(ab*Ni+cd)*Nfree+(n*Ni+i)] += 
		-period[tsh->gR0[g]].Rcart[a]*period[tsh->gR0[g]].Rcart[b]
		*g_i_ginv[g][i][cd];
	// scale final result:
	double scale = (double)tsh->Nelem / (double)Nop;
	for (int cd=0; cd<Ni; ++cd)
	  for (i=0; i<Ni; ++i)
	    Cmat[(ab*Ni+cd)*Nfree+(n*Ni+i)] *= scale;
      }
    }
  }
  // D(0):
  k = Nconst-Ni;
  for (n=0; n<Nsh; ++n) {
    shell_type *tsh = sh+n;
    for (j=0; j<Ni; ++j)
      for (i=0; i<Ni; ++i)
	Cmat[(k +j)*Nfree+(n*Ni+i)] = 0;
    for (int g=0; g<Nop; ++g)
      for (j=0; j<Ni; ++j)
	for (i=0; i<Ni; ++i)
	  Cmat[(k +j)*Nfree+(n*Ni+i)] += g_i_ginv[g][i][j];
    // scale final result:
    double scale = (double)tsh->Nelem / (double)Nop;
    for (j=0; j<Ni; ++j)
      for (i=0; i<Ni; ++i)
	Cmat[(k +j)*Nfree+(n*Ni+i)] *= scale;
  }
  
  // alpha: 
  for (k=0; k<Nkpt; ++k) {
    double tk[3] = {klist_cart[k][0], klist_cart[k][1], klist_cart[k][2]};
    for (n=0; n<Nsh; ++n) {
      shell_type *tsh = sh+n;
      // Begin by clearing out the Ni*Ni block of entries:
      for (j=0; j<Ni; ++j)
	for (i=0; i<Ni; ++i)
	  alpha[(k*Ni+j)*Nfree + (n*Ni+i)] = 0;
      for (int g=0; g<Nop; ++g) {
	double cos_kgR = cos(dot(tk, period[tsh->gR0[g]].Rcart));
	for (j=0; j<Ni; ++j)
	  for (i=0; i<Ni; ++i)
	    alpha[(k*Ni+j)*Nfree + (n*Ni+i)] += cos_kgR*g_i_ginv[g][i][j];
      }
      // scale final result:
      double scale = (double)tsh->Nelem / (double)Nop;
      for (j=0; j<Ni; ++j)
	for (i=0; i<Ni; ++i)
	  alpha[(k*Ni+j)*Nfree + (n*Ni+i)] *= scale;
    }
  }

  // 3.b. Solve
  double* FC = new double[Nfree*Nconst];
  double* FE = new double[Nfree*Neval];
  int fit_return;
  
  fit_return = 
    constrain_fit(Nfree, Nconst, Neval, Cmat, alpha, weight, FC, FE, TESTING);

  if (TESTING) {
    printf("## \"Meaningless\" warnings from constrain_fit:\n");
    if (fit_return & CONST_EXTRAC) {
      printf("## Extraneous constraints...\n");
    }
    if (fit_return & CONST_UNDER) {
      printf("## Underdetermined fit problem...\n");
    }
  }

  // 3.c. determine the cvect and Dk values
  double cvect[Nconst];
  double Dk_d[Neval];
  double DR[Nfree]; // vector to hold the values  

  double Dk[Nkpt][9];
  // fourier transform:
  for (k=0; k<Nkpt; ++k)
    fourier(Np, latt, klist_cart[k], Dk[k], KPT_CART);

  if (!FREE) {
    if (TESTING) {
      printf("## constraints:\n");
    }
    for (int ab=0; ab<Ni; ++ab) {
      int a=i2rank2[ab][0], b=i2rank2[ab][1];
      for (int cd=0; cd<Ni; ++cd) {
	int c=i2rank2[cd][0], d=i2rank2[cd][1];
	cvect[ab*Ni+cd] = 
	  cell_vol*(Cijkl[a*3+c][b*3+d] + Cijkl[a*3+d][b*3+c]);
	if (TESTING) {
	  printf("## d^2D_%c%c/dk_%c dk_%c = %.5le\n",
		 coord[c], coord[d], coord[a], coord[b], cvect[ab*Ni+cd]);
	}
      }
    }
  }
  
  for (i=0; i<Ni; ++i)
    cvect[Nconst-Ni + i] = 0.; // D(0)=0


  for (k=0; k<Nkpt; ++k)
    for (j=0; j<Ni; ++j)
      Dk_d[k*Ni+j] = Dk[k][i6tom9[j]];

  // solve
  for (i=0; i<Nfree; ++i) {
    DR[i]=0;
    for (j=0; j<Nconst; ++j) DR[i] += FC[i*Nconst+j]*cvect[j];
    for (k=0; k<Neval; ++k) DR[i] += FE[i*Neval+k]*Dk_d[k];
  }
  if (TESTING) {
    printf("## Raw variable values:\n");
    for (n=0; n<Nsh; ++n) {
      shell_type *tsh = sh+n;
      point_type *tp = period + tsh->elem[0];
      printf("## (%3d%3d%3d):", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      for (i=0; i<Ni; ++i)
	printf(" %s= %.5le", v_coord[i], DR[n*Ni+i]);
      printf("\n");
    }
  }

  for (n=0; n<Nsh; ++n) {
    shell_type *tsh = sh+n;
    // fill in the representative point first
    point_type *tp = period+tsh->elem[0];
    for (d=0; d<9; ++d) tp->mat[d] = DR[n*Ni+m9toi6[d]];
    // now, work through our other shell members:
    for (i=1; i<(tsh->Nelem);  ++i) {
      m = tsh->elem[i];
      point_type *tnp = period+m;
      // find the group elem that takes n to m:
      for (j=0; (j<Nop) && (m != tsh->gR0[j]); ++j) ;
      // now, multiply tp over to tnp:
      mult(gcart[j], tp->mat, gcart[inv_index[j]], tnp->mat);
    }
  }
  // error estimate:
  double error = 0;
  for (k=0; k<Nkpt; ++k) 
    for (j=0; j<Ni; ++j) {
      double Dinter= -Dk_d[k*Ni+j];
      for (i=0; i<Nfree; ++i) Dinter += alpha[(k*Ni+j)*Nfree+i]*DR[i];
      error += weight[k*Ni+j] * Dinter*Dinter;
    }
  
  // scale error:
  double scale=0;
  for (k=0; k<Nkpt; ++k)
    for (j=0; j<Ni; ++j)
      scale += weight[k*Ni+j];
  error = sqrt(error/(scale));

  // GC:
  delete[] weight;
  delete[] Cmat;
  delete[] alpha;
  delete[] FC;
  delete[] FE;

  if (TESTING) {
    printf("## Fit dynamical matrix:\n");
    for (n=0; n<Nsuper; ++n) {
      point_type *tp = period + n;
      printf("## %3d %3d %3d ", 
             tp->Runit[0], tp->Runit[1], tp->Runit[2]);
      print_mat(tp->mat);
      printf("\n");
    }
  }


  // ****************************** OUTPUT ***************************
  printf("%d 9 # super= ", Nsuper);
  for (d=0; d<9; ++d) printf("%2d ", super_n[d]);
  printf("  estimated error = %.3le\n", error);
  for (np=0; np<Nsuper; ++np) {
    point_type *tp = period + np;
    printf("%3d %3d %3d", tp->Runit[0], tp->Runit[1], tp->Runit[2]);
    for (d=0; d<9; ++d)
      printf(" %.12le", tp->mat[d]);
    printf("\n");
  }

  // ************************* GARBAGE COLLECTION ********************
  delete[] sh;
  delete[] period;
  delete[] latt;
  
  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}
