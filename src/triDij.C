/*
  Program: triDij.C
  Author:  D. Trinkle
  Date:    June 5, 2004
  Purpose: This takes in:
           1. unit cell definition
	   2. lattice function defined at a set of points
	   3. a periodic supercell

	   and we compute the tricubic expansion of the dynamical matrix
	   in Fourier space, using a grid defined by the periodic supercell
	   in real space.  Whew.

	   Basically, we're trying to create an interpolation scheme
	   for direct force, instead of using the usual fourier interpolation.
	   
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

  Algo.:   We construct the reciprocal lattice points, and the grid
           spacing for the unit cube.  We FT our dynamical matrix at those
	   points, and then throw it at our tricubic interpolator.

  Output:  Output the supercell in *reciprocal* space (which is [n]^T)
           followed by the full expansion; suitable for reading by
	   the corresponding routine in tricubic.H
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
#include "tricubic.H"

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
const char* ARGLIST = "<cell> <latt(R)> <supercell>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  cell:      cell file describing our lattice\n\
  Dij(R):    dynamical matrix evaluated at a set of points\n\
  supercell: supercell geometry in real space";

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
  int d, i, j, k, l, m, np; // General counting variables.

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
  point_type *Dij;
  
  infile = myopenr(latt_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", latt_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, Dij, READ_DIJ_MAT);
  myclose(infile);
  //-- ==== latt(R) ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Error encountered with file %s\n", latt_name);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Dynamical matrix:\n");
    for (np=0; np<Np; ++np) {
      point_type *tp = Dij + np;
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

  ERROR = (det(super_n) == 0);

  // If we've thrown an error already, then get out.
  if (ERROR) {
    fprintf(stderr, "Problem with supercell... det=0?\n");
    exit(ERROR);
  }

  // ***************************** ANALYSIS **************************
  // NOTE: We're very "clever" and do everything with our kpoints
  // in UNIT coordinates, so we don't even have to bother computing
  // our RLV.

  int nRL[9], nRLinv[9]; // our [n] in kspace, and its inverse:
  int detn;
  transpose(super_n, nRL);
  detn = inverse(nRL, nRLinv);

  // 1. Determine grid spacing--through the magic of gcd.
  int Ni[3];
  for (d=0; d<3; ++d) {
    Ni[d] = abs( detn / gcd(nRLinv[3*d],nRLinv[3*d+1],nRLinv[3*d+2]) );
  }
  
  if (TESTING) {
    printf("## supercell [n] (real space):     ");
    for (d=0; d<9; ++d) printf("%3d", super_n[d]);
    printf("\n");
    printf("## supercell [n'] (fourier space): ");
    for (d=0; d<9; ++d) printf("%3d", nRL[d]);
    printf("\n");
    printf("## [n']^-1: detn = %3d points      ", detn);
    for (d=0; d<9; ++d) printf("%3d", nRLinv[d]);
    printf("\n");
    printf("## grid spacing: %3d %3d %3d  -> %3d points\n", 
	   Ni[0], Ni[1], Ni[2], Ni[0]*Ni[1]*Ni[2]);
    printf("## ratio of additional points: %3d\n", 
	   abs(Ni[0]*Ni[1]*Ni[2]/detn));
  }
  
  // 2. Allocate a lookup table for our FT:
  // I *think* our grid spacing should actually be lcm(N0,N1,N2),
  // but to be on the safe side, we know that detn is the upper limit...
  coskR_table_type coskR;
  init_coskR_table(abs(detn), coskR);

  // 3. Evaluate our function on our grid spacing:
  int N=Ni[0]*Ni[1]*Ni[2];
  int*** index, n[3], ind;
  int Ndim=6; // we use Voigt notation to cut down our work slightly.
  double** func = new double*[Ndim];
  for (d=0; d<Ndim; ++d) func[d] = new double[N];

  double kpt[3]; // in unit coordinates, in fourier space
  double Dk[9];  // the FT of D(R)
  double n_inv[9];
  for (d=0; d<9; ++d) n_inv[d] = (double)nRLinv[d] / (double)detn;

  index = new int**[Ni[0]];
  for (ind=0, n[0]=0; n[0]<Ni[0]; ++(n[0])) {
    index[n[0]] = new int*[Ni[1]];
    for (n[1]=0; n[1]<Ni[1]; ++(n[1])) {
      index[n[0]][n[1]] = new int[Ni[1]];
      for (n[2]=0; n[2]<Ni[2]; ++(n[2]), ++ind) {
	mult_vect(n_inv, n, kpt); // our kpt.
	fourier_table(Np, Dij, kpt, Dk, KPT_LATT, coskR); // D(k)
        for (d=0; d<Ndim; ++d) {
          func[d][ind] = Dk[i6tom9[d]]; // transform to Voigt notation
        }
	// store this index:
        index[n[0]][n[1]][n[2]] = ind;
      }
    }
  }

  if (TESTING) {
    printf("## D(k) on our grid:\n");
    for (n[0]=0; n[0]<Ni[0]; ++(n[0])) {
      for (n[1]=0; n[1]<Ni[1]; ++(n[1])) {
	for (n[2]=0; n[2]<Ni[2]; ++(n[2]), ++ind) {
	  mult_vect(n_inv, n, kpt); // our kpt.
	  ind = index[n[0]][n[1]][n[2]];
	  printf("## %8.5lf %8.5lf %8.5lf", kpt[0], kpt[1], kpt[2]);
	  for (d=0; d<Ndim; ++d)
	    printf(" D%s= %12.5le", v_coord[d], func[d][ind]);
	  printf("\n");
	}
      }
    }
  }

  // 4. calculate the tricubic expansion:
  tricubic_type* t;
  calc_tricubic(Ni[0], Ni[1], Ni[2], Ndim, index, func, t);


  // ****************************** OUTPUT ***************************
  for (d=0; d<9; ++d) printf("%3d ", nRL[d]);
  printf("# supercell in fourier space\n");
  write_tricubic(stdout, t, Ndim);


  // ************************* GARBAGE COLLECTION ********************
  for (d=0; d<Ndim; ++d)
    free_tricubic(t[d]);
  delete[] t;

  for (n[0]=0; n[0]<Ni[0]; ++(n[0])) {
    for (n[1]=0; n[1]<Ni[1]; ++(n[1])) 
      delete[] index[n[0]][n[1]];
    delete[] index[n[0]];
  }
  delete[] index;
  for (d=0; d<Ndim; ++d) delete[] func[d];
  delete[] func;

  free_coskR_table(coskR);
  delete[] Dij;
  
  return 0;
}
