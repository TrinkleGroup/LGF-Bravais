/*
  Program: symm-fillin.C
  Author:  D. Trinkle
  Date:    February 24, 2004
  Purpose: Read in a lattice GF and fill in the missing pieces using
           symmetry.  Recently changed to use *open* shells.  This is
	   a relatively minor modification, but important for calculating
	   the GF from periodic cells.
           
  Param.:  <cell> <GL(R)> <unknowns>
           cell:   cell file describing our lattice
	   GL(R):  lattice green function
	   unknowns: string containing some combination of xyz

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

           ==== GL(R) ====
           N                           # number of points
           n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GG
           ...
           nN.1 nN.2 nN.3  Gxx .. Gzz
           ==== GL(R) ====


  Flags:   MEMORY:  not used
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Wow... I need to write this all out sometime.  Basically, 
           we put together a bunch of linear equations for each symmetry
	   shell that we solve one by one for the missing LGF entries.

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
#include <gsl/gsl_linalg.h>  // for solving linear systems...
#include "pointgroup.H"
#include "Dij.H"
#include "shell.H"
#include "voigt.H"

//***************************** SUBROUTINES ****************************

int parse_unknown(char* arg, int &Nunk, int unk[3]);

void construct_matrix(point_type* GL, int gunit[MAXop][9], int Nop, 
		      double ginv_i_g[MAXop][Ni][Ni], shell_type *shell_list,
		      int Nunk, int unk[3],
		      int dim, double* Amat, double* bvect);


inline void print_smat (double a[9]) 
{
  printf("xx= %7.4lf yy= %7.4lf zz= %7.4lf  xy= %7.4lf yz= %7.4lf zx= %7.4lf",
         a[0], a[4], a[8],
         0.5*(a[1]+a[3]), 0.5*(a[5]+a[7]), 0.5*(a[2]+a[6]));
}

inline void print_mat (double a[9]) 
{
  printf("%7.4lf %7.4lf %7.4lf  %7.4lf %7.4lf %7.4lf  %7.4lf %7.4lf %7.4lf",
         a[0], a[1], a[2],  a[3], a[4], a[5],  a[6], a[7], a[8]);
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "<cell> <GL(R)> <unknowns>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:   cell file describing our lattice\n\
  GL(R):  lattice greens function\n\
  unknowns: string containing some combination of xyz";

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
==== GL(R) ====\n\
N                           # number of points\n\
n1.1 n1.2 n1.3  Gxx .. Gzz  # unit cell coord, GF\n\
...\n\
nN.1 nN.2 nN.3  Gxx .. Gzz\n\
==== GL(R) ====\n";

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

  // flags:

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters:
  char *cell_name, *GL_name;
  int Nunk, unk[3];
  
  // Let's pull off the args:
  cell_name = args[0];
  GL_name = args[1];
  ERROR = parse_unknown(args[2], Nunk, unk);
  if (ERROR) {
    printf("Bad unknown string: %s\n", args[2]);
    exit(ERROR);
  }
  if (TESTING) {
    printf("## Unknown elements: %s\n", args[2]);
    printf("## Nunk: %d  i6 unknown(s): ", Nunk);
    for (i=0; i<Nunk; ++i) printf(" %d", unk[i]);
    printf("\n");
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

  //++ ==== GL(R) ====
  int Np, Ndim = 9;
  point_type *GL;
  
  infile = myopenr(GL_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", GL_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_Dij(infile, cart, Np, GL, READ_DIJ_MAT);
  if (ERROR) exit(ERROR);
  //-- ==== GL(R) ====


  // ***************************** ANALYSIS **************************

  //++ ==== gen point group ====
  // First, generate our point group operations:
  double gcart[MAXop][9];
  int gunit[MAXop][9];
  int inv_index[MAXop];
  int Nop;
  
  Nop = gen_pointgroup(cart, gcart, gunit, inv_index);

  if (TESTING) {
    printf("## Nops = %d\n", Nop);
    for (n=0; n<Nop; ++n) {
      printf("## %2d: %2d %2d %2d  %2d %2d %2d  %2d %2d %2d  inv= %d\n", n+1, 
	     gunit[n][0], gunit[n][1], gunit[n][2], 
	     gunit[n][3], gunit[n][4], gunit[n][5], 
	     gunit[n][6], gunit[n][7], gunit[n][8], inv_index[n]+1);
      printf("##     %6.3lf %6.3lf %6.3lf  %6.3lf %6.3lf %6.3lf  %6.3lf %6.3lf %6.3lf\n",
	     gcart[n][0], gcart[n][1], gcart[n][2], 
	     gcart[n][3], gcart[n][4], gcart[n][5], 
	     gcart[n][6], gcart[n][7], gcart[n][8]);
    }
  }
  //-- ==== gen point group ====


  //++ ==== gen ginv [i] g elements ====
  double ginv_i_g[MAXop][Ni][Ni];
  
  for (n=0; n<Nop; ++n) 
    for (i=0; i<Ni; ++i)
      gen_ginv_i_g (i, gcart[n], gcart[inv_index[n]], ginv_i_g[n][i]);

  if (TESTING) {
    printf("## ginv [i] g testing:\n");
    for (n=0; n<Nop; ++n) {
      printf("## gcart: ");
      print_mat(gcart[n]);
      printf("\n");
      for (i=0; i<Ni; ++i) {
	printf("##  ginv [%d] g = (", i+1);
	for (j=0; j<Ni; ++j) {
	  printf("%7.4lf ", ginv_i_g[n][i][j]);
	}
	printf(")\n");
      }
    }
  }
  //-- ==== gen ginv [i] g elements ====
  
  
  //++ ==== generate the symmetric shells ====
  shell_type* shell_list;
  int Nshell;
  
  sort_points(GL, Np);
  ERROR = gen_shell_list_open(GL, Np, gunit, Nop, shell_list, Nshell);

  if (TESTING) {
    printf("## Nshells = %d\n", Nshell);
    for (n=0; n<Nshell; ++n) {
      shell_type *tsh = shell_list+n;
      printf("## shell %d -- %d entries\n", n, tsh->Nelem);
      for (i=0; i<(tsh->Nelem); ++i) {
	np = tsh->elem[i];
	printf("##      %11.7lf %11.7lf %11.7lf\n", 
	       GL[np].Rcart[0], GL[np].Rcart[1], GL[np].Rcart[2]);
      }
    }
  }
  
  if (ERROR) {exit(ERROR);}

  // Construct the index list for each point (used for making the matrices)
  int np_index[Np];
  for (n=0; n<Nshell; ++n) {
    shell_type *tsh = shell_list+n;
    for (i=0; i<(tsh->Nelem); ++i)
      np_index[tsh->elem[i]] = i;
  }
  //-- ==== generate the symmetric shells ====


  //++ ==== finally, loop through each shell, and solve the linear problem ====
  if (TESTING) {
    printf("## Matrices!!\n");
  }
  for (n=0; n<Nshell; ++n) {
    shell_type *tsh = shell_list+n;
    int dim = tsh->Nelem*Nunk;
    double Amat[dim*dim], bvect[dim];
    construct_matrix(GL, gunit, Nop, ginv_i_g, tsh, Nunk, unk, dim,
		     Amat, bvect);
    if (TESTING) {
      printf("## shell %d  entries: %d  dimen: %d\n", n, dim/Nunk, dim);
      for (i=0; i<dim; ++i) {
	np = tsh->elem[i/Nunk];
	printf("##  [");
	for (j=0; j<dim; ++j)
	  printf("%6.2lf", Amat[i*dim+j]);
	printf("][G(%6.2lf,%6.2lf,%6.2lf).%d]", 
	       GL[np].Rcart[0], GL[np].Rcart[1], GL[np].Rcart[2],
	       i6tom9[unk[i%Nunk]]);
	if (i == (dim/2)) printf("=");
	else printf(" ");
	printf("[%.4le]\n", bvect[i]);
      }
    }

    // Now, solve it!
    gsl_matrix_view m = gsl_matrix_view_array(Amat, dim, dim);
    gsl_vector_view b = gsl_vector_view_array(bvect, dim);
    gsl_vector *x = gsl_vector_alloc(dim);
    int s;
    gsl_permutation * p = gsl_permutation_alloc (dim);

    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    double detA;
    detA = gsl_linalg_LU_det(&m.matrix, s);
    if (fabs(detA) < 1e-3) {
      fprintf(stderr, "ERROR solving shell %d\n", n);
      fprintf(stderr, "This *most likely* means that you don't have enough entries to solve.\n");
      fprintf(stderr, "Rerun with the -t switch to see more information.\n");
      ERROR = 1;
    }
    else {
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    }

    if (TESTING) {
      printf("## det(Amat) = %.5le\n", detA);
      printf("## x = \n");
      gsl_vector_fprintf (stdout, x, "##   %g");
    }
    // Now, put in the entries!
    for (i=0; i<dim; ++i) {
      np = tsh->elem[i/Nunk];
      int ent1 = i6tom9[unk[i%Nunk]];
      int ent2 = trans[ent1];
      GL[np].mat[ent1] = gsl_vector_get(x, i);
      GL[np].mat[ent2] = gsl_vector_get(x, i);
    }
    // garbage collection...
    gsl_permutation_free(p);
    gsl_vector_free(x);
  }

  //-- ==== loop through each shell ====


  // ****************************** OUTPUT ***************************

  write_Dij(stdout, Np, GL);

  // ************************* GARBAGE COLLECTION ********************
  delete[] shell_list;
  delete[] GL;

  free_cell(Cmn_list, u_atoms, 0);

  return 0;
}


int parse_unknown(char* arg, int &Nunk, int unk[3]) 
{
  int ERROR = 0;
  if (arg[0] == 0) ERROR = 1;
  else {
    if (arg[1] == 0) {
      Nunk = 1;
      switch (arg[0]) {
      case 'x':
      case 'X': unk[0] = m9toi6[0]; break;
      case 'y':
      case 'Y': unk[0] = m9toi6[4]; break;
      case 'z':
      case 'Z': unk[0] = m9toi6[8]; break;
      default: ERROR = 1;
      }
    }
    else {
      if (arg[2] == 0) {
	Nunk = 3;
	int i, j;
	switch (arg[0]) {
	case 'x':
	case 'X': i=0; break;
	case 'y':
	case 'Y': i=1; break;
	case 'z':
	case 'Z': i=2; break;
	default: ERROR = 1;
	}
	switch (arg[1]) {
	case 'x':
	case 'X': j=0; break;
	case 'y':
	case 'Y': j=1; break;
	case 'z':
	case 'Z': j=2; break;
	default: ERROR = 1;
	}
	if (i==j) ERROR = 1;
	if (!ERROR) {
	  unk[0] = m9toi6[4*i];
	  unk[1] = m9toi6[4*j];
	  unk[2] = m9toi6[3*i+j];
	}
      }
      else
	ERROR = 1;
    }
  }
  return ERROR;
}



inline int unknown(int Nunk, int unk[3], int j) 
{
  if (Nunk==1) return (j==unk[0]);
  if (j==unk[0]) return -1;
  if (j==unk[1]) return -1;
  if (j==unk[2]) return -1;
  return 0;
}

void construct_matrix(point_type *GL, int gunit[MAXop][9], int Nop, 
		      double ginv_i_g[MAXop][Ni][Ni], shell_type *shell,
		      int Nunk, int unk[3],
		      int dim, double* Amat, double* bvect)
{
  int i, j, ji;
  int ni, n, npi, np;
  int nops;
  int g;

  // Since we now do open shells, we need to calculate nops by hand:
  for (nops=0, n=0; n<Nop; ++n)
    if (shell->gR0[n] != -1) ++nops;

  for (i=0; i<dim; ++i) {
    bvect[i] = 0.;
    for (j=0; j<dim; ++j)
      Amat[i*dim+j] = 0.;
  }

  int len = shell->Nelem;
  // Loop over R:
  for (ni=0; ni<len; ++ni) {
    n = shell->elem[ni];
    int* runit = GL[n].Runit;
    int gR[3];
    
    // Set the base value:
    for (i=0; i<Nunk; ++i) Amat[(Nunk*ni+i)*(dim+1)] = nops;

    // Loop over R' = g.R:
    for (g=0; g<Nop; ++g) {
      // see if we've got a point that maps there...
      mult_vect(gunit[g], runit, gR);
      for (npi=0; 
	   (npi<shell->Nelem) && (! equal_vect(gR,GL[shell->elem[npi]].Runit));
	   ++npi) ;
      if (npi != shell->Nelem) {
	// We've got a live one...
	np = shell->elem[npi];
      
	// Now, we loop over i and j (i for R, j for R'):
	for (i=0; i<Nunk; ++i) {
	  for (j=0; j<Ni; ++j) {
	    //	  double alpha = ginv_i_g[g][unk[i]][j];
	    // Note: we "switch" j and i here because of the way
	    //       we've definied ginv_i_g.
	    double alpha = ginv_i_g[g][j][unk[i]];
	    if (unknown(Nunk, unk, j)) {
	      for (ji=0; (ji<Nunk) && (j!=unk[ji]); ++ji) ;
	      // subtract from A
	      Amat[(Nunk*ni+i)*dim + (Nunk*npi+ji)] -= alpha;
	    }
	    else {
	      // add to b
	      bvect[Nunk*ni+i] += GL[np].mat[i6tom9[j]] * alpha;
	    }
	  }
	}
      }
    }
  }
  // Now, divide everything by Nop to make the determinant scaling "nicer"
  double scale = 1./nops;
  for (i=0; i<dim; ++i) {
    bvect[i] *= scale;
    for (j=0; j<dim; ++j) 
      Amat[i*dim+j] *= scale;
  }
}

