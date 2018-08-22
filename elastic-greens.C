/*
  Program: elastic-greens.C
  Author:  D. Trinkle
  Date:    February 9, 2004
  Purpose: Calculate the anisotropic elastic greens function for a particular
           direction given:
	   1. elastic constants Cmn, crystal class c

	   Note: Cij are read in in GPa, and converted to meV/A^3.

  Param.:  <cell> <grid>
           cell:  cell file (see below for format)
           grid:  "grid" of direction cosines for which to compute Gij(R)

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
	   
	   ==== grid ====
	   N         # Number of points
	   l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)
	   l2 m2 n2
	   ...
	   lN mN nN
	   ==== grid ====

  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   Read in everything, and just go *nuts* a calculatin'.
           The only thing we'll need from our cell file are the elastic
	   constants.

           For each direction point, we construct a normalized 
	   perpendicular direction m0.  This is done by choosing the
	   largest vector of the set {x-(x.R)R, y-(y.R)R, z-(z.R)R} and
	   normalizing it.  We then construct n0 = t x m0.

	   We then define the vectors z(theta) as:

	   z(theta) =  m0*cos(theta) + n0*sin(theta)

	   and the matrix (ab)_ij as

	   (zz)_ij = sum  z_k C_ikjl z_l
	              kl

           We have to one integrals from 0 to Pi:

	                   1    Pi      -1
	   g_ij(theta) = ----- Int  (zz)   dtheta
	                 4Pi^2  0       ij

	   We integrate using a stepping scheme based on Simpson's
	   extended rule; basically, to integrate f(x) from 0 to x,
	   we calculate f(x) at a grid of points, and if F(N-1) is
	   the integral to x-h, and f[i] = f(x-ih), then:

   	     F(N) = F(N-1) + h*SUM(i=0..3, int_weight[i]*f[i])

	   which works very well.  To get the first two points, we 
	   actually have to evaluate f at -h and -2h, and then
	   start with F(0) = 0.  It works (go figure).

	   We do 16384 integration steps (2^14... woohoo!) to make
	   sure that we have something reasonable :)

  Output:  I figure for now, I'll just output g_ij for each direction
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
#include "elastic.H"
#include "cell.H"
#include "integrate.H"

// 1 GPa = 1e9 J/m^3 = 10/1.60217733 meV/A^3 = 6.241506363 meV/A^3
// or 0.00624 eV/A^3.
// This makes the units of gij(R) A^3 / eV; multiplied by a force
// in eV/A and divided by the distance in A, gives the displacement
// in A.
const double meV_GPa = 10/1.60217733;  // = 10/(ec in coulumbs x 10^19)
const double eV_GPa = 0.01/1.60217733;  // = 0.01/(ec in coulumbs x 10^19)

// This is the permutation matrix; eps[i][j][k] =
//  1: if ijk is an even permutation of (012)
// -1: if ijk is an odd permutation of (012)
//  0: otherwise
const int eps[3][3][3] = {
  {{0,0,0}, {0,0,1}, {0,-1,0}},
  {{0,0,-1}, {0,0,0}, {1,0,0}},
  {{0,1,0}, {-1,0,0}, {0,0,0}}
};

//****************************** SUBROUTINES ****************************

void construct_mn(double lmn[3], double m0[3], double n0[3]);


void z_theta(double theta, double m[3], double n[3], double zt[3]) 
{
  zt[0] = cos(theta)*m[0] + sin(theta)*n[0];
  zt[1] = cos(theta)*m[1] + sin(theta)*n[1];
  zt[2] = cos(theta)*m[2] + sin(theta)*n[2];
}


// Calculate (zz)
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

void print_mat (double a[9]) 
{
  int i, j;
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j)
      printf("  %10.5lf", a[index(i,j)]);
    printf("\n");
  }
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <grid>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'g'}; // Would be the flag characters.

const char* ARGEXPL = 
"  cell:  cell file (see below for format)\n\
  grid:  \"grid\" of direction cosines for which to compute Gij(R)\n\
\n\
  -g               don't convert (default: assume GPa and convert to eV/A^3)";

const char* FILEEXPL =
"==== cell ====\n\
a0                               # Scale factor for unit cell\n\
a1.x a1.y a1.z                   # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.\n\
Natoms                           # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3                   # Atom locations, in direct coord.\n\
...\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== grid ====\n\
N         # Number of points\n\
l1 m1 n1  # Direction cosines (l^2 + m^2 + n^2 = 1)\n\
...\n\
lN mN nN\n\
==== infile ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 16384; // 2^14, default (gets turned into Nsteps for int.)

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
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "Crystal classes:\n%s\n", CRYSTAL_CLASS);
      fprintf(stderr, "\nElastic constants ordering:\n");
      for (k=0; k<NCLASSES; ++k) {
        fprintf(stderr, "  Class %2d (%d):", k, class_len[k]);
        for (i=0; i<class_len[k]; ++i)
          printf(" C_%2d", class_Cij[k][i]);
        fprintf(stderr, "\n");
      }
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  int Nsteps;
  
  // We're going to use the number of steps according to our preferred
  // amount of memory allocation.
  Nsteps = MEMORY;
  if (Nsteps < 4) {
    fprintf(stderr, "Nsteps (MEMORY = %d) must be 4 or larger.\n", Nsteps);
    exit(ERROR_BADFILE);
  }

  int EV_CONV = !(flagon[0]);  // Convert from GPa to eV/A^3?
  
  // ****************************** INPUT ****************************
  // Command line parameters:
  char *cell_name, *infile_name;
  // Let's pull off the args:
  cell_name = args[0];
  infile_name = args[1];

  char dump[512];
  FILE* infile;

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  double Cijkl[9][9];
  int Natoms = NO_ATOMS;
  double** u_atoms = NULL;

  //++ ==== cell ====
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(ERROR_NOFILE);
  }
  ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms);
  myclose(infile);
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }
  if (TESTING)
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
  //-- ==== cell ====

  if (EV_CONV) {
    // Convert Cmn_list from GPa to eV/ang^3:
    for (n=0; n<class_len[crystal]; ++n)
      Cmn_list[n] *= eV_GPa;
    if (TESTING) {
      printf("## Converted elastic constants:\n");
      verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);
    }
  }
  // Calculate elastic constant matrix:
  make_Cijkl(crystal, Cmn_list, Cijkl);


  //++ ==== grid ====
  infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", infile_name);
    exit(ERROR_NOFILE);
  }

  int Npoints;
  double** lmn_list;
  
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d", &Npoints);
  if (Npoints <= 0) {
    fprintf(stderr, "No points listed in %s\n", infile_name);
    ERROR = ERROR_BADFILE;
  }
  else lmn_list = new double* [Npoints];
  for (n=0; (!ERROR) && (n<Npoints) && (!feof(infile)); ++n) {
    double magn;
    lmn_list[n] = new double[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", 
	   &(lmn_list[n][0]), &(lmn_list[n][1]), &(lmn_list[n][2]) );
    // Normalize.
    magn = sqrt(dot(lmn_list[n],lmn_list[n]));
    if (dcomp(magn, 0.)) {
      fprintf(stderr, 
	      "%s at line %d: %.5le %.5le %.5le with magn = %.5le\n",
	      infile_name, i+1,
	      lmn_list[0][0], lmn_list[0][1], lmn_list[0][2], magn);
      ERROR = ERROR_BADFILE;
    }
    else {
      magn = 1./magn;
      for (i=0; i<3; ++i) lmn_list[n][i] *= magn;
    }
  }
  // Make sure we read enough points...
  if (n != Npoints) ERROR = ERROR_BADFILE;
  myclose(infile);
  //-- ==== grid ====

  // If we've thrown an error already, then get out.
  if (ERROR) {
    exit(ERROR);
  }
  

  // ***************************** ANALYSIS **************************

  // Spit out our first bit of info:
  printf("%d 9 # Number of points: l m n  g[xx xy xz  yx yy yz  zx zy zz]",
	 Npoints);
  if (EV_CONV)
    printf(" (A^3/eV)\n");
  else
    printf(" (--)\n");
  // For each point...
  for (n=0; n<Npoints; ++n) {
    double *lmn;
    lmn = lmn_list[n]; // Our specific direction.
    
    // Our integration variables
    double m0[3], n0[3];
    // First, construct m0 and n0:
    construct_mn(lmn, m0, n0);

    if (TESTING) {
      printf("##\n## Normalized vectors:\n");
      printf("## Direction (%.5lf %.5lf %.5lf)\n", lmn[0],lmn[1],lmn[2]);
      printf("## m0        (%.5lf %.5lf %.5lf)\n", m0[0],m0[1],m0[2]);
      printf("## n0        (%.5lf %.5lf %.5lf)\n", n0[0],n0[1],n0[2]);
    }
    
    // Now some evaluating of integrals :)
    double theta;
    double dtheta;
    dtheta = M_PI / Nsteps;
    double zt[3];
    double zz[9], zzi[9];
    double detzz;
    
    // We have to integrate this function:
    double gint[9];

    // Function evaluations, stored for integration purposes.
    double zzi_old[4][9];

    // First, prime the integration pump:
    for (k=1; k<=3; ++k) {
      theta = -(k-1)*dtheta;
      // Eval (zz), (zz)^-1
      z_theta(theta, m0, n0, zt);
      a_mult_a(zt, Cijkl, zz);
      detzz = 1./inverse(zz, zzi);
      for (i=0; i<9; ++i) zzi[i] *= detzz;
      // Now, put into the function evaluations:
      for (i=0; i<9; ++i) zzi_old[k][i] = zzi[i];
      // And HERE's where we'd integrate, if we wanted to... :)
    }
    
    // Now we've got EVERYTHING, let's integrate!
    // theta = 0 is easy...
    for (i=0; i<9; ++i) {
      gint[i] = 0.;
    }
    
    for (k=1; k<=Nsteps; ++k) {
      theta = k*dtheta;
      // Eval (zz), (zz)^-1
      z_theta(theta, m0, n0, zt);
      a_mult_a(zt, Cijkl, zz);
      detzz = 1./inverse(zz, zzi);
      for (i=0; i<9; ++i) zzi[i] *= detzz;

      // Now, put into the function evaluations:
      for (i=0; i<9; ++i) zzi_old[0][i] = zzi[i];
      // Now, we can integrate!
      for (i=0; i<9; ++i) {
	for (j=0; j<4; ++j) {
	  gint[i]    += dtheta*int_weight[j]*zzi_old[j][i];
	}
      }
      // Now, we slide down all of our "old" values:
      for (j=3; j>0; --j)
	for (i=0; i<9; ++i) {
	  zzi_old[j][i] = zzi_old[j-1][i];
	}
      // And do it all again!
    }
    // Scale by 1/(4*Pi^2)
    for (i=0; i<9; ++i) {
      gint[i] *= 0.25*M_1_PI*M_1_PI;
    }      

    // ****************************** OUTPUT ***************************
    // Human readable (sorta) first:
    
    if (VERBOSE) {
      // What, exactly, is "verbose" output...?
    }
    // Spit back out the lmn point, and include the gij matrix values:
    printf("%15.12lf %15.12lf %15.12lf", lmn[0], lmn[1], lmn[2]);
    for (i=0; i<3; ++i) {
      printf(" ");
      for (j=0; j<3; ++j) printf(" %.8le", gint[3*i+j]);
    }
    printf("\n");
  }
  
  // ************************* GARBAGE COLLECTION ********************
  for (n=0; n<Npoints; ++n)
    delete lmn_list[n];
  delete[] lmn_list;
  
  free_cell(Cmn_list, u_atoms, Natoms);
  delete[] Cmn_list;

  return 0;
}
  

const double delta[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  
//=========================== construct_mn ===========================
void construct_mn(double lmn[3], double m0[3], double n0[3]) 
{
  double eproj[3][3];
  double magn[3], maxmagn;
  int i, j, k, maxi;
  
  maxmagn = 0; maxi = 0;
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j)
      eproj[i][j] = delta[i][j] - lmn[i]*lmn[j];
    magn[i] = dot(eproj[i],eproj[i]);
    if (maxmagn < magn[i]) {
      maxmagn = magn[i];
      maxi = i;
    }
  }
  // Now, make m0 = eproj[maxi][j], properly normalized:
  maxmagn = 1./sqrt(maxmagn);
  for (j=0; j<3; ++j)
    m0[j] = eproj[maxi][j] * maxmagn;

  // Finally, n0 = lmn x m0:
  for (i=0; i<3; ++i) {
    n0[i] = 0;
    for (j=0; j<3; ++j)
      for (k=0; k<3; ++k)
	n0[i] += eps[i][j][k] * lmn[j] * m0[k];
  }
  // Normalize:
  maxmagn = 1./sqrt(dot(n0,n0));
  for (i=0; i<3; ++i) n0[i] *= maxmagn;
}
