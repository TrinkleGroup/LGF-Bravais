/*
  Program: cont2posfor.C
  Author:  D. Trinkle
  Date:    October 9, 2004
  Purpose: Utility code to read in a cell, and a CONTCAR/POSCAR file
           and convert it into the positions for a posfor file.  It attempts
	   to make sure that atoms near the edge of the supercell are
	   translated accordingly by assuming they should be in a bulk-like
	   environment.

  Param.:  <cell> <contcar>
           cell:    cell file (see below for format)
           contcar: output positions from VASP

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
	   VERBOSE: ??
	   TESTING: output practically everything as we do it.

  Algo.:   Read in our lattice, and calculate.

  Output:  POSFOR format
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "cell.H"

//***************************** SUBROUTINES ****************************

typedef struct {double u[3], r[3], rmagn;} cart_type;

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "<cell> <contcar> [ucenter.1 ucenter.2]";

const int NFLAGS = 2;
const char USERFLAGLIST[NFLAGS] = {'u', 'r'}; // flag characters.

const char* ARGEXPL = 
"  cell:    cell file (see below for format)\n\
  contcar: output positions from VASP\n\
  ucenter: optional center of cell for determining edges; default: 0,0\n\
  -u       output unit coordinates instead\n\
  -r       output distance from center";

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
==== cell ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, d, n; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // not used

  char* args[NUMARGS+128]; // extra "buffer"
  int flagon[NFLAGS]; // We use this to determine which flags are on.

  for (i=0; i<NFLAGS; ++i) flagon[i] = 0;
  // Read our commandline.
  int Nargs=-NUMARGS; // variable number of arguments...
  ERROR = parse_commandline_var(argc, argv, Nargs, args,
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

  int UNIT = flagon[0];
  int RMAGN = flagon[1];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  // Command line parameters: 
  // Let's pull off the args:
  char* cell_name = args[0];
  char* contcar_name = args[1];
  double ucenter[3] = {0.,0.,0.};
  if (Nargs >= 4) {
    // optional center shift:
    sscanf(args[2], "%lf", ucenter);
    sscanf(args[3], "%lf", ucenter+1);
    for (d=0; d<2; ++d) {
      if (ucenter[d] < 0.) ucenter[d] += 1.;
      if (ucenter[d] >= 1.) ucenter[d] -= 1.;
    }
  }
  double cart[9], cart_inv[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
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
  //-- ==== cell ====
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }
  if (TESTING)
    verbose_output_cell(cart, crystal, Cmn_list, u_atoms, Natoms);

  careful_inverse(cart, cart_inv);

  //++ ==== contcar ====
  int Np;
  double super[9], a0;
  cart_type* p;
  
  infile = myopenr(contcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", contcar_name);
    exit(ERROR_NOFILE);
  }
  
  // comment line
  nextnoncomment(infile, dump, sizeof(dump));
  // a0
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%lf", &a0);
  for (d=0; d<3; ++d) {
    // lattice:
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%lf %lf %lf", super+d, super+d+3, super+d+6);
  }
  mult(super, a0, super);
  // Number of atoms:
  nextnoncomment(infile, dump, sizeof(dump));
  char *startp, *endp = dump;
  Np=0;
  do {
    startp = endp;
    Np += strtol(startp, &endp, 10);
  } while (startp != endp);
  if (Np < 1) {
    fprintf(stderr, "Need at least one atom.\n");
    exit(ERROR_BADFILE);
  }
  p = new cart_type[Np];
  // go till we get the atoms:
  nextnoncomment(infile, dump, sizeof(dump));
  if ( (dump[0] == 's') || (dump[0] == 'S') )
    nextnoncomment(infile, dump, sizeof(dump));
  // now, determine if we have cartesian or direct coord:
  int DIRECT = ( (dump[0] == 'D') || (dump[0] == 'd') );
  double super_inv[9];
  careful_inverse(super, super_inv);
  for (n=0; n<Np; ++n) {
    cart_type *tp = p + n;
    double r[3];
    nextnoncomment(infile, dump, sizeof(dump));
    if (feof(infile)) break;
    sscanf(dump, "%lf %lf %lf", r, r+1, r+2);
    if (!DIRECT)
      mult_vect(super_inv, r, tp->u);
    else
      for (d=0; d<3; ++d) tp->u[d] = r[d];
    // Now, shift by center...
    for (d=0; d<3; ++d) {
      tp->u[d] -= ucenter[d];
      if (tp->u[d]  >  0.5) tp->u[d] -= 1.;
      if (tp->u[d] <= -0.5) tp->u[d] += 1.;
    }
    mult_vect(super, tp->u, tp->r);
  }
  myclose(infile);
  ERROR = (n<Np);
  //-- ==== latt(R) ====

  if (ERROR) {delete p; exit(ERROR);}

  double rshift[3];
  mult_vect(super, ucenter, rshift); // we'll need this to shift back...

  // ***************************** ANALYSIS **************************
  // 0. We're going to tag which atoms are "okay" and which are not.
  int Ncheck = 2; // which axes we check  
  int okay[Np];
  double toler[Ncheck];
  a0 = exp(log(det(cart))/3.); // estimate of our atomic dist...
  for (d=0; d<Ncheck; ++d)
    toler[d] = a0/sqrt(super[d]*super[d]+super[d+3]*super[d+3]
		       +super[d+6]*super[d+6]);
  
  for (n=0; n<Np; ++n) {
    cart_type *tp = p + n;
    // only use the first two as real criteria...
    for (okay[n]=1, d=0; d<Ncheck; ++d)
      okay[n] = okay[n] && ((0.5-fabs(tp->u[d])) > toler[d]);
  }
  
  // 1. Now, let's do our check.
  int du0[Ncheck], du1[Ncheck], du[Ncheck];
  int Ntry = 1<<Ncheck; // 2^Ncheck
  double rtry[Ntry][3], utry[3];
  for (n=0; n<Np; ++n)
    if (! okay[n]) {
      cart_type *tp = p+n;
      // Construct our possiblities...
      for (d=0; d<Ncheck; ++d) {
	if (tp->u[d] < 0) {
	  du0[d] = 0; du1[d] = 1;
	} else {
	  du0[d] = -1; du1[d] = 0;
	}
	du[d] = du0[d];
      }
      for (d=Ncheck; d<3; ++d)
	utry[d] = tp->u[d];
      for (i=0; i<Ntry; ++i) {
	for (d=0; d<Ncheck; ++d)
	  utry[d] = tp->u[d] + du[d];
	mult_vect(super, utry, rtry[i]);
	++(du[0]);
	if (du[0] > du1[0]) {
	  du[0] = du0[0];
	  ++(du[1]);
	  if (du[1] > du1[1]) {
	    du[1] = du0[1];
	    ++(du[2]);
	  }
	}
      }
      // Test them to find the best one.  We do this by finding the
      // closest neighbor, and then seeing how "bulk-like" it is:
      int ibest=-1;
      double dev_min = a0;
      for (i=0; i<Ntry; ++i) {
	int closest=-1;
	double dr[3], Rmin=3*a0;
	for (j=0; j<Np; ++j) 
	  if (okay[j] && (j!=n)) {
	    for (d=0; d<3; ++d) dr[d] = rtry[i][d] - p[j].r[d];
	    if (sqrt(dot(dr,dr)) < Rmin) {
	      Rmin=sqrt(dot(dr,dr));
	      closest=j;
	    }
	  }
	if (closest == -1) continue;
	for (d=0; d<3; ++d) dr[d] = rtry[i][d] - p[closest].r[d];
	mult_vect(cart_inv, dr, utry);
	for (d=0; d<3; ++d) utry[d] = round(utry[d]);
	double t[3];
	mult_vect(cart, utry, t);
	for (d=0; d<3; ++d) dr[d] -= t[d];
	if (sqrt(dot(dr,dr)) < dev_min) {
	  dev_min = sqrt(dot(dr,dr));
	  ibest = i;
	}
      }
      if (ibest == -1) {
	fprintf(stderr, "Could not find a reasonable bulk-like location for:\n");
	fprintf(stderr, "  %d u= %.8lf %.8lf %.8lf  r= %.8lf %.8lf %.8lf\n", 
		n+1, tp->u[0],tp->u[1], tp->u[2], tp->r[0],tp->r[1], tp->r[2]);
	break;
      }
      for (d=0; d<3; ++d)
	tp->r[d] = rtry[ibest][d];
      mult_vect(super_inv, tp->r, tp->u);
      okay[n] = 1;
    }

  // before shifting back, determine distance from center:
  double cart2[9];
  square(super, cart2);
  // zero out the elements from the threading direction:
  cart2[2] = 0; cart2[5] = 0;
  cart2[6] = 0; cart2[7] = 0; cart2[8] = 0;
  for (n=0; n<Np; ++n) {
    cart_type *tp = p+n;
    tp->rmagn = sqrt(magnsq(cart2, tp->u));
  }

  // shift back:
  for (n=0; n<Np; ++n) {
    cart_type *tp = p+n;
    for (d=0; d<3; ++d) tp->r[d] += rshift[d];
    for (d=0; d<3; ++d) tp->u[d] += ucenter[d];
  }

  // determine threading direction:
  double temp[9];
  int t_unit[3];
  mult(cart_inv, super, temp);
  for (d=0; d<3; ++d) t_unit[d] = lround(temp[2+3*d]);

  // ****************************** OUTPUT ***************************
  printf("%d %d %d %d\n", Np, t_unit[0], t_unit[1], t_unit[2]);
  for (n=0; n<Np; ++n) {
    if (UNIT)
      printf("%18.12lf %18.12lf %18.12lf", p[n].u[0], p[n].u[1], p[n].u[2]);
    else
      printf("%18.12lf %18.12lf %18.12lf", p[n].r[0], p[n].r[1], p[n].r[2]);
    if (RMAGN)
      printf("  %18.12lf\n", p[n].rmagn);
    else
      printf("\n");
  }
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] p;
  free_cell(Cmn_list, u_atoms, Natoms);

  return 0;
}
