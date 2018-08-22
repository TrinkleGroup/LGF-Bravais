#ifndef __SHELL_H
#define __SHELL_H

#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "pointgroup.H"
#include "Dij.H"

//***************************** SUBROUTINES ****************************
typedef struct 
{
  int Nelem;
  int* elem;
  int* gR0;
} shell_type;

void free_shell(shell_type* &shell_list, int Nshell);

const int SHELL_LIST_CLOSED = 0;
const int SHELL_LIST_OPEN = 1;

int gen_shell_list(point_type* p, int Np, int gunit[MAXop][9], int Nop,
                   shell_type* &shell_list, int &Nshell, const int SHELL_TYPE);

inline int gen_shell_list(point_type* p, int Np, int gunit[MAXop][9], int Nop,
			  shell_type* &shell_list, int &Nshell) 
{
  return
    gen_shell_list(p, Np, gunit, Nop, shell_list, Nshell, SHELL_LIST_CLOSED);
}

inline int gen_shell_list_open(point_type* p, int Np, int gunit[MAXop][9],
			       int Nop, shell_type* &shell_list, int &Nshell) 
{
  return
    gen_shell_list(p, Np, gunit, Nop, shell_list, Nshell, SHELL_LIST_OPEN);
}



// Comparison routine used for sorting
// Sorts 'em by magnitude, but also results in a unique sorting
// which is useful for merging purposes.
int __COMPARE_POINT_TYPE__ (const void* a, const void* b) 
{
  double amag, bmag;
  amag = ((point_type *)a)->Rmagn; // I know--it looks weird :)
  bmag = ((point_type *)b)->Rmagn;
  if (dcomp(amag,bmag)) {
    int *aunit = ((point_type *)a)->Runit;
    int *bunit = ((point_type *)b)->Runit;
    int better;
    better = (aunit[0] > bunit[0]);
    if (aunit[0] == bunit[0]) {
      better = (aunit[1] > bunit[1]);
      if (aunit[1] == bunit[1]) {
	better = (aunit[2] > bunit[2]);
	if (aunit[2] == bunit[2]) return 0;
      }
    }
    if (better) return 1;
    else        return -1;
  }
  if (amag>bmag) return 1;
  return -1;
}

inline void sort_points(point_type *p, int Np) 
{
  qsort((void*)p, Np, sizeof(point_type), __COMPARE_POINT_TYPE__);
}


void free_shell(shell_type* &shell_list, int Nshell) 
{
  int i;
  for (i=0; i<Nshell; ++i) {
    delete[] shell_list[i].elem;
    delete[] shell_list[i].gR0;
  }
  delete[] shell_list;
  shell_list = NULL;
}


int gen_shell_list(point_type* p, int Np, int gunit[MAXop][9], int Nop,
		   shell_type* &shell_list, int &Nshell, const int SHELL_TYPE) 
{
  int ERROR = 0;
  int i, j, n;
  int shelled[Np];
  int elem_list[Nop]; // use as temporary storage...
  int Nelem;
  int OPEN = SHELL_TYPE & SHELL_LIST_OPEN;
  for (i=0; i<Np; ++i) shelled[i] = 0;
  
  shell_list = new shell_type[Np];
  Nshell = 0;
  for (i=0; i<Np; ++i)
    // Only bother doing the points that haven't been put in a shell already
    if (! shelled[i]) {
      int* runit = p[i].Runit;
      double rmagn = p[i].Rmagn;
      // New shell...
      shell_type *ts = shell_list+Nshell;
      ts->gR0 = new int[Nop];
      Nelem = 1;        // one entry to start!
      elem_list[0] = i;
      shelled[i] = -1; // it's in a shell.
      // Now, let's go through the group operations.
      int gR[9];
      for (n=0; n<Nop; ++n) {
        mult_vect(gunit[n], runit, gR);
        for (j=i; (j<Np) && (! equal_vect(gR, p[j].Runit)) && 
	       (dcomp(rmagn, p[j].Rmagn)); ++j) ;
        if ( (j==Np) || (! dcomp(rmagn, p[j].Rmagn)) ) {
	  if (!OPEN) {
	    fprintf(stderr, "Could not find matching point!\n");
	    fprintf(stderr, "R(%d) = %3d %3d %3d  ", i, 
		    runit[0], runit[1], runit[2]);
	    fprintf(stderr, "g(%d) = %2d %2d %2d %2d %2d %2d %2d %2d %2d\n",
		    n, gunit[n][0], gunit[n][1], gunit[n][2],
		    gunit[n][3], gunit[n][4], gunit[n][5],
		    gunit[n][6], gunit[n][7], gunit[n][8]);
	    ERROR = 1;
	  }
	  ts->gR0[n] = -1;
        }
        else {
	  ts->gR0[n] = j;
          if (! shelled[j]) {
            // Add it to our list, if we haven't already.
            elem_list[Nelem] = j;
	    ++Nelem;
            shelled[j] = -1; // we put it in a shell.
          }
	}
      }
      // Now, finish putting this info into the shell:
      ts->Nelem = Nelem;
      ts->elem = new int[Nelem];
      for (j=0; j<Nelem; ++j) ts->elem[j] = elem_list[j];
      ++Nshell;
    }

  return ERROR;
}


// Probably only used by fitsuper, but here it is nonetheless:
// close_shell_list.  This takes in a (presumably open, though not
// necessarily) shelllist for a set of points, and proceeds to construct
// the list of new points needed to "close" the shells.  It operates
// on all of the points and shells in place.
// If the shells are already closed, then nothing will happen to anything.
// NOTE: This algorithm deallocates memory associated with p if it needs
// to add points, so one must be warned about keeping copies of a (possibly
// bad) pointer.
// We return 0 if the shell list was closed, and 1 if a change was made.

int close_shell_list(point_type* &p, int &Np, 
		     double cart[9], int gunit[MAXop][9], int Nop,
		     shell_type* &shell_list, int &Nshell) 
{
  int n, i, j, k, d;
  int Nnew=0;
  // first, count how many spaces we have:
  for (n=0; n<Nshell; ++n)
    for (i=0; i<Nop; ++i)
      if (shell_list[n].gR0[i] == -1) ++Nnew;
  if (Nnew == 0) // closed shell
    return 0;

  point_type* pnew = new point_type[Np+Nnew]; // overallocate space
  for (Nnew=0; Nnew<Np; ++Nnew) {
    point_type *tp = p+Nnew;
    point_type *tnp = pnew+Nnew;
    for (d=0; d<3; ++d) tnp->Runit[d] = tp->Runit[d];
    for (d=0; d<3; ++d) tnp->Rcart[d] = tp->Rcart[d];
    tnp->Rmagn = tp->Rmagn;
    for (d=0; d<9; ++d) tnp->mat[d] = tp->mat[d];
  }
  int elem_list[Nop], gR[3];
  int Nelem;
  for (n=0; n<Nshell; ++n) {
    shell_type* tsh = shell_list + n;
    for (Nelem=0; Nelem<tsh->Nelem; ++Nelem)
      elem_list[Nelem] = tsh->elem[Nelem];
    point_type *tp = p+tsh->elem[0];
    for (i=0; i<Nop; ++i)
      if (shell_list[n].gR0[i] == -1) {
	// We need to find an entry for this; but we have to be
	// careful, because this point may already be in our list.
	mult_vect(gunit[i], tp->Runit, gR);
	for (j=0; (j<Nnew) && (! equal_vect(pnew[j].Runit, gR)); ++j) ;
	if (j == Nnew) {
	  point_type *tnp = pnew+j;
	  // we've got a new point to add to the list.
	  for (d=0; d<3; ++d) tnp->Runit[d] = gR[d];
	  mult_vect(cart, gR, tnp->Rcart);
	  tnp->Rmagn = p[tsh->elem[0]].Rmagn;
	  for (d=0; d<9; ++d) tnp->mat[d] = 0;
	  ++Nnew; // increment the point count
	}
	// See if j is already in our list of elements:
	for (k=0; (k<Nelem) && (elem_list[k] != j); ++k) ;
	if (k == Nelem) {
	  elem_list[k] = j;
	  ++Nelem; // increment the count
	}
	// Finally, set gR0 accordingly:
	tsh->gR0[i] = j;
      }
    // Now, we need to copy over Nelem and elem_list:
    delete[] tsh->elem;
    tsh->elem = new int[Nelem];
    for (i=0; i<Nelem; ++i) tsh->elem[i] = elem_list[i];
    tsh->Nelem = Nelem;
  }
  // Finally, we need to copy over our final list of points--which is 
  // NO LONGER SORTED BY MAGNITUDE.  There's no easy way to do that,
  // so we don't bother.
  delete[] p;
  p = new point_type[Nnew];
  // copy back the other way:
  for (Np=0; Np<Nnew; ++Np) {
    point_type *tp = p+Np;
    point_type *tnp = pnew+Np;
    for (d=0; d<3; ++d) tp->Runit[d] = tnp->Runit[d];
    for (d=0; d<3; ++d) tp->Rcart[d] = tnp->Rcart[d];
    tp->Rmagn = tnp->Rmagn;
    for (d=0; d<9; ++d) tp->mat[d] = tnp->mat[d];
  }
  delete[] pnew;
  
  return 1;
}


// I think we'll only use this in gf-polish, but still... here it is:
// This takes a set of unit coordinates, and shells it directly into a
// set of representative points.  The initial list needs to be sorted, 
// and it has to be closed.

int gen_representative_points(double cart[9], int Ncut, int **Runit,
			      int Nop, int gunit[MAXop][9],
			      int &Np, point_type* &p) 
{
  int ERROR = 0;
  int i, j, d, n;
  int shelled[Ncut];
  //  int elem_list[MAXsame][Nop]; // use as temporary storage...
  double cartsq[9];
  for (i=0; i<Ncut; ++i) shelled[i] = 0;

  square(cart, cartsq);

  int MAXp;
  MAXp = (2*Ncut)/Nop; // our initial guess for how many points we'll need
  p = new point_type[MAXp];

  Np = 0;
  for (i=0; i<Ncut; ++i)
    // Only bother doing the points that haven't been put in a shell already
    if (! shelled[i]) {
      int* runit = Runit[i];
      // Add the point to our list:
      for (d=0; d<3; ++d) p[Np].Runit[d] = runit[d];
      mult_vect(cart, runit, p[Np].Rcart);
      p[Np].Rmagn = sqrt(dot(p[Np].Rcart, p[Np].Rcart));
      ++Np;
      if (Np==MAXp) {
	// We have to allocate more space...
	point_type *newp = new point_type[MAXp*2];
	for (j=0; j<Np; ++j) {
	  for (d=0; d<3; ++d) newp[j].Runit[d] = p[j].Runit[d];
	  for (d=0; d<3; ++d) newp[j].Rcart[d] = p[j].Rcart[d];
	  newp[j].Rmagn = p[j].Rmagn;
	}
	MAXp *= 2;
	delete[] p;
	p = newp;
      }

      shelled[i] = -1; // it's in a shell.
      // Now, let's go through the group operations.
      int gR[3];
      for (n=0; n<Nop; ++n) {
        mult_vect(gunit[n], runit, gR);
        for (j=i; (j<Ncut) && (! equal_vect(gR, Runit[j])); ++j) ;

        if ( j==Ncut ) {
	  fprintf(stderr, "Could not find matching point!\n");
	  fprintf(stderr, "R(%d)= %3d %3d %3d | gR= %3d %3d %3d  ", i, 
		  runit[0], runit[1], runit[2],
		  gR[0], gR[1], gR[2]);
	  fprintf(stderr, "g(%d) = %2d %2d %2d %2d %2d %2d %2d %2d %2d\n",
		  n, gunit[n][0], gunit[n][1], gunit[n][2],
		  gunit[n][3], gunit[n][4], gunit[n][5],
		  gunit[n][6], gunit[n][7], gunit[n][8]);
	  ERROR = 1;
        }
        else {
	  shelled[j] = -1; // we put it in a shell.
	}
      }
    }

  return ERROR;
}

#endif
