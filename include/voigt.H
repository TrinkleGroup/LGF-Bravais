#ifndef __VOIGT_H
#define __VOIGT_H

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "dcomp.H"
#include "matrix.H"

//***************************** SUBROUTINES ****************************

const char coord[3] = {'x','y','z'};

const int Ni = 6;
const double i_mat[Ni][9] = {
  {1,0,0, 0,0,0, 0,0,0}, // [0]
  {0,1,0, 1,0,0, 0,0,0}, // [1]
  {0,0,1, 0,0,0, 1,0,0}, // [2]
  {0,0,0, 0,1,0, 0,0,0}, // [3]
  {0,0,0, 0,0,1, 0,1,0}, // [4]
  {0,0,0, 0,0,0, 0,0,1}};// [5]
const char* v_coord[Ni] = {"xx", "xy", "xz", "yy", "yz", "zz"};
const int i2rank2[Ni][2] = {{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};
const int m9toi6[9] = {0,1,2,1,3,4,2,4,5};
const int i6tom9[Ni] = {0,1,2,4,5,8};
const int trans[9] = {0,3,6, 1,4,7, 2,5,8};

const int k2i_mat[6] = {0,1,2,1,2,0};
const int k2j_mat[6] = {0,1,2,2,0,1};

inline void k2ij (int k, int& i, int& j) 
{ i = k2i_mat[k]; j = k2j_mat[k]; }

const int ij2k_mat[3][3] = {{0,5,4},
			    {5,1,3},
			    {4,3,2}};

inline int ij2k (int i, int j) 
{ return ij2k_mat[i][j]; }


// For our 4th rank tensors, with Cauchy symmetry:
const int Nii = 15;
const int ii2rank4[Nii][4] = {
  {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2},
  {0,0,2,2}, {0,1,1,1}, {0,1,1,2}, {0,1,2,2}, {0,2,2,2},
  {1,1,1,1}, {1,1,1,2}, {1,1,2,2}, {1,2,2,2}, {2,2,2,2}};

const int rank42ii[3][3][3][3] = {
{{{0,1,2},{1,3,4},{2,4,5}},
 {{1,3,4},{3,6,7},{4,7,8}},
 {{2,4,5},{4,7,8},{5,8,9}}},
{{{1,3,4},{3,6,7},{4,7,8}},
 {{3,6,7},{6,10,11},{7,11,12}},
 {{4,7,8},{7,11,12},{8,12,13}}},
{{{2,4,5},{4,7,8},{5,8,9}},
 {{4,7,8},{7,11,12},{8,12,13}},
 {{5,8,9},{8,12,13},{9,13,14}}}};

const int mmtoii[9][9] = {
  {0,1,2,1,3,4,2,4,5},
  {1,3,4,3,6,7,4,7,8},
  {2,4,5,4,7,8,5,8,9},
  {1,3,4,3,6,7,4,7,8},
  {3,6,7,6,10,11,7,11,12},
  {4,7,8,7,11,12,8,12,13},
  {2,4,5,4,7,8,5,8,9},
  {4,7,8,7,11,12,8,12,13},
  {5,8,9,8,12,13,9,13,14}};

   
void gen_ginv_i_g (int i, double gcart[9], double gcart_inv[9], 
		   double ginv_i_g[Ni]) 
{
  int j;
  double gmult[9], temp[9];
  if ( (i<0) || (i>=Ni) ) {
    fprintf(stderr, "Bad i value (%d) in gen_ginv_i_g\n", i);
    exit(-1);
  }
  mult(gcart_inv, i_mat[i], temp);
  mult(temp, gcart, gmult);
  for (j=0; j<Ni; ++j) 
    ginv_i_g[j] = gmult[i6tom9[j]];
}

void gen_g_i_ginv (int i, double gcart[9], double gcart_inv[9], 
		   double g_i_ginv[Ni]) 
{
  int j;
  double gmult[9], temp[9];
  if ( (i<0) || (i>=Ni) ) {
    fprintf(stderr, "Bad i value (%d) in gen_g_i_ginv\n", i);
    exit(-1);
  }
  mult(gcart, i_mat[i], temp);
  mult(temp, gcart_inv, gmult);
  for (j=0; j<Ni; ++j) 
    g_i_ginv[j] = gmult[i6tom9[j]];
}

#endif
