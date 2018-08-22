#ifndef __EXPR_FUNCTIONS_H
#define __EXPR_FUNCTIONS_H
/*
  Program: expr-functions.H
  Author:  D. Trinkle
  Date:    April 27, 2004
  Purpose: Built in functions and constants for expr-compile.
*/

#include <math.h>
#include "dcomp.H"

// ****************************** STRUCTURES ***************************
// New type: we take an array of *pointers* to the individual arguments
// This facilitates compilation of efficient code, as we can fully
// construct the argument list at compile time, rather than at 
// execution time.
typedef double numerical_function_type (int Nargs, double** arglist);

typedef struct 
{
  char* name;
  int Nargs;
  numerical_function_type *num_func;
  int func_type;
} num_func_type;

typedef struct 
{
  char* name;
  double value;
} builtin_const_type;


// ***************************** SUBROUTINES ***************************
double num_neg (int Nargs, double** arglist) 
{
  return -**arglist;
}

double num_abs (int Nargs, double** arglist) 
{
  if (arglist[0] >= 0) return **arglist;
  else                 return -**arglist;
}

double num_if_gt (int Nargs, double** arglist) 
{
  if (**arglist > 0.) return **(arglist+1);
  else                return **(arglist+2);
}

double num_if_geq (int Nargs, double** arglist) 
{
  if (**arglist > -TOLER) return **(arglist+1);
  else                    return **(arglist+2);
}

double num_if_zero (int Nargs, double** arglist) 
{
  if (zero(**arglist)) return **(arglist+1);
  else                 return **(arglist+2);
}


double num_gt (int Nargs, double** arglist) 
{
  if (**arglist > **(arglist+1)) return 1.;
  else return 0.;
}

double num_geq (int Nargs, double** arglist) 
{
  if (**arglist >= **(arglist+1)) return 1.;
  else return 0.;
}

double num_lt (int Nargs, double** arglist) 
{
  if (**arglist < **(arglist+1)) return 1.;
  else return 0.;
}

double num_leq (int Nargs, double** arglist) 
{
  if (**arglist <= **(arglist+1)) return 1.;
  else return 0.;
}

double num_eq (int Nargs, double** arglist) 
{
  if (**arglist == **(arglist+1)) return 1.;
  else return 0.;
}

double num_neq (int Nargs, double** arglist) 
{
  if (**arglist != **(arglist+1)) return 1.;
  else return 0.;
}


double num_if2 (int Nargs, double** arglist) 
{
  if (! zero(**arglist)) return **(arglist+1);
  else                   return 0.;
}

double num_if3 (int Nargs, double** arglist) 
{
  if (! zero(**arglist)) return **(arglist+1);
  else                   return **(arglist+2);
}

// (if (test) >0 ==0 <0)
double num_if4 (int Nargs, double** arglist) 
{
  if (zero(**arglist)) return **(arglist+2);
  if (**arglist > 0.)  return **(arglist+1);
  else                 return **(arglist+3);
}


double num_add (int Nargs, double** arglist) 
{
  double result = **arglist;
  for (int i=1; i<Nargs; ++i) result += **(arglist+i);
  return result;
}

double num_subt (int Nargs, double** arglist) 
{
  return **arglist - **(arglist+1);
}

double num_mult (int Nargs, double** arglist) 
{
  double result = **arglist;
  for (int i=1; i<Nargs; ++i) result *= **(arglist+i);
  return result;
}

double num_div (int Nargs, double** arglist) 
{
  return **arglist/ **(arglist+1);
}

double num_inv (int Nargs, double** arglist) 
{
  return 1./ **arglist;
}

double num_pow (int Nargs, double** arglist) 
{
  return exp(log(fabs(**arglist))* **(arglist+1));
}

double num_sqrt (int Nargs, double** arglist) 
{
  return sqrt(fabs(**arglist));
}

double num_exp (int Nargs, double** arglist) 
{
  return exp(**arglist);
}

double num_log (int Nargs, double** arglist) 
{
  return log(fabs(**arglist));
}

double num_sin (int Nargs, double** arglist) 
{
  return sin(**arglist);
}

double num_cos (int Nargs, double** arglist) 
{
  return cos(**arglist);
}

double num_tan (int Nargs, double** arglist) 
{
  return tan(**arglist);
}

double num_csc (int Nargs, double** arglist) 
{
  return 1./sin(**arglist);
}

double num_sec (int Nargs, double** arglist) 
{
  return 1./cos(**arglist);
}

double num_cot (int Nargs, double** arglist) 
{
  return 1./tan(**arglist);
}


double num_asin (int Nargs, double** arglist) 
{
  return asin(**arglist);
}

double num_acos (int Nargs, double** arglist) 
{
  return acos(**arglist);
}

double num_atan (int Nargs, double** arglist) 
{
  return atan(**arglist);
}

double num_atan2 (int Nargs, double** arglist) 
{
  return atan2(**arglist, **(arglist+1));
}

double num_sinh (int Nargs, double** arglist) 
{
  return sinh(**arglist);
}

double num_cosh (int Nargs, double** arglist) 
{
  return cosh(**arglist);
}

double num_tanh (int Nargs, double** arglist) 
{
  return tanh(**arglist);
}

// returns a random number between [0,1]
const double __RAND_SCALE__ = (1./(double)( (1<<31) - 1));
double num_rand (int Nargs, double** arglist) 
{
  return __RAND_SCALE__ * (double)(random());
}

// returns a random number between [-1,1]
const double __SRAND_SCALE__ = (2./(double)( (1<<31) - 1));
double num_srand (int Nargs, double** arglist) 
{
  return __SRAND_SCALE__ * (double)(random()) - 1.;
}

// returns a gaussian random number with width given by the first
// argument.  We do this with the Box-Mueller algorithm, and use
// srand as a shortcut to give random numbers between [-1,1]
// We also store the "next" number to make things a little more effecient.
int NUM_GRAND_NEED_NEW = 1;
double NUM_GRAND_NEXT_RAND = 0.;

double num_grand (int Nargs, double** arglist) 
{
  double a = **arglist;
  if (NUM_GRAND_NEED_NEW) {
    double x, y, magn2, scale;
    do {
      x = num_srand(Nargs, arglist);
      y = num_srand(Nargs, arglist);
      magn2 = x*x+y*y;
    } while (magn2 > 1.0);
    scale = a* sqrt(-2 * log(magn2) / magn2);
    NUM_GRAND_NEED_NEW = 0; // we've got another random number on deck.
    NUM_GRAND_NEXT_RAND = (scale*y);
    return (scale*x);
  }
  else {
    NUM_GRAND_NEED_NEW = 1;
    return NUM_GRAND_NEXT_RAND;
  }
}



// ****************************** CONSTANTS ****************************

// different function types--needed for optimization
const int FUNC_ASSOCIATIVE = 1;
const int FUNC_SYMMETRIC = 2;
const int FUNC_STATE = 4;

const int FUNCTION_COUNT = 39;
const num_func_type NUM_FUNCTIONS[FUNCTION_COUNT] = {
  {"-", 1, num_neg, 0},
  {"||", 1, num_abs, 0},
  {"if>0", 3, num_if_gt, 0},
  {"if>=0", 3, num_if_geq, 0},
  {"if==0", 3, num_if_zero, 0},
  {">", 2, num_gt, 0},
  {">=", 2, num_geq, 0},
  {"<", 2, num_lt, 0},
  {"<=", 2, num_leq, 0},
  {"==", 2, num_eq, 0},
  {"!=", 2, num_neq, 0},
  {"if", 2, num_if2, 0},
  {"if", 3, num_if3, 0},
  {"if", 4, num_if4, 0},
  {"+", -2, num_add, FUNC_ASSOCIATIVE | FUNC_SYMMETRIC},
  {"-", 2, num_subt, 0},
  {"*", -2, num_mult, FUNC_ASSOCIATIVE | FUNC_SYMMETRIC},
  {"/", 2, num_div, 0},
  {"inv", 1, num_inv, 0},
  {"^", 2, num_pow, 0},
  {"sqrt", 1, num_sqrt, 0},
  {"exp", 1, num_exp, 0},
  {"log", 1, num_log, 0},
  {"sin", 1, num_sin, 0},
  {"cos", 1, num_cos, 0},
  {"tan", 1, num_tan, 0},
  {"csc", 1, num_csc, 0},
  {"sec", 1, num_sec, 0},
  {"cot", 1, num_cot, 0},
  {"asin", 1, num_asin, 0},
  {"acos", 1, num_acos, 0},
  {"atan", 1, num_atan, 0},
  {"atan2", 2, num_atan2, 0},
  {"sinh", 1, num_sinh, 0},
  {"cosh", 1, num_cosh, 0},
  {"tanh", 1, num_tanh, 0},
  {"rand", 0, num_rand, FUNC_STATE},
  {"srand", 0, num_srand, FUNC_STATE},
  {"grand", 1, num_grand, FUNC_STATE}
};

 
// state tests
inline int is_symm_func(const int &func_index) 
{ if (func_index < 0) return 0;
  return (NUM_FUNCTIONS[func_index].func_type & FUNC_SYMMETRIC); }

inline int is_assoc_func(const int &func_index) 
{ if (func_index < 0) return 0;
  return (NUM_FUNCTIONS[func_index].func_type & FUNC_ASSOCIATIVE); }

inline int is_func_withstate(const int &func_index) 
{ if (func_index < 0) return 0;
  return (NUM_FUNCTIONS[func_index].func_type & FUNC_STATE); }



// three "special" function types:
// const num_func_type symbol_load={"S", 0, NULL};
// const num_func_type const_load={"C", 0, NULL};
const int symbol_load=-1;
const int const_load=-2;
const int BAD_FUNC=-3;

inline int is_symbol_load(const int &func_index) 
{ return (func_index == symbol_load); }

inline int is_const_load(const int &func_index) 
{ return (func_index == const_load); }

inline int is_BAD_FUNC(const int &func_index) 
{ return (func_index == BAD_FUNC); }


inline int argmatch(const int exactargs, const int nargs) 
{
  if (exactargs >= 0) return (nargs == exactargs);
  // else
  return (nargs >= (-exactargs));
}

int get_num_func_index (const char *sym_string, const int &nargs) 
{
  int i;
  // go until we match both the name and the arg. count
  for (i=0;
       (i<FUNCTION_COUNT) && 
	 ((strcmp(sym_string, NUM_FUNCTIONS[i].name) != 0) || 
	  (! argmatch(NUM_FUNCTIONS[i].Nargs, nargs))); 
       ++i) ;
  if (i!=FUNCTION_COUNT) return i;
  return BAD_FUNC; // if we didn't find it.
}

void verbose_printfuncs (FILE* outfile, const char *pre_str) 
{
  for (int i=0; i<FUNCTION_COUNT; ++i) {
    fprintf(outfile, "%s(%s", pre_str, NUM_FUNCTIONS[i].name);
    int nargs = NUM_FUNCTIONS[i].Nargs;
    for (int j=0; j<abs(nargs); ++j)
      fprintf(outfile, " v%d", j+1);
    if (nargs < 0) fprintf(outfile, " ...)");
    else fprintf(outfile, ")");
    if (is_symm_func(i)) fprintf(outfile, " --symm");
    if (is_assoc_func(i)) fprintf(outfile, " --assoc");
    if (is_func_withstate(i)) fprintf(outfile, " --state");
    
    fprintf(outfile, "\n");
  }
}


const int CONST_COUNT = 3;
const builtin_const_type BUILTIN_CONSTS[CONST_COUNT] = {
  {"pi", M_PI},
  {"e", M_E},
  {"sqrt2", M_SQRT2}};


void verbose_printconsts (FILE* outfile, const char *pre_str) 
{
  for (int i=0; i<CONST_COUNT; ++i) {
    fprintf(outfile, "%s%s = %.17le\n", pre_str, BUILTIN_CONSTS[i].name,
	    BUILTIN_CONSTS[i].value);
  }
}

#endif
