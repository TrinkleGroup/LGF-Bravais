#ifndef __EXPR_OPT_FUNC_H
#define __EXPR_OPT_FUNC_H
/*
  Program: expr-opt-func.H
  Author:  D. Trinkle
  Date:    April 27, 2004
  Purpose: Our different in-place expression optimizations
*/

#include <math.h>
#include "expr-functions.H"
#include "expr-compile.H" // do we need this?

// Just pass the pointer, because we do all modifications *in place*
// We return 0 if we didn't touch anything; 1 else.
typedef int optimization_func_type (exprtree_elem_type* &node,
				    const_table_type& ct);


// some of the functions that we'll simplify
int INDEXED = 0;
int FUNC_INDEX_NEG, FUNC_INDEX_ADD, FUNC_INDEX_SUB, FUNC_INDEX_MUL,
  FUNC_INDEX_DIV, FUNC_INDEX_POW, FUNC_INDEX_INV, FUNC_INDEX_SQRT;

void index_funcs_current() 
{
  int ERROR = 0;
  // get 'em all...
  ERROR=( is_BAD_FUNC(FUNC_INDEX_NEG = get_num_func_index("-", 1)) ||
	  is_BAD_FUNC(FUNC_INDEX_ADD = get_num_func_index("+", 2)) ||
	  is_BAD_FUNC(FUNC_INDEX_SUB = get_num_func_index("-", 2)) ||
	  is_BAD_FUNC(FUNC_INDEX_MUL = get_num_func_index("*", 2)) ||
	  is_BAD_FUNC(FUNC_INDEX_DIV = get_num_func_index("/", 2)) ||
	  is_BAD_FUNC(FUNC_INDEX_POW = get_num_func_index("^", 2)) ||
	  is_BAD_FUNC(FUNC_INDEX_INV = get_num_func_index("inv", 1)) ||
	  is_BAD_FUNC(FUNC_INDEX_SQRT = get_num_func_index("sqrt", 1))
	  );
  if (ERROR)
    fprintf(stderr, "ERROR indexing functions--should NEVER happen\n");
}

// call this function for low overhead
inline void index_funcs() 
{
  if (! INDEXED) {
    INDEXED = 0;
    index_funcs_current();
  }
}


// this is our base optimization for associative functions
// By definition, an associative function has multiple arguments.
int optimize_assoc (exprtree_elem_type* &node, const_table_type& ct) 
{
  // we can only work on associative functions
  if (! is_assoc_func(node)) return 0;
  // else, we want to collapse our entries "in place"
  int count=0, i;
  int nargs = node->Nargs;
  int samefunc[nargs];
  for (i=0; i<nargs; ++i)
    if (is_const_load(node->arglist[i]) || is_symbol_load(node->arglist[i])) {
      ++count; samefunc[i] = 0;
    }
    else if (node->func_index == node->arglist[i]->func_index) {
      count += node->arglist[i]->Nargs;
      samefunc[i] = 1;
    }
    else {++count; samefunc[i] = 0; }
  if (count == nargs) // no collapse to do...
    return 0;
  // else, let's do it:
  exprtree_elem_type** newarglist = new exprtree_elem_type *[count];
  node->Nargs = count;
  count = 0;
  // collapse the arguments, in place:
  for (i=0; i<nargs; ++i) 
    if ( ! samefunc[i] )
      newarglist[count++] = node->arglist[i];
    else 
      for (int j=0; j<node->arglist[i]->Nargs; ++j) 
	newarglist[count++] = node->arglist[i]->arglist[j];
  // replace the arglist
  delete[] node->arglist;
  node->arglist = newarglist;
  return 1;
}

// this is our base optimization for symmetric functions
// We "sort" the arguments so that we can identify identical expressions
// more easily (i.e., (+ x 1) == (+ 1 x), but (/ 1 x) != (/ x 1)
// Because our function is symmetric, the order doesn't matter, so long
// as it's unique--we sort using the pointer addresses
int optimize_symm (exprtree_elem_type* &node, const_table_type& ct)  
{
  if (! is_symm_func(node)) return 0;
  int i, j, changed = 0;
  int nargs = node->Nargs;
  exprtree_elem_type *temp; // for swappin'
  exprtree_elem_type **arglist = node->arglist; // for swappin'
  // we do a little insertion sort
  for (i=1; i<nargs; ++i) {
    temp = arglist[i];
    j = i;
    while (arglist[j-1] > temp) {
      changed = 1;
      arglist[j] = arglist[j-1];
      --j;
      if (j==0) break;
    }
    arglist[j] = temp;
  }
  return changed;
}

// See if we're subtracting a constant; if so, negate the constant and 
// make it an addition operation.
int optimize_sub (exprtree_elem_type* &node, const_table_type& ct) 
{
  if (node->func_index != FUNC_INDEX_SUB) return 0;
  if (! is_const_load(node->arglist[1])) return 0;
  // else, negate the constant and convert to addition:
  exprtree_elem_type* newnode = new exprtree_elem_type;
  newnode->func_index = const_load;
  newnode->arglist = NULL;
  newnode->Nargs = get_value_tree( -(ct.table[node->arglist[1]->Nargs]), ct);
  node->func_index = FUNC_INDEX_ADD;
  node->arglist[1] = newnode;
  return 1;
}


// Look to see if the constant we're adding is 0 (at this point,
// we should have already reduced addition to at most one constant).
int optimize_add (exprtree_elem_type* &node, const_table_type& ct) 
{
  if (node->func_index != FUNC_INDEX_ADD) return 0;
  int i, nargs = node->Nargs;
  for (i=0; (i<nargs) && (! is_const_load(node->arglist[i])); ++i) ;
  if (i == nargs) return 0;
  // else, check out the value of that constant:
  if (zero(ct.table[node->arglist[i]->Nargs])) {
    // we can safely remove this single item; two ways to go:
    // a. there were but two arguments; replace with the non-constant one:
    if (nargs == 2) {
      node = node->arglist[(1-i)];
      return 1;
    }
    // else, we just need to snip out this one argument:
    nargs = (--(node->Nargs));
    for (int j=i; j<nargs; ++j) node->arglist[j] = node->arglist[j+1];
    return 1;
  }
  return 0;
}


// NOTE: This code currently *won't* combine like powers for small integer
// cases, because it expands them out first--this has to do with the
// order in which we optimize by depth.

// First, replace (x^a)^b with x^(ab).  Then see if we're taking a small 
// integer power; replace with multiplications or divisions if need be.
//  Also replace 1/2 and -1/2 with sqrt() and sqrt(inv()), respectively.
int optimize_pow (exprtree_elem_type* &node, const_table_type& ct) 
{
  if (node->func_index != FUNC_INDEX_POW) return 0;
  if (! is_const_load(node->arglist[1])) return 0;
  double value = ct.table[node->arglist[1]->Nargs];
  
  // first effort--can we combine like powers?
  if (node->arglist[0]->func_index == FUNC_INDEX_POW) {
    // then, (^ (^ x a) b) == (^ x (* a b))
    exprtree_elem_type* newnode = new exprtree_elem_type;
    newnode->func_index = FUNC_INDEX_MUL;
    newnode->Nargs = 2;
    newnode->arglist = new exprtree_elem_type*[2];
    newnode->arglist[0] = node->arglist[1];
    newnode->arglist[1] = node->arglist[0]->arglist[1];
    // now, fix up node:
    node->arglist[0] = node->arglist[0]->arglist[0];
    node->arglist[1] = newnode;
    return 1;
  }

  // check against our basic roots:
  if (dcomp(fabs(value), 0.5)) {
    exprtree_elem_type* newnode = new exprtree_elem_type;
    newnode->func_index = FUNC_INDEX_SQRT;
    newnode->Nargs = 1;
    newnode->arglist = new exprtree_elem_type*[1];
    if (value > 0)
      newnode->arglist[0] = node->arglist[0];
    else {
      exprtree_elem_type *newnewnode = new exprtree_elem_type;
      newnewnode->func_index = FUNC_INDEX_INV;
      newnewnode->Nargs = 1;
      newnewnode->arglist = new exprtree_elem_type*[1];
      newnewnode->arglist[0] = node->arglist[0];
      newnode->arglist[0] = newnewnode;
    }
    node = newnode;
    return 1;
  }

  int intpow; // see if we have an integer power...
  const int MAXPOW = 6;
  for (intpow = -MAXPOW; (intpow<=MAXPOW) && (! dcomp(value, (double)intpow));
       ++intpow) ;
  if (intpow > MAXPOW) return 0;
  switch (intpow) {
  case 0: // x^0 == 1
    node = new exprtree_elem_type;
    node->func_index = const_load;
    node->arglist = NULL;
    node->Nargs = get_value_tree(1., ct);
    return 1;
  case 1: // x^1 == x
    node = node->arglist[0];
    return 1;
  case -1: // x^-1 == 1/x
    {
      exprtree_elem_type* newnode = new exprtree_elem_type;
      newnode->func_index = FUNC_INDEX_INV;
      newnode->Nargs = 1;
      newnode->arglist = new exprtree_elem_type*[1];
      newnode->arglist[0] = node;
      node = newnode;
    }
    return 1;
  default: 
    // construct the multiple multiplications first:
    // x^n = x*x*...*x
    {
      exprtree_elem_type* newnode = new exprtree_elem_type;
      newnode->func_index = FUNC_INDEX_MUL;
      newnode->Nargs = abs(intpow);
      newnode->arglist = new exprtree_elem_type*[abs(intpow)];
      for (int i=0; i<abs(intpow); ++i)
	newnode->arglist[i] = node->arglist[0];
      node = newnode;
      if (intpow < 0) {
	// x^-n = 1./(x*x*...*x)
	newnode = new exprtree_elem_type;
	newnode->func_index = FUNC_INDEX_INV;
	newnode->Nargs = 1;
	newnode->arglist = new exprtree_elem_type*[1];
	newnode->arglist[0] = node;
	node = newnode;
      }
    }
    return 1;
  }
}



// See if we're dividing by a constant; if so, invert the constant and 
// make it an addition operation.
int optimize_div (exprtree_elem_type* &node, const_table_type& ct) 
{
  if (node->func_index != FUNC_INDEX_DIV) return 0;
  if (! is_const_load(node->arglist[1])) return 0;
  // else, negate the constant and convert to addition:
  exprtree_elem_type* newnode = new exprtree_elem_type;
  newnode->func_index = const_load;
  newnode->arglist = NULL;
  newnode->Nargs = get_value_tree( 1./(ct.table[node->arglist[1]->Nargs]), ct);
  node->func_index = FUNC_INDEX_MUL;
  node->arglist[1] = newnode;
  return 1;
}


// Look to see if the constant we're multiplying is 0, 1 (at 
// this point, we should have already reduced multiplication to at 
// most one constant).
// We might at some point want to put multiplication of powers in here,
// and add up the powers of like terms
int optimize_mul (exprtree_elem_type* &node, const_table_type& ct) 
{
  if (node->func_index != FUNC_INDEX_MUL) return 0;
  int i, nargs = node->Nargs;
  for (i=0; (i<nargs) && (! is_const_load(node->arglist[i])); ++i) ;
  if (i == nargs) return 0;
  // else, check out the value of that constant:
  double value = ct.table[node->arglist[i]->Nargs];
  if (zero(value)) {
    // then the whole shebang is zero:
    node = new exprtree_elem_type;
    node->func_index = const_load;
    node->arglist = NULL;
    node->Nargs = get_value_tree(0., ct);
    return 1;
  }
  if (dcomp(value, 1.)) {
    // then we remove the multiplication by unity
    // we can safely remove this single item; two ways to go:
    // a. there were but two arguments; replace with the non-constant one:
    if (nargs == 2) {
      node = node->arglist[(1-i)];
      return 1;
    }
    // else, we just need to snip out this one argument:
    nargs = (--(node->Nargs));
    for (int j=i; j<nargs; ++j) node->arglist[j] = node->arglist[j+1];
    return 1;
  }
  
  return 0;
}


const int NUM_FUNC_OPTS = 7; // how many different optimizations do we have?

// Note: the most general ones should go first, and the least general
// need to go last
optimization_func_type* num_func_opt[NUM_FUNC_OPTS] = {
  optimize_assoc, optimize_symm, optimize_sub, optimize_add,
  optimize_pow, optimize_div, optimize_mul};


#endif
