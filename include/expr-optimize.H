#ifndef __EXPR_OPTIMIZE_H
#define __EXPR_OPTIMIZE_H
/*
  Program: expr-optimize.H
  Author:  D. Trinkle
  Date:    April 27, 2004
  Purpose: To provide expression tree optimization.
*/

#include <ctype.h>
#include <string.h> // needed for parsing...
#include "dcomp.H"
#include "expr-compile.H"   // contains the tree definitions, parser, compiler
#include "expr-opt-func.H"  // our special optimizations

// ****************************** STRUCTURES ***************************
// First, a table of all the subexpressions we've encountered so far;
// useful for combining subexpressions and garbage collection..
typedef struct 
{
  int Nelem;
  int maxelem;
  exprtree_elem_type** table;
} exprtree_table_type;

// ****************************** CONSTANTS ****************************
const int MAXEXPR = 16384; // pretty high limit...

// ***************************** SUBROUTINES ***************************

// Code to optimize a given expression tree (optimized under for evaluation
// by a register machine).  We do this by folding down constants when
// we can, and by computing only a single instance of each subexpression.
// Some of these optimizations won't be seen until the tree is compiled;
// ie, they won't improve evaluation using the tree with eval_exprtree.

exprtree_table_type* init_exprtree_table (int Nelem) 
{
  exprtree_table_type* etp = new exprtree_table_type;
  etp->Nelem = 0;
  etp->maxelem = Nelem;
  etp->table = new exprtree_elem_type*[Nelem];
  for (int i=0; i<Nelem; ++i) etp->table[i] = NULL;
  return etp;
}

void free_exprtree_table (exprtree_table_type* &etp) 
{
  if (etp != NULL) {
    delete[] etp->table;
    delete etp;
    etp = NULL;
  }
}

int equiv_exprtree(exprtree_elem_type* node1, exprtree_elem_type* node2) 
{
  if (node1 == node2) return 1; // same pointer
  if (node1 == NULL) return 0;  // not the same, and one is NULL
  if (node2 == NULL) return 0;
  if (node1->func_index != node2->func_index) return 0;
  if (node1->Nargs != node2->Nargs) return 0;
  // if we're loading a symbol/const, you only need to check Nargs
  if (is_symbol_load(node1) || is_const_load(node1)) return 1;
  // normally, we'd check subarguments at this point; however, if
  // we have a function that has state, we can *only* match if
  // the pointers are *equal*
  if (is_func_withstate(node1)) return 0;
  // else, look at the subarguments.  We don't *recurse* here, because
  // we assume they've been simplified already and made to equal
  // pointers if they're equivalent.
  for (int i=0; i<node1->Nargs; ++i)
    if (node1->arglist[i] != node2->arglist[i]) return 0;
  // if we got this far, they have identical arguments, so they're
  // equivalent:
  return 1;
}

// Returns a pointer to an equivalent expression in our list, if it
// exists.  If it hasn't seen this *pointer* before, it adds it to the
// end of the list (to help with garbage collection efforts)
exprtree_elem_type* get_exprtree_tree (exprtree_elem_type* node,
				       exprtree_table_type &et) 
{
  int i;
  // first, see if the pointer is there as is:
  for (i=0; (i<et.Nelem) && (node != et.table[i]); ++i) ;
  if (i!=et.Nelem) return node;
  // else, we have to add it:
  if (et.Nelem == et.maxelem) return NULL; // sanity check
  et.table[i] = node;
  ++(et.Nelem);
  // Now, we need to sweep through again, in case it's *equivalent*
  // (we wanted to add it in case we have to garbage collect it later;
  // in this way, we keep track of every node we've looked at).  Plus,
  // we're guaranteed to match at some point.
  exprtree_elem_type** testnode;
  for (testnode=et.table; (! equiv_exprtree(node, *testnode)); ++testnode) ;
  return *testnode; // at worse, we'll end up returning the same pointer.
}

void mark_used_expr(exprtree_elem_type* tree, exprtree_table_type &et,
		    int* usedlist);
void mark_used_expr(exprtree_elem_type* tree, exprtree_table_type &et,
		    int* usedlist)
{
  int i;
  for (i=0; (i<et.Nelem) && (et.table[i] != tree); ++i) ;
  usedlist[i] = 1;
  if (is_const_load(tree) || is_symbol_load(tree)) return;
  // else, we have to check all the arguments:
  for (i=0; i<tree->Nargs; ++i)
    mark_used_expr(tree->arglist[i], et, usedlist);
}

void garbage_collect(exprtree_elem_type* tree, exprtree_table_type &et) 
{
  int* usedlist = new int[et.Nelem];
  for (int i=0; i<et.Nelem; ++i) usedlist[i] = 0;
  mark_used_expr(tree, et, usedlist);
  for (int i=0; i<et.Nelem; ++i)
    if (! usedlist[i]) {
      delete et.table[i];
      et.table[i]=NULL;
    }
  delete[] usedlist;
}


void mark_used_const(exprtree_elem_type* tree, int* usedlist);
void mark_used_const(exprtree_elem_type* tree, int* usedlist)
{
  if (is_const_load(tree)) 
    usedlist[tree->Nargs] = 1;
  if (is_const_load(tree) || is_symbol_load(tree)) return;
  // else, we have to check all the arguments:
  for (int i=0; i<tree->Nargs; ++i)
    mark_used_const(tree->arglist[i], usedlist);
}


void const_collect(exprtree_elem_type* tree, const_table_type &ct,
		   exprtree_table_type &et)
{
  int usedlist[ct.Nelem];
  for (int i=0; i<ct.Nelem; ++i) usedlist[i] = 0;
  mark_used_const(tree, usedlist);
  // now, construct our map
  int map[ct.Nelem]; // map[i] = new constant index
  int current=0;
  for (int i=0; i<ct.Nelem; ++i)
    if (usedlist[i]) {
      // move down the constant
      ct.table[current] = ct.table[i];
      map[i] = current; // store the mapping
      ++current;
    }
  ct.Nelem = current; // our count of how many constants
  // we apply the constant map using the *exprtree* so as to not
  // apply the map multiple times (this would happen if we traversed
  // our tree, and a subexpression, or constant, was referenced more
  // than once).
  for (int i=0; i<et.Nelem; ++i)
    if (et.table[i] != NULL) {
      if (is_const_load(et.table[i]->func_index))
	et.table[i]->Nargs = map[et.table[i]->Nargs];
    }
}


void mark_used_symbols(exprtree_elem_type* tree, int* usedlist);
void mark_used_symbols(exprtree_elem_type* tree, int* usedlist)
{
  if (is_symbol_load(tree)) 
    usedlist[tree->Nargs] = 1;
  if (is_const_load(tree) || is_symbol_load(tree)) return;
  // else, we have to check all the arguments:
  for (int i=0; i<tree->Nargs; ++i)
    mark_used_symbols(tree->arglist[i], usedlist);
}

// same mentality as const collect: we toss out unused symbols.  This
// can be somewhat disconcerting--for example, (* x 0) will produce
// code that has *no* arguments, since as a symbol x will be
// obliterated when (* x 0) is simplified to 0.
void symbol_collect(exprtree_elem_type* tree, symbol_table_type &st,
		    exprtree_table_type &et)
{
  int usedlist[st.Nelem];
  for (int i=0; i<st.Nelem; ++i) usedlist[i] = 0;
  mark_used_symbols(tree, usedlist);
  // now, construct our map
  int map[st.Nelem]; // map[i] = new symbol index
  int current=0;
  for (int i=0; i<st.Nelem; ++i)
    if (usedlist[i]) {
      // move down the constant
      st.table[current] = st.table[i];
      map[i] = current; // store the mapping
      ++current;
    }
  st.Nelem = current; // our count of how many constants
  // we apply the constant map using the *exprtree* so as to not
  // apply the map multiple times (this would happen if we traversed
  // our tree, and a subexpression, or constant, was referenced more
  // than once).
  for (int i=0; i<et.Nelem; ++i)
    if (et.table[i] != NULL) {
      if (is_symbol_load(et.table[i]->func_index))
	et.table[i]->Nargs = map[et.table[i]->Nargs];
    }
}




// the recursive workhorse...
int optimize_exprtree(exprtree_elem_type* &tree, const_table_type &ct,
		      symbol_table_type &st, exprtree_table_type &et);
int optimize_exprtree(exprtree_elem_type* &tree, const_table_type &ct,
		      symbol_table_type &st, exprtree_table_type &et) 
{
  int ERROR = 0;
  int i;
  if (is_const_load(tree) || is_symbol_load(tree)) {
    // there's only one simplification--reduction of symbol/const
    // loads.  This is done so that we recognize that 1==1, etc. 
    // and also so that we manage to "mark" every subtree element,
    // for garbage collection purposes.
    tree = get_exprtree_tree(tree, et);
    return ERROR;
  }
  // now, how to handle a function call:
  int nconst = 0;
  int narg=tree->Nargs;
  int constarg[narg];

  // 1. we have to optimize all of the argument list first:
  for (i=0; (i<narg) && (!ERROR); ++i)
    ERROR = optimize_exprtree(tree->arglist[i], ct, st, et);
  if (ERROR) return ERROR;

  if (! is_func_withstate(tree)) {
    // 2. pool together our number of constant args:
    double* arglist[narg], result; // in case we need to eval. something...
    for (i=0; i<narg; ++i)
      if (is_const_load(tree->arglist[i])) {++nconst; constarg[i] = 1;}
      else constarg[i] = 0;
    if (nconst == narg) {
      // We can just evaluate this piece right now, and even do it
      // in place... because garbage collection will take care of the rest.
      // copy our arguments in...
      for (i=0; i<narg; ++i) {
	arglist[i] = new double;
	*(arglist[i]) = ct.table[tree->arglist[i]->Nargs];
      }
      // evaluate:
      result = 
	(NUM_FUNCTIONS[tree->func_index].num_func) (narg, arglist);
      
      // we change tree to point to our new constant value:
      tree->func_index = const_load;
      tree->Nargs = get_value_tree(result, ct);
      delete[] tree->arglist; // *just* destroy the pointer list...
      tree->arglist = NULL;
      
      // garbage collect arglist:
      for (i=0; i<narg; ++i) delete arglist[i];
    }
    else 
      // NOTE: associative -> variable number of arguments
      if ( (nconst > 1) && (is_assoc_func(tree)) && (is_symm_func(tree)) ) {
	// we've got variable number of args, so we can replace those
	// constants with the partial sum/product:
	int count = 0;
	for (i=0; i<narg; ++i) 
	  if (constarg[i]) {
	    arglist[count] = new double;
	    *(arglist[count]) = ct.table[tree->arglist[i]->Nargs];
	    ++count;
	  }
	// evaluate:
	result = 
	  (NUM_FUNCTIONS[tree->func_index].num_func) (count, arglist);
	// roll up constants:
	count = 0;
	for (i=0; i<narg; ++i)
	  if (! constarg[i]) {
	    tree->arglist[count] = tree->arglist[i];
	    ++count;
	  }
	// create a new constant node:
	exprtree_elem_type *constnode = new exprtree_elem_type;
	constnode->func_index = const_load;
	constnode->Nargs = get_value_tree(result, ct);
	constnode->arglist = NULL;
	// check it in with the expression table:
	constnode = get_exprtree_tree(constnode, et);
	// add it to the arglist:
	tree->arglist[count] = constnode;
	tree->Nargs = count+1;
      }
  }
  // finally, run through possible "in place" optimizations:
  int modified = 0;
  // save at this stage, before we do more modification:
  tree = get_exprtree_tree(tree, et);
  for (i=0; i<NUM_FUNC_OPTS; ++i)
    if ( (num_func_opt[i]) (tree, ct) ) {
      modified = 1; // but keep going...
      // try to save--not always necessary; just in case of garbage collection
      tree = get_exprtree_tree(tree, et);
      // also make sure we index the arguments, in case any of them
      // were changed; but if we've now become just a symbol or
      // constant, break out now.
      if (is_symbol_load(tree) || is_const_load(tree)) break;

      for (int j=0; j<narg; ++j)
	tree->arglist[j] = get_exprtree_tree(tree->arglist[j], et);
    }
  
  // always see if we have something equivalent to a know expression...
  if (modified) // reoptimize it (since we modified it)
    ERROR = optimize_exprtree(tree, ct, st, et);

  return ERROR;
}


// NOTE: we have the code setup to remove unused symbols after optimization.
// This can have the unexpected consequence of *reducing* the number of
// arguments for a given program.  So, for example, (* x 0) initially
// has a single argument.  After optimization, the exprtree is simply
// 0, and now has *no* arguments.
int optimize_exprtree(exprtree_elem_type* &tree, const_table_type &ct,
		      symbol_table_type &st)
{
  int ERROR = 0;
  // create some workspace
  exprtree_table_type* et = init_exprtree_table(MAXEXPR);

  index_funcs(); // allows us to do function optimizations.  
  ERROR = optimize_exprtree(tree, ct, st, *et);
  if (!ERROR) {
    garbage_collect(tree, *et);    // remove unused subexpressions
    const_collect(tree, ct, *et);  // remove unused constants
    symbol_collect(tree, st, *et); // remove unused symbols
  }
  
  free_exprtree_table(et);
  return ERROR;
}


#endif
