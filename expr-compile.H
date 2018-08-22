#ifndef __EXPR_COMPILE_H
#define __EXPR_COMPILE_H
/*
  Program: expr-compile.H
  Author:  D. Trinkle
  Date:    April 27, 2004
  Purpose: To provide expression evaluation in realtime.  We input a
           function expression in lisp-style, and:

	   1. Parse it into an expression tree, symbol table and constant
	      table.
	   2. The exp. tree is optimized, folding down as many constants
	      and sub-expressions as possible 
	   3. The exp. tree is compiled in a register machine implementation.
	   4. The compiled function is optimized for memory usage.

	   Finally, the compiled function can be evaluated for a given
	   symbol list in real time.

  To-do:   It might be worthwhile to add two (or three?) things to this 
           code:

	   1. Macro substitution.  The simplest would be the addition of
	      a "subs" command: (subs u (+ x x) (* u u)) would be parsed,
	      and expanded out to (* (+ x x) (+ x x)).  Other possible
	      implemetations exist, such as parameterized macros:
	      (subs (u y) (+ y y) (* (u x) (u 2))) ->
	      (* (+ x x) (+ 2 2))
	      and independently defined parameterized macros:
	      (defmacro (u y) (+ y y))
	      (* (u x) (u 2)) -> (* (+ x x) (+ 2 2)).
	      This gets very close to the second possibility:
	   2. User-defined functions.  This is probably best implemented
	      via 1.
	   3. Differentiation.  Another take on macro substitution, this
	      becomes a game in subexpression dependendent macro subs.
	      That is, (diff x x) -> 1, (diff x (+ x x)) -> (+ 1 1),
	      (diff x (+ x (^ x 2))) -> (+ 1 (* 2 (^ x 1))).  To really
	      make this efficient, it might be worth revisiting and
	      improving upon some of the expression optimization we do.

	   And there's always room for more and more clever optimization
	   routines to be written (expr-opt-func.H).
*/

#include <ctype.h>
#include <string.h> // needed for parsing...
#include "dcomp.H"
#include "expr-functions.H" // contains all of the functions, constants, etc.
#include "expr-programs.H"  // assembler, as well as constant and symbol tables

// ****************************** STRUCTURES ***************************
// Now defined in expr-functions:
// typedef double numerical_function_type (int Nargs, double* arglist);

typedef struct exprtree_elem_type;
struct exprtree_elem_type 
{
  //  num_func_type *func;
  int func_index;
  int Nargs;
  exprtree_elem_type **arglist; // list of pointers
};

// We need this table when we're compiling expression trees that might
// have crossig branches (i.e., "stored" results from previous evaluations)
typedef struct 
{
  int Nelem;
  int maxelem;
  exprtree_elem_type **table;
  int *reg;
} reg_table_type;


// ****************************** CONSTANTS ****************************
// "Hard" limits on number of unique numerical constants and symbols
// when parsing; may be overridden by allocating own tables before
// parse_expr
const int MAXCONSTS = 1024;
const int MAXSYMBOLS = 128;

const int PARSE_ERROR_EARLYTERM =      1; // \0 appears before we're done
const int PARSE_ERROR_TOOMANYCONST =   2; // constant table too small
const int PARSE_ERROR_TOOMANYSYMBOLS = 4; // symbol table too small
const int PARSE_ERROR_BADSYMBOL =      8; // badly constructed symbol (eg: 1a)
const int PARSE_ERROR_BADFUNC =       16; // bad function name

const int MAXINST = 2048;  // maximum length of any program
  

// ***************************** SUBROUTINES ***************************

inline int is_symbol_load(const exprtree_elem_type* node) 
{ return is_symbol_load(node->func_index); }

inline int is_const_load(const exprtree_elem_type* node) 
{ return is_const_load(node->func_index); }

inline int is_BAD_FUNC(const exprtree_elem_type* node) 
{ if (node == NULL) return 1;
  return is_BAD_FUNC(node->func_index); }

inline int is_symm_func(const exprtree_elem_type* node) 
{ return is_symm_func(node->func_index); }

inline int is_assoc_func(const exprtree_elem_type* node) 
{ return is_assoc_func(node->func_index); }

inline int is_func_withstate(const exprtree_elem_type* node) 
{ return is_func_withstate(node->func_index); }



// string manipulation: move s up to the first non-blank character,
// but not past the end.
void next_nonblank(char* &s) 
{ for ( ; isblank(*s) && (*s != 0); ++s) ; }

// picks out the substring consisting of non-separator characters,
// and increments s just past.
const char* SEPARATOR_CHAR = " \t()\0";

char* get_symbol_string (char* &s) 
{
  int i;
  char* sym_string;
  i = strcspn(s, SEPARATOR_CHAR);
  sym_string = new char[i+1];
  strncpy(sym_string, s, i); // copy to sym_string
  sym_string[i] = 0;         // add null termination
  s += i;                    // increment s past it.
  // convert symbol string to all lowercase:
  for (int j=0; j<i; ++j) sym_string[j] = tolower(sym_string[j]);
  return sym_string;
}

// determine if we have a constant, and if so, extract the value.
int is_constant(char* sym_string, double &value) 
{
  char* endptr;
  // naively try to convert it...
  value = strtod(sym_string, &endptr);
  // see if we got to the end...
  if (*endptr == 0) return 1;
  // see if we got partway through (indicated a bad symbol...)
  if (endptr != sym_string) return 0;
  // now, final check is... see if it's a built-in symbol:
  for(int i=0; i<CONST_COUNT; ++i)
    if (strcmp(sym_string, BUILTIN_CONSTS[i].name) == 0) {
      value = BUILTIN_CONSTS[i].value;
      return 1;
    }
  // if we got this far, it *should* be a symbol:
  return 0;
}

// VERY simple test (assumes we failed on being a builtin constant):
int is_symbol(char* sym_string) 
{
  return isalpha(*sym_string);
}

// all this does is make our tree list, so we can safely free it:
void free_exprtree_recurse (exprtree_elem_type* &tree, int &Nelem,
			    exprtree_elem_type** treelist);
void free_exprtree_recurse (exprtree_elem_type* &tree, int &Nelem,
			    exprtree_elem_type** treelist)
{
  int i;
  for (i=0; (i<Nelem) && (tree != treelist[i]); ++i) ;
  if (i != Nelem) return;
  // else, we have to mark this pointer, and all the children:
  treelist[Nelem] = tree;
  ++Nelem;
  if ( (tree->func_index >= 0) && (tree->arglist != NULL) ) 
    // happens if parsing error, or load command
    for (i=0; i<tree->Nargs; ++i) 
      free_exprtree_recurse(tree->arglist[i], Nelem, treelist);
}

void free_exprtree(exprtree_elem_type* &tree) 
{
  int Nelem = 0;
  exprtree_elem_type** treelist = new exprtree_elem_type*[MAXINST];
  
  if (tree != NULL) {
    free_exprtree_recurse(tree, Nelem, treelist); // list of unique pointers
    for (int i=0; i<Nelem; ++i) {
      delete[] treelist[i]->arglist;
      delete[] treelist[i];
      treelist[i] = NULL;
    }
    delete[] treelist;
    tree = NULL;
  }
}


// Recursive method for parsing a string composed of an expression into
// an expression tree, along with a constant and symbol table.
// Very exciting.

// Our recursive workhorse for going through the string
exprtree_elem_type* parse_expr_recurse (char* &s, const_table_type &ct,
					symbol_table_type &st, int &ERROR);

int parse_expr(char* s, const_table_type* &ctp, symbol_table_type* &stp,
	       exprtree_elem_type* &exprtree) 
{
  int ERROR = 0;
  if (ctp == NULL) ctp = init_const_table(MAXCONSTS);
  if (stp == NULL) stp = init_symbol_table(MAXSYMBOLS);
  exprtree = parse_expr_recurse(s, *ctp, *stp, ERROR);
  return ERROR;
}


exprtree_elem_type* parse_expr_recurse (char* &s, const_table_type &ct,
					symbol_table_type &st, int &ERROR) 
{
  int i;
  char *sym_string;
  exprtree_elem_type *node;  

  // We always assume ERROR = 0 at this point
  next_nonblank(s); // move our string up.
  if (*s == 0) {
    ERROR = PARSE_ERROR_EARLYTERM;
    return NULL;
  }
  if (*s != '(') {
    // we've got a constant or a symbol.
    sym_string = get_symbol_string(s); // also moves s up.
    double value;
    if (is_constant(sym_string, value)) {
      node = new exprtree_elem_type;
      node->func_index = const_load;
      node->Nargs = get_value_tree(value, ct); // may add element
      if (node->Nargs < 0) ERROR = PARSE_ERROR_TOOMANYCONST;
      node->arglist = NULL;
      delete[] sym_string;
      return node;
    }
    if (is_symbol(sym_string)) {
      node = new exprtree_elem_type;
      node->func_index = symbol_load;
      node->Nargs = get_symbol_tree(sym_string, st); // may add element
      if (node->Nargs < 0) ERROR = PARSE_ERROR_TOOMANYSYMBOLS;
      node->arglist = NULL;
      delete[] sym_string;
      return node;
    }
    // if we got to this point, it's not a symbol or a constant...
    ERROR = PARSE_ERROR_BADSYMBOL;
    delete[] sym_string;
    return NULL;
  }
  // Else, we're dealing with a function, so get what kind.
  ++s; // increment past the "("
  next_nonblank(s);
  if (*s == 0) {
    ERROR = PARSE_ERROR_EARLYTERM;
    return NULL;
  }
  sym_string = get_symbol_string(s); // also moves s up.
  // Now we can start making our new exprtree_elem.
  node = new exprtree_elem_type;
  // we don't know the exact function till we count the args

  int nargs=0, maxargs=1;
  // assume we have at least one argument, expand by powers of 2
  node->arglist = new exprtree_elem_type *[maxargs];
  next_nonblank(s);
  while ( (*s != ')') && (*s != 0) && (!ERROR)) {
    if (nargs == maxargs) {
      // expand list by a factor of 2 if need be
      exprtree_elem_type **temp = new exprtree_elem_type *[2*maxargs];
      for (i=0; i<nargs; ++i) temp[i] = node->arglist[i];
      delete[] node->arglist; // just delete the list of *pointers*
      node->arglist = temp;
      maxargs *= 2;
    }
    // and... recursion:
    node->arglist[nargs] = parse_expr_recurse(s, ct, st, ERROR);
    ++nargs;
    next_nonblank(s);
  }
  // if we've got an error already, don't bother with the rest of this...
  if (!ERROR) {
    if (*s == ')') ++s; // advance past the trailing paren.
    node->Nargs = nargs;
    node->func_index = get_num_func_index(sym_string, nargs);
    if (node->func_index < 0) {
      ERROR = PARSE_ERROR_BADFUNC;
    }  
  }
  else node->func_index = BAD_FUNC;
  
  delete[] sym_string;
  return node;
}


// Recursively traverse the tree, and evaluate:
// Note: we *technically* don't need the symbol table at this point...
double eval_exprtree(exprtree_elem_type *tree, const_table_type &ct,
		     symbol_table_type &st, double* arglist);
double eval_exprtree(exprtree_elem_type *tree, const_table_type &ct,
		     symbol_table_type &st, double* arglist)
{
  int i;
  if (is_BAD_FUNC(tree->func_index)) return 0; // should NEVER do this.
  if (is_const_load(tree->func_index))
    return ct.table[tree->Nargs];
  if (is_symbol_load(tree->func_index))
    return arglist[tree->Nargs];
  // else, we've got some function to work with.
  // 1. evaluate all of the arguments
  double* args[tree->Nargs];
  for (i=0; i<tree->Nargs; ++i) {
    args[i] = new double;
    *(args[i]) = eval_exprtree(tree->arglist[i], ct, st, arglist);
  }
  // 2. evaluate, and return, the value from our function
  double value = 
    NUM_FUNCTIONS[tree->func_index].num_func(tree->Nargs, args);
  for (i=0; i<tree->Nargs; ++i) delete args[i];
  return value;
}



// Prints out everything in traversal:
void traverse_exprtree(exprtree_elem_type *tree, const_table_type &ct,
		       symbol_table_type &st, int depth, const char* prestr);
void traverse_exprtree(exprtree_elem_type *tree, const_table_type &ct,
		       symbol_table_type &st, int depth, const char* prestr)
{
  int i;
  if (prestr != NULL) printf("%s", prestr);
  for (i=0; i<depth; ++i) printf("  ");
  if (is_const_load(tree->func_index)) {
    printf("c%d %.8le\n", tree->Nargs, ct.table[tree->Nargs]);
    return;
  }
  if (is_symbol_load(tree->func_index)) {
    printf("s%d %s\n", tree->Nargs, st.table[tree->Nargs]);
    return;
  }
  // else, it's a function:
  printf("%s", NUM_FUNCTIONS[tree->func_index].name); // print our function
  printf(" (%d arguments):\n", tree->Nargs);
  for (i=0; i<tree->Nargs; ++i)
    traverse_exprtree(tree->arglist[i], ct, st, depth+1, prestr);
}

// initial depth = 0
inline void traverse_exprtree(exprtree_elem_type *tree, const_table_type &ct,
			      symbol_table_type& st, const char* prestr) 
{ traverse_exprtree(tree, ct, st, 0, prestr);}

inline void traverse_exprtree(exprtree_elem_type *tree, const_table_type &ct,
			      symbol_table_type& st) 
{ traverse_exprtree(tree, ct, st, 0, NULL);}


// ******************************* COMPILER ****************************

reg_table_type* init_reg_table (int Nreg) 
{
  reg_table_type* rtp = new reg_table_type;
  rtp->Nelem = 0;
  rtp->maxelem = Nreg;
  rtp->table = new exprtree_elem_type *[Nreg];
  for (int i=0; i<Nreg; ++i) rtp->table[i] = NULL;
  rtp->reg = new int[Nreg];
  return rtp;
}

void free_reg_table (reg_table_type* &rtp) 
{
  if (rtp != NULL) {
    delete[] rtp->table;
    delete[] rtp->reg;
    delete rtp;
    rtp = NULL;
  }
}

void add_reg_table (reg_table_type &rt, exprtree_elem_type *node, int reg) 
{
  if (rt.Nelem < rt.maxelem) {
    rt.table[rt.Nelem] = node;
    rt.reg[rt.Nelem] = reg;
    ++(rt.Nelem);
  }
}


// the return value from a call to compile is the register pointing
// to the value of it's part of the subtree
int compile_asm_recurse(exprtree_elem_type* tree, program_asm_type* &asmp,
			reg_table_type& rt);
int compile_asm_recurse(exprtree_elem_type* tree, program_asm_type* &asmp,
			reg_table_type& rt)
{
  int result; // where are we going to put our result?
  // handle constant and symbol loads...
  if (is_symbol_load(tree))
    return tree->Nargs;
  if (is_const_load(tree))
    return -(tree->Nargs+1);
  // we've got a function; check to see that we haven't already compiled
  // this piece of the tree:
  int k;
  for (k=0; (k<rt.Nelem) && (tree != rt.table[k]); ++k) ;
  if (k != rt.Nelem) return rt.reg[k]; // register where compiled value lives

  // else, we've got something new to compile.
  int arglist[tree->Nargs]; // list of registers
  for (int i=0; i<tree->Nargs; ++i)
    arglist[i] = compile_asm_recurse(tree->arglist[i], asmp, rt);
  result = asmp->Nreg; // our result goes in the first available register

  // let create our instruction:
  int inst = asmp->Ninst;
  asmp->inst_index[inst] = tree->func_index;
  asmp->Nargs[inst] = tree->Nargs;
  asmp->reslist[inst] = result;
  asmp->arglist[inst] = new int[tree->Nargs];
  for (int j=0; j<tree->Nargs; ++j)
    asmp->arglist[inst][j] = arglist[j];

  // increment counters...
  ++(asmp->Nreg);
  ++(asmp->Ninst);
  // store in our register table:
  add_reg_table(rt, tree, result);
  
  return result; // return the register holding our result
}


// We go through backwards, counting register usage; we then
// construct a new register map and change all of the registers to
// minimize the amount used
void optimize_asm_registers(program_asm_type* &asmp) 
{
  int i, j, r;
  /*
  // minimal:
  asmp->reslist[asmp->Ninst - 1] = 0; // where we will ALWAYS put the result.
  --(asmp->Nreg);
  */

  // the rest of what we do is just gravy on top of that.
  int ninst = asmp->Ninst;
  int nreg = asmp->Nreg;
  int removable[ninst][nreg];
  int accessed[nreg];

  for (r=0; r<nreg; ++r) accessed[r] = 0;
  // now, go through our instructions and find out what registers are
  // removable and which have been accessed at each stage:
  for (i=(ninst-1); i>=0; --i) {
    // we work *backwards*
    // see if our arguments are removable:
    // assume it's removable, unless it will be accessed in the future
    for (r=0; r<nreg; ++r) 
      if (accessed[r]) removable[i][r] = 0;
      else             removable[i][r] = -1;
    // add our list of args to the accessed list
    for (j=0; j<asmp->Nargs[i]; ++j) {
      r = asmp->arglist[i][j];
      if (r>=0)
	accessed[r] = 1;
    }
  }
  // Now, we use this to construct our register map.
  int map[nreg];
  int needed[nreg];
  for (r=0; r<nreg; ++r) map[r] = r;

  //  for (r=0; r<nreg; ++r) needed[r] = 0;
  for (r=0; r<asmp->Narg; ++r) needed[r] = 1;
  for (; r<nreg; ++r) needed[r] = 0;

  for (i=0; i<ninst; ++i) {
    for (j=0; j<asmp->Nargs[i]; ++j) {
      r = asmp->arglist[i][j];
      if (r>=0) {
	asmp->arglist[i][j] = map[r]; // map to new register
	if (! removable[i][r])
	  needed[map[r]] = 1; // mark this register as needed later
	else
	  needed[map[r]] = 0;
      }
    }
    // result remapping--find the lowest unneeded register:
    int resreg;
    for (resreg=0; needed[resreg]; ++resreg) ;
    map[asmp->reslist[i]] = resreg;
    asmp->reslist[i] = resreg;
    needed[resreg] = 1;
  }

  // Count how many registers we *really* use now:
  nreg = -1;
  for (i=0; i<ninst; ++i) {
    r = asmp->reslist[i];
    if (r>nreg) nreg = r;
    for (j=0; j<asmp->Nargs[i]; ++j) {
      r = asmp->arglist[i][j];
      if (r>nreg) nreg = r;
    }
  }
  asmp->Nreg = nreg+1;
}


program_asm_type* compile_asm(exprtree_elem_type* tree, symbol_table_type& st) 
{
  program_asm_type* asmp = init_program_asm(MAXINST);
  asmp->Nreg = st.Nelem; // we need enough space at *least* for our input
  asmp->Narg = st.Nelem; // number of arguments
  // first compilation pass:
  reg_table_type* rtp = init_reg_table(MAXINST);
  compile_asm_recurse(tree, asmp, *rtp);
  free_reg_table(rtp);
  // second pass (optimize register usage):
  optimize_asm_registers(asmp);
  // that's the best we can hope for:
  return asmp;
}

// do the compilation *and* assembly steps; don't keep the assembler version
// only useful if we don't want to store the assembler externally.
inline program_pointer_type* compile(exprtree_elem_type* tree,
				     symbol_table_type& st,
				     const_table_type& ct) 
{
  program_asm_type* asmp = compile_asm(tree, st);
  program_pointer_type* mach = assemble(asmp, ct);
  free_program_asm(asmp);
  return mach;
}


#endif
