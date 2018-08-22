/*
  Program: compiler.C
  Author:  D. Trinkle
  Date:    April 29, 2004
  Purpose: Inputs a mathematical expression with variables, etc., in lisp
           style (see below), and outputs an optimized "assembler" version
	   of the expression that can be efficiently evaluated by other
	   programs.
	   
	   LISP style follows the lisp definition of an expression; i.e.,
	   parenthetical prefix notation.  Namely, an expression is contained
	   in one pair of parentheses, with the first entry being the
	   function name, and subsequent entries are parameters.

	   For example, 1+x becomes (+ 1 x), and sin(pi*x^2) becomes
	   (sin (* pi (^ x 2))).  At some point, someone *might* want
	   to bother writing a normal mixed-infix/prefix style to LISP
	   style converter.  Today is not that day.

  Param.:  <exprfile>
           exprfile: a file where the first non-comment (#) line is a 
	             LISP expression (parenthetic prefix)
	   -g  do not optimize expression before compilation

  Flags:   MEMORY:  used to change const and symbol table sizes
           VERBOSE: output final optimized expression tree and tables
           TESTING: output initial expression tree and tables

  Algo.:   We read in the first non-comment line from our infile as
           our expression.  We parse it, optimize it (unless -g),
	   compile it to assembler, and then output to stdout the resulting
	   "assembly" language program.

  Output:  "Assembly" language to be linked by a separate program upon
           reading.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"
#include "expr-functions.H" // function definitions
#include "expr-programs.H"  // defines assembler and machine code
#include "expr-compile.H"  // expression parser and compiler
#include "expr-optimize.H" // includes our expression tree optimizations

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<exprfile>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'g'};  // flag characters.

const char* ARGEXPL = 
"  exprfile: a file where the first non-comment (#) line is a\n\
	    LISP expression (parenthetic prefix)\n\
  -g  do not optimize expression before compilation\n\
  -m  can set how large the symbol and constant tables should be";


int main (int argc, char** argv) 
{
  int i;
  
  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // If set, we roll our own tables.
  int NO_OPT = 0;   // decide if we don't want to optimize

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
      fprintf(stderr, "Command list:\n");
      verbose_printfuncs(stderr, "  ");
      fprintf(stderr, "Builtin constants:\n");
      verbose_printconsts(stderr, "  ");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags
  NO_OPT = flagon[0];
  
  // ******************************* INPUT *****************************
  char* exprname = args[0];
  FILE* infile;
  char exprstring[16384]; // a VERY big expression
  
  infile = myopenr(exprname);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", exprname);
    exit(ERROR_NOFILE);
  }

  nextnoncomment(infile, exprstring, sizeof(exprstring));
  if (feof(infile)) ERROR=ERROR_BADFILE;
  myclose(infile);

  if (ERROR) {
    fprintf(stderr, "%s didn't contain an expression...?\n", exprname);
    exit(ERROR);
  }

  i = strlen(exprstring); // remove trailing newline
  if (exprstring[i-1] == '\n') exprstring[i-1] = 0;

  // ****************************** ANALYSIS ***************************
  const_table_type *ctp = NULL;
  symbol_table_type *stp = NULL;
  exprtree_elem_type *tree = NULL;

  if (MEMORY > 0) {
    // roll our own...
    ctp = init_const_table(MEMORY);
    stp = init_symbol_table(MEMORY);
  }
  
  if (VERBOSE) {
    printf("# parsing expression: %s\n", exprstring);
    if (NO_OPT)
      printf("# no optimization will be performed.\n");
  }
  
  ERROR = parse_expr(exprstring, ctp, stp, tree);

  if (ERROR) {
    fprintf(stderr, "Parsing error:\n");
    if (has_error(ERROR, PARSE_ERROR_EARLYTERM))
      fprintf(stderr, "String terminated early--not enough args or )\n");
    if (has_error(ERROR, PARSE_ERROR_TOOMANYCONST))
      fprintf(stderr, "Too many constants needed.  Increase -m value\n");
    if (has_error(ERROR, PARSE_ERROR_TOOMANYSYMBOLS))
      fprintf(stderr, "Too many symbols needed.  Increase -m value\n");
    if (has_error(ERROR, PARSE_ERROR_BADSYMBOL))
      fprintf(stderr, "Read a poorly constructed symbol (eg, 1a)--probably missing space\n");
    if (has_error(ERROR, PARSE_ERROR_BADFUNC))
      fprintf(stderr, "Read a bad function name--probably wrong number of parameters, missing ) or typo\n");
  }
  
  if (TESTING) {
    printf("## ==== pre-optimization ====\n");
    printf("## Constant table: %d constants\n", ctp->Nelem);
    for (i=0; i<ctp->Nelem; ++i)
      printf("##  c%d %.8le\n", i, ctp->table[i]);
    printf("##\n## Symbol table: %d symbols\n", stp->Nelem);
    for (i=0; i<stp->Nelem; ++i)
      printf("##  s%d %s\n", i, stp->table[i]);
    if (!ERROR) {
      // Now, traverse the tree:
      printf("##\n## expression tree:\n");
      traverse_exprtree(tree, *ctp, *stp, "## ");
    }
    printf("##\n");
  }
  
  if ( (!ERROR) && (!NO_OPT) ) {
    // optimize!
    ERROR = optimize_exprtree(tree, *ctp, *stp);
  }
  
  if (ERROR) {
    fprintf(stderr, "Error occured during optimization...\n");
  }
  if (VERBOSE) {
    printf("# ==== final expression tree ====\n");
    printf("# Constant table: %d constants\n", ctp->Nelem);
    for (i=0; i<ctp->Nelem; ++i)
      printf("#  c%d %.8le\n", i, ctp->table[i]);
    printf("#\n# Symbol table: %d symbols\n", stp->Nelem);
    for (i=0; i<stp->Nelem; ++i)
      printf("#  s%d %s\n", i, stp->table[i]);
    if (!ERROR) {
      // Now, traverse the tree:
      printf("#\n# expression tree:\n");
      traverse_exprtree(tree, *ctp, *stp, "# ");
    }
    printf("#\n");
  }
  
  if (!ERROR) {
    // compile!!
    program_asm_type* asmp = compile_asm(tree, *stp);
    if (VERBOSE)
      printf("# assembly version:\n#\n");

    // use the expression as comment
    write_asm(stdout, asmp, *ctp, stp, exprstring);
    
    free_program_asm(asmp);
  }

  // Garbage collection
  free_const_table(ctp);
  free_symbol_table(stp);
  free_exprtree(tree);

  return ERROR;
}
