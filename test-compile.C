#include <stdio.h>
#include <stdlib.h>
#include "expr-compile.H"  // expression parser and compiler
#include "expr-optimize.H" // includes our expression tree optimizations

int main (int argc, char** argv) 
{
  int i;
  
  int ERROR = 0;
  char teststring[512];

  verbose_printfuncs(stdout, "# ");
  verbose_printconsts(stdout, "# ");
  
  printf("Enter expression to parse:\n");
  fgets(teststring, sizeof(teststring), stdin);
  i = strlen(teststring); // remove trailing newline
  if (teststring[i-1] == '\n') teststring[i-1] = 0;
  
  const_table_type *ctp = NULL;
  symbol_table_type *stp = NULL;
  exprtree_elem_type *tree = NULL;
  
  ERROR = parse_expr(teststring, ctp, stp, tree);

  if (ERROR) printf("!!!! ERROR: %d !!!!\n", ERROR);
  printf("\nConstant table: %d constants\n", ctp->Nelem);
  for (i=0; i<ctp->Nelem; ++i)
    printf("  c%d %.8le\n", i, ctp->table[i]);
  printf("Symbol table: %d symbols\n", stp->Nelem);
  for (i=0; i<stp->Nelem; ++i)
    printf("  s%d %s\n", i, stp->table[i]);

  if (!ERROR) {
    // Now, traverse the tree:
    printf("\nexpression: %s\n", teststring);
    printf("preoptimization:\n");
    traverse_exprtree(tree, *ctp, *stp);
    // optimize!
    ERROR = optimize_exprtree(tree, *ctp, *stp);
  }
  
  if (!ERROR) {
    // Now, traverse the tree:
    printf("\npostoptimization:\n");
    traverse_exprtree(tree, *ctp, *stp);

    // compile, and compare results:
    program_asm_type* asmp = compile_asm(tree, *stp);
    printf("\nassembly version:\n");
    write_asm(stdout, asmp, *ctp);
    
    program_pointer_type* mach = assemble(asmp, *ctp);

    // do one evaluation cycle:
    double arglist[stp->Nelem];
    char dump[512];
    printf("\nEnter values to evaluate:\n");
    for (i=0; i<stp->Nelem; ++i) {
      printf("%s: ", stp->table[i]);
      fgets(dump, sizeof(dump), stdin);
      sscanf(dump, "%lf", &(arglist[i]));
    }
    printf("tree value = %.8le\n", 
	   eval_exprtree(tree, *ctp, *stp, arglist));
    printf("prog value = %.8le\n", 
	   evalprogram(mach, arglist));
    
    free_program_pointer(mach);
    free_program_asm(asmp);
  }
  else {
    printf("Error in parsing...\n");
  }

  // Garbage collection
  free_const_table(ctp);
  free_symbol_table(stp);
  free_exprtree(tree);

  return 0;
}
