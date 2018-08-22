#include <stdio.h>
#include <stdlib.h>
#include "io.H"
#include "expr-programs.H"  // our assembler, machine code (contains functions)

// designed to test the assembler...

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<asm>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; 
// flag characters.

const char* ARGEXPL = 
"  asm:  file containing \"assembly\" code\n\
  -m <seed>  for random number generator";


int main (int argc, char** argv) 
{
  int i;

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 0;   // our seed for the random number generator

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
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  srandom(MEMORY); // use our mem. value to seed.  

  char* asmname = args[0];
  FILE* infile;
  
  infile = myopenr(asmname);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", asmname);
    exit(ERROR_NOFILE);
  }
  program_asm_type* asmp = NULL;
  const_table_type* ctp = NULL;
  
  read_asm(infile, asmp, ctp);
  myclose(infile);

  if (VERBOSE) {
    printf("Assembly file contents:\n");
    write_asm(stdout, asmp, *ctp);
  }
  
  program_pointer_type* mach = assemble(asmp, *ctp);

  // do one evaluation cycle:
  double arglist[mach->Narg];
  char dump[512];
  printf("Enter values to evaluate:\n");
  for (i=0; i<mach->Narg; ++i) {
    printf("r%d: ", i);
    fgets(dump, sizeof(dump), stdin);
    sscanf(dump, "%lf", &(arglist[i]));
  }
  printf("function value = %.8le\n", 
	 evalprogram(mach, arglist));

  // Garbage collection
  free_program_asm(asmp);
  free_program_pointer(mach);
  free_const_table(ctp);
  
  return 0;
}
