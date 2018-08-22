#ifndef __EXPR_PROGRAMS_H
#define __EXPR_PROGRAMS_H
/*
  Program: expr-programs.H
  Author:  D. Trinkle
  Date:    April 27, 2004
  Purpose: To provide expression evaluation in realtime.  We store
           the different program types, an assembler, and an execution
	   of the virtual register machine.
*/

#include <ctype.h>
#include <string.h> // needed for parsing...
#include "dcomp.H"
#include "expr-functions.H" // contains all of the functions, constants, etc.

// ****************************** STRUCTURES ***************************
// We have to define the constant table *here* because it's used by the
// assembler.  For completeness, we also define the symbol table here, too.

typedef struct 
{
  int Nelem;
  int maxelem;
  double* table;
} const_table_type;

typedef struct 
{
  int Nelem;
  int maxelem;
  char** table;
} symbol_table_type;



// Two ways of specifying a program; call the first "assembler" and
// the second "machine code."
typedef struct 
{ // note: there is NO register space allocated here--that comes in assembly
  int Nreg; // size of memory space *without* constants
  int Narg; // number of unbound symbols
  int Ninst;       // number of instructions to execute
  int* inst_index; // list of instructions, indexed
  int* Nargs;      // number of arguments for each
  int* reslist;    // result list: where we store each instruction result
  int** arglist;   // argument list: index to registers
} program_asm_type;

typedef struct 
{
  int Nreg; // total number of registers--includes constants
  int Narg; // number of arguments we expect (fill first set of registers)
  double *regs;
  int Ninst;                      // number of instructions to execute
  numerical_function_type** inst; // direct index to the function pointer
  int* Nargs;                     // number of arguments for each
  double** reslist;               // result list: list of pointer to reg.
  double*** arglist;              // argument list: list of pointers to reg.
} program_pointer_type;
  

// ****************************** CONSTANTS ****************************
// "Hard" limits on number of unique numerical constants and symbols
// when parsing; may be overridden by allocating own tables before
// parse_expr


// ***************************** SUBROUTINES **************************
// Table creation and access routines:

const_table_type* init_const_table (int Nelem) 
{
  const_table_type* ctp = new const_table_type;
  ctp->Nelem = 0;
  ctp->maxelem = Nelem;
  ctp->table = new double[Nelem];
  return ctp;
}

void free_const_table (const_table_type* &ctp) 
{
  if (ctp != NULL) {
    delete[] ctp->table;
    delete ctp;
    ctp = NULL;
  }
}

int get_value_tree (double value, const_table_type &ct) 
{
  int i;
  // first, see if value is there:
  for (i=0; (i<ct.Nelem) && (value != ct.table[i]); ++i) ;
  if (i!=ct.Nelem) return i;
  // else, we have to add it:
  if (ct.Nelem == ct.maxelem) return -1; // sanity check
  ct.table[i] = value;
  ++(ct.Nelem);
  return i;
}

symbol_table_type* init_symbol_table (int Nelem) 
{
  int i;
  symbol_table_type* stp = new symbol_table_type;
  stp->Nelem = 0;
  stp->maxelem = Nelem;
  stp->table = new char*[Nelem];
  for (i=0; i<Nelem; ++i) stp->table[i] = NULL;
  return stp;
}

void free_symbol_table (symbol_table_type* &stp) 
{
  int i;
  if (stp != NULL) {
    for (i=0; i<(stp->Nelem); ++i)
      if (stp->table[i] != NULL) delete[] stp->table[i];
    delete[] stp->table;
    delete stp;
    stp = NULL;
  }
}

int get_symbol_tree(char* sym_string, symbol_table_type &st) 
{
  int i;
  // first, see if symbol is there:
  for (i=0; (i<st.Nelem) && (strcmp(sym_string, st.table[i]) != 0); ++i) ;
  if (i!=st.Nelem) return i;
  // else, we have to add it:
  if (st.Nelem == st.maxelem) return -1; // sanity check
  st.table[i] = new char[strlen(sym_string)+1];
  strcpy(st.table[i], sym_string);
  ++(st.Nelem);
  return i;
}


program_asm_type* init_program_asm(int Ninst) 
{
  program_asm_type* asmp = new program_asm_type;
  asmp->Nreg = 0;
  asmp->Narg = 0;
  asmp->Ninst = 0;
  asmp->inst_index = new int[Ninst];
  asmp->Nargs = new int[Ninst];
  asmp->reslist = new int[Ninst];
  asmp->arglist = new int *[Ninst];
  return asmp;
}

void free_program_asm (program_asm_type* &asmp) 
{
  if (asmp != NULL) {
    for (int i=0; i<asmp->Ninst; ++i)
      delete[] asmp->arglist[i];
    delete[] asmp->arglist;
    delete[] asmp->reslist;
    delete[] asmp->Nargs;
    delete[] asmp->inst_index;
    delete asmp;
    asmp = NULL;
  }
}

void free_program_pointer (program_pointer_type* &mach) 
{
  if (mach != NULL) {
    for (int i=0; i<mach->Ninst; ++i)
      delete[] mach->arglist[i]; // we only allocated a list of pointers
    delete[] mach->arglist;
    delete[] mach->reslist;
    delete[] mach->Nargs;
    delete[] mach->inst;
    delete[] mach->regs;
    delete mach;
    mach = NULL;
  }
}


// Our compiler and our assembler.  The assembler takes an "assembly lang."
// version of our program and turns all the indices into pointers for
// efficient evaluation.  It also puts all the constants into registers

// The constants go at the top of the register space; in the assembler,
// a constant has address -(ct.index+1); we can then access each constant
// in register (nreg + address)%nreg.  This doesn't move the registers,
// and puts the constants at the top.
program_pointer_type* assemble(program_asm_type* asmp, const_table_type& ct) 
{
  int ninst = asmp->Ninst, nreg = asmp->Nreg + ct.Nelem;
  program_pointer_type* mach = new program_pointer_type;
  mach->Nreg = nreg; // space for constant and registers
  mach->Narg = asmp->Narg;
  mach->regs = new double[nreg];
  mach->Ninst = ninst;
  mach->inst = new numerical_function_type*[ninst];
  mach->Nargs = new int[ninst];
  mach->reslist = new double*[ninst];
  mach->arglist = new double**[ninst];
  // put constants into place:
  for (int c=0; c<ct.Nelem; ++c)
    mach->regs[(nreg-1-c)%nreg] = ct.table[c];

  // assemble the functions and arguments:
  for (int i=0; i<asmp->Ninst; ++i) {
    mach->inst[i] = NUM_FUNCTIONS[asmp->inst_index[i]].num_func;
    mach->reslist[i] = mach->regs + asmp->reslist[i];
    int nargs = asmp->Nargs[i];
    mach->Nargs[i] = nargs;
    mach->arglist[i] = new double*[nargs];
    for (int j=0; j<nargs; ++j)
      mach->arglist[i][j] = mach->regs + ((nreg+asmp->arglist[i][j])%nreg);
  }
  
  return mach;
}


// try to inline this... dunno if the g++ will actually do that.
inline double evalprogram(program_pointer_type* mach, double *arglist) 
{
  // copy arglist to our registers:
  for (int i=0; i<mach->Narg; ++i) mach->regs[i] = arglist[i];
  // execute commands, and put result in mach->reslist[i]
  for (int i=0; i<mach->Ninst; ++i)
    *(mach->reslist[i]) = mach->inst[i] (mach->Nargs[i], mach->arglist[i]);
  return mach->regs[0]; // should always have our final value
}


// output assembler...
void write_asm(FILE* outfile, program_asm_type* &asmp, const_table_type& ct,
	       symbol_table_type* stp, const char* commentline);

inline void write_asm(FILE* outfile, program_asm_type* &asmp,
		      const_table_type& ct, symbol_table_type* stp)
{
  write_asm(outfile, asmp, ct, stp, NULL);
}

inline void write_asm(FILE* outfile, program_asm_type* &asmp,
		      const_table_type& ct, const char* commentline) 
{
  write_asm(outfile, asmp, ct, NULL, commentline);
}

inline void write_asm(FILE* outfile, program_asm_type* &asmp,
		      const_table_type& ct) 
{
  write_asm(outfile, asmp, ct, NULL, NULL);
}

// we only use the symbol table to output verbosely which symbols
// have remained (in case some are lost), and gives the register to
// symbol map.  Very important...
void write_asm(FILE* outfile, program_asm_type* &asmp, const_table_type& ct,
	       symbol_table_type* stp, const char* commentline)
{
  fprintf(outfile, "%d %d %d # Nreg, Narg ( ",
	  asmp->Nreg, asmp->Narg, ct.Nelem);
  if (stp!=NULL)
    for (int s=0; s<stp->Nelem; ++s)
      fprintf(outfile, "r%d:%s ", s, stp->table[s]);
  fprintf(outfile, "), Nconst");
  if (commentline == NULL) fprintf(outfile, "\n");
  else fprintf(outfile, " | %s\n", commentline);
  for (int c=0; c<ct.Nelem; ++c)
    fprintf(outfile, "  %.17le  # c%d\n", ct.table[c], c);
  fprintf(outfile, "%d # number of instructions\n", asmp->Ninst);
  for (int i=0; i<asmp->Ninst; ++i) {
    fprintf(outfile, "%s %d %d", NUM_FUNCTIONS[asmp->inst_index[i]].name,
	    asmp->reslist[i], asmp->Nargs[i]);
    for (int j=0; j<asmp->Nargs[i]; ++j)
      fprintf(outfile, " %d", asmp->arglist[i][j]);
    fprintf(outfile, " # %s",  NUM_FUNCTIONS[asmp->inst_index[i]].name);
    for (int j=0; j<asmp->Nargs[i]; ++j)
      if (asmp->arglist[i][j] < 0)
	fprintf(outfile, " %.8le", ct.table[-(asmp->arglist[i][j]+1)]);
      else
	fprintf(outfile, " r%d", asmp->arglist[i][j]);
    fprintf(outfile, " -> r%d\n", asmp->reslist[i]);
  }
}

// read in assembler...
void read_asm(FILE* infile, program_asm_type* &asmp, const_table_type* &ctp) 
{
  char dump[512];
  int Nreg, Narg, Nconst, Ninst;

  // number of register, arguments, and constants:
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", &Nreg, &Narg, &Nconst);

  // read constants
  ctp = init_const_table(Nconst);
  ctp->Nelem = Nconst;
  for (int i=0; i<Nconst; ++i) {
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%lf", ctp->table+i);
  }

  // number of instructions:
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d", &Ninst);
  
  // read instructions
  asmp = init_program_asm(Ninst);
  asmp->Ninst = Ninst;
  asmp->Nreg = Nreg;
  asmp->Narg = Narg;
  for (int i=0; i<Ninst; ++i) {
    char instname[512];
    nextnoncomment(infile, dump, sizeof(dump));
    sscanf(dump, "%s %d %d", instname, &(asmp->reslist[i]), &(asmp->Nargs[i]));
    asmp->inst_index[i] = get_num_func_index(instname, asmp->Nargs[i]);
    // Now, read the arguments:
    char *beginp = dump, *endp;
    while (isblank(*beginp)) ++beginp;  // pass up whitespace...
    while (!isblank(*beginp)) ++beginp; // pass up the instruction...
    endp = beginp;
    // pass up the result arg and number of args:
    strtol(beginp, &endp, 0); beginp=endp;
    strtol(beginp, &endp, 0); beginp=endp;
    // NOW, we can read the arguments
    asmp->arglist[i] = new int[asmp->Nargs[i]];
    for (int j=0; j<asmp->Nargs[i]; ++j) {
      asmp->arglist[i][j] = strtol(beginp, &endp, 0);
      beginp=endp;
    }
  }
}


// read in the assembler code, assemble it, and toss it.  Return the
// machine language version:
void read_mach(FILE* infile, program_pointer_type* &mach)
{
  program_asm_type* asmp = NULL;
  const_table_type* ctp = NULL;

  read_asm(infile, asmp, ctp);
  mach = assemble(asmp, *ctp);

  free_const_table(ctp);
  free_program_asm(asmp);
}

#endif
