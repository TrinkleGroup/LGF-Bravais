# This is the Makefile for the library routines used in for the random
# number generator

# setenv MallocHelp   -- to see help about tracking memory

OS ?= macosx
# OS ?= linux

CC := gcc
CCPP := g++

INCLUDE := -I. -I../include/
LIBS := -L.

LIBM := -lm
LIBGSL := -lgsl -lgslcblas
ifeq ($(OS),macosx)
#  the first one protects the -framework vecLib: only needed for xlf/xlc
##  LIBLAPACK := -Wl,-framework -Wl,vecLib
#  LIBLAPACK := -framework vecLib
  LIBLAPACK := 
  OSDEFS := -DMACOSX
  LIBS := -L/opt/local/lib $(LIBS)
  INCLUDE := -I/opt/local/include $(INCLUDE)
endif
ifeq ($(OS),linux)
#  LIBLAPACK := -L/usr/local/intel/mkl61/lib/32 -lmkl_def -lmkl_lapack -lguide -lpthread -lmkl
  LIBLAPACK = $(MKL)
  LIBGSL := $(LIBGSL) -lgslcblas
endif

CFLAGS ?= -O5
# have_inline defined for gsl to use extern inline functions when possible:
CXXFLAGS ?= -O5 -DHAVE_INLINE $(OSDEFS)

# bcc map removed, as well as nnpair.H drawfig.H
TARGET = addlatt anal-phonon check-harm compiler difflatt discDij egf-FT elasDij \
	elastic-greens errorDij eval-greens foldsuper fourthorder-FT \
	func gauss-grid gauss gf-disc gf-period gf-phonon gf-polish \
	harmonic-inter harmonic hybrid4 kptmesh magnlatt make-ball phonon \
	phononDij phonon-DOS point-to-maple randDij real-GF scalelatt semicont sphere-tria \
	super-poscar symm-fillin symm-harm symmetrize test-disc test-force \
	test-gsl test-int

INCLUDES = Dij.H ball.H cell.H cont-int.H dcomp.H eigen.H elastic.H \
	integrate.H io.H kpts.H matrix.H pointgroup.H shell.H sphere-harm.H \
	supercell.H unit.H voigt.H fourier.H

.PHONY: all clean archive

all: $(TARGET)

sextic: sextic.C $(INCLUDES)
	$(CCPP) $(CXXFLAGS) $(INCLUDE) $< -o $@ $(LIBLAPACK) $(LIBM)

.SUFFIXES: .C .c .o

.C: $< $(INCLUDES)
	$(CCPP) $(CXXFLAGS) $(INCLUDE) $< -o $@ $(LIBS) $(LIBGSL) $(LIBLAPACK) $(LIBM)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.C.o:
	$(CCPP) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm -f *.o *.a

SOURCES = *.C *.H
GRACEIT = *graceit accumulate.awk
MAKEFILES = Makefile dij.makefile
SCRIPTS = comp-phonon* make-Dij make-seed make.all plot-Dij plot-gf\
	quad-error* run-rand-Dij test-error
ARCHIVEFILES = $(SOURCES) $(GRACEIT) $(MAKEFILES) $(SCRIPTS)

archive:
	@makearchive -z $(ARCHIVEFILES)

GFBC.tar: $(ARCHIVEFILES)
	tar uf GFBC.tar $(ARCHIVEFILES)
