# we set this so that the makefile doesn't use builtin implicit rules
MAKEFLAGS = -r

# external parameters:
base ?= fcc
kptsize ?= s120
# superind ?= 01 02 03 04 05 06 07 08 09 10 12 14 16 18 20
superind ?= 01 02 03 04 05 06 07 08 09 10 12 14
maxani ?= 3

# this is the directory where every executable lives:
# exec ?= /wight/disk2/trinkle/work/GFBC/EGF
exec ?= $(HOME)/work/GFBC/EGF

# time and echo commands:
time ?= /usr/bin/time
echo ?= echo "====" 

# use lookup table and symmetrization on discDij
DISCOPT ?= -ls

# input files: allow to be overwritten at the command line
cell ?= $(base).cell
potfile ?= $(base).pot
kpttest ?= $(exec)/$(base).kpt.s6
disckpt ?= $(exec)/$(base).kpt.$(kptsize)

# special executables:
graceit := $(exec)/$(base)-graceit
phonongraceit := $(exec)/phonon-graceit
errgraceit := $(exec)/err-graceit
compphonon := $(exec)/comp-phonon.$(base)
quaderror := $(exec)/quad-error.$(base)

# output files
potasmfile := $(base).asm
Dij := $(base).Dij
spectra := $(base).agr
spectralatt := $(base).latt.agr
Ylm := $(base).Ylm
Ylmdc := $(base).Ylm.dc
elas := $(base).elas
disc := $(base).disc
hybrid := $(base).hyb
disckspace := $(base).disc.kpt
discelas := $(base).disc.elas
discelaslong := $(base).disc.elas.long
latt := $(base).latt
diff := $(base).lgf-elas
relerr := $(base).err.rel
phononout := $(base).err.phonon
phononsummary := $(base).err.phonon.agr
quadout := $(base).err.quad
errsummary := $(base).err.rel.agr
## errDij := $(base).err.Dij
## quartic := $(base).err.quartic

# temporary files:
Ylmgrid ?= 256 # for calculating the elastic GF
gridfile := $(exec)/gauss.$(Ylmgrid)

## supercell construction:
superfile := $(foreach index,$(superind),super.$(index))
Dij-fold := $(superfile:super.%=$(Dij).%)
disc-fold := $(superfile:super.%=$(disc).%)
Ylmdc-hyb := $(superfile:super.%=$(Ylmdc).%)
disc-hyb := $(superfile:super.%=$(hybrid).%)

## parse potential file:
potsize  := $(shell head -2 $(potfile) | tail -1)
potname  := $(shell head -3 $(potfile) | tail -1)
gfcutoff ?= $(shell head -4 $(potfile) | tail -1)
longcut  ?= 30.0001
Rmax     := $(shell echo $(gfcutoff) | awk '{print int($$1)}')

## output files:
outfile :=  $(potasmfile) $(Dij) $(spectra) $(spectralatt) \
	$(Ylm) $(Ylmdc) $(elas) $(discelas) $(disc) $(diff) $(latt) \
	$(relerr) $(Dij-fold) \
	$(quadout) $(errsummary)

##	$(disc-fold) $(Ylmdc-hyb) $(disc-hyb) $(phononout) $(phononsummary) 
##	$(errDij) $(quartic) 

.PHONY: clean restart help output super

output: $(outfile)

help:
	@echo "Environment variables to set:"
	@echo "basename:              base     def: $(base)"
	@echo "kpt mesh size for LGF: kptsize  def: $(kptsize)"
	@echo "maximum anisotropy:    maxani   def: $(maxani)"
	@echo "supercells:            supedind def: $(superind)"
	@echo
	@echo "output files:"
	@echo "  $(outfile)"
	@echo
	@echo "Settings:"
	@echo "Executable directory:  exec     def: $(exec)"
	@echo "time command:          time     def: $(time)"
	@echo "echo command:          echo     def: $(echo)"
	@echo "discDij options:       DISCOPT  def: $(DISCOPT)"
	@echo
	@echo "Use restart to restart run from scratch"

super: $(superfile)

######################################################################

## compile our potential into an .asm file
$(potasmfile): $(potfile)
	@$(echo) "Parsing and compiling $(potfile)"
	$(exec)/compiler $(potfile) > $(potasmfile)

## construct a stable dynamical matrix
$(Dij): $(potasmfile) 
	@echo "-----------------------------"
	@echo "Read potential from $(potfile)"
	@echo "Potential function:        `head -1 $(potfile)`"
	@echo "Potential functional form: $(potname)"
	@echo "Potential cutoff:          $(potsize)"
	@echo "GF cutoff:                 $(gfcutoff)"
	@echo "-----------------------------"
	$(time) $(exec)/make-Dij $(cell) $(potsize) $(potasmfile) $(kpttest) $(maxani) $(Dij)

## generate phonon spectra
$(spectra): $(Dij) 
	@$(echo) "Making Dij phonon spectra"
	$(time) $(exec)/phononDij $(cell) $(Dij) \
		$(exec)/kpt.$(base).path -n |\
		$(graceit) "$(potname)" - > $(spectra)

## calculate the elastic GF Ylm, discontinuity correct.
## now including symmetrization, and error checking
$(gridfile):
	@$(exec)/gauss $(Ylmgrid) > $(gridfile)
$(Ylm): $(cell) $(Dij) $(gridfile)
	@$(echo) "Calculating elastic GF"
	$(exec)/gauss-grid $(Ylmgrid) $(gridfile) |\
	$(exec)/egf-FT $(cell) - -gc |\
	$(time) $(exec)/harmonic $(gridfile) - -e -m 92 |\
	$(time) $(exec)/symm-harm $(cell) - > $(Ylm)
	@$(exec)/sphere-tria 8 | $(exec)/egf-FT $(cell) - -gc |\
	$(exec)/check-harm $(Ylm) -


$(Ylmdc): $(cell) $(Dij) $(gridfile) $(Ylm)
	@$(echo) "Calculating the discontinuity correction"
	$(exec)/gauss-grid $(Ylmgrid) $(gridfile) |\
	$(exec)/fourthorder-FT $(cell) $(Ylm) $(Dij) - |\
	$(time) $(exec)/harmonic $(gridfile) - -e -m 92 |\
	$(time) $(exec)/symm-harm $(cell) - > $(Ylmdc)
	@$(exec)/sphere-tria 8 | $(exec)/fourthorder-FT $(cell) $(Ylm) $(Dij) -|\
	$(exec)/check-harm $(Ylmdc) -

## Calculate the EGF in real space
$(elas): $(cell) $(Ylm) 
	@$(echo) "Calculating elastic GF in real space"
	$(exec)/make-ball $(cell) $(gfcutoff) |\
		$(time) $(exec)/real-GF $(cell) $(Ylm) - > $(elas)
$(discelas): $(cell) $(Ylm) $(elas)
	@$(echo) "Calculating elastic GF in real space discretization correc."
	$(time) $(exec)/gf-disc $(cell) $(Ylm) $(elas) > $(discelas)

## calculate the lattice GF: (DISCOPT decides if we use lookup table)
$(disckspace): $(cell) $(Ylm) $(Ylmdc) $(Dij) $(discelas) $(disckpt) 
	@$(echo) "Calculating lattice GF disc. (kspace): $(disckspace)"
	$(exec)/make-ball $(cell) $(gfcutoff) |\
	$(time) $(exec)/discDij $(cell) $(Ylm) $(Ylmdc) $(Dij) \
		$(disckpt) - $(DISCOPT) > $(disckspace)

## polish the lattice GF
$(disc): $(cell) $(Ylm) $(Dij) $(discelas) $(disckpt) $(disckspace)
#	@$(echo) "Polishing lattice GF disc. in real space: $(disc)"
#	$(time) $(exec)/gf-polish $(cell) $(Ylm) $(discelas) $(disckspace) $(Dij) -m 1024 > $(disc)
	@$(echo) "NOT polishing lattice GF disc. in real space: $(disc)"
	cat $(disckspace) > $(disc)

## calculate the difference between LGF and EGF
$(diff): $(disc) $(discelas)
	@$(echo) "Calculating the difference between LGF and EGF"
	$(exec)/scalelatt $(discelas) -1 | $(exec)/addlatt $(disc) - > $(diff)

## lattice GF
$(latt): $(diff)
	@$(echo) "Calculating the LGF"
	$(exec)/addlatt $(diff) $(elas) > $(latt)
$(relerr): $(elas) $(latt)
	$(echo) "Calculating the relative error"
	$(exec)/difflatt $(cell) $(elas) $(latt) -r | sort -n > $(relerr)

## generate phonon spectra:
$(spectralatt): $(Ylm) $(disc) 
	@$(echo) "Making GF phonon spectra: $(spectralatt)"
	$(time) $(exec)/gf-phonon $(cell) $(Ylm) $(disc) \
		$(exec)/kpt.$(base).path -n |\
		$(graceit) "$(potname)" - > $(spectralatt)

## if we need a given supercell...
super.%:
	@echo "-$* $* $*  $* -$* $*  $* $* -$*" > $@
#	@echo "$* 0 0 0 $* 0 0 0 $*" > $@

## now, construct the folded versions:
$(Dij).%: super.% $(cell) $(Dij)
	$(exec)/foldsuper $(cell) $(Dij) $< -s > $@

# note--this is VERY tricky to fix the first line of the discretization correc.
$(disc).%: super.% $(cell) $(disc) $(discelas)
	$(exec)/foldsuper $(cell) $(diff) $< -s |\
		$(exec)/addlatt - $(discelas) |\
		sed "1s/.*/`head -1 $(disc)`/" > $@

$(Ylmdc).%: $(Dij).% $(cell) $(Ylm) $(gridfile)
	$(exec)/gauss-grid $(Ylmgrid) $(gridfile) |\
	$(exec)/fourthorder-FT $(cell) $(Ylm) $< - -y |\
	$(time) $(exec)/harmonic $(gridfile) - -e -m 92 |\
	$(time) $(exec)/symm-harm $(cell) - > $@

$(discelaslong): $(cell) $(Ylm) 
	@$(echo) "Calculating long EGF for hybridization--now empty"
	$(exec)/make-ball $(cell) $(longcut) |\
	$(exec)/randDij $(cell) - $(exec)/empty.asm \
		> $(discelaslong)

#		$(time) $(exec)/real-GF $(cell) $(Ylm) - | \
#		$(time) $(exec)/gf-disc $(cell) $(Ylm) $(elas) \

$(hybrid).%: super.% $(Dij).% $(Ylmdc).% $(cell) $(discelaslong)
	$(time) $(exec)/hybrid4 $(cell) super.$* $(Dij).$* $(Ylmdc).$* -k \
		$(discelaslong) > $@

## phonon analysis:
$(phononout) $(phononsummary): $(Dij-fold) $(disc-fold) $(disc-hyb)
	@$(echo) "Calculating phonon errors"
	$(time) $(compphonon) > $(phononout)
	@$(echo) "Summarizing phonon errors"
	$(phonongraceit) $(phononout) "$(potname)" > $(phononsummary)

## deviation estimates:
$(quadout): $(Dij) $(Ylm) $(Dij-fold)
	@$(echo) "Calculating supercell quadratic errors:"
	$(time) $(quaderror) > $(quadout)

$(errsummary): $(quadout) $(relerr)
	@$(echo) "Summarizing deviations"
	$(errgraceit) $(base) "$(potname)" > $(errsummary)

# Old error estimates (no longer used)
#$(errDij): $(cell) $(Dij) 
#	@$(echo) "Calculating quadratic errors:"
#	$(exec)/errorDij $(cell) $(Dij) $(Rmax) > $(errDij)
#$(quartic): $(disc)
#	@$(echo) "Calculating quartic error:"
## set ld = `head -1 $(disc) | awk '{print $NF}'`
#	echo $(Rmax) `head -1 $(disc) | awk '{print $$NF}'` |\
#	awk 'BEGIN{dx=0.001} {for (x=1.; x<=($$1+0.5*dx); x+=dx) printf("%.5lf %.8le\n", x, $$2/(x*x))}' \
#		> $(quartic)



## designed to remove temporary cruft.
clean:
	-@rm pot.?????? Cij.?????? phonon.?????? cell.??????
	clean-empty -f

restart:
	-rm $(Dij) $(potasmfile) make.dump
