#!/bin/tcsh

# Recently rewritten to calculate _both_ error estimates, and add 'em.
# This all assumes that our cell volume = 1 -- in the future, we'll fix this

set exec = $0:h

set base = sc
set cell = $base.cell

set Dij = $cell:r.Dij
set Dij_comp = ($Dij.*)
set Ylm = $base.Ylm
# pull off the 00 component:
set Ylm00 = `head -2 $Ylm | tail -1 | awk '{print $3}'`

# needed for getting 00 component of the fourth order term:
set gaussgrid = $exec/gauss.16

set C11_2C44 = `$exec/elasDij $cell $Dij | head -1 | awk '{print $2+2*$4}'`

foreach i ($Dij_comp)
    # Real space first:
    set C11_2C44_new = `$exec/elasDij $cell $i | head -1 | awk '{print $2+2*$4}'`
    set D00 = `head -2 $i | tail -1 | awk '{print $4}'`
    set err_R = `echo $C11_2C44 $C11_2C44_new $D00 $i:e | awk '{x = 2.*($1-$2)/($3*$4*4); if (x<0) x=-x; print sqrt(3)*x}'`
    # handle division by zero errors...
    if ("$err_R" == "") set err_R=0

    # Fourier space next:
    set Y4lm00 = `$exec/gauss-grid $gaussgrid:e $gaussgrid | $exec/fourthorder-FT $cell $Ylm $i - -y | $exec/harmonic $gaussgrid - -m 0 | $exec/symm-harm $cell - | head -2 | tail -1 | awk '{print $3}'`
    set err_k = `echo $Ylm00 $Y4lm00 $i:e | awk '{x=4.9*$2/($1*$3*$3); if (x>0) print x; else print -x;}'`
    # handle division by zero errors...
    if ("$err_k" == "") set err_k=0

    set err=`echo $err_R $err_k| awk '{print $1+$2}'`

    echo $i:e $err $err_R $err_k
end
