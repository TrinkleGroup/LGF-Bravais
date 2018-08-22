#!/bin/tcsh
onintr cleanup

set exec = $0:h

set cell = sc.cell
set basephonon = sc.phonon.grid
set Dij = sc.Dij
set disc = sc.disc
set hyb = sc.hyb
set Ylm = sc.Ylm
set kpt = $exec/sc.kpt.20
# set kpt = $exec/sc.kpt.10

set supercells = (super.*)
# set Dij_comp = (sc.Dij.*)
# set LGF_comp = (sc.disc*)

set output = `mktemp "out.XXXXXX"`

echo "n theta0 theta2 rms : errors"
echo "Dij comparison"
$exec/phononDij $cell $Dij $kpt >! $basephonon
foreach s ($supercells)
    set i = $Dij.$s:e
    $exec/phononDij $cell $i $kpt | $exec/anal-phonon $kpt $basephonon - >! $output
    set theta0 = `head -1 $output | awk '{print $1}'`
    set theta2 = `head -2 $output | tail -1 | awk '{print $1}'`
    set rms = `tail -1 $output | awk '{print $1}'`

    echo "$i:e $theta0 $theta2 $rms"    
end

echo "LGF comparison"
# try comparing to the "true" phonons instead...
# $exec/gf-phonon $cell $Ylm $disc $kpt >! $basephonon
foreach s ($supercells)
    set i = $disc.$s:e
    $exec/gf-phonon $cell $Ylm $i $kpt | $exec/anal-phonon $kpt $basephonon - >! $output
    set theta0 = `head -1 $output | awk '{print $1}'`
    set theta2 = `head -2 $output | tail -1 | awk '{print $1}'`
    set rms = `tail -1 $output | awk '{print $1}'`

    echo "$i:e $theta0 $theta2 $rms"
end


echo "Hybrid comparison"
# we don't bother recalculating basephonon...
# $exec/gf-phonon $cell $Ylm $hyb.20 $kpt >! $basephonon
foreach s ($supercells)
    if ( ($s:e != 1) && ($s:e != 01) ) then
	set i = $hyb.$s:e
	$exec/gf-phonon $cell $Ylm $i $kpt | $exec/anal-phonon $kpt $basephonon - >! $output
	set theta0 = `head -1 $output | awk '{print $1}'`
	set theta2 = `head -2 $output | tail -1 | awk '{print $1}'`
	set rms = `tail -1 $output | awk '{print $1}'`

	echo "$i:e $theta0 $theta2 $rms"
    endif
end

cleanup:
\rm $output
