#!/bin/bash

SOURCES="baboon barbara lena peppers spiral"
FLOWS=('tr "1 0"'
'tr "0.5 0.5"'
'tr "0.1 0.1"'
#'tr "0.01 0.01"'
#'hom "1 0.010101 -2.02022 8.78994e-09 1.0202 -4.04042 3.78497e-10 5.0505e-05 0.989899"'
#'hom "1 0.0344828 -2.06897 -4.14107e-08 1.06897 -4.13792 2.6441e-10 0.000574713 0.965517"'
#'hom "1 0.0204082 -2.04081 8.84994e-09 1.04082 -4.08164 -3.43306e-11 0.000204082 0.979592"'
'hom "1 0.015873 -2.03177 0 1.03175 -4.06351 1.38243e-09 0.000124008 0.984127"'
)
for i in $SOURCES; do
	C=0
	echo downsa v 2 $i.png p_$i.png
	for j in "${FLOWS[@]}"; do
		echo synflow $j $i.png PNG:- f_${C}_$i.tiff \| qeasy 0 65535 \| downsa v 2 - po_${C}_$i.png
		echo downsa v 2 f_${C}_$i.tiff - \| faxpb 0.5 0 \| tee pf_${C}_$i.tiff \| viewflow -1 - pv_${C}_$i.png
		C=$[C+1]
	done
done

#usage:
#	synflow {tr|aff|hom|smooth} "params" in out flow
#usage:
#	viewflow satscale flow view
#synflow hom "`^Ctryhom 200 200 200 400 400 200 400 400    200 200  200 400  400 200 398 400 2>/dev/null`" /tmp/clena.png /tmp/oclena.png /tmp/u.tiff
#downsa v 2 f_3_spiral.tiff - | faxpb 0.5 0 | viewflow -1 - PNG:- | d
#1 0.0344828 -2.06897 -4.14107e-08 1.06897 -4.13792 2.6441e-10 0.000574713 0.965517

