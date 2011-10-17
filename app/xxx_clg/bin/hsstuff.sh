#!/bin/bash

# usage:
# /path/to/stuff.sh alpha niter

# IMPLICIT INPUT FILES
#  a.png
#  b.png
#  t.tiff
#
# IMPLICIT OUTPUT FILES
#
#  stuff_hs.{tiff,png}       "F"
#  stuff_hs_div.{tiff,png}   "divF"
#  stuff_hs_abs.{tiff,png}   "|F|"
#  stuff_hs_inv.png          "F*B"
#  stuff_hs_apinv.png        "(A+F*B)/2"
#  stuff_hs_aminv.{tiff,png} "A-F*B"
#  stuff_hs_fmt.{tiff,png}   "F-T"
#  stuff_hs_afmt.{tiff,png}  "|F-T|"
#




DIVRANG=4
FLORANG=10
IDIFRAD=50
FDIFRAD=4



set -e

testex() {
	if which $1 > /dev/null ; then
		echo > /dev/null
	else
		echo "ERROR: executable file $1 not available" >&2
		exit 1
	fi
}

usgexit() {
	echo -e "usage:\n\t `basename $0` [nico|luis|warp|ldof...]"
	rm -rf $TPD
	exit 1
}

echo STUFF ARGC: $#
echo STUFF ARGV: $*

# check input
if [ $# != "2" ]; then
	usgexit
fi

echo STUFF ARGS: $*

test -f a.png || (echo "I need a file \"a.png\"!" ; exit 1)
test -f b.png || (echo "I need a file \"b.png\"!" ; exit 1)
test -f t.tiff || (echo "I need a file \"t.tiff\"!" ; exit 1)

ALPHA=$1
NITER=$2

FTYPE=hs
P=stuff_hs

echo "ALPHA = ${ALPHA}"
echo "NITER = ${NITER}"
echo "FTYPE = ${FTYPE}"
echo "P = ${P}"

#if [ $FTYPE = "truth" ]; then
#	cp t.tiff ${P}.tiff
#else
#	flowany.sh ${FTYPE} a.png b.png > ${P}.tiff
#fi

testex hs

hs $NITER $ALPHA a.png b.png ${P}.tiff

#  stuff_X.{tiff,png}       "F"
#  stuff_X_div.{tiff,png}   "divF"
#  stuff_X_abs.{tiff,png}   "|F|"
#  stuff_X_inv.png          "F*B"
#  stuff_X_apinv.png        "(A+F*B)/2"
#  stuff_X_aminv.{tiff,png} "A-F*B"
#  stuff_X_ofce.png         "|grad(A)*F + dA/dt|"

viewflow -1 ${P}.tiff ${P}.png
flowdiv ${P}.tiff ${P}_div.tiff
qeasy -$DIVRANG $DIVRANG  ${P}_div.tiff ${P}_div.png
fnorm ${P}.tiff ${P}_abs.tiff
qauto ${P}_abs.tiff ${P}_abs.png
backflow ${P}.tiff b.png | qeasy 0 255 - ${P}_inv.png
plambda a.png ${P}_inv.png "x y + 0 512 qe" | pnmtopng > ${P}_apinv.png
plambda a.png ${P}_inv.png "x y - -${IDIFRAD} ${IDIFRAD} qe" | pnmtopng > ${P}_aminv.png
ofc a.png b.png ${P}.tiff | plambda - "x fabs" | qeasy 0 100 | pnmtopng > ${P}_ofce.png


#  stuff_X_fmt.{tiff,png}   "F-T"
#  stuff_X_afmt.{tiff,png}  "|F-T|"
#  stuff_X_aerr.png  "angle(F,T)"

if [ $FTYPE != "truth" ]; then
	plambda ${P}.tiff t.tiff "x y -" > ${P}_fmt.tiff
	viewflow -1 ${P}_fmt.tiff ${P}_fmt.png
	fnorm ${P}_fmt.tiff ${P}_afmt.tiff
	qeasy 0 $FDIFRAD ${P}_afmt.tiff ${P}_afmt.png
	plambda t.tiff ${P}.tiff "x[0] y[0] * x[1] y[1] * + x[0] x[1] hypot y[0] y[1] hypot * / acos 0 pi qe" | pnmtopng > ${P}_aerr.png
fi

#flowany.sh $FTYPE a.png b.png > f.tiff
#viewflow -1 f.tiff f.png
#flowdiv f.2d
