#!/bin/bash

# usage:
# /path/to/flowbrox_warp.sh a.pgm b.pgm > ab.nd


PLLUIS=/home/coco/src/flows/nico
export PATH=$PLLUIS:$PATH

export LD_LIBRARY_PATH="/usr/local/cuda/lib"

NICO_EXEC=optical_flow_narrow_GPU
NICO_OPTIONS="-b 2 -p 10 -mu -5 -Mu 5 -mv -5 -Mv 5 -nu 20 -nv 20 -ru 50 -rv 50 -o 50 -mi 0"


set -e

testex() {
	if which $1 > /dev/null ; then
		echo > /dev/null
	else
		echo "ERROR: executable file $1 not available" >&2
		exit 1
	fi
}

testex qnm
testex ppmtopgm
testex pnmtopng
testex $NICO_EXEC


TPD=`mktemp -d /tmp/flownico.XXXXXX` || exit 1
FILEA=$1
FILEB=$2

usgexit() {
	echo -e "usage:\n\t `basename $0` a.ppm b.ppm > ab.nd"
	rm -rf $TPD
	exit 1
}


# check input
if [ $# != "2" ]; then
	usgexit
fi

#if [ "`file -Lb $IFILE|cut -c1-10`" != "Netpbm PPM" ]; then
#	echo "input file \"$IFILE\" is not in PPM format"
#	exit 1
#fi

echo "TPD = ${TPD}" 1>&2

cat $FILEA > ${TPD}/filea.ppm
cat $FILEB > ${TPD}/fileb.ppm
cd ${TPD}
MIDA=`ppmtopgm < filea.ppm | qnm getsize`
mkdir data
pnmtopng filea.ppm > data/frame10.png
pnmtopng fileb.ppm > data/frame11.png
rm filea.ppm fileb.ppm
$NICO_EXEC $NICO_OPTIONS >oo 2>oe
qnm rawatoasc 1 $MIDA 1 < Results2/Field_Field_00_u > fx
qnm rawatoasc 1 $MIDA 1 < Results2/Field_Field_00_v > fy
qnm vecstack fx fy


rm -rf $TPD
