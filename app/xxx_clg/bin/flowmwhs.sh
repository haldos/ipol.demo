#!/bin/bash

# usage:
# /path/to/flowmwhs.sh a.pgm b.pgm > ab.nd

export LD_LIBRARY_PATH=/home/coco/bak/mw64/lib

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
testex vecstack
testex hs_flow
testex Mkmovie
testex fnorm
testex qeasy


TPD=`mktemp -d /tmp/flowmwhs.XXXXXX` || exit 1
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

cat $FILEA | fnorm | qeasy 0 255 > ${TPD}/a_1
cat $FILEB | fnorm | qeasy 0 255 > ${TPD}/a_2
cd ${TPD}
Mkmovie Cmovie a 1 2
hs_flow -n 2000 -a 100 a x y >/dev/null
vecstack x_01 y_01


#densemotionfield -g filea.pgm fileb.pgm >/dev/null
#FLUNINDEX_HACKLLU=1 qnm flunindex < motion_field_0.data
#cat motion_field_0.data > /tmp/motion_field_0.data


rm -rf $TPD
