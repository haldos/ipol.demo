#!/bin/bash

# usage:
# /path/to/flowbrox_ldof.sh a.pgm b.pgm > ab.nd

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
testex ldof


TPD=`mktemp -d /tmp/flowbrox_ldof.XXXXXX` || exit 1
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

cat $FILEA | qeasy 0 255 > ${TPD}/filea.ppm
cat $FILEB | qeasy 0 255 > ${TPD}/fileb.ppm
cd ${TPD}
ldof filea.ppm fileb.ppm >/dev/null
cat fileaLDOF.flo | qnm flo2nd


rm -rf $TPD
