#!/bin/bash

# usage:
# /path/to/flow_thing.sh a.pgm b.pgm > ab.nd

export LD_LIBRARY_PATH=/home/coco/src/prefix/lib:$LD_LIBRARY_PATH

POCV=/home/coco/src/ocvflows
export PATH=$POCV:$PATH


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
testex hs


FILEA=$1
FILEB=$2

usgexit() {
	echo -e "usage:\n\t `basename $0` a.ppm b.ppm > ab.nd"
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

hs $FILEA $FILEB
