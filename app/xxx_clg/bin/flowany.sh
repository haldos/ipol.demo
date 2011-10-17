#!/bin/bash

# usage:
# /path/to/flowany.sh [nico|luis|warp|ldof] a.pgm b.pgm > ab.nd


set -e

testex() {
	if which $1 > /dev/null ; then
		echo > /dev/null
	else
		echo "ERROR: executable file $1 not available" >&2
		exit 1
	fi
}


FTYPE=$1
FILEA=$2
FILEB=$3

usgexit() {
	echo -e "usage:\n\t `basename $0` [nico|luis|warp|ldof] a b > ab.nd"
	rm -rf $TPD
	exit 1
}


# check input
if [ $# != "3" ]; then
	usgexit
fi

testex flowbrox_warp.sh
testex flowbrox_ldof.sh
testex flowmwhs.sh
testex flownico.sh
testex flowlluis.sh


#if [ "`file -Lb $IFILE|cut -c1-10`" != "Netpbm PPM" ]; then
#	echo "input file \"$IFILE\" is not in PPM format"
#	exit 1
#fi

case "$1" in
"nico")
	flownico.sh $FILEA $FILEB
	;;
"luis")
	flowlluis.sh $FILEA $FILEB
	;;
"warp")
	flowbrox_warp.sh $FILEA $FILEB
	;;
"ldof")
	flowbrox_ldof.sh $FILEA $FILEB
	;;
"mwhs")
	flowmwhs.sh $FILEA $FILEB
	;;
"ocvlk")
	flowocv_lk.sh $FILEA $FILEB
	;;
"ocvff")
	flowocv_ff.sh $FILEA $FILEB
	;;
"ocvbm")
	flowocv_bm.sh $FILEA $FILEB
	;;
"ocvhs")
	flowocv_hs.sh $FILEA $FILEB
	;;
"ipolhs")
	hs 2000 100 $FILEA $FILEB -
	;;
"jhs")
	jhs 100 1000 $FILEA $FILEB -
	;;
"ipollk")
	lk 13 5 $FILEA $FILEB -
	;;
"ipollkaff")
	lk -2 1 $FILEA $FILEB -
	;;
"ipollkp2")
	lk -3 2 $FILEA $FILEB -
	;;
"ipollkp3")
	lk -3 3 $FILEA $FILEB -
	;;
"ipollkp4")
	lk -3 4 $FILEA $FILEB -
	;;
"ipollkp5")
	lk -3 5 $FILEA $FILEB -
	;;
"ipollkp6")
	lk -3 6 $FILEA $FILEB -
	;;
*)
	usgexit
	;;
esac
