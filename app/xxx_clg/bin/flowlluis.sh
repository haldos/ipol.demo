#!/bin/bash

# usage:
# /path/to/flowlluis.sh a.pgm b.pgm > ab.nd

#PLLUIS=/home/coco/src/flows/lluis
#export LD_LIBRARY_PATH=$PLLUIS
#export PATH=$PLLUIS:$PATH

export LD_LIBRARY_PATH=/home/coco/bak/mw64/lib
PFILE=parameters_34k.cfg

cat <<EOF >$PFILE
# 1) Scheme, 2) Resolution levels  
2 6	
# General aspects of functional: 1) type of data term, 2) regularitzation term, 3) border conditions
1 1 0
# Functional parameters 1) weight regul 2) weight geom 3) weight intensity 4) eps regul 5) c normalisation 6) size_k_sigma, 7) flag gauss, 8) sigma_gauss
180.0 1.0 1.0 0.1 0.05 1 0 0
# Functional gradient computation 1) data, 2) regularization, 3) border image, 4) delta
1 1 0 1e-06
# Truncated newton parameters: 1) ftol, 2) xtol, 3) pgtol, 4) internal iterations
1e-15 1e-15 1e-15 50 #-1 -1 -1 50 
# Image gradient computation 1) type, 2) interpolated pixels, 2) stencil size, 3) threshold, 4) delta precision
1 8 5 0.05 1.0
# Max iterations of TN in 1) multiresolution, 2) coarsest scale of multigrid
20 25
# 1) Initialization of flow field, 2) correct flow available, 3) print objective function
1 0 0
# Specific Multigrid: 1) v or w-cycle, 2) check alpha = 1, 3) gamma contraint, 4) a, 5) e 
0 1 8.0 0.1 0.1
# V-cycle iterations 
2 1 1 1 2 2 2 1 1 1 1 1 1 1
# Presmoothing 
1 20 20 20 20 20 20 20 2 2 2 2 2 2 2 2 2
# Postsmoothing
2 1 1 1 1 1 1

### COMMENTS ON THE PREVIOUS VALUES

#define SCHEME_MULTIRESOLUTION    0
#define SCHEME_MULTIGRID          1
#define SCHEME_FMG                2

#define FUNCTIONAL_DATA_NONE      0
#define FUNCTIONAL_DATA_INTENSITY 1   /* Intensity */
#define FUNCTIONAL_DATA_GEOMETRY  2   /* Gradient orientation */
#define FUNCTIONAL_DATA_OFE       3   /* Optical flow equation */
#define FUNCTIONAL_DATA_COMB_IG   4   /* Combined intensity/geomety */
#define FUNCTIONAL_DATA_GRADCONST 5   /* Gradient constancy */
#define FUNCTIONAL_DATA_GRAD      6   /* Gradient */

#define FUNCTIONAL_REGUL_NONE     0
#define FUNCTIONAL_REGUL_HS       1
#define FUNCTIONAL_REGUL_TV       2
#define FUNCTIONAL_REGUL_TV2      3

#define BORDER_NEUMANN            0
#define BORDER_DIRICHLET          1
#define BORDER_TEST               2

#define GRADIENT_DATA_FINDIFF     0
#define GRADIENT_DATA_ANAL        1

#define GRADIENT_REGUL_FINDIFF    0
#define GRADIENT_REGUL_ANAL       1

#define FILTER_FOURIER            0
#define FILTER_OPTIMUM            1
#define FILTER_GAUSS              2
#define FILTER_OFE_HS             3

#define INIT_FLOW_ZERO            0
#define INIT_FLOW_RAND            1
#define INIT_FLOW_FILE            2

#define CORRECT_FLOW_NOT_AVAIL    0
#define CORRECT_FLOW_ZERO         1
#define CORRECT_FLOW_FILE         2
EOF






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
testex densemotionfield
testex fnorm
testex qeasy


TPD=`mktemp -d /tmp/flowlluis.XXXXXX` || exit 1
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

cat $FILEA | fnorm | qeasy 0 255 >  ${TPD}/filea.pgm
cat $FILEB | fnorm | qeasy 0 255 >  ${TPD}/fileb.pgm
cat $PFILE > ${TPD}/parameters_34k.cfg
cd ${TPD}
densemotionfield -g filea.pgm fileb.pgm >/dev/null
FLUNINDEX_NEGATEY=1 FLUNINDEX_HACKLLU=1 qnm flunindex < motion_field_0.data
cat motion_field_0.data > /tmp/motion_field_0.data


rm -rf $TPD


