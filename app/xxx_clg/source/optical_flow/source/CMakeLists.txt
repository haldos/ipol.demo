SET(LIBS iio nrutil)
ADD_LIBRARY(iio SHARED iio.c iio.h)
ADD_LIBRARY(nrutil SHARED nrutil.c nrutil.h)