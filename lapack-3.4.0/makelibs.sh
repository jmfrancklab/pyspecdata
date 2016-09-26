#!/bin/bash
# April 2012
# CREATE BLAS/LAPACK 3.4.0 LIBRARIES 
# START IN LAPACK TOP DIR
# make.inc must also be present
# copyright by anonymous june 2012, all rights to public domain
# for gfortran (mingw - windows)

## MACHINE DEPENDENT TEST SUITE
##echo "machine dependent routines and testing"
##cd install
##make
###rm *.o
### run machine dependent tests
##testlsame > testlsame.out
##testslamch > testslamch.out
##testdlamch > testdlamch.out 
##testsecond > testsecond.out
##testdsecnd > testdsecnd.out
##testieee > testieee.out
##testversion > testversion.out
##cd ../
## CREATE BLAS
#echo "creating blas library"
#echo
#cd blas/src
#make
##rm *.o
#cd ../../
##  CREATE VARIANTS LIBRARIES
#echo "creating variants libraries"
#echo
#cd src/variants
#make
##rm *.o
#cd ../../
# CREATE LAPACK LIBRARY
echo "creating lapack library"
echo
cd src
make
#rm *.o
cd ../
echo "done creating blas-lapack libraries"
