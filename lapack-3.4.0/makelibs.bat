REM April 2012
REM CREATE BLAS/LAPACK 3.4.0 LIBRARIES 
REM START IN LAPACK TOP DIR
REM make.inc must also be present
rem copyright by anonymous june 2012, all rights to public domain
rem for gfortran (mingw - windows)

@echo off
REM MACHINE DEPENDENT TEST SUITE
echo machine dependent routines and testing
cd install
make
del *.o
REM run machine dependent tests
testlsame > testlsame.out
testslamch > testslamch.out
testdlamch > testdlamch.out 
testsecond > testsecond.out
testdsecnd > testdsecnd.out
testieee > testieee.out
testversion > testversion.out
cd ..\
REM CREATE BLAS
echo creating blas library
echo
cd blas\src
make
del *.o
cd ..\..\
REM  CREATE VARIANTS LIBRARIES
echo creating variants libraries
echo
cd src\variants
make
del *.o
cd ..\..\
REM CREATE LAPACK LIBRARY
echo creating lapack library
echo
cd src
make
del *.o
cd ..\
echo done creating blas-lapack libraries
pause
