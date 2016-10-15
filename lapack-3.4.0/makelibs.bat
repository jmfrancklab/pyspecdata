@echo off
echo %PATH%
echo creating blas library
echo
cd blas\src
make
cd ..\..\
echo creating variants libraries
echo
cd src\variants
make
cd ..\..\
echo creating lapack library
echo
cd src
make
cd ..\
echo done creating blas-lapack libraries
pause
