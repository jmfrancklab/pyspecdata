REM April 2012
REM TEST BLAS/LAPACK LIBRARIES (please make them first, eh?;)
REM START IN LAPACK 3.4.1 TOP DIR
REM make.inc, M_Lintst, M_Eigtst,  must also be present in lapack dir
rem copyright by anonymous April 2012, all rights to public domain
rem for gfortran (mingw - windows)
@echo off
rem TEST LAPACK 3.4.0 libraries
REM TEST BLAS
cd blas\testing
make -f Makeblat1
make -f Makeblat2
make -f Makeblat3
cd ..\
REM level 1 blas tests
xblat1s > sblat1.out
xblat1d > dblat1.out
xblat1c > cblat1.out
xblat1z > zblat1.out
REM level 2 blas tests
xblat2s < sblat2.in
xblat2d < dblat2.in
xblat2c < cblat2.in
xblat2z < zblat2.in
REM level 3 blas tests
xblat3s < sblat3.in
xblat3d < dblat3.in
xblat3c < cblat3.in
xblat3z < zblat3.in
cd ..\
REM TEST LAPACK
REM CREATE TEST MATRIX GENERATION LIBRARY
echo creating matrix generation library
echo
cd testing\matgen
make
cd ..\..\
REM CREATE LAPACK TEST SUITE
echo testing blas\lapack
REM REPLACE EIG and LIN makefiles
xCOPY M_LINTST .\TESTING\LIN\MAKEFILE /Y
xCOPY M_EIGTST .\TESTING\EIG\MAKEFILE /Y
echo create lin test programs
cd testing\lin
make
cd ..\
echo create eig test programs
cd eig
make
cd ..\..\
REM RUN LAPACK TEST SUITE
echo run all tests
cd testing
copy .\eig\*.exe
copy .\lin\*.exe
rem Lin tests 
rem ======== SINGLE LIN TESTS ===========================
xlintsts < stest.in > stest.out
rem ======== COMPLEX LIN TESTS ==========================
xlintstc < ctest.in > ctest.out
rem ======== DOUBLE LIN TESTS ===========================
xlintstd < dtest.in > dtest.out
rem ======== COMPLEX16 LIN TESTS ========================
xlintstz < ztest.in > ztest.out
rem ======== SINGLE-DOUBLE PROTO LIN TESTS ==============
xlintstds < dstest.in > dstest.out
rem ======== COMPLEX-COMPLEX16 LIN TESTS ========================
xlintstzc < zctest.in > zctest.out
rem ======== SINGLE RFP LIN TESTS ========================
xlintstrfs < stest_rfp.in > stest_rfp.out
rem ======== COMPLEX16 RFP LIN TESTS ========================
xlintstrfd < dtest_rfp.in > dtest_rfp.out
rem ======== COMPLEX16 RFP LIN TESTS ========================
xlintstrfc < ctest_rfp.in > ctest_rfp.out
rem ======== COMPLEX16 RFP LIN TESTS ========================
xlintstrfz < ztest_rfp.in > ztest_rfp.out

rem  ======== END LIN TESTS ================
rem Eig tests 
rem ======= SINGLE EIG TESTS ===========================
xeigtsts < nep.in > snep.out
xeigtsts < sep.in > ssep.out
xeigtsts < svd.in > ssvd.out
xeigtsts < sec.in > sec.out
xeigtsts < sed.in > sed.out
xeigtsts < sgg.in > sgg.out
xeigtsts < sgd.in > sgd.out
xeigtsts < ssb.in > ssb.out
xeigtsts < ssg.in > ssg.out
xeigtsts < sbal.in > sbal.out
xeigtsts < sbak.in > sbak.out
xeigtsts < sgbal.in > sgbal.out
xeigtsts < sgbak.in > sgbak.out
xeigtsts < sbb.in > sbb.out
xeigtsts < glm.in > sglm.out
xeigtsts < gqr.in > sgqr.out
xeigtsts < gsv.in > sgsv.out
xeigtsts < csd.in > scsd.out
xeigtsts < lse.in > slse.out
rem ======= COMPLEX EIG TESTS ===========================
xeigtstc < nep.in > cnep.out
xeigtstc < sep.in > csep.out
xeigtstc < svd.in > csvd.out
xeigtstc < cec.in > cec.out
xeigtstc < ced.in > ced.out
xeigtstc < cgg.in > cgg.out
xeigtstc < cgd.in > cgd.out
xeigtstc < csb.in > csb.out
xeigtstc < csg.in > csg.out
xeigtstc < cbal.in > cbal.out
xeigtstc < cbak.in > cbak.out
xeigtstc < cgbal.in > cgbal.out
xeigtstc < cgbak.in > cgbak.out
xeigtstc < cbb.in > cbb.out
xeigtstc < glm.in > cglm.out
xeigtstc < gqr.in > cgqr.out
xeigtstc < gsv.in > cgsv.out
xeigtstc < csd.in > ccsd.out
xeigtstc < lse.in > clse.out
rem ======= DOUBLE EIG TESTS ===========================
xeigtstd < nep.in > dnep.out
xeigtstd < sep.in > dsep.out
xeigtstd < svd.in > dsvd.out
xeigtstd < dec.in > dec.out
xeigtstd < ded.in > ded.out
xeigtstd < dgg.in > dgg.out
xeigtstd < dgd.in > dgd.out
xeigtstd < dsb.in > dsb.out
xeigtstd < dsg.in > dsg.out
xeigtstd < dbal.in > dbal.out
xeigtstd < dbak.in > dbak.out
xeigtstd < dgbal.in > dgbal.out
xeigtstd < dgbak.in > dgbak.out
xeigtstd < dbb.in > dbb.out
xeigtstd < glm.in > dglm.out
xeigtstd < gqr.in > dgqr.out
xeigtstd < gsv.in > dgsv.out
xeigtstd < csd.in > dcsd.out
xeigtstd < lse.in > dlse.out
rem ======= COMPLEX16 EIG TESTS ===========================
xeigtstz < nep.in > znep.out
xeigtstz < sep.in > zsep.out
xeigtstz < svd.in > zsvd.out
xeigtstz < zec.in > zec.out
xeigtstz < zed.in > zed.out
xeigtstz < zgg.in > zgg.out
xeigtstz < zgd.in > zgd.out
xeigtstz < zsb.in > zsb.out
xeigtstz < zsg.in > zsg.out
xeigtstz < zbal.in > zbal.out
xeigtstz < zbak.in > zbak.out
xeigtstz < zgbal.in > zgbal.out
xeigtstz < zgbak.in > zgbak.out
xeigtstz < zbb.in > zbb.out
xeigtstz < glm.in > zglm.out
xeigtstz < gqr.in > zgqr.out
xeigtstz < gsv.in > zgsv.out
xeigtstz < csd.in > zcsd.out
xeigtstz < lse.in > zlse.out
cd ..\
echo done
pause
