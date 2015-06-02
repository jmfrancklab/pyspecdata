!compile with
!f2py.py -c --fcompiler=gnu95 --compiler=mingw32 -lmsvcr71 -m foo foo.f90
subroutine otherhello ()
    write(*,*)'This is a different fortran90 test'
end subroutine otherhello
