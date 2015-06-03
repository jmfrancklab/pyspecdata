!compile with
!f2py.py -c --fcompiler=gnu95 --compiler=mingw32 -lmsvcr71 -m foo foo.f90
subroutine hello ()
    write(*,*)'Hello from Fortran90!!! (modified)'
    write(*,*)'and here is something from c:'
    call lprmpt
    write(*,*)'Yay!'
end subroutine hello
