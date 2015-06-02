module propagator
!complex*16, allocatable, dimension(:,:) :: b
double complex, allocatable, dimension(:,:,:) :: operators ! dim(h1,h2,operators)
double precision, allocatable, dimension(:,:) :: time_dependence ! dim(operators,timecounter)
double precision, allocatable, dimension(:,:) :: p ! dim(operators,parameter index)
double complex, allocatable, dimension(:,:,:,:) :: propagators ! dim(h1,h2,timecounter,parameter index)
double complex, allocatable, dimension(:,:,:,:) :: rho ! dim(h1,h2,timecounter,parameter index)
double complex, allocatable, dimension(:,:) :: traceout ! dim(timecounter,parameter index)
double complex, allocatable, dimension(:,:,:) :: C ! dim(h1,h2,parameter index)
double complex, allocatable, dimension(:,:) :: gradient ! dim(timecounter,parameter)
double complex, allocatable, dimension(:,:) :: H_c ! dim(h1,h2)
contains

subroutine clear
implicit none
write(*,*) '----------------------------------------------------------'
if (allocated(operators)) then
	deallocate(operators)
endif
if (allocated(time_dependence)) then
	deallocate(time_dependence)
endif
if (allocated(propagators)) then
	deallocate(propagators)
endif
if (allocated(rho)) then
	deallocate(rho)
endif
if (allocated(p)) then
	deallocate(p)
endif
write(*,*) '----------------------------------------------------------'
return
end subroutine clear

subroutine printvar(input,n_j,n_k)
implicit none
real*8, intent(in), dimension(:,:):: input
integer, intent(in):: n_j,n_k
integer:: j,k
write(*,*) '----------------------------------------------------------'
write (*,*) 'size is: ',n_j,n_k
do j = 1,n_j
	do k = 1,n_k
		write(*,*) input(j,k)
	enddo
enddo
write(*,*) '----------------------------------------------------------'
return
end subroutine printvar

subroutine nothing
implicit none
write(*,*) '----------------------------------------------------------'
write(*,*) '----------------------------------------------------------'
return
end subroutine nothing

subroutine sandwich
implicit none
integer p_counter,t_counter
double complex tempmat(size(propagators,1),size(propagators,2)) 
write(*,*) '------------------------------------------------------------'
write (*,*) 'sandwiching rho of size ',size(rho,1),size(rho,2),size(rho,3),size(rho,4),&
    ' with propagator of size ',size(propagators,1),size(propagators,2),&
    size(propagators,3),size(propagators,4)
if (allocated(propagators)) then
	write (*,*) 'propagators is allocated for size ',size(propagators)
else
	write (*,*) 'propagators is NOT allocated'
endif
if (size(rho,4).eq.size(propagators,4).and.size(rho,3).eq.size(rho,3)) then
	write (*,*) 'case 1'
	do p_counter = 1,size(propagators,4)
		do t_counter = 1,size(propagators,3)
			call sandwich_internal(rho(:,:,t_counter,p_counter),propagators(:,:,t_counter,p_counter),propagators(:,:,t_counter,p_counter))
		enddo
	enddo
elseif (size(rho,4).eq.1.and.size(rho,3).eq.1) then
	write (*,*) 'case 2'
	do p_counter = 1,size(propagators,4)
		do t_counter = 1,size(propagators,3)
			!write (*,*) 'multiply',t_counter,p_counter
			propagators(:,:,t_counter,p_counter) =&
                matmul(propagators(:,:,t_counter,p_counter),&
                matmul(rho(:,:,1,1) ,&
                conjg(transpose(propagators(:,:,t_counter,p_counter)))))
		enddo
	enddo
else
	write(*,*) 'Error! for sandwiching, the t and p dimensions must either match the propagator or be 1!'
endif
write(*,*) '------------------------------------------------------------'
return
end subroutine sandwich

subroutine lvdot
implicit none
integer p_counter,t_counter
double complex tempmat(size(rho,1),size(rho,2)) 
write(*,*) '------------------------------------------------------------'
write (*,*) 'lvdot rho of size ',size(rho,1),size(rho,2),size(rho,3),size(rho,4)
write (*,*) 'with C of size ',size(C,1),size(C,2),size(C,3)
write (*,*) 'putting into traceout of size ',size(traceout,1),size(traceout,2)
if (size(rho,3).eq.size(traceout,1).and.size(rho,4).eq.size(traceout,2)) then
	write (*,*) 'case 2'
	do t_counter = 1,size(rho,3)
		do p_counter = 1,size(rho,4)
			!write (*,*) 'multiply',t_counter,p_counter
			if (size(C,3).gt.1) then
				tempmat = matmul(rho(:,:,t_counter,p_counter),conjg(transpose(C(:,:,p_counter))))
			else
				tempmat = matmul(rho(:,:,t_counter,p_counter),conjg(transpose(C(:,:,1))))
			endif
			traceout(t_counter,p_counter) = 0
			call trace_internal(tempmat,traceout(t_counter,p_counter))
		enddo
	enddo
else
	write(*,*) 'Error! for lvdot, traceout isn''t properly allocated'
endif
write(*,*) '------------------------------------------------------------'
return
end subroutine lvdot

subroutine commutator_internal(mat1,mat2)
double complex, intent(inout), dimension(:,:):: mat1
double complex, intent(in), dimension(:,:):: mat2
mat1 = matmul(mat1,mat2) - matmul(mat2,mat1)
return
end subroutine commutator_internal

subroutine trace_internal(mat,trace)
double complex, intent(out):: trace
double complex, intent(in), dimension(:,:):: mat
integer i_counter, j_counter
trace = 0
do i_counter = 1,size(mat,1)
	trace = trace + mat(i_counter,i_counter)
enddo
return
end subroutine trace_internal

subroutine propagate
implicit none
double complex temp_exp(size(operators,1),size(operators,2)) 
integer t_counter,op_counter,p_counter
!allocate(propagators(size(operators,1),size(operators,2)),size(time_dependence))
write(*,*) '------------------------------------------------------------'
write (*,*) 'shape of operators',size(operators,1),size(operators,2),size(operators,3)
write (*,*) 'shape of p',size(p,1),size(p,2)
write (*,*) 'shape of time_dependence',size(time_dependence,1),size(time_dependence,2)
if (size(operators,1) .ne. size(operators,2)) then
	write (*,*) 'Operators are not square!!'
	return
endif
if (size(operators,3) .ne. size(time_dependence,1)) then
	write (*,*) 'Number of operators specified by array of operators and array of time dependence is not the same!!'
	return
endif
if (size(operators,3) .ne. size(p,1)) then
	write (*,*) 'Number of operators specified by array of operators and array of parameters is not the same!!'
	return
endif
if (size(propagators,3) .gt. 1) then
	do p_counter = 1,size(p,2)
		!write(*,*) 'propagating parameter ',p_counter,' out of ',size(p,2)
		do t_counter = 1,size(time_dependence,2)
			propagators(:,:,t_counter,p_counter) = 0
			do op_counter = 1,size(time_dependence,1)
				propagators(:,:,t_counter,p_counter) =&
                    propagators(:,:,t_counter,p_counter) +&
                    p(op_counter,p_counter) *&
                    time_dependence(op_counter,t_counter) *&
                    operators(:,:,op_counter)
			end do
			call complexexp(propagators(:,:,t_counter,p_counter)) ! this passes a HERMITIAN argument, and gets back the exponent of the 2*pi*i* that argument
			!call exact(propagators(:,:,t_counter,p_counter)) ! this passes a HERMITIAN argument, and gets back the exponent
			! for now, i'm not actually propagating, just exponentiating
			if (t_counter.gt.1) then
				propagators(:,:,t_counter,p_counter) = matmul( propagators(:,:,t_counter,p_counter), propagators(:,:,t_counter-1,p_counter))
			endif
		enddo
		!write(*,*) 'done propagating parameter ',p_counter,' out of ',size(p,2)
	enddo
else
	write (*,*) 'not saving all time points'
	do p_counter = 1,size(p,2)
		do t_counter = 1,size(time_dependence,2)
			temp_exp = 0
			do op_counter = 1,size(time_dependence,1)
				temp_exp = temp_exp + p(op_counter,p_counter) * time_dependence(op_counter,t_counter) * operators(:,:,op_counter)
			end do
			call complexexp(temp_exp) ! this passes a HERMITIAN argument, and gets back the exponent of the 2*pi*i* that argument
			! for now, i'm not actually propagating, just exponentiating
			if (t_counter.gt.1) then
				propagators(:,:,1,p_counter) = matmul( temp_exp, propagators(:,:,1,p_counter))
			else
				propagators(:,:,1,p_counter) = temp_exp
			endif
		enddo
	enddo
endif
write(*,*) '------------------------------------------------------------'
return
end subroutine propagate

subroutine grape
implicit none
! this gives the derivative of the inner product
integer t_steps,t_counter,p_counter
double complex tempmat(size(propagators,1),size(propagators,2)) 
double complex end_sandwich(size(propagators,1),size(propagators,2)) 
write(*,*) '------------------------------------------------------------'
if (size(rho,3).gt.1) then
	write(*,*) 'GRAPE error!  You aren''t allowed to pass a rho w/ a time dimension'
else
	write(*,*) 'rho_0 of shape',size(rho,1),size(rho,2),size(rho,3),size(rho,4)
	write(*,*) 'propagators of shape',size(propagators,1),size(propagators,2),size(propagators,3),size(propagators,4)
	write(*,*) 'gradient of shape',size(gradient,1),size(gradient,2)
	t_steps = size(propagators,3)
	do p_counter = 1,size(propagators,4)
		if (size(C,3).gt.1) then
			call sandwich_internal( conjg(transpose(C(:,:,p_counter))),&
                conjg(transpose(propagators(:,:,size(propagators,3),&
                p_counter))),end_sandwich)
				!propagators(:,:,size(propagators,3),p_counter),&
		else
			call sandwich_internal(conjg(transpose(C(:,:,1))),&
                conjg(transpose(propagators(:,:,size(propagators,3),p_counter))),&
                end_sandwich)
				!propagators(:,:,size(propagators,3),p_counter),&
		endif
		do t_counter = 1,t_steps
			call sandwich_internal(H_c,conjg(transpose(propagators(:,:,t_counter,p_counter))),tempmat)
            !propagators(:,:,t_counter,p_counter),&
			if (size(rho,3).gt.1) then
				call commutator_internal(tempmat,rho(:,:,1,p_counter))
			else
				call commutator_internal(tempmat,rho(:,:,1,1))
			endif
			tempmat = matmul(end_sandwich,tempmat)
			call trace_internal(tempmat,gradient(t_counter,p_counter))
		enddo
	enddo
	gradient = gradient * (0,1)
endif
write(*,*) '------------------------------------------------------------'
end subroutine grape

subroutine complexexp(input)
implicit none
! finding the eigenvalues of a complex matrix using LAPACK
!	Implicit none
! declarations, notice double precision
double complex, dimension(:,:):: input
double precision eigenvalues(size(input,1))
double complex neweig(size(input,1))
double complex DUMMY(1,1), WORK(size(input,1)*4)
double complex RWORK(3*size(input,1)-2)
double complex a(size(input,1),size(input,2))
double precision :: pi = 3.14159265358979323846264338327950288419716939937510
integer i, ok
!write (*,*) 'input to fortran'
!write (*,*) input
call zheev('V','L',size(input,1), input, size(input,1), eigenvalues, WORK, 4*size(input,1), RWORK, ok) ! diagonalize hermitian matrix with lapack routine zheev
if (ok .ne. 0) then
	write (*,*) "An error occured with the diagonalization"
endif
a = input
neweig = exp((0.,1.D0) * (2.D0,0.) * pi * eigenvalues)
do i=1,size(input,1) ! multiply the rows of the eigenvalues by the diagonals
	a(:,i) = input(:,i)*neweig(i)
end do
input = matmul(a,transpose(conjg(input)))
!write (*,*) 'output from fortran'
!write (*,*) input
return
end subroutine complexexp

subroutine exact(input)
implicit none
complex*16 input(2,2),output(2,2),a,b,c,d,norm
complex*16, parameter:: sigma_z(2,2)=reshape((/1,0,0,-1/),(/2,2/))
complex*16, parameter:: sigma_y(2,2)=reshape((/(0,0),(0,1),(0,-1),(0,0)/),(/2,2/))
complex*16, parameter:: sigma_x(2,2)=reshape((/0,1,1,0/),(/2,2/))
complex*16, parameter:: identity(2,2)=reshape((/1,0,0,1/),(/2,2/))
double precision :: pi = 3.14159265358979323846264338327950288419716939937510
input=input*(0,1)*2.D0*pi
a=(input(1,1)-input(2,2))/2
b=(input(2,1)-input(1,2))/(0,2)
c=(input(2,1)+input(1,2))/2
d=(input(1,1)+input(2,2))/2
norm=sqrt(a**2+b**2+c**2)
if(norm.ne.0) then
	output=sigma_z*a/norm+sigma_y*b/norm+sigma_x*c/norm
else
	output=0
endif
output=(identity*cos(norm)-(0,1)*output*sin(norm))*exp((0,-1)*d)
input=output
end subroutine exact

subroutine sandwich_internal(mat,with,output)
double complex, intent(in), dimension(:,:):: mat
double complex, intent(in), dimension(:,:):: with
double complex, intent(out), dimension(:,:):: output
output = matmul(with,matmul(mat , conjg(transpose(with))))
return
end subroutine sandwich_internal

end module propagator
