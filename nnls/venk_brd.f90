      subroutine venk_brd(initial_alpha,k0,m,n,mr,f,alpha_out,tol,maxiter)
!f2py threadsafe
! implement the Butler–Reeds–Dawson solver in Fortran
      integer,intent(in) :: m,n,maxiter
      double precision,intent(in) :: initial_alpha,tol
      double precision,intent(in) :: k0(m,n),mr(m)
      double precision,intent(out) :: f(n),alpha_out
      double precision,allocatable :: c(:),c_new(:),grad(:),newgrad(:)
      double precision,allocatable :: delta(:)
      double precision,allocatable :: G(:,:),H(:,:),K0t(:)
      double precision,allocatable :: tempvec(:)
      integer,allocatable :: piv(:)
      double precision :: alpha,alpha_new,sqrt_n
      double precision :: chi_old,chi_new,s,denom
      integer :: iter,j,k,info
      external dgesv
      double precision :: norm_grad,norm_mr
      sqrt_n = sqrt(dble(m))
      alpha = initial_alpha
      allocate(c(m),c_new(m),grad(m),newgrad(m),delta(m),G(m,m),H(m,m),K0t(n),tempvec(m))
      allocate(piv(m))
      c = 1.0d0
      norm_mr = norm2(mr)
      do iter=1,maxiter
        do j=1,100
            ! IN: k0,m,n,c OUT: G,K0t
            call compute_g(k0,m,n,c,G,K0t)
            ! IN: G,c,m,alpha,mr OUT: grad
            call compute_grad(G,c,m,alpha,mr,grad)
            ! IN: G,m,alpha OUT: H
            call add_diag(G,m,alpha,H)
            delta = grad
            ! INOUT: H,delta; IN: m,piv; OUT: info
            call dgesv(m,1,H,m,piv,delta,m,info)
            tempvec = matmul(H,delta)
            denom = dot_product(delta,tempvec)
            if (denom == 0.d0) denom = 1.d-12
            s = dot_product(delta,grad) / denom
            c_new = c - s*delta
            ! IN: G,c_new,m,alpha,mr OUT: newgrad
            call compute_grad(G,c_new,m,alpha,mr,newgrad)
            ! IN: c,G,alpha,mr OUT: chi_old
            call chi_func(c,G,alpha,mr,chi_old)
            ! IN: c_new,G,alpha,mr OUT: chi_new
            call chi_func(c_new,G,alpha,mr,chi_new)
            k=0
            do while (chi_new >= chi_old .and. k < 20)
               k = k + 1
               c_new = c - s*(0.5d0**k)*delta
               ! IN: G,c_new,m,alpha,mr OUT: newgrad
               call compute_grad(G,c_new,m,alpha,mr,newgrad)
               ! IN: c_new,G,alpha,mr OUT: chi_new
               call chi_func(c_new,G,alpha,mr,chi_new)
            end do
            c = c_new
            grad = newgrad
            norm_grad = norm2(grad)
            if (norm_grad/norm_mr < 1.d-8) exit
         end do
         alpha_new = sqrt_n / norm2(c)
         write(*,*) 'alpha iteration', iter, 'value', alpha_new
         if (abs(alpha_new-alpha)/alpha < tol) exit
         alpha = alpha_new
      end do
      tempvec = matmul(transpose(k0),c)
      do j=1,n
         if (tempvec(j) > 0.d0) then
            f(j) = tempvec(j)
         else
            f(j) = 0.d0
         end if
      end do
      deallocate(c,c_new,grad,newgrad,delta,G,H,K0t,tempvec,piv)
      alpha_out = alpha_new
      return
      contains
         subroutine compute_g(k0,m,n,c,G,k0t)
            integer,intent(in)::m,n
            double precision,intent(in)::k0(m,n),c(m)
            double precision,intent(out)::G(m,m),k0t(n)
            integer :: p
            G = 0.d0
            k0t = 0.d0
            do p=1,n
               k0t(p) = dot_product(k0(:,p),c)
               if (k0t(p) > 0.d0) then
                  ! IN: k0(:,p) INOUT: G
                  call add_outer(k0(:,p),G)
               end if
            end do
         end subroutine compute_g
         subroutine add_outer(col,Mat)
            double precision,intent(in)::col(:)
            double precision,intent(inout)::Mat(:,:)
            integer :: i,j
            do i=1,size(col)
               do j=1,size(col)
                  Mat(i,j) = Mat(i,j) + col(i)*col(j)
               end do
            end do
         end subroutine add_outer
         subroutine compute_grad(G,c,m,alpha,mr,grad)
            integer,intent(in)::m
            double precision,intent(in)::G(m,m),c(m),alpha,mr(m)
            double precision,intent(out)::grad(m)
            grad = matmul(G,c) + alpha*c - mr
         end subroutine compute_grad
        subroutine add_diag(G,m,alpha,H)
            integer,intent(in)::m
            double precision,intent(in)::G(m,m),alpha
            double precision,intent(out)::H(m,m)
            integer::i
            H = G
            do i=1,m
               H(i,i) = H(i,i) + alpha
            end do
        end subroutine add_diag

         subroutine chi_func(c,G,alpha,mr,val)
            double precision,intent(in)::c(:),G(size(c),size(c)),alpha,mr(:)
            double precision,intent(out)::val
            val = dot_product(c,matmul(G,c)) + alpha*dot_product(c,c) - dot_product(c,mr)
         end subroutine chi_func
      end subroutine venk_brd
