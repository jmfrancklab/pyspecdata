      subroutine venk_brd(initial_alpha,k0,m,n,mr,f,alpha_out,tol,maxiter)
!f2py threadsafe
      integer,intent(in) :: m,n,maxiter
      double precision,intent(in) :: initial_alpha,tol
      double precision,intent(in) :: k0(m,n),mr(m)
      double precision,intent(out) :: f(n),alpha_out
      double precision,allocatable :: c(:),c_new(:),grad(:),newgrad(:)
      double precision,allocatable :: gmat(:,:),h(:,:),hd(:,:),k0t(:)
      double precision,allocatable :: tempvec(:)
      double precision :: alpha,alpha_new,sqrt_n
      double precision :: chi_old,chi_new,s,denom
      integer :: iter,j,k,info
      double precision :: norm_grad,norm_mr
      sqrt_n = sqrt(dble(m))
      alpha = 1.0d-3
      allocate(c(m),c_new(m),grad(m),newgrad(m),gmat(m,m),h(m,m),hd(m,1),k0t(n),tempvec(m))
      c = 1.0d0
      norm_mr = sqrt(sum(mr*mr))
      do iter=1,maxiter
         do j=1,100
            call compute_g(k0,m,n,c,gmat,k0t)
            call compute_grad(gmat,c,m,alpha,mr,grad)
            call compute_h(gmat,m,alpha,h)
            call solve_linear(h,grad,hd(:,1),m)
            tempvec = matmul(h,hd(:,1))
            denom = dot_product(hd(:,1),tempvec)
            if (denom == 0.d0) denom = 1.d-12
            s = dot_product(hd(:,1),grad) / denom
            c_new = c - s*hd(:,1)
            call compute_grad(gmat,c_new,m,alpha,mr,newgrad)
            chi_old = chi_func(c,gmat,alpha,mr)
            chi_new = chi_func(c_new,gmat,alpha,mr)
            k=0
            do while (chi_new >= chi_old .and. k < 20)
               k = k + 1
               c_new = c - s*(0.5d0**k)*hd(:,1)
               call compute_grad(gmat,c_new,m,alpha,mr,newgrad)
               chi_new = chi_func(c_new,gmat,alpha,mr)
            end do
            c = c_new
            grad = newgrad
            norm_grad = sqrt(sum(grad*grad))
            if (norm_grad/norm_mr < 1.d-8) exit
         end do
         alpha_new = sqrt_n / sqrt(sum(c*c))
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
      deallocate(c,c_new,grad,newgrad,gmat,h,hd,k0t,tempvec)
      alpha_out = alpha_new
      return
      contains
         subroutine compute_g(k0,m,n,c,gmat,k0t)
            integer,intent(in)::m,n
            double precision,intent(in)::k0(m,n),c(m)
            double precision,intent(out)::gmat(m,m),k0t(n)
            integer :: p
            gmat = 0.d0
            k0t = 0.d0
            do p=1,n
               k0t(p) = dot_product(k0(:,p),c)
               if (k0t(p) > 0.d0) then
                  gmat = gmat + outer_prod(k0(:,p))
               end if
            end do
         end subroutine compute_g
         function outer_prod(col) result(mat)
            double precision,intent(in)::col(:)
            double precision :: mat(size(col),size(col))
            integer :: i,j
            do i=1,size(col)
               do j=1,size(col)
                  mat(i,j) = col(i)*col(j)
               end do
            end do
         end function outer_prod
         subroutine compute_grad(gmat,c,m,alpha,mr,grad)
            integer,intent(in)::m
            double precision,intent(in)::gmat(m,m),c(m),alpha,mr(m)
            double precision,intent(out)::grad(m)
            grad = matmul(gmat,c) + alpha*c - mr
         end subroutine compute_grad
        subroutine compute_h(gmat,m,alpha,h)
            integer,intent(in)::m
            double precision,intent(in)::gmat(m,m),alpha
            double precision,intent(out)::h(m,m)
            integer::i
            h = gmat
            do i=1,m
               h(i,i) = h(i,i) + alpha
            end do
        end subroutine compute_h

        subroutine solve_linear(a,b,x,n)
            integer,intent(in) :: n
            double precision,intent(in) :: a(n,n),b(n)
            double precision,intent(out) :: x(n)
            double precision :: aa(n,n), bb(n)
            integer :: i,j,k
            aa = a
            bb = b
            do k=1,n-1
                do i=k+1,n
                    if (aa(k,k) /= 0.d0) then
                        bb(i) = bb(i) - aa(i,k)/aa(k,k)*bb(k)
                        do j=k+1,n
                            aa(i,j) = aa(i,j) - aa(i,k)/aa(k,k)*aa(k,j)
                        end do
                    end if
                end do
            end do
            do i=n,1,-1
                x(i) = bb(i)
                do j=i+1,n
                    x(i) = x(i) - aa(i,j)*x(j)
                end do
                if (aa(i,i) /= 0.d0) then
                    x(i) = x(i)/aa(i,i)
                else
                    x(i) = 0.d0
                end if
            end do
        end subroutine solve_linear
         function chi_func(c,gmat,alpha,mr) result(val)
            double precision,intent(in)::c(:),gmat(size(c),size(c)),alpha,mr(:)
            double precision::val
            val = dot_product(c,matmul(gmat,c)) + alpha*dot_product(c,c) - dot_product(c,mr)
         end function chi_func
      end subroutine venk_brd
